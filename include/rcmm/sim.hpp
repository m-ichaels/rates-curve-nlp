#pragma once
// Event-driven RFQ simulator: factor-driven curve, tiered RFQ flow with informed pre-event trading,
// event jumps, post-trade drift, hedging in the futures at a cost, P&L decomposition, mark-outs.
#include "mm.hpp"

namespace rcmm {

struct SimDay {
    std::vector<Event> events;               // intraday events with their jumps (bp per tenor) and signal
    std::vector<double> close_change;        // optional realised close-to-close change per tenor (bp); empty = free diffusion
    std::string date;
};

struct SimResult {
    std::string quoter; double pnl = 0, spread = 0, inventory = 0, hedge_cost = 0;
    std::vector<double> pnl_5min; double pnl_5min_std = 0; double risk_var_mean = 0; double risk_penalty = 0;   // mean gamma-free q'Sigma q ($^2/day)
    size_t n_rfq = 0, n_fills = 0; std::vector<double> hits_by_tier, rfq_by_tier;
    double hedged_dv01 = 0, max_abs_q = 0;
    std::vector<Fill> fills;
    std::vector<double> factor_exposure_std;   // per factor, $ per unit factor
    Json markouts;                             // by tier x {pre_event, normal}
    Json quote_structure;                      // mean half-spread by tier x {risk-reducing, risk-adding side} x {normal, pre-event}: the policy shape
    double pnl_per_var() const { return pnl_5min_std > 0 ? pnl / (pnl_5min_std * pnl_5min_std) : 0; }
};

class MMSim {
public:
    MMSim(const LadderSpec& ls, const LadderRisk& risk, const std::vector<ClientTier>& tiers, MMParams p, uint64_t seed) : ls_(ls), risk_(risk), tiers_(tiers), p_(p), rng_(seed) {}
    // signal(day, t_day, tenor) -> expected yield change in bp (0 if none); event_window_min: informed window before events
    SimResult run(MMQuoter& quoter, const std::vector<SimDay>& days, std::function<double(int, double, int)> signal = nullptr, double event_window_min = 30) {
        SimResult R; R.quoter = quoter.name(); size_t n = ls_.T.size(); size_t nt = tiers_.size();
        quoter_hook_ = [&](const Fill& f) { quoter.on_markout(f); };
        double risk_pen = 0; cursor_ = 0; pending_markout_.clear();
        R.hits_by_tier.assign(nt, 0); R.rfq_by_tier.assign(nt, 0); std::vector<double> qs_sum(nt * 4, 0.0), qs_n(nt * 4, 0.0);
        std::vector<double> q(n, 0.0), n_pos(n, 0.0);   // q in $/bp, n in $M
        std::vector<double> y(n, 0.0), yi(n, 0.0);        // yield change since start, bp; yi = information component (client drift + event jumps)
        const double DAY = p_.day_hours * 3600;
        std::vector<double> esum(risk_.L.size(), 0.0), esq(risk_.L.size(), 0.0); size_t esamples = 0; double var_acc = 0;
        struct Drift { double t_end, rate; int tenor; }; std::vector<Drift> drifts;
        double last_pnl = 0;
        for (size_t d = 0; d < days.size(); ++d) {
            const SimDay& day = days[d];
            double t = 0; size_t ev_i = 0; double next_review = p_.hedge_review_min * 60, next_5 = 300;
            std::vector<double> y_open = y;
            // factor path on a one-minute grid drawn up front, so the curve is common across quoters (paired comparisons)
            int M = (int)std::ceil(DAY / 60.0); std::vector<std::vector<double>> cum(M + 1, std::vector<double>(n, 0.0));
            for (int m = 1; m <= M; ++m) { cum[m] = cum[m - 1]; for (size_t k = 0; k < risk_.L.size(); ++k) { double e = std::sqrt(risk_.lam[k] * 60.0 / DAY) * rng_.normal(); for (size_t j = 0; j < n; ++j) cum[m][j] += risk_.L[k][j] * e; } }
            // Brownian bridge to a realised close if given: extra drift each step
            auto step_curve = [&](double dt) {
                int m0 = std::min(M, (int)(t / 60.0)), m1 = std::min(M, (int)((t + dt) / 60.0));
                for (size_t j = 0; j < n; ++j) y[j] += cum[m1][j] - cum[m0][j];
                if (!day.close_change.empty()) { double rem = DAY - t; for (size_t j = 0; j < n; ++j) { double target = y_open[j] + day.close_change[j] - jumps_remaining(day, ev_i, j); y[j] += (target - y[j]) * dt / std::max(rem, dt); } }
                for (auto& dr : drifts) { double a = std::max(0.0, std::min(dt, dr.t_end - t)); y[dr.tenor] += dr.rate * a; yi[dr.tenor] += dr.rate * a; }
                drifts.erase(std::remove_if(drifts.begin(), drifts.end(), [&](const Drift& x) { return x.t_end <= t + dt; }), drifts.end());
            };
            auto pnl_now = [&]() { return R.spread + R.inventory - R.hedge_cost; };   // inventory is marked to the curve at every step
            auto mark_inventory = [&](const std::vector<double>& y_before) { for (size_t j = 0; j < n; ++j) R.inventory += -q[j] * (y[j] - y_before[j]); };
            while (t < DAY) {
                // intensities (per second) for each (tier, tenor, side), with the informed pre-event window
                const Event* next_ev = ev_i < day.events.size() ? &day.events[ev_i] : nullptr;
                bool pre = next_ev && (next_ev->t_day * DAY - t) <= event_window_min * 60 && next_ev->t_day * DAY > t;
                std::vector<double> lam; lam.reserve(nt * n * 2); double tot = 0;
                for (size_t i = 0; i < nt; ++i) for (size_t j = 0; j < n; ++j) for (int s = -1; s <= 1; s += 2) {
                    double l = tiers_[i].lambda / 3600.0;
                    if (pre && tiers_[i].p_info > 0) { double sgn = next_ev->jump_bp[j] < 0 ? +1 : -1; l *= next_ev->info_mult * ((s == sgn) ? (0.5 + tiers_[i].p_info / 2) * 2 : (0.5 - tiers_[i].p_info / 2) * 2); }   // client buys bonds when yields will fall
                    lam.push_back(l); tot += l;
                }
                double tau = rng_.exponential(tot);
                double t_ev = next_ev ? next_ev->t_day * DAY : 1e18;
                double t_next = std::min({t + tau, next_review, next_5, t_ev, DAY});
                risk_pen += p_.gamma * risk_.variance(q) * (t_next - t) / DAY;
                std::vector<double> yb = y; step_curve(t_next - t); mark_inventory(yb); t = t_next;
                process_markouts(R, t + d * DAY, y, yi);
                if (t >= DAY) break;
                if (t_next == t_ev) { std::vector<double> y0 = y; double W = p_.release_window_min * 60; for (size_t j = 0; j < n; ++j) { double J = next_ev->jump_bp[j]; y[j] += p_.jump_instant_frac * J; yi[j] += p_.jump_instant_frac * J; if (W > 0 && p_.jump_instant_frac < 1) drifts.push_back({t + W, (1 - p_.jump_instant_frac) * J / W, (int)j}); } mark_inventory(y0); ++ev_i; continue; }   // part of the move at the release, the rest over the discovery window
                if (t_next == next_5) { double pn = pnl_now(); R.pnl_5min.push_back(pn - last_pnl); last_pnl = pn; next_5 += 300; var_acc += risk_.variance(q); ++esamples; auto e = risk_.factor_exposure(q); for (size_t k = 0; k < e.size(); ++k) { esum[k] += e[k]; esq[k] += e[k] * e[k]; } continue; }
                if (t_next == next_review) { next_review += p_.hedge_review_min * 60; for (size_t j = 0; j < n; ++j) if (std::fabs(q[j]) > p_.hedge_band) { double target = q[j] > 0 ? p_.hedge_band : -p_.hedge_band; double dq = target - q[j]; R.hedge_cost += std::fabs(dq) * p_.hedge_cost_bp; R.hedged_dv01 += std::fabs(dq); q[j] = target; n_pos[j] = q[j] / ls_.dv01[j]; } continue; }
                // RFQ
                size_t k = rng_.weighted(lam); size_t i = k / (n * 2), j = (k / 2) % n; int side = (k % 2 == 0) ? -1 : +1;
                size_t si = rng_.weighted(ls_.size_w); double z = ls_.sizes[si];
                double sig = signal ? signal((int)d, t / DAY, (int)j) : 0.0;
                quoter.on_context(pre && use_calendar);
                double delta = quoter.quote((int)i, (int)j, side, z, q, sig, t / DAY);
                R.n_rfq++; R.rfq_by_tier[i] += 1;
                { double sq = 0; for (size_t b = 0; b < n; ++b) sq += risk_.Sigma[j][b] * q[b]; int red = (side * sq > 0) ? 0 : 1; size_t cell = i * 4 + (pre ? 2 : 0) + red; qs_sum[cell] += delta; qs_n[cell] += 1; }   // client buy (side +1) reduces risk when the dealer is long: side * (Sigma q)_j > 0
                bool hit = rng_.uniform() < tiers_[i].f(delta);
                if (auto* bc = dynamic_cast<BarzykinCiceriQuoter*>(&quoter)) bc->on_request_outcome((int)i, hit);
                if (!hit) continue;
                R.n_fills++; R.hits_by_tier[i] += 1;
                double dq = side * z * ls_.dv01[j];   // client buys => dealer short DV01
                q[j] -= dq; n_pos[j] -= side * z;
                R.spread += z * ls_.dv01[j] * delta;
                R.max_abs_q = std::max(R.max_abs_q, std::fabs(q[j]));
                if (tiers_[i].mu > 0) drifts.push_back({t + 1800.0, -side * tiers_[i].mu * (z / 10.0) / 1800.0, (int)j});   // client buys => yields fall
                Fill f; f.t = t + d * DAY; f.tier = (int)i; f.tenor = (int)j; f.side = side; f.z = z; f.delta = delta; f.yield_mid = y[j]; f.yinfo_mid = yi[j]; f.pre_event = pre;
                R.fills.push_back(f); quoter.on_fill(f);
            }
            { double pn = pnl_now(); if (pn != last_pnl) { R.pnl_5min.push_back(pn - last_pnl); last_pnl = pn; } }   // close of day: remainder after the last 5-minute mark
        }
        process_markouts(R, 1e18, y, yi, true);
        R.pnl = R.spread + R.inventory - R.hedge_cost; R.risk_penalty = risk_pen; R.pnl_5min_std = stdev(R.pnl_5min); R.risk_var_mean = esamples ? var_acc / esamples : 0;
        for (size_t k = 0; k < esum.size(); ++k) R.factor_exposure_std.push_back(esamples ? std::sqrt(std::max(0.0, esq[k] / esamples - std::pow(esum[k] / esamples, 2))) : 0);
        R.markouts = markout_table(R);
        { Json qs = Json::array(); for (size_t i = 0; i < nt; ++i) for (int c = 0; c < 4; ++c) { Json e = Json::object(); e["tier"] = tiers_[i].name; e["pre_event"] = c >= 2; e["risk_adding"] = (c % 2) == 1; e["n"] = (long long)qs_n[i * 4 + c]; e["mean_delta_bp"] = qs_n[i * 4 + c] > 0 ? qs_sum[i * 4 + c] / qs_n[i * 4 + c] : 0.0; qs.push(e); } R.quote_structure = qs; }
        return R;
    }
private:
    double jumps_remaining(const SimDay& day, size_t ev_i, size_t j) const { double s = 0; for (size_t e = ev_i; e < day.events.size(); ++e) s += day.events[e].jump_bp[j]; return s; }
    void process_markouts(SimResult& R, double t_now, const std::vector<double>& y, const std::vector<double>& yi, bool flush = false) {
        const double H[3] = {300, 1800, 3600};
        for (size_t k = cursor_; k < R.fills.size(); ++k) {
            Fill& f = R.fills[k]; size_t h = f.markouts.size();
            while (h < 3 && (t_now >= f.t + H[h] || flush)) { f.markouts.push_back(-f.side * (y[f.tenor] - f.yield_mid)); f.markouts_info.push_back(-f.side * (yi[f.tenor] - f.yinfo_mid)); ++h; }   // client buys: gains if yields fall
            if (h == 3) { if (k == cursor_) ++cursor_; pending_markout_.push_back(k); } else break;
        }
        // notify the quoter for completed fills (quality-adjusted hit ratio)
        for (size_t k : pending_markout_) if (quoter_hook_) quoter_hook_(R.fills[k]);
        pending_markout_.clear();
    }
    Json markout_table(const SimResult& R) const {
        Json a = Json::array();
        for (size_t i = 0; i < tiers_.size(); ++i) for (int pre = 0; pre < 2; ++pre) {
            double m[3] = {0, 0, 0}; size_t cnt = 0; double dsum = 0;
            for (auto& f : R.fills) if (f.tier == (int)i && (int)f.pre_event == pre && f.markouts.size() == 3) { ++cnt; for (int h = 0; h < 3; ++h) m[h] += f.markouts[h]; dsum += f.delta; }
            Json e = Json::object(); e["tier"] = tiers_[i].name; e["pre_event"] = pre == 1; e["n"] = (long long)cnt; e["mean_delta_bp"] = cnt ? dsum / cnt : 0;
            e["markout_5m_bp"] = cnt ? m[0] / cnt : 0; e["markout_30m_bp"] = cnt ? m[1] / cnt : 0; e["markout_60m_bp"] = cnt ? m[2] / cnt : 0; a.push(e);
        }
        return a;
    }
public:
    std::function<void(const Fill&)> quoter_hook_; bool use_calendar = true;   // pass the scheduled-event window to the quoter
private:
    LadderSpec ls_; LadderRisk risk_; std::vector<ClientTier> tiers_; MMParams p_; Rng rng_;
    size_t cursor_ = 0; std::vector<size_t> pending_markout_;
};

// Synthetic days: n_events per day on average, jump size from the factor model scaled by `event_sigma_mult`,
// signal = alpha = (signal_r2 correlation) * jump on the 10y + noise.
inline std::vector<SimDay> synthetic_days(int n_days, const LadderRisk& risk, Rng& rng, double p_event = 0.25, double event_sigma_mult = 1.0, double signal_corr = 0.5) {
    std::vector<SimDay> days; size_t n = risk.Sigma.size();
    for (int d = 0; d < n_days; ++d) {
        SimDay day; day.date = "synthetic-" + std::to_string(d);
        if (rng.uniform() < p_event) {
            Event e; e.day = d; e.t_day = 0.2 + 0.7 * rng.uniform(); e.type = "FOMC"; e.jump_bp.assign(n, 0.0);
            for (size_t k = 0; k < risk.L.size(); ++k) { double f = event_sigma_mult * std::sqrt(risk.lam[k]) * rng.normal(); for (size_t j = 0; j < n; ++j) e.jump_bp[j] += risk.L[k][j] * f; }
            double j10 = e.jump_bp[std::min<size_t>(2, n - 1)]; double sd = event_sigma_mult * std::sqrt(risk.Sigma[2][2]);
            e.alpha_bp = signal_corr * j10 + std::sqrt(1 - signal_corr * signal_corr) * sd * rng.normal();   // a signal with the stated correlation to the jump
            day.events.push_back(e);
        }
        days.push_back(day);
    }
    return days;
}

} // namespace rcmm
