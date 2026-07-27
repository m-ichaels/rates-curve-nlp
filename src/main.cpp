// rcmm - yield-curve market making: curve, factors, RFQ simulator, closed-form quoters, events.
//   rcmm curve     --date D                         bootstrap vs kernel ridge, LOO errors, key-rate DV01s, step front end
//   rcmm curvehist --from D0 --to D1                daily repricing errors, LOO 5y, 10y DV01 stability, both methods
//   rcmm factors   --windows W1,W2                  PCA level/slope/curvature per window
//   rcmm mm        --train D0:D1 [--gammas ...]     closed-form quoters on synthetic days: frontier, decomposition
//   rcmm events    --features F --from D0 --to D1   replay real days with event/text features, adverse selection, placebo
#include "rcmm/sim.hpp"
#include "rcmm/json.hpp"
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
using namespace rcmm;

struct Args { std::map<std::string, std::string> kv; std::vector<std::string> pos;
    std::string get(const std::string& k, const std::string& d = "") const { auto it = kv.find(k); return it == kv.end() ? d : it->second; }
    double num(const std::string& k, double d) const { auto it = kv.find(k); return it == kv.end() ? d : std::atof(it->second.c_str()); }
    bool has(const std::string& k) const { return kv.count(k) > 0; }
    std::vector<std::string> list(const std::string& k, const std::vector<std::string>& d) const { if (!has(k)) return d; std::vector<std::string> v; std::string cur; for (char c : get(k)) { if (c == ',') { if (!cur.empty()) v.push_back(cur); cur.clear(); } else cur += c; } if (!cur.empty()) v.push_back(cur); return v; }
    std::vector<double> nums(const std::string& k, const std::vector<double>& d) const { if (!has(k)) return d; std::vector<double> v; for (auto& s : list(k, {})) v.push_back(std::atof(s.c_str())); return v; } };
static Args parse_args(int argc, char** argv) { Args a; for (int i = 2; i < argc; ++i) { std::string s = argv[i]; if (s.rfind("--", 0) == 0) { std::string k = s.substr(2); if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) a.kv[k] = argv[++i]; else a.kv[k] = "1"; } else a.pos.push_back(s); } return a; }
static void logmsg(const std::string& s) { std::cerr << "[rcmm] " << s << std::endl; }
static const CmtRow* find_row(const std::vector<CmtRow>& rows, const std::string& d) { const CmtRow* best = nullptr; for (auto& r : rows) if (r.date <= d) best = &r; return best; }
static std::pair<std::string, std::string> split_window(const std::string& w) { auto c = w.find(':'); if (c == std::string::npos) throw Error("window must be D0:D1"); return {w.substr(0, c), w.substr(c + 1)}; }

// ladder DV01s ($ per bp per $M) from a parallel 1 bp bump of the bootstrapped curve on a date
static LadderSpec ladder_from(const CmtRow& r) {
    auto inst = cmt_instruments(r); Bootstrap b; b.fit(inst); std::vector<Instrument> up = inst; for (auto& i : up) i.rate += 1e-4; Bootstrap bu; bu.fit(up); LadderSpec ls;
    for (double T : ls.T) { Instrument bond{InstKind::ParBond, T, b.par_yield(T), 2}; ls.dv01.push_back((price_instrument(bond, [&](double t) { return b.P(t); }) - price_instrument(bond, [&](double t) { return bu.P(t); })) * 1e6); }
    return ls;
}

// ---- curve --------------------------------------------------------------------------------------------
static int cmd_curve(const Args& a) {
    auto rows = read_cmt(a.get("cmt", "data/raw/treasury_cmt.csv"));
    const CmtRow* r = find_row(rows, a.get("date", rows.back().date)); if (!r) { std::cerr << "no such date\n"; return 1; }
    auto inst = cmt_instruments(*r); double lam = a.num("lambda", 1e-4);
    Bootstrap b; b.fit(inst); KernelRidge kr; kr.fit(inst, lam);
    Json out = Json::object(); out["date"] = r->date; out["bootstrap_iterations"] = b.iterations; out["bootstrap_max_reprice_error_bp"] = b.max_error_bp; out["kernel_ridge_max_reprice_error_bp"] = kr.max_error_bp(inst); out["kernel_ridge_lambda"] = lam;
    Json ins = Json::array(); for (auto& i : inst) { Json e = Json::object(); e["T"] = i.T; e["par"] = i.rate; e["zero_boot"] = b.zero(i.T); e["fwd_boot"] = b.mc.forward(i.T); e["zero_kr"] = -std::log(kr.P(i.T)) / i.T; e["fwd_kr"] = kr.forward(i.T); ins.push(e); } out["instruments"] = ins;
    Json grid = Json::array(); for (double t = 0.05; t <= 30.0001; t += 0.05) { Json e = Json::object(); e["t"] = t; e["fwd_boot"] = b.mc.forward(t); e["fwd_kr"] = kr.forward(t); e["zero_boot"] = b.zero(t); e["zero_kr"] = -std::log(kr.P(t)) / t; grid.push(e); } out["grid"] = grid;
    out["roughness_boot"] = forward_roughness([&](double t) { return b.mc.forward(t); }, 30.0); out["roughness_kr"] = forward_roughness([&](double t) { return kr.forward(t); }, 30.0);
    Json loo = Json::array();
    for (size_t k = 0; k < inst.size(); ++k) { std::vector<Instrument> sub; for (size_t j = 0; j < inst.size(); ++j) if (j != k) sub.push_back(inst[j]); Bootstrap bb; bb.fit(sub); KernelRidge kk; kk.fit(sub, lam);
        Json e = Json::object(); e["T"] = inst[k].T; e["loo_err_bp_boot"] = (price_instrument(inst[k], [&](double t) { return bb.P(t); }) - 1) * 1e4; e["loo_err_bp_kr"] = (price_instrument(inst[k], [&](double t) { return kk.P(t); }) - 1) * 1e4; loo.push(e); }
    out["leave_one_out"] = loo;
    Json kr01 = Json::array(); for (double T : {2.0, 5.0, 10.0, 30.0}) { double c = b.par_yield(T); auto v = key_rate_dv01(inst, T, c); Json e = Json::object(); e["T"] = T; e["coupon"] = c; e["krdv01_per_1"] = Json(v); double tot = 0; for (double x : v) tot += x; e["dv01_per_1"] = tot; kr01.push(e); } out["ladder_krdv01"] = kr01;
    // meeting-date step front end: SOFR futures when given (--futures file: T_start,T_end,rate per line), else bills as the proxy
    { std::vector<double> meetings = a.nums("meetings", {0.12, 0.25, 0.37, 0.5, 0.62, 0.75, 0.87, 1.0}); std::vector<StepFrontEnd::Obs> obs;
      if (a.has("futures")) { for (auto& line : split_lines(read_file(a.get("futures")))) { auto c = split_csv(line); if (c.size() >= 3) obs.push_back({std::atof(c[0].c_str()), std::atof(c[1].c_str()), std::atof(c[2].c_str())}); } }
      else for (auto& i : inst) if (i.kind == InstKind::Bill) obs.push_back({0.0, i.T, std::log(1 + i.rate * i.T) / i.T});
      double sofr = a.num("sofr", inst[0].rate); StepFrontEnd fe; fe.fit(meetings, obs, sofr);
      Json f = Json::object(); f["meetings"] = Json(meetings); f["levels"] = Json(fe.level); f["source"] = a.has("futures") ? "SOFR futures" : "bills + fixing (proxy)"; out["front_end"] = f; }
    write_file(a.get("out", "results/curve_" + r->date + ".json"), out.dump(1));
    std::cout << "date " << r->date << " bootstrap max reprice " << b.max_error_bp << " bp (" << b.iterations << " it); KR max " << kr.max_error_bp(inst) << " bp; roughness boot " << out["roughness_boot"].num() << " kr " << out["roughness_kr"].num() << std::endl;
    for (auto& e : loo.arr()) std::cout << "  LOO T=" << e["T"].num() << " boot " << e["loo_err_bp_boot"].num() << " bp, kr " << e["loo_err_bp_kr"].num() << " bp\n";
    return 0;
}

// Daily curve builds: out-of-sample comparison of the bootstrap and the kernel ridge (lambda grid).
//   leave-one-out pricing error on the interior tenors (6m..10y), in price bp and yield bp
//   10y par DV01 day-to-day stability, forward roughness, max repricing error
//   hedge test: a 7y (3y) par bond priced off the curve built WITHOUT the 7y (3y) point, hedged with the remaining inputs
//   using that curve's key-rate DV01s; the realised next-day P&L of the hedged package (yield bp) measures how well the
//   interpolation captures the co-movement of an off-curve bond with its neighbours.
static int cmd_curvehist(const Args& a) {
    auto rows = read_cmt(a.get("cmt", "data/raw/treasury_cmt.csv")); std::string d0 = a.get("from", "2015-01-01"), d1 = a.get("to", "2030-01-01");
    std::vector<double> lams = a.nums("lambdas", {1e-2, 1e-3, 1e-4, 1e-5, 1e-6}); std::vector<double> hedge_T = a.nums("hedge", {3.0, 7.0});
    struct Curve { std::function<double(double)> P; std::function<double(double)> fwd; };
    auto fit_boot = [](const std::vector<Instrument>& in) { auto b = std::make_shared<Bootstrap>(); b->fit(in); return Curve{[b](double t) { return b->P(t); }, [b](double t) { return b->mc.forward(t); }}; };
    auto fit_kr = [](const std::vector<Instrument>& in, double lam) { auto k = std::make_shared<KernelRidge>(); k->fit(in, lam); return Curve{[k](double t) { return k->P(t); }, [k](double t) { return k->forward(t); }}; };
    size_t M = 1 + lams.size();   // method 0 = bootstrap, 1.. = kernel ridge per lambda
    auto fit = [&](size_t m, const std::vector<Instrument>& in) { return m == 0 ? fit_boot(in) : fit_kr(in, lams[m - 1]); };
    auto annuity = [](const Instrument& ins, const Curve& c) { double s = 0; if (ins.kind == InstKind::ParBond) { int n = (int)std::lround(ins.T * ins.freq); for (int k = 1; k <= n; ++k) s += c.P((double)k / ins.freq) / ins.freq; } else s = ins.T * c.P(ins.T); return s; };   // dP/dy per unit of yield
    Json days = Json::array(); int n = 0;
    std::vector<std::vector<double>> loo_px(M), loo_y(M), dv01(M), rough(M), maxerr(M);
    std::map<double, std::vector<std::vector<double>>> hedged;
    for (double T : hedge_T) hedged[T].assign(M, {});
    std::map<double, std::vector<double>> unhedged; std::map<double, std::vector<std::vector<double>>> hw;   // hw[T][m] = hedge weights on the ex-T inputs, kept for the next day
    std::map<double, std::vector<double>> prev_y; std::vector<double> prev_rates; std::vector<double> prev_T;
    for (auto& r : rows) { if (r.date < d0 || r.date > d1) continue; auto inst = cmt_instruments(r); if (inst.size() < 10) continue;
        std::sort(inst.begin(), inst.end(), [](const Instrument& x, const Instrument& y) { return x.T < y.T; });
        // realised next-day hedge P&L from yesterday's weights
        if (!prev_rates.empty() && prev_rates.size() == inst.size()) for (double T : hedge_T) if (hw.count(T)) { size_t iT = 0; for (size_t i = 0; i < inst.size(); ++i) if (std::fabs(inst[i].T - T) < 1e-9) iT = i;
            double dyT = (inst[iT].rate - prev_rates[iT]) * 1e4; unhedged[T].push_back(dyT);
            for (size_t m = 0; m < M; ++m) { double h = 0; size_t k = 0; for (size_t i = 0; i < inst.size(); ++i) if (i != iT) { h += hw[T][m][k] * (inst[i].rate - prev_rates[i]) * 1e4; ++k; } hedged[T][m].push_back(dyT - h); } }
        Json e = Json::object(); e["date"] = r.date;
        for (size_t m = 0; m < M; ++m) {
            Curve c = fit(m, inst); double me = 0; for (auto& ins : inst) me = std::max(me, std::fabs(price_instrument(ins, c.P) - 1) * 1e4); maxerr[m].push_back(me);
            Instrument bond{InstKind::ParBond, 10.0, 0.0, 2}; { double lo = 0, hi = 0.2; for (int it = 0; it < 60; ++it) { bond.rate = 0.5 * (lo + hi); (price_instrument(bond, c.P) > 1 ? hi : lo) = bond.rate; } }
            std::vector<Instrument> up = inst; for (auto& i : up) i.rate += 1e-4; Curve cu = fit(m, up); dv01[m].push_back((price_instrument(bond, c.P) - price_instrument(bond, cu.P)) * 1e6);
            rough[m].push_back(forward_roughness(c.fwd, 30.0));
            double lp = 0, ly = 0; int nl = 0;
            for (size_t k = 0; k < inst.size(); ++k) { if (inst[k].T < 0.4 || inst[k].T > 10.0) continue; std::vector<Instrument> sub; for (size_t j = 0; j < inst.size(); ++j) if (j != k) sub.push_back(inst[j]); Curve cs = fit(m, sub);
                double pe = (price_instrument(inst[k], cs.P) - 1) * 1e4; lp += std::fabs(pe); ly += std::fabs(pe) / std::max(1e-9, annuity(inst[k], cs)); ++nl; }   // price bp / (dP/dy per unit yield) = yield bp
            loo_px[m].push_back(lp / std::max(nl, 1)); loo_y[m].push_back(ly / std::max(nl, 1));
            // hedge weights for the next day: key-rate DV01s of the T bond off the ex-T curve, normalised by the bond's own DV01
            for (double T : hedge_T) { size_t iT = 0; for (size_t i = 0; i < inst.size(); ++i) if (std::fabs(inst[i].T - T) < 1e-9) iT = i; std::vector<Instrument> sub; for (size_t i = 0; i < inst.size(); ++i) if (i != iT) sub.push_back(inst[i]);
                Curve cs = fit(m, sub); Instrument bT = inst[iT]; double p0 = price_instrument(bT, cs.P), own = annuity(bT, cs) * 1e-4; std::vector<double> w;
                for (size_t k = 0; k < sub.size(); ++k) { std::vector<Instrument> b2 = sub; b2[k].rate += 1e-4; Curve cb = fit(m, b2); w.push_back((p0 - price_instrument(bT, cb.P)) / own); }   // yield-bp of the T bond per bp of input k
                if (!hw.count(T)) hw[T].assign(M, {}); hw[T][m] = w; if (m <= 3) e["hw_" + std::to_string((int)T) + "y_" + (m == 0 ? std::string("boot") : "kr_" + std::to_string(m - 1))] = Json(w); }
            std::string tag = m == 0 ? "boot" : "kr_" + std::to_string(m - 1); e["dv01_10y_" + tag] = dv01[m].back(); e["loo_px_" + tag] = loo_px[m].back(); e["loo_y_" + tag] = loo_y[m].back(); e["reprice_" + tag] = me;
        }
        e["y10"] = r.y[10]; Json rates = Json::array(); for (auto& i : inst) rates.push(i.rate); e["rates"] = rates; days.push(e); ++n; prev_rates.clear(); for (auto& i : inst) prev_rates.push_back(i.rate);
        if (n % 200 == 0) logmsg("curvehist " + r.date); }
    auto dstd = [](const std::vector<double>& v) { std::vector<double> d; for (size_t i = 1; i < v.size(); ++i) d.push_back(v[i] - v[i - 1]); return stdev(d); };
    Json out = Json::object(); out["days"] = days; out["n"] = n; out["lambdas"] = Json(lams); Json methods = Json::array();
    for (size_t m = 0; m < M; ++m) { Json j = Json::object(); j["method"] = m == 0 ? "bootstrap" : "kernel_ridge"; j["lambda"] = m == 0 ? 0.0 : lams[m - 1]; j["max_reprice_bp"] = *std::max_element(maxerr[m].begin(), maxerr[m].end()); j["loo_price_bp_mean"] = mean(loo_px[m]); j["loo_yield_bp_mean"] = mean(loo_y[m]); j["loo_yield_bp_median"] = median(loo_y[m]);
        j["dv01_10y_daily_change_std"] = dstd(dv01[m]); j["roughness_mean"] = mean(rough[m]); Json hj = Json::object(); for (double T : hedge_T) { Json t = Json::object(); t["hedged_pnl_bp_std"] = stdev(hedged[T][m]); t["unhedged_bp_std"] = stdev(unhedged[T]); t["n"] = (long long)hedged[T][m].size(); hj[std::to_string((int)T) + "y"] = t; } j["hedge_test"] = hj;
        // per-year LOO (yield bp) for the time-series cross-validation in Python
        std::map<std::string, std::pair<double, int>> by; size_t di = 0; for (auto& d : days.arr()) { std::string y = d["date"].str().substr(0, 4); by[y].first += loo_y[m][di]; by[y].second++; ++di; } Json yj = Json::object(); for (auto& [y, v] : by) yj[y] = v.first / v.second; j["loo_yield_bp_by_year"] = yj;
        methods.push(j); }
    out["methods"] = methods; write_file(a.get("out", "results/curvehist.json"), out.dump(1));
    for (auto& j : methods.arr()) std::cout << j["method"].str() << " lambda " << j["lambda"].num() << ": LOO |err| " << j["loo_price_bp_mean"].num() << " price bp = " << j["loo_yield_bp_mean"].num() << " yield bp; max reprice " << j["max_reprice_bp"].num() << "; 10y DV01 day-to-day std " << j["dv01_10y_daily_change_std"].num() << "; roughness " << j["roughness_mean"].num() << "; hedged 7y std " << j["hedge_test"]["7y"]["hedged_pnl_bp_std"].num() << " (unhedged " << j["hedge_test"]["7y"]["unhedged_bp_std"].num() << ")\n";
    return 0;
}

static int cmd_factors(const Args& a) {
    auto rows = read_cmt(a.get("cmt", "data/raw/treasury_cmt.csv")); Json out = Json::object();
    for (auto& w : a.list("windows", {"2015-01-01:2026-12-31", "2023-01-01:2026-12-31", "2022-01-01:2022-12-31", "2020-01-01:2022-12-31"})) { auto [d0, d1] = split_window(w); FactorModel fm; fm.fit(rows, d0, d1); out[w] = fm.to_json();
        std::cout << w << ": " << fm.n_days << " days, explained " << fm.explained(1) << " / " << fm.explained(2) << " / " << fm.explained(3) << ", factor sd (bp/day) " << std::sqrt(fm.evals[0]) << " " << std::sqrt(fm.evals[1]) << " " << std::sqrt(fm.evals[2]) << std::endl; }
    write_file(a.get("out", "results/factors.json"), out.dump(1)); return 0;
}

// ---- market making -----------------------------------------------------------------------------------
static std::unique_ptr<MMQuoter> make_quoter(const std::string& n, const std::vector<ClientTier>& tiers, const LadderSpec& ls, const LadderRisk& risk, MMParams p) {
    if (n == "Symmetric") return std::make_unique<SymmetricQuoter>(tiers, p);
    if (n == "Bergault") return std::make_unique<BergaultQuoter>(tiers, ls, risk, p, false);
    if (n == "CarteaWang") return std::make_unique<BergaultQuoter>(tiers, ls, risk, p, true);
    if (n == "Bergault_oracle") return std::make_unique<BergaultQuoter>(tiers, ls, risk, p, false, true);
    if (n == "BarzykinCiceri") return std::make_unique<BarzykinCiceriQuoter>(tiers, ls, risk, p, false);
    if (n == "BarzykinCiceri_QAHR") return std::make_unique<BarzykinCiceriQuoter>(tiers, ls, risk, p, true);
    throw Error("unknown quoter " + n);
}
struct Setup { LadderSpec ls; FactorModel fm; LadderRisk risk; std::vector<ClientTier> tiers; };
static Setup setup_from(const Args& a, const std::string& w) {
    auto rows = read_cmt(a.get("cmt", "data/raw/treasury_cmt.csv")); Setup s; auto [d0, d1] = split_window(w);
    s.fm.fit(rows, d0, d1); const CmtRow* r = find_row(rows, d1); if (!r) throw Error("no CMT row before " + d1); s.ls = ladder_from(*r); s.risk.build(s.fm, s.ls); s.tiers = default_tiers();
    if (a.has("toxicity")) { double m = a.num("toxicity", 1); for (auto& t : s.tiers) { t.p_info = std::min(1.0, t.p_info * m); t.mu *= m; } }
    return s;
}
// the text nowcast is available from the release for release_window_min (default 30 min): never before the statement exists
static std::function<double(int, double, int)> release_signal(const std::vector<SimDay>& days, double window_min = 30, double day_hours = 6.5) { double w = window_min / 60.0 / day_hours; return [&days, w](int d, double t, int j) { for (auto& e : days[d].events) if (t >= e.t_day && t < e.t_day + w) return e.alpha_at((size_t)j); return 0.0; }; }
static Json hit_json(const SimResult& r) { Json h = Json::array(); for (size_t i = 0; i < r.hits_by_tier.size(); ++i) h.push(r.rfq_by_tier[i] > 0 ? r.hits_by_tier[i] / r.rfq_by_tier[i] : 0); return h; }

// frontier of the closed-form quoters on synthetic days, common random numbers across quoters
static int cmd_mm(const Args& a) {
    Setup S = setup_from(a, a.get("train", "2023-01-01:2026-12-31"));
    std::vector<double> gammas = a.nums("gammas", {5e-7, 1e-6, 2e-6, 5e-6, 1e-5}); int n_days = (int)a.num("days", 20), seeds = (int)a.num("seeds", 4);
    std::vector<std::string> names = a.list("quoters", {"Symmetric", "Bergault", "BarzykinCiceri", "BarzykinCiceri_QAHR", "CarteaWang", "Bergault_oracle"});
    Json cells = Json::array();
    for (auto& n : names) for (double g : gammas) {
        std::vector<double> pnl, sd, hv, resid; double sp = 0, inv = 0, hc = 0; std::vector<double> hits(S.tiers.size(), 0); Json mk, qsj;
        for (int s = 0; s < seeds; ++s) {
            MMParams p; p.gamma = g; Rng rng(100 + s); auto days = synthetic_days(n_days, S.risk, rng); auto q = make_quoter(n, S.tiers, S.ls, S.risk, p);
            MMSim sim(S.ls, S.risk, S.tiers, p, 1000 + s);
            SimResult r = sim.run(*q, days, n == "CarteaWang" ? release_signal(days) : nullptr);
            pnl.push_back(r.pnl / n_days); sd.push_back(r.pnl_5min_std); hv.push_back(r.hedged_dv01 / n_days); resid.push_back(r.pnl - (r.spread + r.inventory - r.hedge_cost)); sp += r.spread / n_days / seeds; inv += r.inventory / n_days / seeds; hc += r.hedge_cost / n_days / seeds;
            for (size_t i = 0; i < hits.size(); ++i) hits[i] += (r.rfq_by_tier[i] > 0 ? r.hits_by_tier[i] / r.rfq_by_tier[i] : 0) / seeds; if (s == 0) { mk = r.markouts; qsj = r.quote_structure; }
        }
        CI ci = bootstrap_mean_ci(pnl); Json e = Json::object(); e["quoter"] = n; e["gamma"] = g; e["pnl_per_day"] = ci.mean; e["pnl_lo"] = ci.lo; e["pnl_hi"] = ci.hi; e["pnl_5min_std"] = mean(sd); e["pnl_per_var"] = mean(sd) > 0 ? ci.mean / (mean(sd) * mean(sd)) : 0; e["hedged_dv01_per_day"] = mean(hv);
        Json dec = Json::object(); dec["spread"] = sp; dec["inventory"] = inv; dec["hedge_cost"] = hc; dec["max_abs_residual"] = *std::max_element(resid.begin(), resid.end(), [](double x, double y) { return std::fabs(x) < std::fabs(y); }); e["decomposition"] = dec; e["hit_ratio_by_tier"] = Json(hits); e["markouts"] = mk; e["quote_structure"] = qsj; cells.push(e);
        logmsg(n + " gamma " + std::to_string(g) + ": $" + std::to_string((int)ci.mean) + "/day, 5-min std " + std::to_string((int)mean(sd)) + ", P&L/var " + std::to_string(e["pnl_per_var"].num()));
    }
    Json out = Json::object(); out["cells"] = cells; out["ladder_dv01"] = Json(S.ls.dv01); out["days"] = n_days; out["seeds"] = seeds; Json sg = Json::array(); for (auto& r : S.risk.Sigma) sg.push(Json(r)); out["Sigma_ladder"] = sg;
    write_file(a.get("out", "results/mm_frontier.json"), out.dump(1)); return 0;
}

// ---- event replay on real days ------------------------------------------------------------------------
// features CSV: date,is_event,type,alpha_bp,info_mult (alpha = signal forecast of the 10y event move, from pre-timestamp documents only)
static int cmd_events(const Args& a) {
    auto rows = read_cmt(a.get("cmt", "data/raw/treasury_cmt.csv")); std::string d0 = a.get("from", "2024-01-01"), d1 = a.get("to", "2026-12-31");
    Setup S = setup_from(a, a.get("train", "2015-01-01:" + d0)); bool ex_covid = a.has("exclude-2020-22");
    if (ex_covid) { FactorModel fm; std::vector<CmtRow> sub; for (auto& r : rows) if (r.date < "2020-01-01" || r.date > "2022-12-31") sub.push_back(r); fm.fit(sub, "2015-01-01", d0); S.fm = fm; S.risk.build(S.fm, S.ls); }
    struct Feat { int is_event; std::string type; double alpha, info; std::vector<double> jumps, alphas; }; std::map<std::string, Feat> feat;
    if (a.has("features")) { std::map<std::string, size_t> col; bool first = true; for (auto& line : split_lines(read_file(a.get("features")))) { auto c = split_csv(line); if (first) { first = false; for (size_t i = 0; i < c.size(); ++i) col[c[i]] = i; continue; }
        auto g = [&](const std::string& k) { auto it = col.find(k); return it != col.end() && it->second < c.size() ? c[it->second] : std::string(); };
        Feat f; f.is_event = std::atoi(g("is_event").c_str()); f.type = g("type"); f.alpha = std::atof(g("alpha_bp").c_str()); f.info = std::atof(g("info_mult").c_str());
        for (auto k : {"alpha2y", "alpha5y", "alpha10y", "alpha30y"}) { std::string v = g(k); if (!v.empty()) f.alphas.push_back(std::atof(v.c_str())); } if (f.alphas.size() != 4) f.alphas.clear();
        for (auto k : {"jump2y", "jump5y", "jump10y", "jump30y"}) { std::string v = g(k); if (!v.empty()) f.jumps.push_back(std::atof(v.c_str())); }
        if (f.jumps.size() != 4) f.jumps.clear(); feat[g("date")] = f; } }
    std::vector<SimDay> days; const CmtRow* prev = nullptr; std::vector<int> ladder_idx = {6, 8, 10, 12}; double jump_frac = a.num("jump-frac", 0.6); size_t n_measured = 0;
    for (auto& r : rows) { if (r.date < d0 || r.date > d1) { prev = &r; continue; } if (!prev) { prev = &r; continue; } bool ok = true; for (int i : ladder_idx) if (std::isnan(r.y[i]) || std::isnan(prev->y[i])) ok = false; if (!ok) { prev = &r; continue; }
        SimDay day; day.date = r.date; for (int i : ladder_idx) day.close_change.push_back((r.y[i] - prev->y[i]) * 1e4);
        auto it = feat.find(r.date);
        if (it != feat.end() && it->second.is_event) { const Feat& f = it->second; Event e; e.day = (int)days.size(); e.type = f.type; e.t_day = (e.type == "FOMC") ? 0.7 : (e.type == "ECB" ? 0.15 : (e.type == "AUCTION" ? 0.6 : 0.05)); e.info_mult = f.info; e.alpha_bp = f.alpha; e.alpha = f.alphas;
            if (!f.jumps.empty()) { e.jump_bp = f.jumps; ++n_measured; } else for (double c : day.close_change) e.jump_bp.push_back(jump_frac * c);   // measured (USMPD window) jump when available, else a fraction of the daily move
            day.events.push_back(e); }
        days.push_back(day); prev = &r; }
    size_t ne = 0; for (auto& d : days) ne += d.events.size(); logmsg("event replay: " + std::to_string(days.size()) + " days, " + std::to_string(ne) + " event days (" + std::to_string(n_measured) + " with measured jumps)");
    std::vector<std::string> names = a.list("quoters", {"Symmetric", "Bergault", "BarzykinCiceri_QAHR", "CarteaWang", "Bergault_oracle"}); int seeds = (int)a.num("seeds", 3); double g = a.num("gamma", 2e-6);
    Json cells = Json::array();
    for (auto& n : names) for (int variant = 0; variant < 4; ++variant) {   // 0 text signal + calendar; 1 placebo: text permuted across events + calendar; 2 calendar only (alpha = 0); 3 neither
        bool uses_signal = n == "CarteaWang"; if (!uses_signal && variant > 0) break;
        std::vector<double> pnl, sd; Json mk; double sp = 0, inv = 0, hc = 0; std::vector<double> hits(S.tiers.size(), 0);
        for (int s = 0; s < seeds; ++s) {
            MMParams p; p.gamma = g; auto cq = make_quoter(n, S.tiers, S.ls, S.risk, p); MMQuoter* qp = cq.get(); MMSim sim(S.ls, S.risk, S.tiers, p, 9000 + s); sim.use_calendar = variant < 3;
            std::vector<std::pair<double, std::vector<double>>> alphas; for (auto& d : days) for (auto& e : d.events) alphas.push_back({e.alpha_bp, e.alpha}); Rng prng(31 + s); if (variant == 1) for (size_t i = alphas.size(); i > 1; --i) std::swap(alphas[i - 1], alphas[prng.index(i)]);   // placebo: the nowcasts permuted across events
            std::vector<SimDay> dd = days; size_t ai = 0; for (auto& d : dd) for (auto& e : d.events) { if (variant >= 2) { e.alpha_bp = 0.0; e.alpha.clear(); } else { e.alpha_bp = alphas[ai].first; e.alpha = alphas[ai].second; ++ai; } }
            SimResult r = sim.run(*qp, dd, uses_signal ? release_signal(dd) : nullptr);
            pnl.push_back(r.pnl / days.size()); sd.push_back(r.pnl_5min_std); sp += r.spread / days.size() / seeds; inv += r.inventory / days.size() / seeds; hc += r.hedge_cost / days.size() / seeds; for (size_t i = 0; i < hits.size(); ++i) hits[i] += (r.rfq_by_tier[i] > 0 ? r.hits_by_tier[i] / r.rfq_by_tier[i] : 0) / seeds; if (s == 0) mk = r.markouts;
        }
        CI ci = bootstrap_mean_ci(pnl); Json e = Json::object(); e["quoter"] = n; e["variant"] = variant == 0 ? "text+calendar" : (variant == 1 ? "placebo_text+calendar" : (variant == 2 ? "calendar_only" : "none")); e["pnl_per_day"] = ci.mean; e["pnl_lo"] = ci.lo; e["pnl_hi"] = ci.hi; e["pnl_5min_std"] = mean(sd); e["pnl_per_var"] = mean(sd) > 0 ? ci.mean / (mean(sd) * mean(sd)) : 0; Json dec = Json::object(); dec["spread"] = sp; dec["inventory"] = inv; dec["hedge_cost"] = hc; e["decomposition"] = dec; e["hit_ratio_by_tier"] = Json(hits); e["markouts"] = mk; cells.push(e);
        logmsg(n + " " + e["variant"].str() + ": $" + std::to_string((int)ci.mean) + "/day, 5-min std " + std::to_string((int)mean(sd)));
    }
    Json out = Json::object(); out["cells"] = cells; out["n_days"] = (long long)days.size(); out["n_event_days"] = (long long)ne; out["n_measured_jumps"] = (long long)n_measured; out["from"] = d0; out["to"] = d1; out["exclude_2020_22"] = ex_covid; out["gamma"] = g; out["features"] = a.get("features", ""); out["ladder_dv01"] = Json(S.ls.dv01);
    std::map<std::string, int> et; for (auto& d : days) for (auto& e : d.events) et[e.type]++; Json etj = Json::object(); for (auto& [k, v] : et) etj[k] = (long long)v; out["events_by_type"] = etj;
    write_file(a.get("out", "results/events.json"), out.dump(1)); return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: rcmm <curve|curvehist|factors|mm|events> [--key value ...]\n"; return 1; }
    std::string cmd = argv[1]; Args a = parse_args(argc, argv);
    try {
        if (cmd == "curve") return cmd_curve(a); if (cmd == "curvehist") return cmd_curvehist(a); if (cmd == "factors") return cmd_factors(a);
        if (cmd == "mm") return cmd_mm(a); if (cmd == "events") return cmd_events(a);
        std::cerr << "unknown command " << cmd << "\n"; return 1;
    } catch (const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 1; }
}
