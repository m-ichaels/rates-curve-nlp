#pragma once
// RFQ market making on a 2y / 5y / 10y / 30y ladder against the curve's fair value.
//
// Units: yields in bp, notional in $M, DV01_j = $ per bp per $M, inventory q_j = n_j DV01_j ($/bp),
// P&L in $.  Yield changes at the ladder come from the K-factor model (loadings interpolated to the
// ladder maturities); risk = gamma q' Sigma q per day with Sigma in bp^2/day.
//
// Quoters (half-spread in bp of yield per RFQ):
//   Symmetric      delta_i solving f_i(delta) = h* (target hit ratio), no skew
//   Bergault       Bergault-Evangelista-Gueant-Vieira (2021) closed-form multi-asset quotes with the
//                  factor reduction: Gamma D Gamma = (gamma/2) Sigma_K, skew_j = (Gamma q)_j,
//                  delta0_j = 1/k_j + Gamma_jj z_j / 2
//   BarzykinCiceri Barzykin-Ciceri (2026) linearisation: riskless spread + inventory correction +
//                  hit-ratio correction  kappa (h_hat - h*), with Niang's quality-adjusted hit ratio
//                  (a toxic hit counts 1 + toxicity, so the target is not met by subsidising toxic flow)
//   CarteaWang     Cartea-Wang (2020): Bergault quotes shifted by the alpha signal (expected yield change)
#include "factors.hpp"
#include <functional>
#include <deque>

namespace rcmm {

struct LadderSpec {
    std::vector<double> T = {2, 5, 10, 30};
    std::vector<double> dv01;            // $ per bp per $M, per tenor
    std::vector<double> sizes = {5, 10, 25};    // RFQ sizes, $M
    std::vector<double> size_w = {0.5, 0.35, 0.15};
};

struct ClientTier {
    std::string name; double lambda;      // RFQs per hour per tenor per side
    double alpha, beta;                   // hit prob f(delta) = 1/(1+exp(alpha + beta delta)), delta in bp
    double p_info;                        // probability the client knows the sign of the next event jump (pre-event window)
    double mu;                            // post-trade yield drift in the client's favour, bp per $10M DV01-equivalent
    double f(double d) const { return 1.0 / (1.0 + std::exp(alpha + beta * d)); }
    double d_riskless() const { double best = 0, bv = -1; for (int i = 1; i <= 300; ++i) { double d = 0.01 * i, v = d * f(d); if (v > bv) { bv = v; best = d; } } return best; }   // argmax delta f(delta)
    // exponential intensity A e^{-k delta} matching f and f' at the riskless optimum (so 1/k = d_riskless), for the closed forms
    void exp_fit(double& A, double& k) const { double d0 = d_riskless(); k = beta * (1 - f(d0)); A = f(d0) * std::exp(k * d0); }
};
inline std::vector<ClientTier> default_tiers() {
    return {{"informed", 0.4, -1.2, 4.0, 0.6, 0.20}, {"real_money", 0.8, -1.0, 3.0, 0.15, 0.05}, {"other", 1.2, -0.8, 2.5, 0.0, 0.0}};
}

struct Event { double t_day; int day; std::string type; std::vector<double> jump_bp; double alpha_bp = 0; std::vector<double> alpha; double info_mult = 3.0;   // jump per ladder tenor at time t_day (fraction of the day)
    double alpha_at(size_t j) const { return j < alpha.size() ? alpha[j] : alpha_bp; } };   // text nowcast of the window move per tenor (alpha_bp = same for all tenors), known at release

struct MMParams {
    double gamma = 2e-6;        // risk aversion, 1/$
    double hit_target = 0.35;   // target hit ratio
    double hedge_cost_bp = 0.15;// futures half bid-ask, bp of yield
    double hedge_band = 1e4;    // $/bp per tenor before hedging
    double hedge_review_min = 5;
    double kappa_hit = 1.0;     // Barzykin-Ciceri hit-ratio controller gain (bp per unit of hit-ratio miss)
    double signal_horizon = 0.5;// Cartea-Wang: fraction of the predicted move to skew by
    double signal_clip_bp = 1.0;// cap on the signal skew
    double event_widen_bp = 0.5;// calendar feature: extra half-spread for the informed tier inside the pre-event window (others get a quarter)
    double release_window_min = 30; // price discovery after a release: the text signal is available during this window (never before the release)
    double jump_instant_frac = 0.5; // share of the event-window move that happens at the release; the rest accrues linearly over the release window
    double day_hours = 6.5;
    int K = 3;
};

struct Fill { double t; int tier, tenor, side; double z, delta, yield_mid, yinfo_mid = 0; std::vector<double> markouts, markouts_info; bool pre_event = false; };   // mark-outs in bp in the client's favour at 5/30/60 min; _info = information component only (drift + event jumps)

class MMQuoter {
public:
    virtual ~MMQuoter() = default;
    virtual std::string name() const = 0;
    // half-spread (bp) for an RFQ; q: inventory in $/bp per tenor; side +1 client buys (dealer sells), -1 client sells
    virtual double quote(int tier, int tenor, int side, double z, const std::vector<double>& q, double signal_bp, double t_day) = 0;
    virtual void on_context(bool /*pre_event*/) {}   // called before quote(): scheduled-event window flag (calendar feature)
    virtual void on_fill(const Fill&) {}
    virtual void on_markout(const Fill&) {}
};

// shared machinery: factor covariance at the ladder, closed-form Gamma
struct LadderRisk {
    std::vector<std::vector<double>> Sigma;   // 4x4 bp^2/day (K-factor)
    std::vector<std::vector<double>> L;       // [K][4] loadings at ladder tenors
    std::vector<double> lam;                  // K eigenvalues
    void build(const FactorModel& fm, const LadderSpec& ls) {
        size_t n = ls.T.size(); L.assign(fm.K, std::vector<double>(n)); lam = std::vector<double>(fm.evals.begin(), fm.evals.begin() + fm.K);
        for (int k = 0; k < fm.K; ++k) for (size_t j = 0; j < n; ++j) L[k][j] = fm.loading_at(k, ls.T[j])[0];
        Sigma.assign(n, std::vector<double>(n, 0.0)); for (int k = 0; k < fm.K; ++k) for (size_t a = 0; a < n; ++a) for (size_t b = 0; b < n; ++b) Sigma[a][b] += lam[k] * L[k][a] * L[k][b];
    }
    double variance(const std::vector<double>& q) const { double v = 0; for (size_t a = 0; a < q.size(); ++a) for (size_t b = 0; b < q.size(); ++b) v += q[a] * Sigma[a][b] * q[b]; return v; }   // $^2 per day
    std::vector<double> factor_exposure(const std::vector<double>& q) const { std::vector<double> e(L.size(), 0.0); for (size_t k = 0; k < L.size(); ++k) for (size_t j = 0; j < q.size(); ++j) e[k] += L[k][j] * q[j]; return e; }   // $ per unit factor
};

class SymmetricQuoter : public MMQuoter {
public:
    SymmetricQuoter(const std::vector<ClientTier>& tiers, MMParams p) { for (auto& t : tiers) { double lo = 0, hi = 50; for (int i = 0; i < 100; ++i) { double m = 0.5 * (lo + hi); (t.f(m) > p.hit_target ? lo : hi) = m; } d0_.push_back(0.5 * (lo + hi)); } }
    std::string name() const override { return "Symmetric"; }
    double quote(int tier, int, int, double, const std::vector<double>&, double, double) override { return d0_[tier]; }
protected: std::vector<double> d0_;
};

class BergaultQuoter : public MMQuoter {
public:
    // oracle = true: the quote also charges each tier its expected adverse selection (post-trade drift mu_i z/10, and in the
    // pre-event window p_i E|jump_j| with E|jump| = event_sigma sqrt(Sigma_jj) sqrt(2/pi)); this is the closed form given the
    // true client model: the benchmark for every quoter that has to infer it
    BergaultQuoter(const std::vector<ClientTier>& tiers, const LadderSpec& ls, const LadderRisk& risk, MMParams p, bool use_signal = false, bool oracle = false, double event_sigma = 1.0) : ls_(ls), p_(p), signal_(use_signal), oracle_(oracle), tiers_(tiers) {
        for (size_t j = 0; j < ls.T.size(); ++j) jump_abs_.push_back(event_sigma * std::sqrt(std::max(0.0, risk.Sigma[j][j])) * std::sqrt(2.0 / 3.141592653589793));
        size_t n = ls.T.size(); A_.resize(tiers.size()); k_.resize(tiers.size());
        for (size_t i = 0; i < tiers.size(); ++i) { tiers[i].exp_fit(A_[i], k_[i]); A_[i] *= tiers[i].lambda * p.day_hours; }   // per day per tenor per side
        // quadratic (Riccati) approximation of the ergodic HJB: Gamma D Gamma = (gamma/2) Sigma with D_j = A_j z_j k_j / e,
        // A_j the total intensity per side per day, z_j the mean size in $/bp; Gamma in bp^2/$ so that the skew Gamma q is in bp
        double zbar = 0; for (size_t s = 0; s < ls.sizes.size(); ++s) zbar += ls.sizes[s] * ls.size_w[s];
        std::vector<double> D(n); double Asum = 0, ksum = 0; for (size_t i = 0; i < tiers.size(); ++i) { Asum += A_[i]; ksum += k_[i] * A_[i]; } double kbar = ksum / Asum;
        for (size_t j = 0; j < n; ++j) D[j] = zbar * ls.dv01[j] * Asum * kbar * std::exp(-1.0);
        // Gamma = D^-1/2 (D^1/2 (gamma/2) Sigma D^1/2)^1/2 D^-1/2  (Bergault et al. 2021)
        std::vector<std::vector<double>> M(n, std::vector<double>(n)); for (size_t a = 0; a < n; ++a) for (size_t b = 0; b < n; ++b) M[a][b] = std::sqrt(D[a]) * 0.5 * p.gamma * risk.Sigma[a][b] * std::sqrt(D[b]);
        std::vector<double> ev; std::vector<std::vector<double>> V; eig_sym(M, ev, V);
        Gamma_.assign(n, std::vector<double>(n, 0.0));
        for (size_t a = 0; a < n; ++a) for (size_t b = 0; b < n; ++b) { double s = 0; for (size_t k = 0; k < n; ++k) s += V[a][k] * std::sqrt(std::max(ev[k], 0.0)) * V[b][k]; Gamma_[a][b] = s / std::sqrt(D[a] * D[b]); }
        kbar_ = kbar; zbar_ = zbar;
    }
    std::string name() const override { return oracle_ ? "Bergault_oracle" : (signal_ ? "CarteaWang" : "Bergault"); }
    void on_context(bool pre) override { pre_event_ = pre; }
    double skew(int tenor, const std::vector<double>& q) const { double s = 0; for (size_t b = 0; b < q.size(); ++b) s += Gamma_[tenor][b] * q[b]; return s; }   // bp; > 0 when long DV01 at that tenor => wants to sell
    double quote(int tier, int tenor, int side, double z, const std::vector<double>& q, double signal_bp, double) override {
        double d0 = 1.0 / k_[tier] + Gamma_[tenor][tenor] * z * ls_.dv01[tenor] / 2;
        double sk = skew(tenor, q);
        // client buys (side +1): dealer sells => wants to sell when long (sk > 0) => tighter; client sells: wider when long
        double d = d0 - side * sk;
        if (oracle_) d += tiers_[tier].mu * z / 10.0 + (pre_event_ ? tiers_[tier].p_info * jump_abs_[tenor] : 0.0);
        if (signal_) { d -= side * std::max(-p_.signal_clip_bp, std::min(p_.signal_clip_bp, p_.signal_horizon * signal_bp));   // yields expected to rise (prices fall): sell cheaper (client buys, side +1), buy dearer
            if (pre_event_) d += p_.event_widen_bp * (tier == 0 ? 1.0 : 0.25); }
        return std::max(0.0, d);
    }
    const std::vector<std::vector<double>>& Gamma() const { return Gamma_; }
protected:
    LadderSpec ls_; MMParams p_; bool signal_; bool pre_event_ = false; bool oracle_ = false; std::vector<ClientTier> tiers_; std::vector<double> jump_abs_; std::vector<double> A_, k_; std::vector<std::vector<double>> Gamma_; double kbar_ = 1, zbar_ = 1;
};

class BarzykinCiceriQuoter : public BergaultQuoter {
public:
    BarzykinCiceriQuoter(const std::vector<ClientTier>& tiers, const LadderSpec& ls, const LadderRisk& risk, MMParams p, bool quality_adjusted = true)
        : BergaultQuoter(tiers, ls, risk, p, false), qa_(quality_adjusted), hits_(tiers.size(), p.hit_target), reqs_(tiers.size(), 1.0), whits_(tiers.size(), p.hit_target) {
        for (auto& t : tiers) { double lo = 0, hi = 50; for (int i = 0; i < 100; ++i) { double m = 0.5 * (lo + hi); (t.f(m) > p.hit_target ? lo : hi) = m; } d0_.push_back(0.5 * (lo + hi)); }
    }
    std::string name() const override { return qa_ ? "BarzykinCiceri_QAHR" : "BarzykinCiceri"; }
    double quote(int tier, int tenor, int side, double, const std::vector<double>& q, double, double) override {
        double h_hat = whits_[tier] / reqs_[tier];
        double d = d0_[tier] - side * skew(tenor, q) + p_.kappa_hit * d0_[tier] * (h_hat - p_.hit_target) / p_.hit_target;
        last_tier_ = tier; reqs_[tier] = 0.98 * reqs_[tier] + 1;   // EWMA of requests
        return std::max(0.0, d);
    }
    void on_request_outcome(int tier, bool hit) { hits_[tier] = 0.98 * hits_[tier] + (hit ? 1 : 0); if (!hit) whits_[tier] = 0.98 * whits_[tier]; }
    void on_fill(const Fill&) override { /* quality-weighted count added when the mark-out is known */ }
    // quality-adjusted count: a hit whose 30-minute mark-out exceeds the spread earned counts 1 + (loss beyond the spread) / spread, capped at 3
    void on_markout(const Fill& f) override { double w = 1.0; if (qa_) { double adverse = f.markouts.size() > 1 ? f.markouts[1] : 0; w = 1.0 + std::min(2.0, std::max(0.0, adverse - f.delta) / std::max(f.delta, 0.1)); } whits_[f.tier] = 0.98 * whits_[f.tier] + w; }
    double hit_ratio(int tier) const { return whits_[tier] / reqs_[tier]; }
private:
    bool qa_; std::vector<double> d0_, hits_, reqs_, whits_; int last_tier_ = 0;
};

} // namespace rcmm
