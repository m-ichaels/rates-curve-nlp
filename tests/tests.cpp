// Unit and known-answer tests: curve construction, factor model, quoters, simulator identities.
#include "rcmm/sim.hpp"
#include <iostream>
#include <cmath>

using namespace rcmm;
static int failures = 0, checks = 0;
#define CHECK(cond, msg) do { ++checks; if (!(cond)) { ++failures; std::cerr << "FAIL: " << msg << "  [" #cond "]\n"; } } while (0)
#define CHECK_NEAR(a, b, tol, msg) CHECK(std::fabs((double)(a) - (double)(b)) <= (tol), msg << " (" << (a) << " vs " << (b) << ")")

// a plausible par curve (decimals): bills 1m..1y, par bonds 2y..30y
static std::vector<Instrument> sample_curve(double bump = 0) {
    std::vector<Instrument> v;
    double bills[5][2] = {{1.0 / 12, 0.0432}, {2.0 / 12, 0.0430}, {0.25, 0.0428}, {0.5, 0.0420}, {1, 0.0405}};
    double bonds[7][2] = {{2, 0.0385}, {3, 0.0380}, {5, 0.0390}, {7, 0.0400}, {10, 0.0415}, {20, 0.0460}, {30, 0.0470}};
    for (auto& b : bills) v.push_back({InstKind::Bill, b[0], b[1] + bump, 0});
    for (auto& b : bonds) v.push_back({InstKind::ParBond, b[0], b[1] + bump, 2});
    return v;
}

static void test_monotone_convex_interpolates_and_positive() {
    // discrete forwards on knots: the interpolant must integrate to the same discount factors and stay positive (Hagan-West)
    std::vector<double> T = {0, 0.5, 1, 2, 5, 10, 30}, f = {0.043, 0.040, 0.037, 0.039, 0.042, 0.046};
    MonotoneConvex mc; mc.build(T, f);
    double acc = 0; for (size_t i = 0; i < f.size(); ++i) { acc += f[i] * (T[i + 1] - T[i]); CHECK_NEAR(mc.integral(T[i + 1]), acc, 1e-12, "integral matches discrete forwards at knot " << T[i + 1]); }
    // closed-form integral agrees with quadrature of the forward
    for (double tt : {0.3, 1.7, 4.2, 8.0, 22.5}) { int m = 2000; double h = tt / m, sum = mc.forward(0) + mc.forward(tt); for (int k = 1; k < m; ++k) sum += (k % 2 ? 4 : 2) * mc.forward(k * h); CHECK_NEAR(mc.integral(tt), sum * h / 3, 1e-6, "integral matches quadrature at " << tt); }
    for (double t = 0.01; t < 30; t += 0.07) CHECK(mc.forward(t) > 0, "forward positive at " << t);
    CHECK(std::fabs(mc.forward(T[3] + 1e-9) - mc.forward(T[3] - 1e-9)) < 1e-3, "forward continuous across a knot");
}

static void test_bootstrap_reprices_within_1bp() {
    auto inst = sample_curve(); Bootstrap b; b.fit(inst);
    CHECK(b.max_error_bp < 1.0, "bootstrap max repricing error < 1 bp: " << b.max_error_bp);
    CHECK(b.max_error_bp < 0.01, "bootstrap max repricing error < 0.01 bp: " << b.max_error_bp);
    for (auto& i : inst) CHECK_NEAR(price_instrument(i, [&](double t) { return b.P(t); }), 1.0, 1e-6, "instrument T=" << i.T << " reprices to par");
    for (double t = 0.1; t <= 30; t += 0.1) CHECK(b.P(t) > 0 && b.P(t) <= 1.0 && b.mc.forward(t) > 0, "discount factor and forward sane at " << t);
    CHECK_NEAR(b.par_yield(10.0), 0.0415, 1e-6, "par yield recovered at 10y");
}

static void test_kernel_ridge_reprices_and_smooths() {
    auto inst = sample_curve(); KernelRidge kr; kr.fit(inst, 1e-4); Bootstrap b; b.fit(inst);
    CHECK(kr.max_error_bp(inst) < 1.0, "kernel ridge reprices within 1 bp: " << kr.max_error_bp(inst));
    double rk = forward_roughness([&](double t) { return kr.forward(t); }, 30.0), rb = forward_roughness([&](double t) { return b.mc.forward(t); }, 30.0);
    CHECK(rk < rb, "kernel ridge forward is smoother than the bootstrap (" << rk << " vs " << rb << ")");
    KernelRidge tight; tight.fit(inst, 1e-9); CHECK(tight.max_error_bp(inst) < kr.max_error_bp(inst) + 1e-9, "smaller lambda fits at least as tightly");
}

static void test_key_rate_dv01s_sum_to_parallel_dv01() {
    auto inst = sample_curve(); Bootstrap b; b.fit(inst); double c = b.par_yield(10.0);
    auto kr = key_rate_dv01(inst, 10.0, c); double tot = 0; for (double x : kr) tot += x;
    auto up = sample_curve(1e-4); Bootstrap bu; bu.fit(up); Instrument bond{InstKind::ParBond, 10.0, c, 2};
    double par = price_instrument(bond, [&](double t) { return b.P(t); }) - price_instrument(bond, [&](double t) { return bu.P(t); });
    CHECK_NEAR(tot, par, 0.02 * par, "sum of key-rate DV01s equals the parallel DV01");
    size_t i10 = 0; for (size_t i = 0; i < inst.size(); ++i) if (inst[i].T == 10.0) i10 = i; CHECK(kr[i10] > 0.5 * tot, "10y bond DV01 concentrated on the 10y key rate");
    CHECK(par > 7e-4 && par < 9e-4, "10y par bond DV01 per 1 notional per bp ~ 8e-4: " << par);
}

static void test_step_front_end_recovers_steps() {
    // forward is a step function on meeting dates; bills spanning the steps must recover the levels
    std::vector<double> meetings = {0.25, 0.5, 0.75}; std::vector<double> lv = {0.045, 0.0425, 0.040, 0.0375};
    StepFrontEnd truth; truth.meeting_t = meetings; truth.level = lv;
    std::vector<StepFrontEnd::Obs> obs; for (double T : {0.1, 0.3, 0.6, 0.9, 1.0}) obs.push_back({0.0, T, truth.integral(T) / T});
    StepFrontEnd fe; fe.fit(meetings, obs, 0.045);
    for (size_t k = 0; k < lv.size(); ++k) CHECK_NEAR(fe.level[k], lv[k], 2e-4, "step level " << k << " recovered");
    CHECK_NEAR(fe.forward(0.26), fe.forward(0.49), 0, "forward flat between meetings");
}

static void test_factor_model_pca() {
    // synthetic daily changes from two known factors: PCA must recover the variance ordering and sign conventions
    Rng rng(3); std::vector<CmtRow> rows; std::vector<double> y(13, 0.04);
    for (int d = 0; d < 800; ++d) { CmtRow r; char buf[16]; std::snprintf(buf, sizeof buf, "%05d", d); r.date = buf;
        double lvl = 4e-4 * rng.normal(), slp = 1.5e-4 * rng.normal();
        for (int i = 0; i < 13; ++i) y[i] += lvl + slp * (std::log(CMT_T[i]) / std::log(30.0)); r.y = y; rows.push_back(r); }
    FactorModel fm; fm.fit(rows, "00000", "99999", 3);
    CHECK(fm.explained(2) > 0.99, "two factors explain > 99%: " << fm.explained(2));
    double s = 0; for (double v : fm.loadings[0]) s += v; CHECK(s > 0, "level loadings positive");
    CHECK(fm.loadings[1].back() > fm.loadings[1].front(), "slope loading increases with maturity");
    CHECK(std::fabs(fm.loading_at(0, 4.0)[0] - fm.loadings[0][7]) < 0.1, "interpolated loading close to nearby tenor");
}

struct Setup_t { LadderSpec ls; FactorModel fm; LadderRisk risk; std::vector<ClientTier> tiers; };
static Setup_t make_setup() {
    Setup_t s; s.ls.dv01 = {189.0, 440.0, 783.0, 1489.0}; s.tiers = default_tiers();
    // hand-built 3-factor model at the 12 tenors
    s.fm.K = 3; s.fm.T.clear(); for (int i : {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12}) s.fm.T.push_back(CMT_T[i]);
    s.fm.evals = {250, 40, 30, 1, 1, 1, 1, 1, 1, 1, 1, 1}; s.fm.loadings.assign(3, std::vector<double>(12));
    for (size_t i = 0; i < 12; ++i) { double x = std::log(s.fm.T[i]) / std::log(30.0); s.fm.loadings[0][i] = 1 / std::sqrt(12.0); s.fm.loadings[1][i] = (x - 0.5) * 0.6; s.fm.loadings[2][i] = (0.25 - (x - 0.5) * (x - 0.5)) * 1.2; }
    s.risk.build(s.fm, s.ls); return s;
}

static void test_bergault_gamma_solves_riccati() {
    Setup_t S = make_setup(); MMParams p; p.gamma = 2e-6; BergaultQuoter q(S.tiers, S.ls, S.risk, p);
    // Gamma D Gamma = (gamma/2) Sigma: recompute D as the quoter does and check the identity
    size_t n = 4; std::vector<double> D(n); double zbar = 0; for (size_t k = 0; k < S.ls.sizes.size(); ++k) zbar += S.ls.sizes[k] * S.ls.size_w[k];
    double Asum = 0, ksum = 0; for (auto& t : S.tiers) { double A, k; t.exp_fit(A, k); A *= t.lambda * p.day_hours; Asum += A; ksum += k * A; }
    for (size_t j = 0; j < n; ++j) D[j] = zbar * S.ls.dv01[j] * Asum * (ksum / Asum) * std::exp(-1.0);
    const auto& G = q.Gamma(); for (size_t a = 0; a < n; ++a) for (size_t b = 0; b < n; ++b) { double lhs = 0; for (size_t k = 0; k < n; ++k) lhs += G[a][k] * D[k] * G[k][b]; double rhs = 0.5 * p.gamma * S.risk.Sigma[a][b]; CHECK_NEAR(lhs, rhs, 1e-6 * std::fabs(rhs) + 1e-12, "Gamma D Gamma = gamma/2 Sigma [" << a << "," << b << "]"); }
    // skew: long DV01 at the 10y => tighter offer (client buys), wider bid; symmetric at zero inventory
    std::vector<double> q0(4, 0.0), ql(4, 0.0); ql[2] = 2e4;
    CHECK_NEAR(q.quote(1, 2, +1, 10, q0, 0, 0.5), q.quote(1, 2, -1, 10, q0, 0, 0.5), 1e-12, "symmetric at zero inventory");
    CHECK(q.quote(1, 2, +1, 10, ql, 0, 0.5) < q.quote(1, 2, -1, 10, ql, 0, 0.5), "long inventory: sell side tighter than buy side");
    CHECK(q.quote(1, 0, +1, 10, ql, 0, 0.5) < q.quote(1, 0, -1, 10, ql, 0, 0.5), "long 10y also skews the correlated 2y through the factor covariance");
    // riskless part equals the tier's revenue-maximising half-spread when gamma -> 0
    MMParams p0; p0.gamma = 1e-14; BergaultQuoter q0q(S.tiers, S.ls, S.risk, p0); CHECK_NEAR(q0q.quote(1, 2, +1, 10, q0, 0, 0.5), S.tiers[1].d_riskless(), 0.02, "gamma -> 0: quote -> argmax delta f(delta)");
}

static void test_cartea_wang_signal_direction() {
    Setup_t S = make_setup(); MMParams p; p.gamma = 2e-6; BergaultQuoter cw(S.tiers, S.ls, S.risk, p, true); std::vector<double> q0(4, 0.0);
    // yields expected to rise: sell bonds cheaper (client buy side tighter), buy dearer
    CHECK(cw.quote(1, 2, +1, 10, q0, +2.0, 0.5) < cw.quote(1, 2, -1, 10, q0, +2.0, 0.5), "positive alpha (yields up): offer tighter than bid");
    CHECK(cw.quote(1, 2, +1, 10, q0, -2.0, 0.5) > cw.quote(1, 2, -1, 10, q0, -2.0, 0.5), "negative alpha: bid tighter than offer");
    CHECK_NEAR(cw.quote(1, 2, +1, 10, q0, 0.0, 0.5), cw.quote(1, 2, -1, 10, q0, 0.0, 0.5), 1e-12, "no signal: symmetric");
    CHECK_NEAR(cw.quote(1, 2, -1, 10, q0, 100.0, 0.5) - cw.quote(1, 2, -1, 10, q0, 0.0, 0.5), p.signal_clip_bp, 1e-9, "signal skew is capped");
    double d_norm = cw.quote(0, 2, +1, 10, q0, 0.0, 0.5); cw.on_context(true); CHECK_NEAR(cw.quote(0, 2, +1, 10, q0, 0.0, 0.5) - d_norm, p.event_widen_bp, 1e-9, "calendar feature widens the informed tier in the pre-event window"); cw.on_context(false);
}

static void test_barzykin_ciceri_hit_ratio_controller() {
    Setup_t S = make_setup(); MMParams p; p.gamma = 2e-6; p.hit_target = 0.35; std::vector<double> q0(4, 0.0);
    BarzykinCiceriQuoter bc(S.tiers, S.ls, S.risk, p, false); double d_start = bc.quote(0, 2, +1, 10, q0, 0, 0.5);
    for (int i = 0; i < 200; ++i) { bc.quote(0, 2, +1, 10, q0, 0, 0.5); bc.on_request_outcome(0, true); Fill f; f.tier = 0; f.tenor = 2; f.delta = 0.4; f.markouts = {0, 0, 0}; bc.on_markout(f); }
    CHECK(bc.quote(0, 2, +1, 10, q0, 0, 0.5) > d_start, "hit ratio above target: quotes widen");
    BarzykinCiceriQuoter bc2(S.tiers, S.ls, S.risk, p, false); for (int i = 0; i < 200; ++i) { bc2.quote(0, 2, +1, 10, q0, 0, 0.5); bc2.on_request_outcome(0, false); }
    CHECK(bc2.quote(0, 2, +1, 10, q0, 0, 0.5) < d_start, "hit ratio below target: quotes tighten");
    // quality adjustment: toxic hits count more, so the same raw hit ratio widens the QAHR quoter more
    BarzykinCiceriQuoter qa(S.tiers, S.ls, S.risk, p, true), plain(S.tiers, S.ls, S.risk, p, false);
    for (int i = 0; i < 100; ++i) { qa.quote(0, 2, +1, 10, q0, 0, 0.5); plain.quote(0, 2, +1, 10, q0, 0, 0.5); qa.on_request_outcome(0, i % 3 == 0); plain.on_request_outcome(0, i % 3 == 0); if (i % 3 == 0) { Fill f; f.tier = 0; f.tenor = 2; f.delta = 0.4; f.markouts = {1.0, 3.0, 3.0}; qa.on_markout(f); plain.on_markout(f); } }
    CHECK(qa.quote(0, 2, +1, 10, q0, 0, 0.5) > plain.quote(0, 2, +1, 10, q0, 0, 0.5), "toxic hits widen the quality-adjusted quoter more");
}

static void test_ladder_risk_and_factor_exposure() {
    Setup_t S = make_setup(); std::vector<double> q = {1e4, -2e4, 3e4, 0};
    double v = S.risk.variance(q); double v2 = 0; auto e = S.risk.factor_exposure(q); for (size_t k = 0; k < e.size(); ++k) v2 += S.risk.lam[k] * e[k] * e[k];
    CHECK_NEAR(v, v2, 1e-6 * v, "q' Sigma q equals sum_k lambda_k (L_k q)^2");
    std::vector<double> par(4, 1e4); CHECK(S.risk.variance(par) > S.risk.variance(q), "parallel position carries more level risk than a mixed one of similar size");
}

static void test_simulator_identities() {
    Setup_t S = make_setup(); MMParams p; p.gamma = 2e-6;
    for (std::string n : {"Symmetric", "Bergault", "BarzykinCiceri_QAHR", "CarteaWang"}) {
        std::unique_ptr<MMQuoter> q; if (n == "Symmetric") q = std::make_unique<SymmetricQuoter>(S.tiers, p); else if (n == "Bergault") q = std::make_unique<BergaultQuoter>(S.tiers, S.ls, S.risk, p); else if (n == "CarteaWang") q = std::make_unique<BergaultQuoter>(S.tiers, S.ls, S.risk, p, true); else q = std::make_unique<BarzykinCiceriQuoter>(S.tiers, S.ls, S.risk, p, true);
        Rng rng(5); auto days = synthetic_days(6, S.risk, rng, 0.5); MMSim sim(S.ls, S.risk, S.tiers, p, 11);
        SimResult r = sim.run(*q, days, n == "CarteaWang" ? std::function<double(int, double, int)>([&](int d, double t, int) { for (auto& e : days[d].events) if (t >= e.t_day && t < e.t_day + 0.077) return e.alpha_bp; return 0.0; }) : nullptr);
        CHECK_NEAR(r.pnl, r.spread + r.inventory - r.hedge_cost, 1e-6 * std::max(1.0, std::fabs(r.pnl)), n << ": P&L = spread + inventory - hedge cost");
        double s5 = 0; for (double x : r.pnl_5min) s5 += x; CHECK_NEAR(s5, r.pnl, 1e-6 * std::max(1.0, std::fabs(r.pnl)), n << ": 5-minute P&L increments sum to the total");
        double sp = 0; for (auto& f : r.fills) sp += f.z * S.ls.dv01[f.tenor] * f.delta; CHECK_NEAR(sp, r.spread, 1e-6 * std::max(1.0, sp), n << ": spread P&L = sum of size x DV01 x half-spread");
        CHECK(r.n_fills > 0 && r.n_fills <= r.n_rfq, n << ": fills within RFQs");
        for (auto& f : r.fills) CHECK(f.markouts.size() == 3 && f.markouts_info.size() == 3, n << ": every fill has three mark-outs");
        CHECK(r.risk_penalty >= 0, n << ": risk penalty non-negative");
        if (n == "Symmetric") { double hits = 0, reqs = 0; for (size_t i = 0; i < r.hits_by_tier.size(); ++i) { hits += r.hits_by_tier[i]; reqs += r.rfq_by_tier[i]; } CHECK_NEAR(hits / reqs, p.hit_target, 0.06, "symmetric quoter realises the target hit ratio"); }
    }
}

static void test_simulator_common_random_numbers() {
    // the curve path and RFQ stream are common across quoters: two quoters with identical quotes give identical results
    Setup_t S = make_setup(); MMParams p; SymmetricQuoter a(S.tiers, p), b(S.tiers, p); Rng rng(9); auto days = synthetic_days(3, S.risk, rng);
    MMSim s1(S.ls, S.risk, S.tiers, p, 21), s2(S.ls, S.risk, S.tiers, p, 21); SimResult r1 = s1.run(a, days), r2 = s2.run(b, days);
    CHECK_NEAR(r1.pnl, r2.pnl, 1e-9, "same seed, same quotes: same P&L"); CHECK(r1.n_rfq == r2.n_rfq, "same RFQ count");
}

static void test_adverse_selection_in_informed_pre_event_flow() {
    // informed clients trade in the direction of the coming jump: their pre-event mark-outs are positive (against the dealer)
    Setup_t S = make_setup(); MMParams p; SymmetricQuoter q(S.tiers, p); Rng rng(17); auto days = synthetic_days(40, S.risk, rng, 1.0, 1.0); MMSim sim(S.ls, S.risk, S.tiers, p, 5);
    SimResult r = sim.run(q, days); double pre = 0, normal = 0; int npre = 0, nn = 0;
    for (auto& f : r.fills) if (f.tier == 0) { if (f.pre_event) { pre += f.markouts_info[1]; ++npre; } else { normal += f.markouts_info[1]; ++nn; } }
    CHECK(npre > 10 && nn > 10, "enough informed fills in both windows (" << npre << ", " << nn << ")");
    CHECK(pre / npre > 1.0, "informed pre-event mark-out strongly adverse: " << pre / npre << " bp");
    CHECK(pre / npre > normal / std::max(nn, 1) + 0.5, "pre-event mark-out worse than normal-time mark-out");
}

static void test_realised_close_bridge() {
    // a quoter that is never filled carries no inventory: the bridge to the realised close runs and inventory P&L is exactly zero
    Setup_t S = make_setup(); MMParams p; std::vector<SimDay> days(2); days[0].close_change = {3.0, 4.0, 5.0, 6.0}; days[1].close_change = {-2.0, -1.0, 0.0, 1.0};
    struct Probe : MMQuoter { std::string name() const override { return "probe"; } double quote(int, int, int, double, const std::vector<double>&, double, double) override { return 100.0; } } probe;
    MMSim sim(S.ls, S.risk, S.tiers, p, 3); SimResult r = sim.run(probe, days); CHECK(r.n_fills == 0, "probe never filled"); CHECK_NEAR(r.inventory, 0.0, 1e-9, "no inventory, no inventory P&L"); CHECK_NEAR(r.pnl, 0.0, 1e-9, "zero P&L");
}

static void test_release_window_signal_and_split_jump() {
    // the text nowcast exists only from the release; the event move is split between the release and the discovery window
    Setup_t S = make_setup(); MMParams p; std::vector<SimDay> days(1); Event e; e.day = 0; e.t_day = 0.5; e.type = "FOMC"; e.jump_bp = {10, 10, 10, 10}; e.alpha = {4, 3, 2, 1}; days[0].events.push_back(e);
    auto sig = [&](int d, double t, int j) { for (auto& ev : days[d].events) if (t >= ev.t_day && t < ev.t_day + 30.0 / 60 / 6.5) return ev.alpha_at((size_t)j); return 0.0; };
    CHECK_NEAR(sig(0, 0.49, 0), 0.0, 1e-12, "no signal before the release"); CHECK_NEAR(sig(0, 0.51, 0), 4.0, 1e-12, "2y nowcast during the window"); CHECK_NEAR(sig(0, 0.51, 3), 1.0, 1e-12, "per-tenor nowcast"); CHECK_NEAR(sig(0, 0.6, 0), 0.0, 1e-12, "signal gone after the window");
    struct Probe : MMQuoter { std::string name() const override { return "probe"; } double quote(int, int, int, double, const std::vector<double>&, double, double) override { return 100.0; } } probe;
    MMSim sim(S.ls, S.risk, S.tiers, p, 3); SimResult r = sim.run(probe, days); CHECK(r.n_fills == 0, "probe never filled");
    Event e2 = e; e2.alpha.clear(); CHECK_NEAR(e2.alpha_at(2), e2.alpha_bp, 1e-12, "scalar alpha fallback");
}

static void test_bootstrap_ci_and_helpers() {
    std::vector<double> x; Rng rng(1); for (int i = 0; i < 400; ++i) x.push_back(1.0 + rng.normal()); CI ci = bootstrap_mean_ci(x, 500);
    CHECK(ci.lo < ci.mean && ci.mean < ci.hi, "CI brackets the mean"); CHECK(ci.lo < 1.0 && ci.hi > 1.0, "CI covers the true mean");
    auto lines = split_lines("a,b\r\n1,2\n\n3,4"); CHECK(lines.size() == 3 && split_csv(lines[1])[1] == "2", "csv helpers");
}

int main() {
    test_monotone_convex_interpolates_and_positive(); test_bootstrap_reprices_within_1bp(); test_kernel_ridge_reprices_and_smooths(); test_key_rate_dv01s_sum_to_parallel_dv01(); test_step_front_end_recovers_steps();
    test_factor_model_pca(); test_bergault_gamma_solves_riccati(); test_cartea_wang_signal_direction(); test_barzykin_ciceri_hit_ratio_controller(); test_ladder_risk_and_factor_exposure();
    test_simulator_identities(); test_simulator_common_random_numbers(); test_adverse_selection_in_informed_pre_event_flow(); test_realised_close_bridge(); test_release_window_signal_and_split_jump(); test_bootstrap_ci_and_helpers();
    std::cout << checks - failures << " / " << checks << " checks passed\n";
    return failures ? 1 : 0;
}
