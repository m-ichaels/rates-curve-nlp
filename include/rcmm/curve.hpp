#pragma once
// Yield curves: Hagan-West monotone-convex forward interpolation, iterative bootstrap of par
// instruments (Treasury par bonds / bills, OIS swaps), a meeting-date step-function front end,
// key-rate DV01s, and a kernel-ridge discount-curve fit (Filipovic, Pelger & Ye) as the alternative.
// Times are years (ACT/365.25); rates are decimals; DV01s are per 1bp on unit notional, in price units.
#include "types.hpp"
#include "stats.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

namespace rcmm {

// ---- Hagan & West (2006) monotone convex interpolation of discrete forwards --------------------------
// Input: knot times t_0 = 0 < t_1 < ... < t_n and discrete forwards fd_i on (t_{i-1}, t_i].
struct MonotoneConvex {
    std::vector<double> t, fd, f;   // knots, discrete forwards (size n), knot forwards (size n+1)
    void build(const std::vector<double>& knots, const std::vector<double>& disc_fwd) {
        t = knots; fd = disc_fwd; size_t n = fd.size(); f.assign(n + 1, 0.0);
        if (n == 1) { f[0] = f[1] = fd[0]; return; }
        for (size_t i = 1; i < n; ++i) f[i] = ((t[i] - t[i - 1]) * fd[i] + (t[i + 1] - t[i]) * fd[i - 1]) / (t[i + 1] - t[i - 1]);
        f[0] = fd[0] - 0.5 * (f[1] - fd[0]); f[n] = fd[n - 1] - 0.5 * (f[n - 1] - fd[n - 1]);
    }
    // g(x) on interval i for x in [0,1]: the four sectors of Hagan & West (2006), with the degenerate cases guarded
    static double g(double g0, double g1, double x) {
        const double eps = 1e-14;
        if (std::fabs(g0) < eps && std::fabs(g1) < eps) return 0;
        if ((g0 > 0 && -0.5 * g0 >= g1 && g1 >= -2 * g0) || (g0 < 0 && -0.5 * g0 <= g1 && g1 <= -2 * g0)) return g0 * (1 - 4 * x + 3 * x * x) + g1 * (-2 * x + 3 * x * x);
        if ((g0 < 0 && g1 > -2 * g0) || (g0 > 0 && g1 < -2 * g0)) { double eta = (g1 + 2 * g0) / (g1 - g0); if (x <= eta || 1 - eta < eps) return g0; return g0 + (g1 - g0) * std::pow((x - eta) / (1 - eta), 2); }
        if ((g0 > 0 && 0 > g1 && g1 > -0.5 * g0) || (g0 < 0 && 0 < g1 && g1 < -0.5 * g0)) { double eta = 3 * g1 / (g1 - g0); if (x >= eta || eta < eps) return g1; return g1 + (g0 - g1) * std::pow((eta - x) / eta, 2); }
        // same sign (or one of them zero)
        if (std::fabs(g0 + g1) < eps) return 0;
        double eta = g1 / (g0 + g1), A = -g0 * g1 / (g0 + g1);
        if (x < eta) return eta < eps ? A : A + (g0 - A) * std::pow((eta - x) / eta, 2);
        return 1 - eta < eps ? A : A + (g1 - A) * std::pow((x - eta) / (1 - eta), 2);
    }
    double forward(double tt) const {
        size_t n = fd.size(); if (n == 0) return 0;
        if (tt <= 0) return f[0]; if (tt >= t[n]) return f[n];
        size_t i = 1; while (i < n && t[i] < tt) ++i;
        double x = (tt - t[i - 1]) / (t[i] - t[i - 1]);
        return fd[i - 1] + g(f[i - 1] - fd[i - 1], f[i] - fd[i - 1], x);
    }
    // integral of the forward from 0 to tt (Simpson on each knot interval)
    // int_0^x g, sector by sector (closed form; int_0^1 g = 0 so each interval reproduces its discrete forward)
    static double G(double g0, double g1, double x) {
        const double eps = 1e-14; auto cube = [](double v) { return v * v * v; };
        if (std::fabs(g0) < eps && std::fabs(g1) < eps) return 0;
        if ((g0 > 0 && -0.5 * g0 >= g1 && g1 >= -2 * g0) || (g0 < 0 && -0.5 * g0 <= g1 && g1 <= -2 * g0)) return g0 * (x - 2 * x * x + x * x * x) + g1 * (-x * x + x * x * x);
        if ((g0 < 0 && g1 > -2 * g0) || (g0 > 0 && g1 < -2 * g0)) { double eta = (g1 + 2 * g0) / (g1 - g0); if (1 - eta < eps) return g0 * x; return g0 * x + (x > eta ? (g1 - g0) * (1 - eta) / 3 * cube((x - eta) / (1 - eta)) : 0); }
        if ((g0 > 0 && 0 > g1 && g1 > -0.5 * g0) || (g0 < 0 && 0 < g1 && g1 < -0.5 * g0)) { double eta = 3 * g1 / (g1 - g0); if (eta < eps) return g1 * x; return g1 * x + (g0 - g1) * eta / 3 * (1 - cube((eta - std::min(x, eta)) / eta)); }
        if (std::fabs(g0 + g1) < eps) return 0;
        double eta = g1 / (g0 + g1), A = -g0 * g1 / (g0 + g1); double acc = A * x;
        if (eta >= eps) acc += (g0 - A) * eta / 3 * (1 - cube((eta - std::min(x, eta)) / eta));
        if (1 - eta >= eps && x > eta) acc += (g1 - A) * (1 - eta) / 3 * cube((x - eta) / (1 - eta));
        return acc;
    }
    double integral(double tt) const {
        double acc = 0; size_t n = fd.size(); if (n == 0) return 0;
        for (size_t i = 1; i <= n && t[i - 1] < tt; ++i) {
            double w = t[i] - t[i - 1], x = std::min(1.0, (tt - t[i - 1]) / w);
            acc += w * (fd[i - 1] * x + G(f[i - 1] - fd[i - 1], f[i] - fd[i - 1], x));
        }
        if (tt > t[n]) acc += f[n] * (tt - t[n]);
        return acc;
    }
    double discount(double tt) const { return std::exp(-integral(tt)); }
};

// ---- instruments ----------------------------------------------------------------------------------------
enum class InstKind { Bill, ParBond, OIS };
struct Instrument { InstKind kind; double T; double rate; int freq = 2; };   // rate as decimal; ParBond: semi-annual coupon; OIS: annual fixed
struct DiscountFn { std::function<double(double)> P; };

inline double price_instrument(const Instrument& ins, const std::function<double(double)>& P) {
    switch (ins.kind) {
        case InstKind::Bill: return P(ins.T) * (1 + ins.rate * ins.T);                 // bond-equivalent simple yield, should equal 1
        case InstKind::ParBond: { double v = 0; int n = (int)std::lround(ins.T * ins.freq); for (int k = 1; k <= n; ++k) v += ins.rate / ins.freq * P((double)k / ins.freq); return v + P(ins.T); }
        case InstKind::OIS: { double v = 0; int n = (int)std::lround(ins.T); if (n < 1) return P(ins.T) * (1 + ins.rate * ins.T); for (int k = 1; k <= n; ++k) v += ins.rate * P((double)k); return v + P(ins.T); }
    }
    return 1;
}

// Bootstrap: discrete forwards on the instrument maturities so that every instrument prices at par
// under the monotone-convex interpolant.  Because the interpolant couples neighbouring intervals,
// the sequential solves are iterated to a global fixed point.
struct Bootstrap {
    MonotoneConvex mc; std::vector<Instrument> inst; int iterations = 0; double max_error_bp = 0;
    void fit(std::vector<Instrument> instruments, double tol = 1e-10, int max_iter = 60) {   // tol on the max forward change per sweep (1e-10 = 1e-6 bp)
        inst = instruments; std::sort(inst.begin(), inst.end(), [](const Instrument& a, const Instrument& b) { return a.T < b.T; });
        size_t n = inst.size(); std::vector<double> knots(n + 1, 0.0); for (size_t i = 0; i < n; ++i) knots[i + 1] = inst[i].T;
        std::vector<double> fd(n, 0.0);
        for (size_t i = 0; i < n; ++i) fd[i] = inst[i].rate;   // start from the par rates
        mc.build(knots, fd);
        for (iterations = 0; iterations < max_iter; ++iterations) {
            double maxchg = 0;
            for (size_t i = 0; i < n; ++i) {
                // solve fd[i] by secant on the pricing error of instrument i
                auto err = [&](double x) { fd[i] = x; mc.build(knots, fd); return price_instrument(inst[i], [&](double tt) { return mc.discount(tt); }) - 1.0; };
                double start = fd[i], x0 = start, e0 = err(x0); if (std::fabs(e0) <= 1e-14) continue; double x1 = x0 + 1e-4, e1 = err(x1);
                for (int k = 0; k < 30 && std::fabs(e1) > 1e-14; ++k) { if (std::fabs(e1 - e0) < 1e-18) break; double x2 = x1 - e1 * (x1 - x0) / (e1 - e0); x0 = x1; e0 = e1; x1 = x2; e1 = err(x1); }
                maxchg = std::max(maxchg, std::fabs(x1 - start));
            }
            if (maxchg < tol) { ++iterations; break; }
        }
        max_error_bp = 0; for (auto& ins : inst) max_error_bp = std::max(max_error_bp, std::fabs(price_instrument(ins, [&](double tt) { return mc.discount(tt); }) - 1.0) * 1e4);
    }
    double P(double tt) const { return mc.discount(tt); }
    double zero(double tt) const { return tt > 0 ? -std::log(P(tt)) / tt : mc.forward(0); }
    double par_yield(double T, int freq = 2) const { double a = 0; int n = (int)std::lround(T * freq); for (int k = 1; k <= n; ++k) a += P((double)k / freq) / freq; return (1 - P(T)) / a; }
};

// Key-rate DV01 of a par bond of maturity T (unit notional, price per 1) w.r.t. each input par rate (1bp bumps).
inline std::vector<double> key_rate_dv01(const std::vector<Instrument>& inst, double T, double coupon, int freq = 2) {
    Bootstrap base; base.fit(inst);
    Instrument bond{InstKind::ParBond, T, coupon, freq};
    double p0 = price_instrument(bond, [&](double tt) { return base.P(tt); });
    std::vector<double> kr;
    for (size_t k = 0; k < inst.size(); ++k) {
        std::vector<Instrument> up = inst; up[k].rate += 1e-4; Bootstrap b; b.fit(up);
        kr.push_back(p0 - price_instrument(bond, [&](double tt) { return b.P(tt); }));   // price fall per 1 bp rise of that input, per 1 notional
    }
    return kr;   // price change per 1bp bump of input k (negative for rising rates)
}

// ---- meeting-date step-function front end -------------------------------------------------------------
// Instantaneous forward piecewise constant between policy meeting dates (times in years), calibrated by
// least squares to a set of observed short rates.  Native input is SOFR futures settlements (SR1: monthly
// arithmetic average of the daily rate, SR3: quarterly compounded); with `--futures` absent the same
// solver is fed the bill yields (a proxy, stated in the README).
struct StepFrontEnd {
    std::vector<double> meeting_t;   // sorted, > 0
    std::vector<double> level;       // forward on [meeting_{k-1}, meeting_k), level[0] applies from 0
    double forward(double tt) const { size_t k = 0; while (k < meeting_t.size() && tt >= meeting_t[k]) ++k; return level[std::min(k, level.size() - 1)]; }
    double integral(double tt) const { double acc = 0, a = 0; for (size_t k = 0; k < level.size(); ++k) { double b = k < meeting_t.size() ? meeting_t[k] : 1e9; if (tt <= a) break; acc += level[k] * (std::min(b, tt) - a); a = b; } return acc; }
    double P(double tt) const { return std::exp(-integral(tt)); }
    // observations: (start, end, average-rate) - a bill from 0 to T with simple yield y gives average forward ln(1+yT)/T;
    // a 3m futures contract [a,b] with price p gives average ~ (100-p)/100 (compounding ignored at these rates)
    struct Obs { double a, b, avg; };
    void fit(const std::vector<double>& meetings, const std::vector<Obs>& obs, double sofr_today) {
        meeting_t = meetings; size_t K = meetings.size() + 1; level.assign(K, sofr_today);
        // least squares with a small ridge on step changes (keeps unobserved steps continuous)
        std::vector<std::vector<double>> A; std::vector<double> y;
        for (auto& o : obs) { std::vector<double> row(K, 0.0); double a = 0; for (size_t k = 0; k < K; ++k) { double b = k < meetings.size() ? meetings[k] : 1e9; double lo = std::max(a, o.a), hi = std::min(b, o.b); if (hi > lo) row[k] += (hi - lo) / (o.b - o.a); a = b; } A.push_back(row); y.push_back(o.avg); }
        { std::vector<double> row(K, 0.0); row[0] = 1; A.push_back(row); y.push_back(sofr_today); }   // today's fixing anchors the first step
        double lam = 1e-3; for (size_t k = 1; k < K; ++k) { std::vector<double> row(K, 0.0); row[k] = lam; row[k - 1] = -lam; A.push_back(row); y.push_back(0); }
        // normal equations
        std::vector<std::vector<double>> M(K, std::vector<double>(K, 0.0)); std::vector<double> v(K, 0.0);
        for (size_t r = 0; r < A.size(); ++r) for (size_t i = 0; i < K; ++i) { v[i] += A[r][i] * y[r]; for (size_t j = 0; j < K; ++j) M[i][j] += A[r][i] * A[r][j]; }
        for (size_t c = 0; c < K; ++c) { size_t piv = c; for (size_t r = c + 1; r < K; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r; std::swap(M[c], M[piv]); std::swap(v[c], v[piv]); for (size_t r = c + 1; r < K; ++r) { double f = M[r][c] / M[c][c]; for (size_t j = c; j < K; ++j) M[r][j] -= f * M[c][j]; v[r] -= f * v[c]; } }
        for (size_t c = K; c-- > 0;) { double s = v[c]; for (size_t j = c + 1; j < K; ++j) s -= M[c][j] * level[j]; level[c] = s / M[c][c]; }
    }
};

// ---- kernel ridge discount curve (Filipovic, Pelger & Ye) ------------------------------------------------
// d(t) = 1 - sum_j w_j k(t, t_j) with the cubic-spline (second-derivative Sobolev) kernel
// k(s,t) = min(s,t)^2 (3 max(s,t) - min(s,t)) / 6 on the maturity-weighted cash-flow grid; weights minimise
// sum_i (price_i - 1)^2 / T_i + lambda w'Kw.  Pricing is linear in d, so this is one linear solve.
struct KernelRidge {
    std::vector<double> tj, w; double lambda = 1e-6, beta = 0;   // d(t) = 1 - beta t - sum_j w_j k(t, t_j); the linear term is the unpenalised null space of the curvature penalty
    static double kern(double s, double t) { double a = std::min(s, t), b = std::max(s, t); return a * a * (3 * b - a) / 6.0; }   // reproducing kernel of int d''(t)^2 dt with d(0) = 1
    double P(double tt) const { double d = 1 - beta * tt; for (size_t j = 0; j < tj.size(); ++j) d -= w[j] * kern(tt, tj[j]); return d; }
    double forward(double tt, double h = 1e-4) const { return -(std::log(P(tt + h)) - std::log(P(std::max(tt - h, 0.0)))) / (tt + h - std::max(tt - h, 0.0)); }
    void fit(const std::vector<Instrument>& inst, double lam = 1e-6) {
        lambda = lam;
        // cash-flow times = knot set
        std::vector<double> times; for (auto& ins : inst) { if (ins.kind == InstKind::ParBond) { int n = (int)std::lround(ins.T * ins.freq); for (int k = 1; k <= n; ++k) times.push_back((double)k / ins.freq); } else if (ins.kind == InstKind::OIS) { int n = (int)std::lround(ins.T); for (int k = 1; k <= n; ++k) times.push_back((double)k); times.push_back(ins.T); } else times.push_back(ins.T); }
        std::sort(times.begin(), times.end()); times.erase(std::unique(times.begin(), times.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), times.end());
        tj = times; size_t m = tj.size(), n = inst.size();
        // price_i = c_i . d(t) where c_i are cash flows on the grid; d = 1 - K w  =>  price_i - 1 = sum_l c_il (1 - (K w)_l) - 1
        std::vector<std::vector<double>> C(n, std::vector<double>(m, 0.0)); std::vector<double> b(n, 0.0), wt(n, 1.0);
        auto idx = [&](double tt) { return (size_t)(std::lower_bound(tj.begin(), tj.end(), tt - 1e-9) - tj.begin()); };
        for (size_t i = 0; i < n; ++i) {
            const Instrument& ins = inst[i]; wt[i] = 1.0 / std::max(ins.T, 0.05);
            if (ins.kind == InstKind::ParBond) { int nn = (int)std::lround(ins.T * ins.freq); for (int k = 1; k <= nn; ++k) C[i][idx((double)k / ins.freq)] += ins.rate / ins.freq; C[i][idx(ins.T)] += 1; }
            else if (ins.kind == InstKind::OIS) { int nn = (int)std::lround(ins.T); for (int k = 1; k <= nn; ++k) C[i][idx((double)k)] += ins.rate; C[i][idx(ins.T)] += 1; }
            else C[i][idx(ins.T)] += 1 + ins.rate * ins.T;
            double s = 0; for (double x : C[i]) s += x; b[i] = s - 1.0;   // residual when d == 1
        }
        std::vector<std::vector<double>> K(m, std::vector<double>(m)); for (size_t a = 0; a < m; ++a) for (size_t c = 0; c < m; ++c) K[a][c] = kern(tj[a], tj[c]);
        // G = [C K | C t] (n x (m+1)), unknowns u = (w, beta); minimise sum_i wt_i (b_i - (G u)_i)^2 + lambda w'Kw  =>  (G' W G + lambda K_aug) u = G' W b
        size_t m1 = m + 1; std::vector<std::vector<double>> G(n, std::vector<double>(m1, 0.0)); for (size_t i = 0; i < n; ++i) { for (size_t c = 0; c < m; ++c) for (size_t l = 0; l < m; ++l) G[i][c] += C[i][l] * K[l][c]; for (size_t l = 0; l < m; ++l) G[i][m] += C[i][l] * tj[l]; }
        std::vector<std::vector<double>> M(m1, std::vector<double>(m1, 0.0)); std::vector<double> v(m1, 0.0);
        for (size_t a = 0; a < m1; ++a) for (size_t c = 0; c < m1; ++c) { double s = (a < m && c < m) ? lambda * K[a][c] : 0.0; for (size_t i = 0; i < n; ++i) s += wt[i] * G[i][a] * G[i][c]; M[a][c] = s; }
        for (size_t a = 0; a < m1; ++a) { double s = 0; for (size_t i = 0; i < n; ++i) s += wt[i] * G[i][a] * b[i]; v[a] = s; }
        w.assign(m1, 0.0); m = m1;
        for (size_t c = 0; c < m; ++c) { size_t piv = c; for (size_t r = c + 1; r < m; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r; std::swap(M[c], M[piv]); std::swap(v[c], v[piv]); if (std::fabs(M[c][c]) < 1e-300) continue; for (size_t r = c + 1; r < m; ++r) { double f = M[r][c] / M[c][c]; for (size_t j = c; j < m; ++j) M[r][j] -= f * M[c][j]; v[r] -= f * v[c]; } }
        for (size_t c = m; c-- > 0;) { double s = v[c]; for (size_t j = c + 1; j < m; ++j) s -= M[c][j] * w[j]; w[c] = std::fabs(M[c][c]) < 1e-300 ? 0 : s / M[c][c]; }
        beta = w.back(); w.pop_back();
    }
    double max_error_bp(const std::vector<Instrument>& inst) const { double e = 0; for (auto& ins : inst) e = std::max(e, std::fabs(price_instrument(ins, [&](double tt) { return P(tt); }) - 1.0) * 1e4); return e; }
};

// forward smoothness: integral of (f'(t))^2 over [0, T]
template <class F> double forward_roughness(F&& fwd, double T, int steps = 600) { double h = T / steps, acc = 0; for (int k = 1; k < steps; ++k) { double d = (fwd((k + 1) * h) - fwd((k - 1) * h)) / (2 * h); acc += d * d * h; } return acc; }

} // namespace rcmm
