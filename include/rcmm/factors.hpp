#pragma once
// Factor risk: PCA of daily par-yield changes (level / slope / curvature), the factor covariance, and
// a curve simulator on the factors.  Inventory risk is measured in factor DV01s.
#include "data.hpp"
#include <array>

namespace rcmm {

// Jacobi eigen-decomposition of a symmetric matrix (small n).
inline void eig_sym(std::vector<std::vector<double>> A, std::vector<double>& evals, std::vector<std::vector<double>>& evecs) {
    size_t n = A.size(); evecs.assign(n, std::vector<double>(n, 0.0)); for (size_t i = 0; i < n; ++i) evecs[i][i] = 1;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0; for (size_t i = 0; i < n; ++i) for (size_t j = i + 1; j < n; ++j) off += A[i][j] * A[i][j];
        if (off < 1e-22) break;
        for (size_t p = 0; p < n; ++p) for (size_t q = p + 1; q < n; ++q) {
            if (std::fabs(A[p][q]) < 1e-300) continue;
            double theta = (A[q][q] - A[p][p]) / (2 * A[p][q]); double t = (theta >= 0 ? 1 : -1) / (std::fabs(theta) + std::sqrt(theta * theta + 1)); double c = 1 / std::sqrt(t * t + 1), s = t * c;
            for (size_t k = 0; k < n; ++k) { double akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq; }
            for (size_t k = 0; k < n; ++k) { double apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk; }
            for (size_t k = 0; k < n; ++k) { double vkp = evecs[k][p], vkq = evecs[k][q]; evecs[k][p] = c * vkp - s * vkq; evecs[k][q] = s * vkp + c * vkq; }
        }
    }
    evals.resize(n); for (size_t i = 0; i < n; ++i) evals[i] = A[i][i];
    // sort descending
    std::vector<size_t> idx(n); for (size_t i = 0; i < n; ++i) idx[i] = i; std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return evals[a] > evals[b]; });
    std::vector<double> ev2(n); std::vector<std::vector<double>> V2(n, std::vector<double>(n));
    for (size_t j = 0; j < n; ++j) { ev2[j] = evals[idx[j]]; for (size_t i = 0; i < n; ++i) V2[i][j] = evecs[i][idx[j]]; }
    evals = ev2; evecs = V2;
}

struct FactorModel {
    std::vector<int> tenor_idx;                    // CMT columns used (default: all but 4M)
    std::vector<double> T;                          // maturities
    std::vector<std::vector<double>> cov;           // daily covariance of yield changes, bp^2
    std::vector<double> evals;                      // bp^2 per day, descending
    std::vector<std::vector<double>> loadings;      // [factor][tenor], unit norm, sign-fixed (level positive, slope rising with maturity)
    int K = 3; size_t n_days = 0; std::string from, to;
    double explained(int k) const { double tot = 0; for (double e : evals) tot += e; double s = 0; for (int i = 0; i < k; ++i) s += evals[i]; return tot > 0 ? s / tot : 0; }

    void fit(const std::vector<CmtRow>& rows, const std::string& d0, const std::string& d1, int k = 3) {
        K = k; tenor_idx = {0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12}; T.clear(); for (int i : tenor_idx) T.push_back(CMT_T[i]);
        std::vector<std::vector<double>> X; const CmtRow* prev = nullptr;
        for (auto& r : rows) {
            if (r.date < d0 || r.date > d1) { prev = &r; continue; }
            bool ok = true; for (int i : tenor_idx) if (std::isnan(r.y[i])) ok = false;
            if (prev && ok) { bool okp = true; for (int i : tenor_idx) if (std::isnan(prev->y[i])) okp = false; if (okp) { std::vector<double> d; for (int i : tenor_idx) d.push_back((r.y[i] - prev->y[i]) * 1e4); X.push_back(d); } }
            prev = &r;
        }
        n_days = X.size(); from = d0; to = d1; size_t n = T.size();
        std::vector<double> mu(n, 0.0); for (auto& x : X) for (size_t i = 0; i < n; ++i) mu[i] += x[i] / X.size();
        cov.assign(n, std::vector<double>(n, 0.0));
        for (auto& x : X) for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) cov[i][j] += (x[i] - mu[i]) * (x[j] - mu[j]) / std::max<size_t>(1, X.size() - 1);
        std::vector<std::vector<double>> V; eig_sym(cov, evals, V);
        loadings.assign(K, std::vector<double>(n));
        for (int f = 0; f < K; ++f) { for (size_t i = 0; i < n; ++i) loadings[f][i] = V[i][f];
            // sign conventions: level positive; slope increasing in maturity; curvature positive in the belly
            double s = 0; if (f == 0) for (double v : loadings[f]) s += v; else if (f == 1) s = loadings[f][n - 1] - loadings[f][0]; else s = loadings[f][n / 2] - 0.5 * (loadings[f][0] + loadings[f][n - 1]);
            if (s < 0) for (double& v : loadings[f]) v = -v; }
    }
    // covariance of yield changes at arbitrary maturities via linear interpolation of the loadings (bp^2/day), K factors
    std::vector<double> loading_at(int f, double Tq) const {
        std::vector<double> l(1); size_t n = T.size();
        if (Tq <= T[0]) l[0] = loadings[f][0]; else if (Tq >= T[n - 1]) l[0] = loadings[f][n - 1];
        else { size_t i = 1; while (T[i] < Tq) ++i; double w = (Tq - T[i - 1]) / (T[i] - T[i - 1]); l[0] = loadings[f][i - 1] * (1 - w) + loadings[f][i] * w; }
        return l;
    }
    Json to_json() const {
        Json j = Json::object(); j["from"] = from; j["to"] = to; j["n_days"] = (long long)n_days; j["T"] = Json(T);
        j["evals_bp2_per_day"] = Json(evals); Json L = Json::array(); for (auto& l : loadings) L.push(Json(l)); j["loadings"] = L;
        Json ex = Json::array(); for (int k = 1; k <= K; ++k) ex.push(explained(k)); j["explained"] = ex;
        Json c = Json::array(); for (auto& r : cov) c.push(Json(r)); j["cov"] = c; return j;
    }
};

} // namespace rcmm
