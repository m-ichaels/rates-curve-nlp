#pragma once
// Loaders for the downloaded data: Treasury CMT par yields (bootstrap inputs), NY Fed SOFR / EFFR.
#include "curve.hpp"
#include "json.hpp"
#include <cstdio>
#include <map>
#include <string>

namespace rcmm {

struct CmtRow { std::string date; std::vector<double> y; };   // 13 tenors in decimals, NaN when missing
static const double CMT_T[13] = {1.0 / 12, 2.0 / 12, 0.25, 4.0 / 12, 0.5, 1, 2, 3, 5, 7, 10, 20, 30};
static const char* CMT_NAMES[13] = {"1M", "2M", "3M", "4M", "6M", "1Y", "2Y", "3Y", "5Y", "7Y", "10Y", "20Y", "30Y"};

inline std::vector<CmtRow> read_cmt(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) throw Error("cannot open " + path);
    std::vector<CmtRow> rows; char line[512]; bool first = true;
    while (std::fgets(line, sizeof line, f)) {
        if (first) { first = false; continue; }
        CmtRow r; std::string s(line); size_t p = 0; int col = 0;
        while (p <= s.size()) { size_t q = s.find(',', p); if (q == std::string::npos) q = s.size(); std::string tok = s.substr(p, q - p); while (!tok.empty() && (tok.back() == '\n' || tok.back() == '\r')) tok.pop_back(); if (col == 0) r.date = tok; else r.y.push_back(tok.empty() ? NAN : std::atof(tok.c_str()) / 100.0); ++col; p = q + 1; }
        if (r.y.size() == 13) rows.push_back(r);
    }
    std::fclose(f); return rows;
}

// Instruments from a CMT row: bills (<= 1y) as bond-equivalent simple yields, the rest as semi-annual par bonds.
inline std::vector<Instrument> cmt_instruments(const CmtRow& r, bool drop_4m = true) {
    std::vector<Instrument> v;
    for (int k = 0; k < 13; ++k) { if (std::isnan(r.y[k])) continue; if (drop_4m && k == 3) continue; v.push_back({CMT_T[k] <= 1.0 + 1e-9 ? InstKind::Bill : InstKind::ParBond, CMT_T[k], r.y[k], 2}); }
    return v;
}

struct RateRow { std::string date; double sofr = NAN, effr = NAN, sofr_p99 = NAN, sofr_vol = NAN; };
inline std::map<std::string, RateRow> read_nyfed(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb"); if (!f) return {};
    std::map<std::string, RateRow> out; char line[1024]; bool first = true;
    while (std::fgets(line, sizeof line, f)) {
        if (first) { first = false; continue; }
        std::vector<std::string> c; std::string s(line); size_t p = 0;
        while (p <= s.size()) { size_t q = s.find(',', p); if (q == std::string::npos) q = s.size(); c.push_back(s.substr(p, q - p)); p = q + 1; }
        if (c.size() < 8) continue;
        // Effective Date (MM/DD/YYYY), Rate Type, Rate, 1st, 25th, 75th, 99th, Volume
        std::string d = c[0]; if (d.size() == 10 && d[2] == '/') d = d.substr(6, 4) + "-" + d.substr(0, 2) + "-" + d.substr(3, 2);
        RateRow& r = out[d]; r.date = d;
        if (c[1] == "SOFR") { r.sofr = std::atof(c[2].c_str()) / 100; r.sofr_p99 = std::atof(c[6].c_str()) / 100; r.sofr_vol = std::atof(c[7].c_str()); }
        else if (c[1] == "EFFR") r.effr = std::atof(c[2].c_str()) / 100;
    }
    std::fclose(f); return out;
}

} // namespace rcmm
