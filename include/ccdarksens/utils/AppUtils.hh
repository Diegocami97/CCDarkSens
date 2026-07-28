#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TH1D.h>
#include <nlohmann/json.hpp>

namespace ccdarksens::utils {

/// Flat ε(n_e) histogram filled with eps (clamped to [0,1]) for n_e in [ne_min, ne_max].
inline std::unique_ptr<TH1D> MakeFlatEfficiency(int ne_min, int ne_max, double eps = 1.0) {
    const int nbin = ne_max - ne_min + 1;
    std::vector<double> edges(nbin + 1);
    for (int i = 0; i <= nbin; ++i) edges[i] = (ne_min - 0.5) + i;
    auto h = std::make_unique<TH1D>("eps_ne", "Pattern efficiency; n_{e}; #epsilon",
                                    nbin, edges.data());
    const double v = std::clamp(eps, 0.0, 1.0);
    for (int b = 1; b <= nbin; ++b) h->SetBinContent(b, v);
    return h;
}

/// Sum histogram bins for n_e values in roi (1-based ROOT bin indexing relative to ne_min).
inline double SumROI(const TH1D& h, const std::vector<int>& roi, int ne_min) {
    double s = 0.0;
    for (int ne : roi) {
        const int b = ne - ne_min + 1;
        if (b >= 1 && b <= h.GetNbinsX()) s += h.GetBinContent(b);
    }
    return s;
}

/// Expand a JSON axis specification into a vector of values.
/// Accepts: plain array, {"values": [...]}, {"linspace": {start, stop, num}},
///          {"logspace": {start_exp, stop_exp, num}} or {"logspace": {start, stop, num}}.
inline std::vector<double> ExpandAxis(const nlohmann::json& spec,
                                      const std::string& axis_name = "") {
    std::vector<double> vals;

    auto append = [&](const nlohmann::json& arr) {
        for (const auto& v : arr)
            if (v.is_number()) vals.push_back(v.get<double>());
    };

    if (spec.is_array()) {
        append(spec);
        return vals;
    }

    if (spec.contains("values")) append(spec.at("values"));

    if (spec.contains("linspace")) {
        const auto& p = spec.at("linspace");
        const double a = p.at("start").get<double>();
        const double b = p.at("stop").get<double>();
        const int    n = p.at("num").get<int>();
        const bool endp = p.value("endpoint", true);
        if (n <= 0) throw std::runtime_error("linspace.num must be > 0" +
                                              (axis_name.empty() ? "" : " for " + axis_name));
        if (n == 1) { vals.push_back(a); }
        else {
            const double step = endp ? (b - a) / (n - 1) : (b - a) / n;
            for (int i = 0; i < n; ++i) vals.push_back(a + i * step);
        }
    }

    if (spec.contains("logspace")) {
        const auto& p = spec.at("logspace");
        double aexp, bexp;
        if (p.contains("start_exp") && p.contains("stop_exp")) {
            aexp = p.at("start_exp").get<double>();
            bexp = p.at("stop_exp").get<double>();
        } else {
            const double start = p.at("start").get<double>();
            const double stop  = p.at("stop").get<double>();
            if (start <= 0.0 || stop <= 0.0)
                throw std::runtime_error("logspace start/stop must be > 0" +
                                          (axis_name.empty() ? "" : " for " + axis_name));
            aexp = std::log10(start);
            bexp = std::log10(stop);
        }
        const int  n    = p.at("num").get<int>();
        const bool endp = p.value("endpoint", true);
        if (n <= 0) throw std::runtime_error("logspace.num must be > 0" +
                                              (axis_name.empty() ? "" : " for " + axis_name));
        if (n == 1) { vals.push_back(std::pow(10.0, aexp)); }
        else {
            const double step = endp ? (bexp - aexp) / (n - 1) : (bexp - aexp) / n;
            for (int i = 0; i < n; ++i) vals.push_back(std::pow(10.0, aexp + i * step));
        }
    }

    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    return vals;
}

} // namespace ccdarksens::utils
