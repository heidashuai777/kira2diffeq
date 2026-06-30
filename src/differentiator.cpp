#include "kira2diffeq/differentiator.h"
#include <stdexcept>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>

namespace kira2diffeq {

namespace {
bool debug_enabled() {
    static bool initialized = false;
    static bool enabled = false;
    if (!initialized) {
        const char* raw = std::getenv("KIRA2DIFFEQ_DEBUG");
        std::string raw_str = raw ? std::string(raw) : "";
        enabled = !raw_str.empty() && raw_str != "0" && raw_str != "false" && raw_str != "False";
        initialized = true;
    }
    return enabled;
}

} // namespace

#define KIRA2DIFFEQ_DEBUG(block) do { if (debug_enabled()) { block; } } while (0)

static std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && (s[start] == ' ' || s[start] == '\t')) start++;
    size_t end = s.size();
    while (end > start && (s[end-1] == ' ' || s[end-1] == '\t')) end--;
    return s.substr(start, end - start);
}

static std::vector<std::pair<int, std::string>> parse_momentum_terms(const std::string& s) {
    std::vector<std::pair<int, std::string>> result;
    std::string expr = s;
    if (expr.empty()) return result;
    if (expr[0] != '-' && expr[0] != '+') expr = "+" + expr;
    size_t j = 0;
    while (j < expr.size()) {
        if (expr[j] == '+' || expr[j] == '-') {
            int sign = (expr[j] == '+') ? 1 : -1;
            j++;
            size_t start = j;
            while (j < expr.size() && expr[j] != '+' && expr[j] != '-') j++;
            std::string vec = trim(expr.substr(start, j - start));
            if (!vec.empty()) result.push_back({sign, vec});
        }
    }
    return result;
}

struct SPDecomp {
    std::vector<std::pair<int, GiNaC::ex>> prop_terms; // (prop_index, coefficient)
    GiNaC::ex constant;
};

// Cycle detection: set of (a,b) pairs currently being decomposed
static std::set<std::pair<std::string, std::string>> decomposing;

static SPDecomp decompose_sp_impl(
    const std::string& mom_a, const std::string& mom_b,
    const Topology& topo, const Kinematics& kin,
    const std::set<std::string>& loop_momenta);

static SPDecomp decompose_sp(
    const std::string& mom_a, const std::string& mom_b,
    const Topology& topo, const Kinematics& kin,
    const std::set<std::string>& loop_momenta)
{
    std::string key_a = mom_a < mom_b ? mom_a : mom_b;
    std::string key_b = mom_a < mom_b ? mom_b : mom_a;
    auto key = std::make_pair(key_a, key_b);
    if (decomposing.count(key)) {
        // Cycle detected - this SP can't be expressed in terms of propagators alone
        // Return it as an unresolved symbol
        SPDecomp result;
        return result;
    }
    decomposing.insert(key);
    SPDecomp result = decompose_sp_impl(mom_a, mom_b, topo, kin, loop_momenta);
    decomposing.erase(key);
    return result;
}

static SPDecomp decompose_sp_impl(
    const std::string& mom_a, const std::string& mom_b,
    const Topology& topo, const Kinematics& kin,
    const std::set<std::string>& loop_momenta)
{
    SPDecomp result;
    result.constant = 0;

    int np = topo.num_propagators;

    // Check for bilinear ISP propagator that directly represents this SP
    for (int ip = 0; ip < np; ip++) {
        if (static_cast<int>(topo.bilinear_pairs.size()) > ip) {
            const auto& bp = topo.bilinear_pairs[ip];
            if ((bp.first == mom_a && bp.second == mom_b) ||
                (bp.first == mom_b && bp.second == mom_a)) {
                result.prop_terms.push_back({ip, 1});
                return result;
            }
        }
    }

    bool a_is_loop = loop_momenta.count(mom_a);
    bool b_is_loop = loop_momenta.count(mom_b);

    // External-external: look up from SP rules directly
    if (!a_is_loop && !b_is_loop) {
        std::string sp_name;
        if (mom_a < mom_b)
            sp_name = "sp_" + mom_a + "_" + mom_b;
        else
            sp_name = "sp_" + mom_b + "_" + mom_a;

        for (size_t ri = 0; ri < kin.sp_rules.nops(); ri++) {
            GiNaC::ex rule = kin.sp_rules.op(ri);
            if (GiNaC::is_a<GiNaC::symbol>(rule.op(0))) {
                if (GiNaC::ex_to<GiNaC::symbol>(rule.op(0)).get_name() == sp_name) {
                    result.constant = rule.op(1);
                    return result;
                }
            }
        }
        return result; // Not in rules → 0
    }

    // Loop-loop or loop-external
    std::string loop_m = a_is_loop ? mom_a : mom_b;
    std::string other_m = a_is_loop ? mom_b : mom_a;

    // Case 1: loop_m == other_m (loop-loop SP)
    // First, check for a propagator that is exactly this loop momentum
    if (loop_m == other_m) {
        for (int ip = 0; ip < np; ip++) {
            auto terms = parse_momentum_terms(topo.propagator_strings[ip]);
            if (terms.size() == 1 && terms[0].second == loop_m) {
                // P_ip = (c * loop_m)^2 = c^2 * sp(loop_m, loop_m)
                // sp(loop_m, loop_m) = P_ip / c^2
                int c = terms[0].first;
                if (c * c == 1) {
                    result.prop_terms.push_back({ip, 1});
                } else {
                    result.prop_terms.push_back({ip,
                        GiNaC::numeric(1) / GiNaC::numeric(c * c)});
                }
                return result;
            }
        }
    }

    // Case 2: Find a propagator containing BOTH loop_m and other_m
    for (int ip = 0; ip < np; ip++) {
        auto terms = parse_momentum_terms(topo.propagator_strings[ip]);
        bool has_loop = false, has_other = false;
        int c_loop = 0, c_other = 0;
        for (const auto& [c, m] : terms) {
            if (m == loop_m) { has_loop = true; c_loop = c; }
            if (m == other_m) { has_other = true; c_other = c; }
        }

        if (has_loop && has_other && c_loop != 0 && c_other != 0) {
            // P = c_loop^2 sp(loop,loop) + c_other^2 sp(other,other)
            //   + 2*c_loop*c_other sp(loop,other) + cross_terms + self_terms - mass^2
            // Solve for sp(loop, other)
            int denom = 2 * c_loop * c_other;
            GiNaC::ex coeff = GiNaC::numeric(1) / GiNaC::numeric(denom);
            result.prop_terms.push_back({ip, coeff});

            // Subtract all contributions except the target SP
            for (const auto& [c_i, m_i] : terms) {
                if (c_i == 0) continue;
                if ((m_i == loop_m && m_i == other_m) ||  // loop==other case
                    (m_i == loop_m) || (m_i == other_m)) {
                    // For m_i == loop_m: subtract c_loop^2 sp(loop,loop)
                    if (m_i == loop_m) {
                        GiNaC::ex sub_c = GiNaC::numeric(-c_i * c_i) * coeff;
                        if (loop_m != other_m) {  // only if not already handling loop==other
                            auto sub = decompose_sp(loop_m, loop_m, topo, kin, loop_momenta);
                            for (const auto& [idx, c] : sub.prop_terms)
                                result.prop_terms.push_back({idx, sub_c * c});
                            result.constant += sub_c * sub.constant;
                        }
                    }
                    // For m_i == other_m: subtract c_other^2 sp(other,other)
                    if (m_i == other_m && loop_m != other_m) {
                        GiNaC::ex sub_c = GiNaC::numeric(-c_i * c_i) * coeff;
                        auto sub = decompose_sp(other_m, other_m, topo, kin, loop_momenta);
                        for (const auto& [idx, c] : sub.prop_terms)
                            result.prop_terms.push_back({idx, sub_c * c});
                        result.constant += sub_c * sub.constant;
                    }
                    continue;
                }

                // Other self-terms
                GiNaC::ex sub_c = GiNaC::numeric(-c_i * c_i) * coeff;
                auto sub = decompose_sp(m_i, m_i, topo, kin, loop_momenta);
                for (const auto& [idx, c] : sub.prop_terms)
                    result.prop_terms.push_back({idx, sub_c * c});
                result.constant += sub_c * sub.constant;
            }

            // Cross-terms (i < j) excluding the target (loop, other) pair
            for (size_t i = 0; i < terms.size(); i++) {
                for (size_t j = i + 1; j < terms.size(); j++) {
                    bool is_target = (terms[i].second == loop_m && terms[j].second == other_m) ||
                                     (terms[i].second == other_m && terms[j].second == loop_m);
                    if (is_target) continue;

                    int c_prod = 2 * terms[i].first * terms[j].first;
                    if (c_prod == 0) continue;
                    GiNaC::ex sub_c = GiNaC::numeric(-c_prod) * coeff;
                    auto sub = decompose_sp(terms[i].second, terms[j].second,
                                            topo, kin, loop_momenta);
                    for (const auto& [idx, c] : sub.prop_terms)
                        result.prop_terms.push_back({idx, sub_c * c});
                    result.constant += sub_c * sub.constant;
                }
            }

            // Add mass^2 / denom
            std::string mass_str = topo.mass_strings[ip];
            if (!mass_str.empty() && mass_str != "0") {
                GiNaC::symtab table;
                for (const auto& [name, sym] : kin.symbols) table[name] = sym;
                GiNaC::parser reader(table);
                result.constant += coeff * GiNaC::pow(reader(mass_str), 2);
            }
            return result;
        }
    }

    // Case 3: loop-loop SP but no single-momentum propagator found.
    // Use any propagator containing the loop momentum.
    if (loop_m == other_m) {
        for (int ip = 0; ip < np; ip++) {
            auto terms = parse_momentum_terms(topo.propagator_strings[ip]);
            int c_loop = 0;
            for (const auto& [c, m] : terms) {
                if (m == loop_m) { c_loop = c; break; }
            }
            if (c_loop == 0) continue;

            // P = c_loop^2 sp(loop,loop) + ... - mass^2
            // sp(loop,loop) = (P + mass^2 - all_other_terms) / c_loop^2
            GiNaC::ex coeff = GiNaC::numeric(1) / GiNaC::numeric(c_loop * c_loop);
            result.prop_terms.push_back({ip, coeff});

            for (const auto& [c_i, m_i] : terms) {
                if (m_i == loop_m) continue;
                GiNaC::ex sub_c = GiNaC::numeric(-c_i * c_i) * coeff;
                auto sub = decompose_sp(m_i, m_i, topo, kin, loop_momenta);
                for (const auto& [idx, c] : sub.prop_terms)
                    result.prop_terms.push_back({idx, sub_c * c});
                result.constant += sub_c * sub.constant;
            }
            for (size_t i = 0; i < terms.size(); i++) {
                for (size_t j = i + 1; j < terms.size(); j++) {
                    if (terms[i].second == loop_m || terms[j].second == loop_m) continue;
                    int c_prod = 2 * terms[i].first * terms[j].first;
                    if (c_prod == 0) continue;
                    GiNaC::ex sub_c = GiNaC::numeric(-c_prod) * coeff;
                    auto sub = decompose_sp(terms[i].second, terms[j].second,
                                            topo, kin, loop_momenta);
                    for (const auto& [idx, c] : sub.prop_terms)
                        result.prop_terms.push_back({idx, sub_c * c});
                    result.constant += sub_c * sub.constant;
                }
            }

            std::string mass_str = topo.mass_strings[ip];
            if (!mass_str.empty() && mass_str != "0") {
                GiNaC::symtab table;
                for (const auto& [name, sym] : kin.symbols) table[name] = sym;
                GiNaC::parser reader(table);
                result.constant += coeff * GiNaC::pow(reader(mass_str), 2);
            }
            return result;
        }
    }

    return result;
}

// --- Gram matrix construction ---

struct GramInfo {
    std::vector<std::string> basis_momenta;
    GiNaC::matrix G;
    GiNaC::matrix G_inv;
    GiNaC::matrix dG_dx;
};

static GramInfo build_gram(const Kinematics& kin, const std::string& kin_var) {
    GramInfo info;

    std::set<std::string> dependent;
    for (const auto& [dep, expr] : kin.momentum_conservation) {
        dependent.insert(dep);
    }

    for (const auto& m : kin.all_external) {
        if (!dependent.count(m)) {
            info.basis_momenta.push_back(m);
        }
    }

    int nb = info.basis_momenta.size();
    if (nb == 0) return info;

    info.G = GiNaC::matrix(nb, nb);
    info.dG_dx = GiNaC::matrix(nb, nb);

    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < nb; j++) {
            std::string sp_name;
            if (info.basis_momenta[i] < info.basis_momenta[j])
                sp_name = "sp_" + info.basis_momenta[i] + "_" + info.basis_momenta[j];
            else
                sp_name = "sp_" + info.basis_momenta[j] + "_" + info.basis_momenta[i];

            GiNaC::ex val = 0;
            for (size_t ri = 0; ri < kin.sp_rules.nops(); ri++) {
                GiNaC::ex rule = kin.sp_rules.op(ri);
                GiNaC::ex lhs = rule.op(0);
                if (GiNaC::is_a<GiNaC::symbol>(lhs)) {
                    if (GiNaC::ex_to<GiNaC::symbol>(lhs).get_name() == sp_name) {
                        val = rule.op(1);
                        break;
                    }
                }
            }
            info.G(i, j) = val;

            GiNaC::symbol x = kin.symbols.at(kin_var);
            info.dG_dx(i, j) = val.diff(x);
        }
    }

    GiNaC::ex det = info.G.determinant();
    if (det != 0) {
        info.G_inv = info.G.inverse();
    }

    return info;
}

// --- Main differentiation ---

IntegralSum differentiate_master(
    const std::string& topology_name,
    const std::vector<int>& master_powers,
    const std::string& kin_var,
    const Topology& topo,
    const Kinematics& kin)
{
    IntegralSum result;
    result.topology = topology_name;

    int n = topo.num_propagators;
    if (static_cast<int>(master_powers.size()) != n) {
        throw std::runtime_error("Master powers size doesn't match number of propagators");
    }

    GiNaC::symbol x = kin.symbols.at(kin_var);
    std::set<std::string> loop_mom(topo.loop_momenta.begin(), topo.loop_momenta.end());

    GramInfo gram = build_gram(kin, kin_var);
    int nb = gram.basis_momenta.size();

    KIRA2DIFFEQ_DEBUG({
        std::cerr << "DEBUG differentiate_master: topology=" << topology_name
                  << " var=" << kin_var << " powers=[";
        for (size_t i = 0; i < master_powers.size(); i++) {
            if (i > 0) std::cerr << ",";
            std::cerr << master_powers[i];
        }
        std::cerr << "]" << std::endl;
        std::cerr << "DEBUG gram: nb=" << nb << " G_inv.nops=" << gram.G_inv.nops() << std::endl;
    });

    for (int ip = 0; ip < n; ip++) {
        int power = master_powers[ip];
        if (power == 0) continue;

        auto terms = parse_momentum_terms(topo.propagator_strings[ip]);
        KIRA2DIFFEQ_DEBUG({
            std::cerr << "DEBUG ip=" << ip << " power=" << power
                      << " prop_str=" << topo.propagator_strings[ip]
                      << " n_terms=" << terms.size() << std::endl;
        });

        // Contribution 1: external-external SPs (direct)
        GiNaC::ex prop_with_rules = topo.propagators[ip].subs(kin.sp_rules);
        GiNaC::ex dP_dx = prop_with_rules.diff(x);
        dP_dx = dP_dx.subs(kin.sp_rules);
        KIRA2DIFFEQ_DEBUG({
            std::cerr << "DEBUG ip=" << ip << " Contribution1 dP_dx=" << dP_dx << std::endl;
        });

        // Contribution 2: mass term
        std::string mass_str = topo.mass_strings[ip];
        if (!mass_str.empty() && mass_str != "0") {
            GiNaC::symtab table;
            for (const auto& [name, sym] : kin.symbols) table[name] = sym;
            table["d"] = GiNaC::symbol("d");
            GiNaC::parser reader(table);
            GiNaC::ex mass = reader(mass_str);
            dP_dx -= 2 * mass * mass.diff(x);
            dP_dx = dP_dx.subs(kin.sp_rules);
        }

        // Contribution 3: loop-external SPs through Gram matrix chain rule
        if (nb > 0 && gram.G_inv.nops() > 0) {
            KIRA2DIFFEQ_DEBUG({
                std::cerr << "DEBUG ip=" << ip << " Entering Contribution 3" << std::endl;
            });
            for (const auto& [c_l, mom_l] : terms) {
                if (!loop_mom.count(mom_l)) continue;
                if (c_l == 0) continue;
                KIRA2DIFFEQ_DEBUG({
                    std::cerr << "DEBUG   c_l=" << c_l << " mom_l=" << mom_l << std::endl;
                });

                for (const auto& [c_e, mom_e] : terms) {
                    if (loop_mom.count(mom_e)) continue;
                    if (c_e == 0) continue;
                    KIRA2DIFFEQ_DEBUG({
                        std::cerr << "DEBUG     c_e=" << c_e << " mom_e=" << mom_e << std::endl;
                    });

                    int coeff_in_P = 2 * c_l * c_e;
                    KIRA2DIFFEQ_DEBUG({
                        std::cerr << "DEBUG     coeff_in_P=" << coeff_in_P << std::endl;
                    });
                    if (coeff_in_P == 0) continue;

                    int e_idx = -1;
                    for (int bi = 0; bi < nb; bi++) {
                        if (gram.basis_momenta[bi] == mom_e) {
                            e_idx = bi; break;
                        }
                    }
                    KIRA2DIFFEQ_DEBUG({
                        std::cerr << "DEBUG     e_idx=" << e_idx << std::endl;
                    });
                    if (e_idx < 0) continue;

                    for (int c = 0; c < nb; c++) {
                        GiNaC::ex sum_over_b = 0;
                        for (int b = 0; b < nb; b++) {
                            sum_over_b += gram.G_inv(e_idx, b) * gram.dG_dx(b, c);
                        }

                        KIRA2DIFFEQ_DEBUG({
                            std::cerr << "DEBUG     c=" << c << " basis=" << gram.basis_momenta[c]
                                      << " sum_over_b=" << sum_over_b << std::endl;
                        });
                        if (sum_over_b == 0) continue;

                        GiNaC::ex contrib = GiNaC::numeric(coeff_in_P)
                            * sum_over_b / GiNaC::numeric(2);
                        KIRA2DIFFEQ_DEBUG({
                            std::cerr << "DEBUG     contrib=" << contrib << std::endl;
                        });

                        auto decomp = decompose_sp(mom_l, gram.basis_momenta[c],
                                                   topo, kin, loop_mom);
                        KIRA2DIFFEQ_DEBUG({
                            std::cerr << "DEBUG     decomp: n_prop_terms=" << decomp.prop_terms.size()
                                      << " constant=" << decomp.constant << std::endl;
                        });

                        // Constant part
                        GiNaC::ex const_part = contrib * decomp.constant;
                        if (const_part != 0) {
                            dP_dx += const_part;
                        }

                        // Propagator parts → integral shifts
                        for (const auto& [p_idx, p_coeff] : decomp.prop_terms) {
                            GiNaC::ex total_coeff = -power * contrib * p_coeff;

                            std::vector<int> shifted = master_powers;
                            shifted[ip]++;         // from original ∂P/∂x
                            shifted[p_idx]--;      // from SP decomposition

                            KIRA2DIFFEQ_DEBUG({
                                std::cerr << "DEBUG     term: p_idx=" << p_idx << " p_coeff=" << p_coeff
                                          << " total_coeff=" << total_coeff << " shifted=[";
                                for (size_t si = 0; si < shifted.size(); si++) {
                                    if (si > 0) std::cerr << ",";
                                    std::cerr << shifted[si];
                                }
                                std::cerr << "]" << std::endl;
                            });

                            result.terms.push_back({shifted, total_coeff});
                        }
                    }
                }
            }
        }

        dP_dx = dP_dx.subs(kin.sp_rules);
        KIRA2DIFFEQ_DEBUG({
            std::cerr << "DEBUG ip=" << ip << " final dP_dx=" << dP_dx << std::endl;
        });
        if (dP_dx != 0) {
            std::vector<int> shifted_powers = master_powers;
            shifted_powers[ip]++;
            KIRA2DIFFEQ_DEBUG({
                std::cerr << "DEBUG   adding dP_dx term with shifted[ip]++ shifted=[";
                for (size_t si = 0; si < shifted_powers.size(); si++) {
                    if (si > 0) std::cerr << ",";
                    std::cerr << shifted_powers[si];
                }
                std::cerr << "] coeff=" << -power * dP_dx << std::endl;
            });
            result.terms.push_back({shifted_powers, -power * dP_dx});
        }
    }

    KIRA2DIFFEQ_DEBUG({
        std::cerr << "DEBUG result.terms.size=" << result.terms.size() << std::endl;
    });
    return result;
}

} // namespace kira2diffeq
