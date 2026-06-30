#include "kira2diffeq/config_reader.h"
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <set>

namespace kira2diffeq {

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
            if (!vec.empty()) {
                result.push_back({sign, vec});
            }
        }
    }
    return result;
}

class SymbolRegistry {
    std::map<std::string, GiNaC::symbol> sp_syms_;
    std::map<std::string, GiNaC::symbol> kin_syms_;
    std::map<std::string, GiNaC::symbol> mom_syms_;
public:
    const GiNaC::symbol& get_sp(const std::string& a, const std::string& b) {
        std::string key;
        if (a < b) key = "sp_" + a + "_" + b;
        else key = "sp_" + b + "_" + a;
        auto it = sp_syms_.find(key);
        if (it == sp_syms_.end()) {
            it = sp_syms_.emplace(key, GiNaC::symbol(key)).first;
        }
        return it->second;
    }
    const GiNaC::symbol& get_kin(const std::string& name) {
        auto it = kin_syms_.find(name);
        if (it == kin_syms_.end()) {
            it = kin_syms_.emplace(name, GiNaC::symbol(name)).first;
        }
        return it->second;
    }
    const GiNaC::symbol& get_mom(const std::string& name) {
        auto it = mom_syms_.find(name);
        if (it == mom_syms_.end()) {
            it = mom_syms_.emplace(name, GiNaC::symbol(name)).first;
        }
        return it->second;
    }
    std::map<std::string, GiNaC::symbol> kin_symbols() const { return kin_syms_; }
    std::map<std::string, GiNaC::symbol> mom_symbols() const { return mom_syms_; }
};

static GiNaC::ex build_propagator_expanded(
    const std::string& mom_expr,
    const std::string& mass_expr,
    SymbolRegistry& syms)
{
    auto terms = parse_momentum_terms(mom_expr);
    GiNaC::ex result = 0;

    for (const auto& [s_i, v_i] : terms) {
        int coeff = s_i * s_i;
        if (coeff != 0) result += coeff * syms.get_sp(v_i, v_i);
    }
    for (size_t i = 0; i < terms.size(); i++) {
        for (size_t k = i + 1; k < terms.size(); k++) {
            int coeff = 2 * terms[i].first * terms[k].first;
            if (coeff != 0) result += coeff * syms.get_sp(terms[i].second, terms[k].second);
        }
    }

    if (!mass_expr.empty() && mass_expr != "0") {
        GiNaC::symtab table;
        table["d"] = GiNaC::symbol("d");
        for (const auto& [name, sym] : syms.kin_symbols()) {
            table[name] = sym;
        }
        GiNaC::parser reader(table);
        result -= reader(mass_expr);
    }
    return result;
}

// Build the propagator as pow(linear_momentum_combo, 2) - mass^2
// using momentum symbols (not SP symbols)
static GiNaC::ex build_propagator_squared(
    const std::string& mom_expr,
    const std::string& mass_expr,
    SymbolRegistry& syms)
{
    auto terms = parse_momentum_terms(mom_expr);
    GiNaC::ex lin = 0;
    for (const auto& [coeff, mom_name] : terms) {
        if (coeff == 1) lin += syms.get_mom(mom_name);
        else if (coeff == -1) lin -= syms.get_mom(mom_name);
        else if (coeff > 0) lin += coeff * syms.get_mom(mom_name);
        else lin += coeff * syms.get_mom(mom_name);
    }

    GiNaC::ex prop = GiNaC::pow(lin, 2);

    if (!mass_expr.empty() && mass_expr != "0") {
        GiNaC::symtab table;
        table["d"] = GiNaC::symbol("d");
        for (const auto& [name, sym] : syms.kin_symbols()) {
            table[name] = sym;
        }
        GiNaC::parser reader(table);
        prop -= GiNaC::pow(reader(mass_expr), 2);
    }
    return prop;
}

// Walk an expression and replace products of momentum symbols with SP symbols.
// Also replace momentum^2 with sp(mom, mom).
static GiNaC::ex momenta_to_sp(const GiNaC::ex& expr,
                                const std::set<std::string>& all_momenta,
                                SymbolRegistry& syms)
{
    if (GiNaC::is_a<GiNaC::numeric>(expr)) return expr;
    if (GiNaC::is_a<GiNaC::symbol>(expr)) {
        // Could be a lone momentum symbol from incomplete expansion
        std::string name = GiNaC::ex_to<GiNaC::symbol>(expr).get_name();
        if (all_momenta.count(name)) return 0; // linear term should not appear after proper expansion
        return expr;
    }
    if (GiNaC::is_a<GiNaC::mul>(expr)) {
        // Check if this is a product of two momentum symbols
        GiNaC::symbol s1, s2;
        GiNaC::numeric coeff = 1;
        int n_syms = 0;
        for (size_t i = 0; i < expr.nops(); i++) {
            if (GiNaC::is_a<GiNaC::symbol>(expr.op(i))) {
                if (n_syms == 0) s1 = GiNaC::ex_to<GiNaC::symbol>(expr.op(i));
                else if (n_syms == 1) s2 = GiNaC::ex_to<GiNaC::symbol>(expr.op(i));
                n_syms++;
            } else if (GiNaC::is_a<GiNaC::numeric>(expr.op(i))) {
                coeff = coeff * GiNaC::ex_to<GiNaC::numeric>(expr.op(i));
            }
        }
        if (n_syms == 2 &&
            all_momenta.count(s1.get_name()) &&
            all_momenta.count(s2.get_name())) {
            return coeff * syms.get_sp(s1.get_name(), s2.get_name());
        }
        // Otherwise recurse
        GiNaC::ex result = coeff;
        for (size_t i = 0; i < expr.nops(); i++) {
            if (!GiNaC::is_a<GiNaC::numeric>(expr.op(i))) {
                result = result * momenta_to_sp(expr.op(i), all_momenta, syms);
            }
        }
        return result;
    }
    if (GiNaC::is_a<GiNaC::power>(expr)) {
        GiNaC::ex base = expr.op(0);
        int exp = GiNaC::ex_to<GiNaC::numeric>(expr.op(1)).to_int();
        if (exp == 2 && GiNaC::is_a<GiNaC::symbol>(base)) {
            std::string name = GiNaC::ex_to<GiNaC::symbol>(base).get_name();
            if (all_momenta.count(name)) {
                return syms.get_sp(name, name);
            }
        }
        return GiNaC::pow(momenta_to_sp(base, all_momenta, syms), exp);
    }
    if (GiNaC::is_a<GiNaC::add>(expr)) {
        GiNaC::ex result = 0;
        for (size_t i = 0; i < expr.nops(); i++) {
            result += momenta_to_sp(expr.op(i), all_momenta, syms);
        }
        return result;
    }
    return expr;
}

// Solve c * sp + rest = rhs for sp, return sp == (rhs - rest_without_sp)/c
static std::pair<GiNaC::symbol, GiNaC::ex> solve_for_sp(
    const GiNaC::ex& lhs, const GiNaC::ex& rhs)
{
    std::map<GiNaC::ex, GiNaC::ex, GiNaC::ex_is_less> coeffs;
    GiNaC::ex rest = lhs;

    // Extract SP symbols and their coefficients from the expression
    std::function<void(const GiNaC::ex&, GiNaC::ex)> extract;
    extract = [&](const GiNaC::ex& term, GiNaC::ex coeff) {
        if (GiNaC::is_a<GiNaC::symbol>(term)) {
            coeffs[term] = coeffs[term] + coeff;
        } else if (GiNaC::is_a<GiNaC::mul>(term)) {
            GiNaC::numeric num = 1;
            GiNaC::ex sym_part;
            for (size_t j = 0; j < term.nops(); j++) {
                if (GiNaC::is_a<GiNaC::numeric>(term.op(j))) {
                    num = num * GiNaC::ex_to<GiNaC::numeric>(term.op(j));
                } else {
                    sym_part = term.op(j);
                }
            }
            if (sym_part.is_zero()) {
                // Pure number, no SP symbol
            } else {
                coeffs[sym_part] = coeffs[sym_part] + coeff * num;
            }
        } else if (GiNaC::is_a<GiNaC::add>(term)) {
            for (size_t i = 0; i < term.nops(); i++) {
                extract(term.op(i), coeff);
            }
        }
    };
    extract(lhs, 1);

    for (const auto& [sp, coeff] : coeffs) {
        if (GiNaC::is_a<GiNaC::symbol>(sp) && GiNaC::is_a<GiNaC::numeric>(coeff)) {
            GiNaC::symbol sp_sym = GiNaC::ex_to<GiNaC::symbol>(sp);
            // sp_sym * coeff + (lhs - coeff * sp_sym) = rhs
            // => sp_sym = (rhs - (lhs - coeff * sp_sym)) / coeff
            GiNaC::ex rhs_part = rhs - (lhs - coeff * sp_sym);
            GiNaC::ex result = rhs_part / coeff;
            return {sp_sym, result.expand()};
        }
    }
    return {GiNaC::symbol(), 0};
}

KiraConfig read_config(const std::string& config_dir) {
    KiraConfig config;
    SymbolRegistry syms;

    // --- Read kinematics.yaml ---
    std::string kin_path = config_dir + "/kinematics.yaml";
    YAML::Node kin_file = YAML::LoadFile(kin_path);
    YAML::Node kin_yaml = kin_file["kinematics"];
    if (!kin_yaml || kin_yaml.IsNull()) {
        kin_yaml = kin_file;
    }

    GiNaC::symtab parser_table;
    for (const auto& inv : kin_yaml["kinematic_invariants"]) {
        std::string name = inv[0].as<std::string>();
        int mass_dim = inv[1].as<int>();
        config.kinematics.kinematic_invariants.push_back({name, mass_dim});
        const auto& sym = syms.get_kin(name);
        config.kinematics.symbols[name] = sym;
        parser_table[name] = sym;
    }

    std::set<std::string> all_momenta;
    if (kin_yaml["incoming_momenta"]) {
        for (const auto& m : kin_yaml["incoming_momenta"]) {
            config.kinematics.incoming_momenta.push_back(m.as<std::string>());
            all_momenta.insert(m.as<std::string>());
        }
    }
    if (kin_yaml["outgoing_momenta"]) {
        for (const auto& m : kin_yaml["outgoing_momenta"]) {
            config.kinematics.outgoing_momenta.push_back(m.as<std::string>());
            all_momenta.insert(m.as<std::string>());
        }
    }
    if (kin_yaml["momentum_conservation"]) {
        auto mc = kin_yaml["momentum_conservation"];
        for (size_t i = 0; i < mc.size(); i += 2) {
            config.kinematics.momentum_conservation.push_back({
                mc[i].as<std::string>(), mc[i+1].as<std::string>()});
        }
    }

    for (const auto& m : all_momenta) {
        config.kinematics.all_external.push_back(m);
    }

    if (kin_yaml["symbol_to_replace_by_one"]) {
        config.kinematics.symbol_to_replace_by_one =
            kin_yaml["symbol_to_replace_by_one"].as<std::string>();
    } else {
        config.kinematics.symbol_to_replace_by_one = "";
    }

    // Parse scalar product rules
    if (kin_yaml["scalarproduct_rules"]) {
        struct SPRule {
            GiNaC::ex lhs;
            GiNaC::ex rhs;
            bool is_simple;
        };
        std::vector<SPRule> all_rules;

        for (const auto& rule : kin_yaml["scalarproduct_rules"]) {
            auto lhs = rule[0];
            std::string v1_str = lhs[0].as<std::string>();
            std::string v2_str = lhs[1].as<std::string>();

            // RHS can be integer (0) or string ("s2/2")
            GiNaC::ex rhs;
            auto rhs_node = rule[1];
            if (rhs_node.IsScalar()) {
                // Try as string first
                try {
                    std::string rhs_str = rhs_node.as<std::string>();
                    GiNaC::parser reader(parser_table);
                    rhs = reader(rhs_str);
                } catch (...) {
                    // Integer or other numeric
                    rhs = GiNaC::numeric(rhs_node.as<int>());
                }
            } else {
                rhs = GiNaC::numeric(rhs_node.as<int>());
            }

            auto v1_terms = parse_momentum_terms(v1_str);
            auto v2_terms = parse_momentum_terms(v2_str);

            GiNaC::ex sp_expr = 0;
            bool is_simple = (v1_terms.size() == 1 && v2_terms.size() == 1);

            for (const auto& [a_i, m_i] : v1_terms) {
                for (const auto& [b_j, m_j] : v2_terms) {
                    int coeff = a_i * b_j;
                    if (coeff != 0) {
                        sp_expr += coeff * syms.get_sp(m_i, m_j);
                    }
                }
            }

            all_rules.push_back({sp_expr, rhs, is_simple});
        }

        // Process rules iteratively
        GiNaC::lst& rules_lst = config.kinematics.sp_rules;
        bool changed = true;
        for (int pass = 0; pass < 20 && changed; pass++) {
            changed = false;
            for (auto& spr : all_rules) {
                GiNaC::ex simplified = spr.lhs.subs(rules_lst);
                if (simplified == spr.lhs && spr.is_simple) {
                    rules_lst.append(spr.lhs == spr.rhs);
                    changed = true;
                    continue;
                }
                if (simplified != spr.lhs) {
                    auto [sp, sol] = solve_for_sp(simplified, spr.rhs);
                    if (sp.get_name() != "" && sol != 0) {
                        rules_lst.append(sp == sol);
                        spr.lhs = simplified;
                        changed = true;
                    }
                }
            }
        }

        for (const auto& spr : all_rules) {
            if (spr.is_simple) {
                GiNaC::ex simplified = spr.lhs.subs(rules_lst);
                if (simplified != GiNaC::ex(0)) {
                    rules_lst.append(simplified == spr.rhs);
                }
            }
        }
    }

    // Derive additional SP rules from momentum conservation.
    // For each dependent momentum d = Σ c_i m_i:
    //   sp(d, v) = Σ c_i sp(m_i, v) for any momentum v.
    // If sp(d, v) is known from SP rules, we can solve for unknown sp(m_i, v).
    if (kin_yaml["momentum_conservation"]) {
        auto mc = kin_yaml["momentum_conservation"];
        GiNaC::lst& rules_lst = config.kinematics.sp_rules;

        // Helper: look up known value of an SP symbol
        auto lookup_sp = [&](const GiNaC::symbol& sp_sym) -> GiNaC::ex {
            for (size_t ri = 0; ri < rules_lst.nops(); ri++) {
                GiNaC::ex rule = rules_lst.op(ri);
                if (GiNaC::is_a<GiNaC::symbol>(rule.op(0))) {
                    if (GiNaC::ex_to<GiNaC::symbol>(rule.op(0)).get_name() == sp_sym.get_name()) {
                        return rule.op(1);
                    }
                }
            }
            return sp_sym; // not found, return the symbol itself (unknown)
        };

        // Iterate to derive all possible relations
        for (int iter = 0; iter < 10; iter++) {
            bool any_added = false;

            for (size_t ci = 0; ci < mc.size(); ci += 2) {
                std::string dep_mom = mc[ci].as<std::string>();
                std::string cons_expr = mc[ci+1].as<std::string>();
                auto cons_terms = parse_momentum_terms(cons_expr);

                // For each external momentum v (including dependent ones):
                // sp(dep, v) = Σ c_i * sp(m_i, v)
                // If sp(dep, v) is known and some sp(m_i, v) are unknown,
                // we can solve for one unknown sp(m_i, v).
                // If sp(dep, v) is unknown, we can add it as a derived relation.

                for (const auto& v_name : config.kinematics.all_external) {
                    if (v_name == dep_mom) continue;

                    GiNaC::symbol lhs_sp = syms.get_sp(dep_mom, v_name);
                    GiNaC::ex lhs_val = lookup_sp(lhs_sp);

                    // Build RHS: Σ c_i * sp(m_i, v)
                    GiNaC::ex rhs_sum = 0;
                    std::vector<std::pair<GiNaC::symbol, int>> unknown_terms;
                    for (const auto& [coeff, mom] : cons_terms) {
                        if (coeff == 0) continue;
                        GiNaC::symbol term_sp = syms.get_sp(mom, v_name);
                        GiNaC::ex term_val = lookup_sp(term_sp);
                        if (GiNaC::is_a<GiNaC::symbol>(term_val) &&
                            GiNaC::ex_to<GiNaC::symbol>(term_val).get_name() == term_sp.get_name()) {
                            // Still unknown
                            unknown_terms.push_back({term_sp, coeff});
                            rhs_sum += coeff * term_sp;
                        } else {
                            rhs_sum += coeff * term_val;
                        }
                    }

                    if (GiNaC::is_a<GiNaC::symbol>(lhs_val) &&
                        GiNaC::ex_to<GiNaC::symbol>(lhs_val).get_name() == lhs_sp.get_name()) {
                        // LHS is unknown, RHS might have known terms → add derived rule
                        GiNaC::ex simplified = rhs_sum.subs(rules_lst);
                        if (simplified != 0) {
                            rules_lst.append(lhs_sp == simplified);
                            any_added = true;
                        }
                    } else if (unknown_terms.size() == 1) {
                        // LHS is known, exactly one RHS term is unknown → solve for it
                        auto [unk_sp, unk_coeff] = unknown_terms[0];
                        GiNaC::ex rest = rhs_sum - unk_coeff * unk_sp;
                        GiNaC::ex solved = (lhs_val - rest) / GiNaC::numeric(unk_coeff);
                        solved = solved.subs(rules_lst);
                        if (solved != 0 && !GiNaC::is_a<GiNaC::symbol>(solved)) {
                            rules_lst.append(unk_sp == solved);
                            any_added = true;
                        }
                    }
                }
            }

            if (!any_added) break;
        }
    }

    // --- Read integralfamilies.yaml ---
    std::string fam_path = config_dir + "/integralfamilies.yaml";
    YAML::Node fam_file = YAML::LoadFile(fam_path);
    YAML::Node fam_yaml = fam_file["integralfamilies"];
    if (!fam_yaml || fam_yaml.IsNull()) {
        fam_yaml = fam_file;
    }

    // Build set of all momentum names (external + loop) for SP conversion
    std::set<std::string> all_mom_names = all_momenta;

    for (const auto& fam : fam_yaml) {
        Topology topo;
        topo.name = fam["name"].as<std::string>();

        for (const auto& lm : fam["loop_momenta"]) {
            topo.loop_momenta.push_back(lm.as<std::string>());
            all_mom_names.insert(lm.as<std::string>());
        }

        if (fam["top_level_sectors"].IsSequence()) {
            auto tls = fam["top_level_sectors"][0];
            if (tls.IsScalar()) {
                std::string tls_str = tls.as<std::string>();
                if (!tls_str.empty() && tls_str[0] == 'b') {
                    topo.top_level_sector = std::stoi(tls_str.substr(1), nullptr, 2);
                } else {
                    topo.top_level_sector = std::stoi(tls_str);
                }
            } else {
                topo.top_level_sector = tls.as<int>();
            }
        }

        for (const auto& prop : fam["propagators"]) {
            std::string mom, mass;
            bool is_bilinear = false;
            std::string bilinear_a, bilinear_b;

            if (prop.IsMap() && prop["bilinear"]) {
                // Bilinear (ISP) propagator: { bilinear: [ ["k2", "l1"], 0] }
                is_bilinear = true;
                auto bl = prop["bilinear"];
                bilinear_a = bl[0][0].as<std::string>();
                bilinear_b = bl[0][1].as<std::string>();
                // Mass (second element) is always 0 for ISPs
                if (bl[1].IsScalar()) {
                    mass = bl[1].as<std::string>();
                } else {
                    mass = "0";
                }
                mom = bilinear_a + "*" + bilinear_b;  // placeholder string
            } else {
                mom = prop[0].as<std::string>();
                mass = prop[1].as<std::string>();
            }

            topo.propagator_strings.push_back(mom);
            topo.mass_strings.push_back(mass);

            if (is_bilinear) {
                // ISP: the propagator is sp(a,b) itself (linear in SP), not squared
                GiNaC::ex sp_expr = syms.get_sp(bilinear_a, bilinear_b);
                topo.propagators.append(sp_expr);

                // For propagator_squares: store the SP symbol itself
                topo.propagator_squares.append(sp_expr);

                // Linear combination: 0 (ISP is not a momentum combination)
                topo.propagator_momentum_exprs.append(0);

                // Store bilinear pair for decompose_sp
                topo.bilinear_pairs.push_back({bilinear_a, bilinear_b});
            } else {
                // Standard propagator
                topo.bilinear_pairs.push_back({"", ""});

                // SP-expanded form (for backward compatibility and reduction table)
                GiNaC::ex prop_expr = build_propagator_expanded(mom, mass, syms);
                topo.propagators.append(prop_expr);

                // Momentum-squared form: pow(linear_combo, 2) - mass^2
                GiNaC::ex prop_sq = build_propagator_squared(mom, mass, syms);
                topo.propagator_squares.append(prop_sq);

                // Linear combination: just the momentum sum
                auto terms = parse_momentum_terms(mom);
                GiNaC::ex lin = 0;
                for (const auto& [coeff, mom_name] : terms) {
                    if (coeff == 1) lin += syms.get_mom(mom_name);
                    else if (coeff == -1) lin -= syms.get_mom(mom_name);
                    else if (coeff > 0) lin += coeff * syms.get_mom(mom_name);
                    else lin += coeff * syms.get_mom(mom_name);
                }
                topo.propagator_momentum_exprs.append(lin);
            }
        }

        topo.num_propagators = topo.propagators.nops();
        config.topologies.push_back(topo);
    }

    config.kinematics.symbols = syms.kin_symbols();
    return config;
}

} // namespace kira2diffeq
