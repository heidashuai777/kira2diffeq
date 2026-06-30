#ifndef KIRA2DIFFEQ_CONFIG_READER_H
#define KIRA2DIFFEQ_CONFIG_READER_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <ginac/ginac.h>

namespace kira2diffeq {

// Scalar product symbol: SP(a,b) = a·b where a,b are momentum vectors
// Stored as sp_a_b in GiNaC
inline GiNaC::symbol make_sp_sym(const std::string& a, const std::string& b) {
    // Ensure canonical ordering: shorter name first, then lexicographic
    if (a < b) return GiNaC::symbol("sp_" + a + "_" + b);
    else return GiNaC::symbol("sp_" + b + "_" + a);
}

struct Kinematics {
    std::vector<std::string> incoming_momenta;
    std::vector<std::string> outgoing_momenta;
    std::vector<std::pair<std::string, std::string>> momentum_conservation;
    std::vector<std::pair<std::string, int>> kinematic_invariants;
    // Substitution rules: p_i * p_j -> expression in invariants
    GiNaC::lst sp_rules;
    std::string symbol_to_replace_by_one;
    // GiNaC symbols for kinematic invariants
    std::map<std::string, GiNaC::symbol> symbols;
    // List of all external momenta names
    std::vector<std::string> all_external;
};

struct MomentumTerm {
    int coeff;
    std::string momentum;  // momentum name (e.g., "k1", "q2", "p1")
};

struct Topology {
    std::string name;
    std::vector<std::string> loop_momenta;
    int top_level_sector;
    // Propagators as GiNaC expressions expanded in scalar products:
    // P_i = sum c_{jk} * sp_vj_vk - mass
    GiNaC::lst propagators;
    int num_propagators;
    // Original propagator strings for reference
    std::vector<std::string> propagator_strings;
    std::vector<std::string> mass_strings;
    // Momentum-linear representation: propagators_momentum[i] = list of (coeff, momentum_name)
    // The i-th propagator is (Σ coeff * momentum)² - mass
    // Stored as GiNaC expression using momentum symbols (for diff wrt momenta)
    GiNaC::lst propagator_squares;  // pow(linear_combo, 2) - mass, with momentum symbols
    GiNaC::lst propagator_momentum_exprs;  // just the linear combination (Σ coeff * momentum), for manual diff
    // Pairs (mom_a, mom_b) for bilinear ISP propagators. Empty for standard propagators.
    std::vector<std::pair<std::string, std::string>> bilinear_pairs;
};

struct KiraConfig {
    std::vector<Topology> topologies;
    Kinematics kinematics;
};

KiraConfig read_config(const std::string& config_dir);

} // namespace kira2diffeq

#endif
