#ifndef KIRA2DIFFEQ_MATH_PARSER_H
#define KIRA2DIFFEQ_MATH_PARSER_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <ginac/ginac.h>
#include <functional>

namespace kira2diffeq {

// Integral key: (topology, powers_vector)
using IntegralKey = std::pair<std::string, std::vector<int>>;

// A reduction rule stored with lazy coefficient parsing.
// Coefficients are stored as raw strings and parsed into GiNaC::ex
// only when accessed.
struct LazyReductionRule {
    std::string topology;
    std::vector<int> target_powers;
    // (master_powers, coefficient_string)
    std::vector<std::pair<std::vector<int>, std::string>> rhs_strings;

    // Parse all coefficients using the provided symbol table
    std::vector<std::pair<std::vector<int>, GiNaC::ex>> parsed_rhs(
        const GiNaC::symtab& table) const;
};

// Hash for IntegralKey
struct IntegralKeyHash {
    size_t operator()(const IntegralKey& k) const {
        size_t h = std::hash<std::string>{}(k.first);
        for (int v : k.second) {
            h ^= std::hash<int>{}(v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

// Reduction table using unordered_map for O(1) lookup with lazy parsing
using ReductionTable = std::unordered_map<IntegralKey, LazyReductionRule, IntegralKeyHash>;

struct MathOutput {
    ReductionTable reduction_table;
    // Masters: topology -> list of power vectors
    std::map<std::string, std::vector<std::vector<int>>> masters;
    // Symbol table for parsing coefficients
    GiNaC::symtab sym_table;
};

// Parse kira2math .m output file + masters.final
MathOutput parse_math_output(const std::string& math_file,
                              const std::string& masters_file,
                              const GiNaC::symtab& sym_table = {});

// Parse a single integral string like "topo[1,2,-1,0]" into name + powers
std::pair<std::string, std::vector<int>> parse_integral(const std::string& s);

} // namespace kira2diffeq

#endif
