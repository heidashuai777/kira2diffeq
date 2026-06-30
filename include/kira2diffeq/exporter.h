#ifndef KIRA2DIFFEQ_EXPORTER_H
#define KIRA2DIFFEQ_EXPORTER_H

#include <string>
#include <vector>
#include <ginac/ginac.h>

namespace kira2diffeq {

// Replace d -> 4-2*eps in the matrix, returning a new matrix
std::vector<std::vector<GiNaC::ex>> to_epsilon_form(
    const std::vector<std::vector<GiNaC::ex>>& matrix);

// Expand matrix entries in epsilon up to given order
// (series expansion of rational functions around eps=0)
std::vector<std::vector<GiNaC::ex>> expand_in_epsilon(
    const std::vector<std::vector<GiNaC::ex>>& matrix,
    int order);

// Export DE matrix to Mathematica .m file (importable with Get[])
// Replaces d with 4-2*eps if use_epsilon is true
std::string export_to_mathematica(
    const std::vector<std::vector<GiNaC::ex>>& matrix,
    const std::string& topology_name,
    const std::string& kin_var,
    const std::vector<std::vector<int>>& master_list,
    bool use_epsilon = false);

// Export DE matrices for ALL kinematic variables to a single Mathematica file
std::string export_all_mathematica(
    const std::map<std::string, std::vector<std::vector<GiNaC::ex>>>& matrices,
    const std::string& topology_name,
    const std::vector<std::vector<int>>& master_list,
    bool use_epsilon = false);

// Export DE matrix to plain text
std::string export_to_text(
    const std::vector<std::vector<GiNaC::ex>>& matrix,
    const std::string& topology_name,
    const std::string& kin_var,
    const std::vector<std::vector<int>>& master_list,
    bool use_epsilon = false);

} // namespace kira2diffeq

#endif
