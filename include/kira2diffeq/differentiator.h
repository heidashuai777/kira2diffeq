#ifndef KIRA2DIFFEQ_DIFFERENTIATOR_H
#define KIRA2DIFFEQ_DIFFERENTIATOR_H

#include <string>
#include <vector>
#include <ginac/ginac.h>
#include "kira2diffeq/config_reader.h"

namespace kira2diffeq {

// Represents a linear combination of integrals: sum_i coefficient_i * integral_i
struct IntegralSum {
    // Each term: coefficient * topology[powers]
    std::vector<std::pair<std::vector<int>, GiNaC::ex>> terms;
    std::string topology;
};

// Differentiate a master integral wrt a kinematic invariant
// Returns a sum of shifted integrals with coefficients
IntegralSum differentiate_master(
    const std::string& topology_name,
    const std::vector<int>& master_powers,
    const std::string& kin_var,
    const Topology& topo,
    const Kinematics& kin);

} // namespace kira2diffeq

#endif
