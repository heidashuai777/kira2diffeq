#ifndef KIRA2DIFFEQ_REDUCER_H
#define KIRA2DIFFEQ_REDUCER_H

#include <string>
#include <vector>
#include <ginac/ginac.h>
#include "kira2diffeq/math_parser.h"
#include "kira2diffeq/differentiator.h"

namespace kira2diffeq {

std::vector<std::pair<std::vector<int>, GiNaC::ex>> reduce_to_masters(
    const IntegralSum& sum,
    const ReductionTable& table,
    const std::vector<std::vector<int>>& master_list,
    const GiNaC::symtab& sym_table);

std::vector<std::vector<GiNaC::ex>> assemble_de_matrix(
    const std::string& topology_name,
    const std::vector<std::vector<int>>& masters,
    const std::string& kin_var,
    const Topology& topo,
    const Kinematics& kin,
    const ReductionTable& table,
    const GiNaC::symtab& sym_table);

} // namespace kira2diffeq
#endif
