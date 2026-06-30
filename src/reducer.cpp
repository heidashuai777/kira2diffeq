#include "kira2diffeq/reducer.h"
#include <iostream>

namespace kira2diffeq {

static bool is_master(const std::vector<int>& powers,
                       const std::vector<std::vector<int>>& masters) {
    for (const auto& m : masters) {
        if (m == powers) return true;
    }
    return false;
}

std::vector<std::pair<std::vector<int>, GiNaC::ex>> reduce_to_masters(
    const IntegralSum& sum,
    const ReductionTable& table,
    const std::vector<std::vector<int>>& master_list,
    const GiNaC::symtab& sym_table)
{
    std::vector<std::pair<std::vector<int>, GiNaC::ex>> result;

    for (const auto& [shifted_powers, coeff] : sum.terms) {
        IntegralKey key = {sum.topology, shifted_powers};
        auto it = table.find(key);

        if (it != table.end()) {
            auto parsed = it->second.parsed_rhs(sym_table);
            for (const auto& [master_powers, rhs_coeff] : parsed) {
                GiNaC::ex total_coeff = coeff * rhs_coeff;
                bool found = false;
                for (auto& [mp, c] : result) {
                    if (mp == master_powers) {
                        c += total_coeff;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.push_back({master_powers, total_coeff});
                }
            }
        } else if (is_master(shifted_powers, master_list)) {
            GiNaC::ex total_coeff = coeff;
            bool found = false;
            for (auto& [mp, c] : result) {
                if (mp == shifted_powers) {
                    c += total_coeff;
                    found = true;
                    break;
                }
            }
            if (!found) {
                result.push_back({shifted_powers, total_coeff});
            }
        } else {
            std::cerr << "Warning: integral " << sum.topology << "[";
            for (size_t i = 0; i < shifted_powers.size(); i++) {
                if (i > 0) std::cerr << ",";
                std::cerr << shifted_powers[i];
            }
            std::cerr << "] not found in reduction table" << std::endl;
        }
    }

    return result;
}

std::vector<std::vector<GiNaC::ex>> assemble_de_matrix(
    const std::string& topology_name,
    const std::vector<std::vector<int>>& masters,
    const std::string& kin_var,
    const Topology& topo,
    const Kinematics& kin,
    const ReductionTable& table,
    const GiNaC::symtab& sym_table)
{
    int n = masters.size();
    std::vector<std::vector<GiNaC::ex>> matrix(n, std::vector<GiNaC::ex>(n, 0));

    for (int row = 0; row < n; row++) {
        IntegralSum dMI = differentiate_master(topology_name, masters[row],
                                                kin_var, topo, kin);
        auto reduced = reduce_to_masters(dMI, table, masters, sym_table);
        for (const auto& [master_p, coeff] : reduced) {
            int col = -1;
            for (int j = 0; j < n; j++) {
                if (masters[j] == master_p) { col = j; break; }
            }
            if (col >= 0) matrix[row][col] = coeff;
        }
    }

    return matrix;
}

} // namespace kira2diffeq
