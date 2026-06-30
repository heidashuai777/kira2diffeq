#include "kira2diffeq/exporter.h"
#include <sstream>
#include <fstream>
#include <set>

namespace kira2diffeq {

static GiNaC::symbol eps_sym("eps");

// GiNaC symbols with the same name from different compilation units are
// different objects. Find all symbols named "d" in the expression and
// substitute each one with 4 - 2*eps.
static GiNaC::ex sub_d_to_eps(const GiNaC::ex& expr) {
    if (expr == 0) return expr;
    GiNaC::ex result = expr;

    // Get all symbols recursively via GiNaC's internal traversal
    // Use a simple recursive approach
    std::function<GiNaC::ex(const GiNaC::ex&)> replacer;
    replacer = [&](const GiNaC::ex& e) -> GiNaC::ex {
        if (GiNaC::is_a<GiNaC::numeric>(e)) return e;
        if (GiNaC::is_a<GiNaC::symbol>(e)) {
            if (GiNaC::ex_to<GiNaC::symbol>(e).get_name() == "d") {
                return 4 - 2 * eps_sym;
            }
            return e;
        }
        // Recurse into subexpressions
        int n = e.nops();
        if (n == 0) return e;
        GiNaC::exvector new_ops;
        new_ops.reserve(n);
        for (int i = 0; i < n; i++) {
            new_ops.push_back(replacer(e.op(i)));
        }
        // Rebuild expression from same type + new operands
        // Use ex::subs approach instead
        return e;
    };

    // Use a simpler approach: collect all unique 'd' symbols, then subs each
    std::set<GiNaC::ex, GiNaC::ex_is_less> all_syms;
    std::function<void(const GiNaC::ex&)> collect;
    collect = [&](const GiNaC::ex& e) {
        if (GiNaC::is_a<GiNaC::symbol>(e)) {
            if (GiNaC::ex_to<GiNaC::symbol>(e).get_name() == "d") {
                all_syms.insert(e);
            }
            return;
        }
        for (int i = 0; i < e.nops(); i++) {
            collect(e.op(i));
        }
    };
    collect(expr);

    for (const auto& d_sym : all_syms) {
        result = result.subs(d_sym == 4 - 2 * eps_sym);
    }
    return result;
}

static std::string format_master(const std::string& topo,
                                  const std::vector<int>& powers) {
    std::ostringstream oss;
    oss << topo << "[";
    for (size_t i = 0; i < powers.size(); i++) {
        if (i > 0) oss << ",";
        oss << powers[i];
    }
    oss << "]";
    return oss.str();
}

std::vector<std::vector<GiNaC::ex>> to_epsilon_form(
    const std::vector<std::vector<GiNaC::ex>>& matrix)
{
    int n = matrix.size();
    std::vector<std::vector<GiNaC::ex>> result(n, std::vector<GiNaC::ex>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = sub_d_to_eps(matrix[i][j]);
        }
    }
    return result;
}

std::vector<std::vector<GiNaC::ex>> expand_in_epsilon(
    const std::vector<std::vector<GiNaC::ex>>& matrix,
    int order)
{
    int n = matrix.size();
    auto eps_matrix = to_epsilon_form(matrix);
    std::vector<std::vector<GiNaC::ex>> result(n, std::vector<GiNaC::ex>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = eps_matrix[i][j].series(eps_sym, order);
        }
    }
    return result;
}

std::string export_to_mathematica(
    const std::vector<std::vector<GiNaC::ex>>& matrix,
    const std::string& topology_name,
    const std::string& kin_var,
    const std::vector<std::vector<int>>& master_list,
    bool use_epsilon)
{
    std::ostringstream oss;
    int n = matrix.size();

    oss << "(* DE matrix for " << topology_name << " wrt " << kin_var << " *)\n";
    oss << "(* d/d" << kin_var << " MI_i = Sum[M" << kin_var << "[[i+1,j+1]] * MI_j, {j,1," << n << "}] *)\n\n";

    oss << "M" << kin_var << " = {\n";
    for (int i = 0; i < n; i++) {
        oss << "  {";
        for (int j = 0; j < n; j++) {
            if (j > 0) oss << ", ";
            GiNaC::ex entry = use_epsilon ? sub_d_to_eps(matrix[i][j]) : matrix[i][j];
            oss << entry;
        }
        if (i < n - 1) oss << "},\n";
        else oss << "}\n";
    }
    oss << "};\n\n";

    oss << "MImasters = {\n";
    for (int i = 0; i < n; i++) {
        oss << "  \"" << format_master(topology_name, master_list[i]) << "\"";
        if (i < n - 1) oss << ",\n";
        else oss << "\n";
    }
    oss << "};\n";

    return oss.str();
}

std::string export_all_mathematica(
    const std::map<std::string, std::vector<std::vector<GiNaC::ex>>>& matrices,
    const std::string& topology_name,
    const std::vector<std::vector<int>>& master_list,
    bool use_epsilon)
{
    std::ostringstream oss;
    int n = master_list.size();

    oss << "(* ======================================================== *)\n";
    oss << "(* DE system for topology: " << topology_name << " *)\n";
    if (use_epsilon) oss << "(* Dimension: d = 4 - 2*eps *)\n";
    oss << "(* ======================================================== *)\n\n";

    oss << "masters = {\n";
    for (int i = 0; i < n; i++) {
        oss << "  \"" << format_master(topology_name, master_list[i]) << "\"";
        if (i < n - 1) oss << ",";
        oss << "\n";
    }
    oss << "};\n\n";
    oss << "nMasters = Length[masters];\n\n";

    for (const auto& [kin_var, matrix] : matrices) {
        oss << "(* --- d/d" << kin_var << " --- *)\n";
        oss << "M" << kin_var << " = {\n";
        for (int i = 0; i < n; i++) {
            oss << "  {";
            for (int j = 0; j < n; j++) {
                if (j > 0) oss << ", ";
                GiNaC::ex entry = use_epsilon ? sub_d_to_eps(matrix[i][j]) : matrix[i][j];
                oss << entry;
            }
            if (i < n - 1) oss << "},\n";
            else oss << "}\n";
        }
        oss << "};\n\n";
    }

    oss << "(* Use with: DiffExp, CANONICA, epsilon, Libra, SeaSyde, AMFlow *)\n";
    return oss.str();
}

std::string export_to_text(
    const std::vector<std::vector<GiNaC::ex>>& matrix,
    const std::string& topology_name,
    const std::string& kin_var,
    const std::vector<std::vector<int>>& master_list,
    bool use_epsilon)
{
    std::ostringstream oss;
    int n = matrix.size();

    oss << "d/d" << kin_var << " MI_i = sum_j M_" << kin_var << "[i][j] * MI_j\n\n";
    oss << "Masters:\n";
    for (int i = 0; i < n; i++) {
        oss << "  MI_" << i << " = "
            << format_master(topology_name, master_list[i]) << "\n";
    }
    oss << "\n";

    for (int i = 0; i < n; i++) {
        oss << "dMI_" << i << "/d" << kin_var << " = ";
        bool first = true;
        for (int j = 0; j < n; j++) {
            GiNaC::ex entry = use_epsilon ? sub_d_to_eps(matrix[i][j]) : matrix[i][j];
            if (entry == 0) continue;
            if (!first) oss << " + ";
            oss << "(" << entry << ")*MI_" << j;
            first = false;
        }
        if (first) oss << "0";
        oss << "\n";
    }

    return oss.str();
}

} // namespace kira2diffeq
