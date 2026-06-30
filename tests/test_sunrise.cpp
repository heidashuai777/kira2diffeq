#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include "kira2diffeq/config_reader.h"
#include "kira2diffeq/math_parser.h"
#include "kira2diffeq/differentiator.h"
#include "kira2diffeq/reducer.h"
#include "kira2diffeq/exporter.h"

int main() {
    const std::string fixture_dir = std::string(TEST_FIXTURE_DIR);
    const std::string config_dir = fixture_dir;
    const std::string math_file = fixture_dir + "/reduction.m";
    const std::string masters_file = fixture_dir + "/masters.final";

    auto config = kira2diffeq::read_config(config_dir);
    assert(!config.topologies.empty());
    const auto& topo = config.topologies[0];
    assert(topo.name == "cli_smoke");
    assert(topo.num_propagators == 1);

    GiNaC::symtab sym_table;
    for (const auto& [name, sym] : config.kinematics.symbols) {
        sym_table[name] = sym;
    }
    sym_table["d"] = GiNaC::symbol("d");

    auto math_output = kira2diffeq::parse_math_output(math_file, masters_file, sym_table);
    assert(math_output.masters.count("cli_smoke") == 1);
    assert(math_output.masters["cli_smoke"].size() == 1);
    const auto& masters = math_output.masters["cli_smoke"];
    assert(!math_output.reduction_table.empty());

    auto deriv = kira2diffeq::differentiate_master(topo.name, masters[0], "t", topo, config.kinematics);
    auto reduced = kira2diffeq::reduce_to_masters(deriv, math_output.reduction_table, masters, sym_table);
    assert(!reduced.empty());
    assert(reduced[0].first == masters[0]);

    auto matrix = kira2diffeq::assemble_de_matrix(topo.name, masters, "t", topo, config.kinematics,
                                                 math_output.reduction_table, sym_table);
    assert(matrix.size() == 1 && matrix[0].size() == 1);
    assert(matrix[0][0] != 0);

    auto text = kira2diffeq::export_to_text(matrix, topo.name, "t", masters);
    assert(text.find("d/dt MI_i =") != std::string::npos);

    std::cout << "test_sunrise passed\n";
    return 0;
}
