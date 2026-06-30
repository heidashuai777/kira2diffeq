#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include "kira2diffeq/config_reader.h"
#include "kira2diffeq/math_parser.h"
#include "kira2diffeq/differentiator.h"
#include "kira2diffeq/reducer.h"
#include "kira2diffeq/exporter.h"

static void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  -c <dir>   Config directory (integralfamilies.yaml, kinematics.yaml)\n"
              << "  -m <file>  kira2math .m output file\n"
              << "  -f <file>  masters.final file\n"
              << "  -t <name>  Topology name\n"
              << "  -v <var>   Kinematic variable (default: all)\n"
              << "  -o <fmt>   Output: text | math | all (default: text)\n"
              << "  -e         Use epsilon (d = 4 - 2*eps) in output\n"
              << "  -O <file>  Write output to file instead of stdout\n"
              << "  -h         Show this help\n";
}

int main(int argc, char* argv[]) {
    std::string config_dir, math_file, masters_file, output_format = "text";
    std::string kin_var, topology_name, output_file;
    bool use_epsilon = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h") { print_usage(argv[0]); return 0; }
        else if (arg == "-c" && i + 1 < argc) config_dir = argv[++i];
        else if (arg == "-m" && i + 1 < argc) math_file = argv[++i];
        else if (arg == "-f" && i + 1 < argc) masters_file = argv[++i];
        else if (arg == "-o" && i + 1 < argc) output_format = argv[++i];
        else if (arg == "-v" && i + 1 < argc) kin_var = argv[++i];
        else if (arg == "-t" && i + 1 < argc) topology_name = argv[++i];
        else if (arg == "-O" && i + 1 < argc) output_file = argv[++i];
        else if (arg == "-e") use_epsilon = true;
        else { std::cerr << "Unknown option: " << arg << "\n"; print_usage(argv[0]); return 1; }
    }

    if (config_dir.empty() || math_file.empty() || masters_file.empty()) {
        std::cerr << "Error: -c, -m, -f are required\n";
        print_usage(argv[0]); return 1;
    }

    try {
        // 1. Read config
        auto config = kira2diffeq::read_config(config_dir);

        if (topology_name.empty()) {
            if (!config.topologies.empty()) topology_name = config.topologies[0].name;
            else throw std::runtime_error("No topologies found");
        }

        const kira2diffeq::Topology* topo = nullptr;
        for (const auto& t : config.topologies) {
            if (t.name == topology_name) { topo = &t; break; }
        }
        if (!topo) throw std::runtime_error("Topology '" + topology_name + "' not found");

        // 2. Build symbol table
        GiNaC::symtab sym_table;
        for (const auto& [name, sym] : config.kinematics.symbols) sym_table[name] = sym;
        sym_table["d"] = GiNaC::symbol("d");

        // 3. Parse math output
        auto math_output = kira2diffeq::parse_math_output(math_file, masters_file, sym_table);

        if (math_output.masters.find(topology_name) == math_output.masters.end()) {
            throw std::runtime_error("No masters for topology '" + topology_name + "'");
        }
        const auto& masters = math_output.masters[topology_name];

        std::cerr << "Topology: " << topology_name
                  << ", masters: " << masters.size()
                  << ", rules: " << math_output.reduction_table.size() << "\n";

        // 4. Determine kinematic variables
        std::vector<std::string> kin_vars;
        if (!kin_var.empty()) kin_vars.push_back(kin_var);
        else for (const auto& [name, dim] : config.kinematics.kinematic_invariants) kin_vars.push_back(name);

        // 5. Compute DE matrices for all variables
        std::map<std::string, std::vector<std::vector<GiNaC::ex>>> all_matrices;
        for (const auto& var : kin_vars) {
            std::cerr << "Computing d/d" << var << "...\n";
            auto matrix = kira2diffeq::assemble_de_matrix(
                topology_name, masters, var, *topo,
                config.kinematics, math_output.reduction_table,
                math_output.sym_table);
            all_matrices[var] = matrix;
        }

        // 6. Output
        std::string output;
        if (output_format == "all") {
            output = kira2diffeq::export_all_mathematica(
                all_matrices, topology_name, masters, use_epsilon);
        } else if (output_format == "math" || output_format == "mathematica") {
            for (const auto& [var, matrix] : all_matrices) {
                output += kira2diffeq::export_to_mathematica(
                    matrix, topology_name, var, masters, use_epsilon);
                output += "\n";
            }
        } else {
            for (const auto& [var, matrix] : all_matrices) {
                output += kira2diffeq::export_to_text(
                    matrix, topology_name, var, masters, use_epsilon);
                output += "\n";
            }
        }

        if (!output_file.empty()) {
            std::ofstream fout(output_file);
            if (!fout) throw std::runtime_error("Cannot write to " + output_file);
            fout << output;
            std::cerr << "Written to " << output_file << "\n";
        } else {
            std::cout << output;
        }

        std::cerr << "Done.\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
