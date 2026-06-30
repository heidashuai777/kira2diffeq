#include "kira2diffeq/math_parser.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace kira2diffeq {

// ---- LazyReductionRule ----

std::vector<std::pair<std::vector<int>, GiNaC::ex>> LazyReductionRule::parsed_rhs(
    const GiNaC::symtab& table) const
{
    std::vector<std::pair<std::vector<int>, GiNaC::ex>> result;
    GiNaC::parser reader(table);
    for (const auto& [mp, coeff_str] : rhs_strings) {
        result.push_back({mp, reader(coeff_str)});
    }
    return result;
}

// ---- Parsing utilities ----

std::pair<std::string, std::vector<int>> parse_integral(const std::string& s) {
    auto bracket_pos = s.find('[');
    if (bracket_pos == std::string::npos) {
        throw std::runtime_error("Invalid integral format: " + s);
    }
    std::string name = s.substr(0, bracket_pos);
    auto end_pos = s.find(']', bracket_pos);
    if (end_pos == std::string::npos) {
        throw std::runtime_error("Invalid integral format: " + s);
    }
    std::string powers_str = s.substr(bracket_pos + 1, end_pos - bracket_pos - 1);
    std::vector<int> powers;
    std::istringstream iss(powers_str);
    std::string token;
    while (std::getline(iss, token, ',')) {
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        powers.push_back(std::stoi(token));
    }
    return {name, powers};
}

static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

static size_t find_matching_paren(const std::string& s, size_t open_pos) {
    int depth = 0;
    for (size_t i = open_pos; i < s.size(); i++) {
        if (s[i] == '(') depth++;
        else if (s[i] == ')') {
            depth--;
            if (depth == 0) return i;
        }
    }
    return std::string::npos;
}

static bool extract_term(const std::string& s, std::string& master_str,
                          std::string& coeff_str) {
    auto star_paren = s.find("*(");
    if (star_paren == std::string::npos) return false;
    master_str = trim(s.substr(0, star_paren));
    auto close_paren = find_matching_paren(s, star_paren + 1);
    if (close_paren == std::string::npos) return false;
    coeff_str = s.substr(star_paren + 2, close_paren - star_paren - 2);
    return true;
}

// ---- Main parser ----

MathOutput parse_math_output(const std::string& math_file,
                              const std::string& masters_file,
                              const GiNaC::symtab& sym_table) {
    MathOutput output;
    output.sym_table = sym_table;
    if (output.sym_table.find("d") == output.sym_table.end()) {
        output.sym_table["d"] = GiNaC::symbol("d");
    }

    std::ifstream fin(math_file);
    if (!fin.is_open()) {
        throw std::runtime_error("Cannot open math file: " + math_file);
    }

    std::string line;
    enum State { WAIT_LHS, WAIT_RHS };
    State state = WAIT_LHS;
    LazyReductionRule current_rule;

    while (std::getline(fin, line)) {
        std::string trimmed = trim(line);
        if (trimmed.empty()) continue;
        if (trimmed == "{" || trimmed == "}") continue;

        bool has_comma = false;
        if (trimmed.back() == ',') {
            trimmed = trimmed.substr(0, trimmed.size() - 1);
            has_comma = true;
        }
        trimmed = trim(trimmed);
        if (trimmed.empty()) {
            if (!current_rule.topology.empty()) {
                IntegralKey key = {current_rule.topology, current_rule.target_powers};
                output.reduction_table[key] = current_rule;
                current_rule = LazyReductionRule{};
                state = WAIT_LHS;
            }
            continue;
        }

        if (state == WAIT_LHS) {
            auto arrow_pos = trimmed.find("->");
            if (arrow_pos == std::string::npos) {
                throw std::runtime_error("Expected LHS with '->', got: " + trimmed);
            }

            std::string lhs_str = trim(trimmed.substr(0, arrow_pos));
            auto [topo_name, powers] = parse_integral(lhs_str);
            current_rule.topology = topo_name;
            current_rule.target_powers = powers;
            current_rule.rhs_strings.clear();

            std::string rhs_part = trim(trimmed.substr(arrow_pos + 2));

            if (rhs_part.empty()) {
                state = WAIT_RHS;
            } else if (rhs_part == "0") {
                state = WAIT_LHS;
                if (has_comma) {
                    IntegralKey key = {current_rule.topology, current_rule.target_powers};
                    output.reduction_table[key] = current_rule;
                    current_rule = LazyReductionRule{};
                }
            } else if (rhs_part[0] == '+') {
                std::string term = trim(rhs_part.substr(1));
                if (term != "0") {
                    std::string master_str, coeff_str;
                    if (extract_term(term, master_str, coeff_str)) {
                        auto [m_topo, m_powers] = parse_integral(master_str);
                        current_rule.rhs_strings.push_back({m_powers, coeff_str});
                    }
                }
                if (has_comma) {
                    IntegralKey key = {current_rule.topology, current_rule.target_powers};
                    output.reduction_table[key] = current_rule;
                    current_rule = LazyReductionRule{};
                    state = WAIT_LHS;
                } else {
                    state = WAIT_RHS;
                }
            }
        } else if (state == WAIT_RHS) {
            if (trimmed[0] != '+') {
                throw std::runtime_error("Expected RHS term starting with '+', got: " + trimmed);
            }
            std::string term = trim(trimmed.substr(1));
            if (term != "0") {
                std::string master_str, coeff_str;
                if (extract_term(term, master_str, coeff_str)) {
                    auto [m_topo, m_powers] = parse_integral(master_str);
                    current_rule.rhs_strings.push_back({m_powers, coeff_str});
                }
            }

            if (has_comma) {
                IntegralKey key = {current_rule.topology, current_rule.target_powers};
                output.reduction_table[key] = current_rule;
                current_rule = LazyReductionRule{};
                state = WAIT_LHS;
            }
        }
    }

    // Last rule (no trailing comma)
    if (!current_rule.topology.empty()) {
        IntegralKey key = {current_rule.topology, current_rule.target_powers};
        if (output.reduction_table.find(key) == output.reduction_table.end()) {
            output.reduction_table[key] = current_rule;
        }
    }

    // Parse masters.final
    std::ifstream mfin(masters_file);
    if (mfin.is_open()) {
        std::string mline;
        while (std::getline(mfin, mline)) {
            mline = trim(mline);
            if (mline.empty() || mline[0] == '#') continue;
            auto comment_pos = mline.find('#');
            if (comment_pos != std::string::npos) {
                mline = trim(mline.substr(0, comment_pos));
            }
            if (mline.empty()) continue;
            auto [topo, powers] = parse_integral(mline);
            output.masters[topo].push_back(powers);
        }
    }

    return output;
}

} // namespace kira2diffeq
