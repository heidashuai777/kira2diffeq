# kira2diffeq Manual

## 1) Project goal

`kira2diffeq` is a local C++ utility that turns **Kira reduction outputs** into differential equations for master integrals.

It does not run Feynman-diagram generation or IBP reduction itself.
It expects:
- a topology/kinematics description,
- a Kira-reduction file in `kira2math` style,
- and a `masters.final` master list.

The output is the DE system in the form:

```text
d/dx MI_i = sum_j M_x[i][j] MI_j
```

where `x` is a kinematic invariant from your config and `M_x` is the matrix collected for that variable.

## 2) What changed by the recent fixes

Recent reliability updates in this tree:

- Added CTest registration for the packaged CLI workflow.
- Added a repo-local CLI smoke test with fixtures.
- Gated verbose differentiation debug output behind env var `KIRA2DIFFEQ_DEBUG`.

You can now verify the package with `ctest` instead of only compiling.

## 3) Build and test workflow

### 3.1 Requirements

- C++17 compiler
- CMake >= 3.16
- GiNaC + yaml-cpp installed and discoverable by pkg-config (`GINAC`, `YAMLCPP`)
- POSIX shell for running commands

### 3.2 Configure and build

```bash
cmake -S . -B build
cmake --build build -j4
```

### 3.3 Run tests

```bash
ctest --test-dir build --output-on-failure
```

Expected tests:
- `test_sunrise`
- `cli_smoke`

If you see dependency or path issues in `test_sunrise`, run:

```bash
ctest --test-dir build --output-on-failure -R test_sunrise
```

### 3.4 Typical output

You should see all tests as `Passed`.
`cli_smoke` writes output to:

- [build/tests/cli_smoke_de.txt](/home/Heidashuai/claude/kira2diffeq/build/tests/cli_smoke_de.txt)

## 4) Command line interface

Executable:

- [build/kira2diffeq-cli](/home/Heidashuai/claude/kira2diffeq/build/kira2diffeq-cli)

Usage:

```text
Usage: kira2diffeq-cli [options]
Options:
  -c <dir>   Config directory (integralfamilies.yaml, kinematics.yaml)
  -m <file>  kira2math .m output file
  -f <file>  masters.final file
  -t <name>  Topology name
  -v <var>   Kinematic variable (default: all)
  -o <fmt>   Output: text | math | all (default: text)
  -e         Use epsilon (d = 4 - 2*eps) in output
  -O <file>  Write output to file instead of stdout
  -h         Show this help
```

Minimal invocation:

```bash
./build/kira2diffeq-cli \
  -c /path/to/config \
  -m /path/to/reduction.m \
  -f /path/to/masters.final \
  -t topology_name \
  -v t \
  -O /tmp/de_t.txt
```

If `-v` is omitted, the tool iterates all invariants defined in `kinematics.yaml`.

### 4.1 Output formats

`-o text` (default):

- Human-readable equations with `MI_i` notation.

`-o math`:

- Mathematica syntax for direct import.

`-o all`:

- All supported variables in Mathematica format.

`-e`:

- Replaces `d` with `4 - 2*eps` in expression output.

### 4.2 Debug logging

By default, differentiation no longer emits verbose traces.
To enable traces for troubleshooting:

```bash
export KIRA2DIFFEQ_DEBUG=1
./build/kira2diffeq-cli ...
```

Any other values (`0`, `false`, `False`, empty) keep debug output off.

## 5) Configuration format in detail

The config directory passed by `-c` must contain:

- `kinematics.yaml`
- `integralfamilies.yaml`

Both files are parsed by:

- [include/kira2diffeq/config_reader.h](/home/Heidashuai/claude/kira2diffeq/include/kira2diffeq/config_reader.h)
- [src/config_reader.cpp](/home/Heidashuai/claude/kira2diffeq/src/config_reader.cpp)

### 5.1 `kinematics.yaml`

Required/typical keys:

- `incoming_momenta` and/or `outgoing_momenta`
- `kinematic_invariants`
- `scalarproduct_rules`
- optional `momentum_conservation`
- optional `symbol_to_replace_by_one`

Example:

```yaml
kinematics:
  incoming_momenta: [p1, p2]
  outgoing_momenta: [k2]
  kinematic_invariants:
    - [s, 2]
    - [p12, 2]
  scalarproduct_rules:
    - [ [p1, p1], 0 ]
    - [ [p2, p2], 0 ]
    - [ [k2, k2], 0 ]
    - [ [p1, p2], p12 ]
    - [ [p1, k2], s/2 ]
    - [ [p2, k2], (s-p12) ]
  momentum_conservation:
    - [ q1, "-p1-p2-k2" ]
  symbol_to_replace_by_one: s
```

The parser handles canonicalisation and derived SP substitutions, including from momentum conservation.

### 5.2 `integralfamilies.yaml`

Required keys:

- `name`
- `loop_momenta`
- `top_level_sectors`
- `propagators`

Standard propagator entry:

```yaml
- [ "l1+p1-k2", 0 ]
```

ISP/bilinear entry:

```yaml
- [bilinear: [ ["k1", "l1"], 0 ]]
```

Complete fixture used by smoke test:

- [tests/fixtures/cli_smoke/integralfamilies.yaml](/home/Heidashuai/claude/kira2diffeq/tests/fixtures/cli_smoke/integralfamilies.yaml)

## 6) Reduction inputs

### 6.1 `masters.final`

One master per line:

```text
topology[1,1,0,0,0]
topology[1,0,1,0,0]
```

This is parsed by:

- [src/math_parser.cpp](/home/Heidashuai/claude/kira2diffeq/src/math_parser.cpp)

### 6.2 `kira2math` output file

Expected shape:

```text
topo[1,1,2,0] -> +topo[1,0,2,0]*(1), +topo[1,1,1,0]*(2*t), 
```

Each rule maps a target integral to a signed sum of master integrals and coefficients.
The parser is whitespace-robust and tolerant of line breaks with trailing commas.

Smoke fixtures:

- [tests/fixtures/cli_smoke/reduction.m](/home/Heidashuai/claude/kira2diffeq/tests/fixtures/cli_smoke/reduction.m)
- [tests/fixtures/cli_smoke/masters.final](/home/Heidashuai/claude/kira2diffeq/tests/fixtures/cli_smoke/masters.final)

## 7) Internal pipeline (recommended mental model)

1. Parse config into `KiraConfig` containing:
   - topologies: loop momenta, propagators, sectors, bilinear metadata
   - kinematics: SP rules and invariants
2. Parse reduction tables and masters into `MathOutput`.
3. For each requested variable `x`:
   - differentiate each master via `differentiate_master`
   - reduce every shifted integral in each derivative with `reduce_to_masters`
   - build dense matrix `M_x`
4. Export matrix in requested format.

Relevant source files:

- [src/reducer.cpp](/home/Heidashuai/claude/kira2diffeq/src/reducer.cpp)
- [src/differentiator.cpp](/home/Heidashuai/claude/kira2diffeq/src/differentiator.cpp)
- [src/exporter.cpp](/home/Heidashuai/claude/kira2diffeq/src/exporter.cpp)

## 8) Data structures at a glance

- `Kinematics` (source of invariants and SP rules)
- `Topology` (propagator definitions and loop momenta)
- `IntegralSum` (sum of shifted integrals + coefficients)
- `ReductionTable` (hash map from integral key to lazy coefficient strings)
- `MathOutput` (`reduction_table`, `masters`, and symbol table)

These are declared in:

- [include/kira2diffeq/config_reader.h](/home/Heidashuai/claude/kira2diffeq/include/kira2diffeq/config_reader.h)
- [include/kira2diffeq/math_parser.h](/home/Heidashuai/claude/kira2diffeq/include/kira2diffeq/math_parser.h)
- [include/kira2diffeq/reducer.h](/home/Heidashuai/claude/kira2diffeq/include/kira2diffeq/reducer.h)

## 9) Example runs

### 9.1 Run bundled smoke fixture

```bash
cmake --build build -j4
./build/kira2diffeq-cli \
  -c tests/fixtures/cli_smoke \
  -m tests/fixtures/cli_smoke/reduction.m \
  -f tests/fixtures/cli_smoke/masters.final \
  -t cli_smoke \
  -v t \
  -o text \
  -O /tmp/cli_smoke_de.txt
```

Then verify the file has `dMI_0/dt = ...` lines.

### 9.2 Full integration in scripts

Use a loop to generate all-variable DEs:

```bash
./build/kira2diffeq-cli \
  -c path/to/config \
  -m path/to/reduction.m \
  -f path/to/masters.final \
  -t topo_name \
  -o all \
  -e \
  -O /tmp/topo_all_de.m
```

## 10) CLI test strategy in this repository

Current CTest suite in:

- [tests/CMakeLists.txt](/home/Heidashuai/claude/kira2diffeq/tests/CMakeLists.txt)

Included:

- `test_sunrise`: exercises parsing, differentiation, reduction, matrix assembly.
- `cli_smoke`: repo-local end-to-end CLI execution using local fixtures.

## 11) Error diagnostics and warnings

### Warning pattern

- `Warning: integral X[... ] not found in reduction table`

Interpretation:
- derivative produced a shifted integral not reducible from your table and not listed as master.

Actions:
- Extend `kira2math` file to include that reduction rule.
- Check whether the topology/invariant combination in config matches the table.
- Confirm propagator strings/indices match exactly between topology and reduction outputs.

### Parse failures

- If config fails to parse, check YAML syntax and scalarproduct expression tokenization.
- If parsing `.m` fails, ensure lines follow `lhs -> +rhs*(coeff)` pattern with commas between terms.
- If topologies mismatch, confirm `-t` exists and exact spelling matches.

## 12) Collinear-workflow alignment notes (for the project context)

This project has explicit physics constraints for collinear two-loop processes.
Keep that workflow external to this package but ensure your files are consistent:

- use regulated substitution `p1.p2 = p12` in intermediate algebra,
- apply final `p12 -> 0` only after cancellation checks,
- verify kinematic identities like

```text
(p1 + p2 + k2)^2 == MH^2
```

and Ward checks in your upstream symbolic layer.

Do not silently apply final constraints before reduction simplification.

## 13) Performance notes

- This implementation currently applies reductions term-by-term with hash lookup; large topologies can produce large term lists.
- For very large reduction tables, the dominant cost is symbolic coefficient algebra.
- Use targeted variable `-v` to avoid building all matrices when only one variable is needed.

## 14) Developer extension points

To add new behavior:

- Parser support:
  - `parse_math_output` / `parse_integral` in [src/math_parser.cpp](/home/Heidashuai/claude/kira2diffeq/src/math_parser.cpp)
- New topology syntax or decomposition logic:
  - [src/config_reader.cpp](/home/Heidashuai/claude/kira2diffeq/src/config_reader.cpp)
- New derivative strategies:
  - [src/differentiator.cpp](/home/Heidashuai/claude/kira2diffeq/src/differentiator.cpp)
- New output layouts:
  - [src/exporter.cpp](/home/Heidashuai/claude/kira2diffeq/src/exporter.cpp)

## 15) Quick project map

- [CMakeLists.txt](/home/Heidashuai/claude/kira2diffeq/CMakeLists.txt)
- [src/main.cpp](/home/Heidashuai/claude/kira2diffeq/src/main.cpp)
- [src/config_reader.cpp](/home/Heidashuai/claude/kira2diffeq/src/config_reader.cpp)
- [src/math_parser.cpp](/home/Heidashuai/claude/kira2diffeq/src/math_parser.cpp)
- [src/differentiator.cpp](/home/Heidashuai/claude/kira2diffeq/src/differentiator.cpp)
- [src/reducer.cpp](/home/Heidashuai/claude/kira2diffeq/src/reducer.cpp)
- [src/exporter.cpp](/home/Heidashuai/claude/kira2diffeq/src/exporter.cpp)
- [tests/CMakeLists.txt](/home/Heidashuai/claude/kira2diffeq/tests/CMakeLists.txt)

## 16) PDF generation

A PDF is already checked in as:
- [MANUAL.pdf](/home/Heidashuai/claude/kira2diffeq/MANUAL.pdf)

To regenerate it:

```bash
cd /home/Heidashuai/claude/kira2diffeq
pdflatex -interaction=nonstopmode MANUAL.tex
```

This repository intentionally avoids external conversion tooling to keep the workflow minimal.
`MANUAL.tex` is the source used for this PDF.
