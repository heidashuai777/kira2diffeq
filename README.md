# kira2diffeq

`kira2diffeq` is a tool to convert results from the Kira integral reduction program into systems of differential equations. It provides a C++ library and command-line interface (CLI) to read kinematic/topology configuration files, differentiate master integrals, export differential equation systems, and reduce integrals using GiNaC for symbolic computations and YAML-CPP for configuration parsing.

## Features

- **Configuration parser**: Reads YAML files describing kinematics and topologies for multi-loop Feynman integrals.
- **Differentiator**: Builds differential equations for master integrals and uses GiNaC to manipulate symbolic algebra.
- **Reducer and exporter**: Supports algebraic reduction of integrals and exports results in formats compatible with the `Kira` program and GiNaC.
- **Command-line interface**: Provides a CLI to run conversions from Kira outputs to differential equation systems.

## Build

This package uses CMake. Ensure you have GiNaC and YAML-CPP installed. To build the project and tests:

```bash
cmake -B build
cmake --build build
ctest --test-dir build
```

After building, you can install the library and CLI with:

```bash
cmake --install build
```

## Usage

Refer to the manual (`MANUAL.md`) for detailed usage instructions and examples. Typical usage involves preparing a configuration YAML file with kinematic definitions and running the CLI to produce a system of differential equations.

## License

This project is provided as-is; see the repository for license information.
