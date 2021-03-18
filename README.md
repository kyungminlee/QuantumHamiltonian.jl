# QuantumHamiltonian.jl

| **Documentation** | **Build Status** | **Code Coverage** |
|:-----------------:|:----------------:|:-----------------:|
| [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] | [![Build][githubaction-img]][githubaction-url] | [![Code Coverage][codecov-img]][codecov-url] |

`QuantumHamiltonian.jl` is a library for constructing quantum many-body Hamiltonians. It aims to provide
- convenient and efficient representation of a generic lattice Hamiltonian and wave function
- reduction of the Hilbert space dimension using symmetry

## Installation

To install, type the following in Julia's package Pkg REPL-mode:
```julia-repl
(v1.5) pkg> registry add https://github.com/kyungminlee/KyungminLeeRegistry.jl.git
(v1.5) pkg> add QuantumHamiltonian
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://kyungminlee.org/QuantumHamiltonian.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: http://kyungminlee.org/QuantumHamiltonian.jl/dev

[githubaction-img]: https://github.com/kyungminlee/QuantumHamiltonian.jl/workflows/Build/badge.svg
[githubaction-url]: https://github.com/kyungminlee/QuantumHamiltonian.jl/actions?query=workflow%3ABuild

[codecov-img]: https://codecov.io/gh/kyungminlee/QuantumHamiltonian.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/QuantumHamiltonian.jl
