# QuantumHamiltonian

## Overview

`QuantumHamiltonian.jl` is a tool for constructing quantum many-body Hamiltonians.
It uses Abelian quantum numbers as well as translation symmetry to reduce dimensions
of the Hilbert space and the corresponding matrix representation of the Hamiltonian.

A schematic for the structure of the package is the following:
```
                State
                  ↓
                Site
                  ↓
                HilbertSpace → HilbertSpaceSector   Operator
                  ↓              ↓                    ↓
                HilbertSpaceRepresentation        → OperatorRepresentation
                  ↓                                   ↓
SymmetryGroup → ReducedHilbertSpaceRepresentation → ReducedOperatorRepresentation
```

The `HilbertSpace`, `HilbertSpaceSector` and `Operator` implement the abstract
Hilbert spaces and operators, while the `...Representation`s implement the representations
of the Hilbert spaces as 𝐂ⁿ (or 𝐑ⁿ), and of operators as n×n matrices.

## Workflow

First you need to create a Hilbert space representation:
1. Define `State`s, and `Site`s
1. Define the `HilbertSpace`
1. If there are quantum numbers, use them to define `HilbertSpaceSector`
1. Define `HilbertSpaceRepresentation` and construct basis set
1. If there is space symmetry, translation or point or both, use that to define `ReducedHilbertSpaceRepresentation`

And then you can create operator representation using the Hilbert space representation from above:
1. Define `Operator`s
1. Create `OperatorRepresentation` or `ReducedOperatorRepresentation` using `HilbertSpaceRepresentation` or `ReducedHilbertSpaceRepresentation`
1. Depending on what is more efficient,

## Installation

`QuantumHamiltonian.jl` is not yet registered on the Julia package registry. You
can install it using its URL as
```julia
]add https://github.com/kyungminlee/QuantumHamiltonian.jl.git
```
Since, however, `QuantumHamiltonian.jl` depends on other packages including [`LatticeTools.jl`](https://github.com/kyungminlee/LatticeTools.jl), it is convenient to add a custom registry.
In Julia, type
```sh
]registry add https://github.com/kyungminlee/KyungminLeeRegistry.jl.git
```
After this, you can
```julia
]add QuantumHamiltonian
```
