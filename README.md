# SoEwald2D

[![Build Status](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml?query=branch%3Amain)

# SoEwald2D

[![Build Status](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml?query=branch%3Amain)

`SoEwald2D.jl` is a package based on `Julia` for rapid calculation of electrostatic interaction for particles quasi-2D charged systems.
It combined the SOE method and the Ewald2D method so that the complexity of its pairwise summation in the k space is of $O(N)$, much better than the original one, which is of $O(N^2)$.

This package can work as a force field of the MD package `ExTinyMD.jl`.

## Getting Started

Install the package in Julia REPL simply by
```julia
pkg> add ExTinyMD, SoEwald2D
```
