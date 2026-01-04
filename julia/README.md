# MOPS.jl

Julia port of the Maple MOPS (Multivariate Orthogonal Polynomials) library.

## Installation

```julia
using Pkg
Pkg.add(url="path/to/MOPS/julia")
```

## Usage

```julia
using MOPS
using SymPy

# Validate partitions
parvalid([3, 2, 1])  # true
parvalid([2, 3, 1])  # false (not non-increasing)

# Shifted factorial (Pochhammer symbol)
SFact(3, 5)  # 5040
@vars n
SFact(n, 3)   # n*(n+1)*(n+2)
```

## Development

Run tests:
```julia
using Pkg
Pkg.test("MOPS")
```
