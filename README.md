# A Note on Numerical Fluxes Conserving Harten's Entropies for the Compressible Euler Equations

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO)

This repository contains information and code to reproduce the results presented in the article
```
@online{ranocha2022note,
  title={A Note on Numerical Fluxes Conserving {H}arten's Entropies for the
         Compressible {E}uler Equations},
  author={Ranocha, Hendrik},
  year={2022},
  month={01},
  eprint={TODO},
  eprinttype={arXiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```
@misc{ranocha2022noteRepro,
  title={Reproducibility repository for
         "A Note on Numerical Fluxes Conserving {H}arten's Entropies for the
         Compressible {E}uler Equations"},
  author={Ranocha, Hendrik},
  year={2022},
  howpublished={\url{https://github.com/ranocha/paper-2022-Euler_Harten_EC}},
  doi={TODO}
}
```


## Abstract

Entropy-conserving numerical fluxes are a cornerstone of modern high-order
entropy-dissipative discretizations of conservation laws. In addition to entropy
conservation, other structural properties mimicking the continuous level such as
pressure equilibrium and kinetic energy preservation are important. This note
proves that there are no numerical fluxes conserving (one of) Harten's entropies
for the compressible Euler equations that also preserve pressure equilibria and
have a density flux independent of the pressure. This is in contrast to fluxes
based on the physical entropy, where even kinetic energy preservation can be
achieved in addition.


## Numerical experiments

The numerical experiments presented in the paper use [Trixi.jl](https://github.com/trixi-framework/Trixi.jl),
providing adaptive high-order numerical simulations of hyperbolic PDEs in Julia.
To reproduce the numerical experiments, you need to install [Julia](https://julialang.org/).

The `code` directory contains Julia code and instructions to reproduce the numerical
experiments, including postprocessing.

The numerical experiments were carried out using Julia v1.7.1.

Moreover, the `code` directory contains a Mathematica notebook and its PDF
version verifying some calculations presented in the paper.


## Authors

* [Hendrik Ranocha](https://ranocha.de) (University of Münster, Germany)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
