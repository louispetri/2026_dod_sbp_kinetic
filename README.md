# Domain-of-dependence-stabilized cut-cell discretizations of linear kinetic models with summation-by-parts properties

# ToDo
- References

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18196533.svg)](https://doi.org/10.5281/zenodo.18196533)

This repository contains information and code to reproduce the results
presented in the article
```bibtex
@online{ToDo
}
```

If you find these results useful, please cite the article mentioned above.
If you use the implementations provided here, please **also** cite this
repository as
```bibtex
@misc{petri2026kineticRepro,
  title={Reproducibility repository for
         "{D}omain-of-dependence-stabilized cut-cell discretizations
         of linear kinetic models with summation-by-parts properties"},
  author={Petri, Louis and Ortleb, Sigrun and Birke, Gunnar and
          Engwer, Christian and Ranocha, Hendrik},
  year={2026},
  howpublished={\url{https://github.com/louispetri/2026_dod_sbp_kinetic}},
  doi={10.5281/zenodo.18196533}
}
```

## Abstract

We employ the summation-by-parts (SBP) framework to extend
the recent domain-of-dependence (DoD) stabilization for cut cells to linear
kinetic models in diffusion scaling. Numerical methods for these models are
challenged by increased stiffness for small scaling parameters and the
necessity of asymptotics preservation regarding a parabolic limit equation.
As a prototype model, we consider the telegraph equation in one spatial
dimension subject to periodic boundary conditions with an asymptotic limit
given by the linear heat equation.
We provide a general semidiscrete stability result for this model when
spatially discretized by arbitrary periodic (upwind) SBP operators and
formally prove that the fully discrete scheme is asymptotic preserving.
Moreover, we prove that DoD with central numerical
fluxes leads to periodic SBP operators. Furthermore, we show that adapting
the upwind DoD scheme yields periodic upwind SBP operators. Consequently,
DoD stabilization possess the desired properties considered in the first part
of this work and thus leads to a stable and asymptotic preserving scheme for
the telegraph equation. We back our theoretical results with numerical
simulations and demonstrate the applicability of this cut-cell
stabilization for implicit time integration in the heat equation limit.


## Numerical experiments

In order to generate the results from this repository, you need to install [Julia](https://julialang.org).
We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org).


The numerical results have been generated using Julia version 1.10.2, and we recommend installing the same.
Once you have installed Julia, you can clone this repository, enter this directory and start the executable
`julia` with the following steps

```shell
git clone https://github.com/louispetri/2026_dod_sbp_kinetic.git
cd 2026_dod_sbp_kinetic
cd Julia_code
julia --project=.
```
Then enter the following commands to generate the data for the one-dimensional case and all the plots of the paper.

```julia
julia> import Pkg; Pkg.instantiate() # Does not need to be re-run the next time you enter the REPL
include("./run_all_scripts.jl")
```

  

All the figures are now ready and available in the respective folders/subfolders.

## Authors

- Gunnar Birke
- [Christian Engwer](https://www.uni-muenster.de/FB10srvi/persdb/MM-member.php?id=724)
- [Sigrun Ortleb](https://acom.rwth-aachen.de/the-lab/submenu/team-people/name:sigrun_ortleb)
- Louis Petri
- [Hendrik Ranocha](https://ranocha.de)

## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
