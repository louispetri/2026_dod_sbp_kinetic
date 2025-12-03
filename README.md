# Domain-of-dependence-stabilized cut-cell discretizations of linear kinetic models with summation-by-parts properties

# ToDo
- Ist Titel final?
- Korrektheit der finalen Shell/Julia commands überprüfen
- DOI
- References
- Abstract

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](ToDo)

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
@misc{ToDo
}
```

## Abstract

ToDo


## Numerical experiments

In order to generate the results from this repository, you need to install [Julia](https://julialang.org).
We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org).


The numerical results have been generated using Julia version 1.10.2, and we recommend installing the same.
Once you have installed Julia, you can clone this repository, enter this directory and start the executable
`julia` with the following steps

```shell
git clone https://github.com/louispetri/2025_dod_sbp_kinetic.git
cd 2025_dod_sbp_kinetic
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

- Gunnar Birke (gunnar.birke@uni-muenster.de)
- [Christian Engwer](https://www.uni-muenster.de/FB10srvi/persdb/MM-member.php?id=724) (christian.engwer@uni-muenster.de)
- [Sigrun Ortleb](https://acom.rwth-aachen.de/the-lab/submenu/team-people/name:sigrun_ortleb) (sigrun.ortleb@rwth-aachen.de)
- Louis Petri (lpetri01@uni-mainz.de)
- [Hendrik Ranocha](https://ranocha.de) (hendrik.ranocha@uni-mainz.de)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
