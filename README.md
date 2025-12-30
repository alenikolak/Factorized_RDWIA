# Factorized RDWIA calculations for nucleon knockout

This is a small example code accompanying [arxiv:xxxx.xxxx](https://arxiv.org/blabla) that puts into practice the factorization of the nucleon knockout cross section in terms of overlap matrices.

The overlap matrices computed using the ED-RMF potential of [arxiv:1904.10696](https://arxiv.org/abs/1904.10696) can be obtained: [here](https://drive.google.com/drive/folders/1wFL1jxN16xFpuRadd8U2he3_be9ANQ5a?usp=drive_link).



## Ingredients for nucleon knockout calculation

The nucleon knockout cross section is factorized into three pieces

1. The overlap matrices which are stored in  `struct S_level` and dealt with in `S_tables.cpp`

2. Couplings to the nucleon stored in `struct FormFactors` and dealt with in `Formfactors.cpp`.

3. The lepton tensor stored in `struct LeptonTensor` and dealt with in `Leptonic.cpp`

The overlap matrices and nucleon form factors are combined to yield the hadron tensor as implemented in `Hadronic.cpp`.
The hadron and lepton tensor are contracted to yield the nucleon knockout cross sections as implemented in `Crosssections.cpp`.

The function `double cross_section(double E_i, double E_f, double costheta_f, double pm, double phi_N, LeptonTensor *Lep, FormFactors *FF, S_level *S)` serves mostly as an illustration. In practice one makes use of the fact that the dependence on azimuth angle factorizes to only compute the relevant terms.

## Running the code

The code can simply be compiled using the `compile.sh` script. 

Running `bash compile.sh` produces two executables:
1. `demo` which reads overlap matrices for the 1s shell in carbon and produces some cross sections.

2. `output_benchmarks` which prints certain response functions to standard output, used for confirming that everything is correct.

These executables attempt to read overlap matrices from the local path `./Overlaps/EDRMF/C12/`.
To save space only the overlaps required for `demo` are included here, the others can be obtained [here](https://drive.google.com/drive/folders/1wFL1jxN16xFpuRadd8U2he3_be9ANQ5a?usp=drive_link).

## The overlap matrices
Overlap matrices for carbon and oxygen are provided. 
These are provided for neutral current reactions (p -> p or n -> n), and charged current interaction (p -> n or n -> p). The `readme.txt` file in `./Overlaps/EDRMF` explains the designation of levels and final-state nucleon.  
The maximal nucleon kinetic energy included in the tables is currently 260 MeV. 
These can readily be extended to higher energy if necessary.

Overlap matrices obtained with other models can be computed upon reasonable request to the authors.

## References

If you use this code or modifications thereof please cite

```
@article{Nikolakopoulos:2025bla,
    author = {Nikolakopoulos, Alexis and Gonzalez-Jimenez, Raul},
    title = {Factorized distorted wave calculations for electron neutrino and BSM processes},
    eprint = "xxxx.xxxx",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    reportNumber = "NT@UW-25-20",
    year = "2025"
}
```

If you use the overlap matrices obtained with the ED-RMF potential cite:

```
@article{Gonzalez-Jimenez:2019qhq,
    author = "Gonz{\'a}lez-Jim{\'e}nez, R. and Nikolakopoulos, A. and Jachowicz, N. and Ud{\'\i}as, J. M.",
    title = "{Nuclear effects in electron-nucleus and neutrino-nucleus scattering within a relativistic quantum mechanical framework}",
    eprint = "1904.10696",
    archivePrefix = "arXiv",
    primaryClass = "nucl-th",
    doi = "10.1103/PhysRevC.100.045501",
    journal = "Phys. Rev. C",
    volume = "100",
    number = "4",
    pages = "045501",
    year = "2019"
```

and 

```
@article{Sharma:1993it,
    author = "Sharma, M. M. and Nagarajan, M. A. and Ring, P.",
    title = "{rho meson coupling in the relativistic mean field theory and description of exotic nuclei}",
    reportNumber = "DL-NUC-P341T",
    doi = "10.1016/0370-2693(93)90970-S",
    journal = "Phys. Lett. B",
    volume = "312",
    pages = "377--381",
    year = "1993"
}
```

