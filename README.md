# Solving Weighted Maximum Satisfiability with Branch and Bound and Clause Learning
## Project description

This project contains the software involved in the research paper *Solving Weighted Maximum Satisfiability with Branch and Bound and Clause Learning*. The content is the following:
- **WMaxCDCL**: This folder contains the code and binary files of the WMaxCDCL solvers that were submitted to MaxSAT Evaluation 2023. The content is exactly the one that is obtained if the solver is downloaded from the (evaluation's website)[https://maxsat-evaluations.github.io/2023/descriptions.html].
- **WMaxCDCL-flags**: This folder contains the above source code of WMaxCDCL with slight modifications introduced to enable and disable some of WMaxCDCL's functionalities for experimentation purposes. Details are given below.
- **LS_UB.txt**: This file contains the upper bounds of the objective function used in the experimental section of the paper.

## WMaxCDCL-flags
In order to compile the solver, being placed at the folder `WMaxCDCL-flags/code/simp`, run:

```sh
make rs
```
In order to obtain the different versions of the solver defined in the paper, the solver bust be compiled with:
  ```sh
# wmcdcl\LA
make NO_LA=1 rs

# wmcdcl\harden
make NO_HARDEN=1 rs

# wmcdcl\initUB
make UPDOWN=1 rs

# wmcdcl+alwaysPB
make ALWAYS_ADDPB=1 rs

# wmcdcl+loLA
make IOLA=1 rs

# wmcdcl\PB
make NO_ADDPB=1 rs

# wmcdcl+alwaysLA
make ALWAYS_LA=1 rs

# wmcdcl\simp
make NO_SIMPSC=1 rs

# wmcdcl and wmcdcl\LS
make rs
#The binary file is the same, but UB must not be provided to wmcdcl\LS
```

## Running the solver
The solver is simply run as:
```sh
./wmaxcdcl [instancefilepath] [UB]
```
where the UB parameter is optional. 
