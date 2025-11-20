# Code for A finite element implementation of a large deformation gradient-damage theory for fracture with Abaqus user material subroutines
This repository contains the input files, along with the user material subroutines that were used in the manuscript:

>**Alkhoury, K.**, Chester, S. A., & Srivastava, V. (2026). A finite element implementation of a large deformation gradient-damage theory for fracture with Abaqus user material subroutines. Engineering Fracture Mechanics, 331, 111677. https://doi.org/10.1016/j.engfracmech.2025.111677
>
If you use this code, or build upon any part of it, **please cite the publication above**.

Specifically, the code provided here is used to produce the results in Section 4: Applications to fracture problems.

1. Penny-shaped specimen in tension: large deformation rubber
2. Single-edge notched tension (SENT) test: small deformation linear elasticity
3. Single-edge notched beam (SENB) test: large deformation rate-dependent plasticity

## Running the codes
❗**Disclaimer**❗  
> These example codes require a supported Fortran and C/C++ compiler setup to run user subroutines in Abaqus.  
> Please ensure that you have installed the appropriate compilers compatible with your Abaqus version  
> (e.g., Intel OneAPI Classic Fortran and Microsoft Visual Studio on Windows, or Intel/GNU compilers on Linux).  
> Without the correct compiler configuration, the subroutines will not compile or run.

1. Copy the Fortran subroutine (.for file) and the Abaqus input file (.inp) into the same directory, e.g., Test_Directory.
2. Open the Abaqus command window.
3. Change the working directory to Test_Directory:
    ```
    cd path/to/Test_Directory
    ```
5. Run Abaqus with the input file and link it to the Fortran subroutine:  
    ```
   abaqus job=JobName inp=InputName.inp user=SubroutineName.for
    ```
- To run example 4.1: Penny-shaped specimen in tension: large deformation rubber, set  
`inp = Penny-Shaped_input_4_1.inp` and `user = Penny-Shaped_UMAT_4_1.for`

- To run example 4.2:  Single-edge notched tension (SENT) test: small deformation linear elasticity, set  
`inp = SENT_input_4_2.inp` and `user = SENT_UMAT_4_2.for`

- To run example 4.3:  Single-edge notched beam (SENB) test: large deformation rate-dependent plasticity, set  
`inp = SENB_input_4_3.inp` and `user = SENB_UMAT_4_3.for`

>For any questions regarding the code or the implementation, please contact **keven_alkhoury@brown.edu**.

