# RuSseL3D
This is the repository of the RuSseL3D code.

The code is named after the British philosopher and mathematician, Bertrand Russell.

It is an open-source code, distributed under the terms of the accompanying license.

The authors of the code are:
Constantinos J. Revelas (cjrevelas@gmail.com)
Aristotelis P. Sgouros (arissgouros@gmail.com)
Apostolos T. Lakkas (tolis1981@gmail.com)

RuSseL3D is a code developed in Fortran which applies the Finite Element method to run three-
dimensional calculations on heterogeneous polymer systems, based on Self-Consistent Field Theory (SCFT).
At the moment, the code can address homopolymer melts in contact with solid surfaces and provide with
useful results regarding the thermodynamics and the structural properties of the system. The solid
surfaces can be either bare or grafted with polymer chains of the same chemical identity as the matrix
chains.

The RuSseL3D repository includes the following files and directories:

README        -> this file
LICENCE       -> an MIT license for the open-source usage of the code
LICENCE-FHASH -> a license allowing the usage of the fhash.f90 subroutine
Makefile      -> file to build the binary executable
.fortls       -> a file pairing neovim with the fortran language server
.gitignore    -> files to be ignored by git version control system
.git/         -> directory containing the settings of the git version control system
obj/          -> directory where all object files are stored
run/          -> directory where the compiled executable file is redirected
src/          -> directory containing all the source files of the code
tools/        -> directory containing all the pre- and post- processing files
