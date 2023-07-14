# Membrane-Simulation

INTRODUCTION:

Programs that perform numerical simulations of vesicles, either approximated as a plane or a sphere, coupled to the diffusion curvature dependent lipids. There are two versions of the program a Python program and a C program:
vesicle_sim.py
vesicle_sim.C 

REQUIREMENTS:

vesicle_sim.py is for use with Python3 and requires modules: NumPy, sys, SciPy (optimize and special).
vesicle_sim.C requires GSL libraries are required for random number generation and spherical harmonics.

INSTALLATION (C program):

autoreconf --install
configure
make

If the automake files are changed you may need to run: autoreconf

program is built in "./src/", to change build location update the automake files: configure.ac and Makefile.am

USAGE:

PYTHON: 
Python vesicle_sim.py PARAM_FILE
See example parameters file ("parametrs.py") for syntax in the parameter input file.

C:
./src/vesicle_sim PARAM_FILE
See example parameters file ("./src/parametrs.inp") for syntax in the parameter input file
