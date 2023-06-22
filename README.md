# Multiple Stellar Evolution (MSE) -- A Population Synthesis Code for Multiple-Star Systems #

MSE is code that models the long-term evolution of hierarchical multiple-star systems (binaries, triples, quadruples, and higher-order systems) from the main sequence until remnant phase. It takes into account gravitational dynamical evolution, stellar evolution (using the `sse` tracks), and binary interactions (such as mass transfer and common-envelope evolution).  It includes routines for external perturbations from flybys in the field, or (to limited extent) encounters in dense stellar systems such as galactic nuclei. 

C++ and Fortran compilers are required, as well as Python (2/3) for the Python interface. Make sure to first compile the code using `make`. Please modify the Makefile according to your installation (`CC`, `CXX`, and `FC` should be correctly assigned).  

The script `test_mse.py` can be used to test the installation. The script `run_system.py` is useful for quickly running a system. 

**See the user guide (doc/doc.pdf) for more detailed information.**

#Updates to stellar evolution routines src/sse/evolv1.f

1. mlwind.f has been replaced with the Cosmic version.
2. hrdiag.f has been replaced with the Cosmic version, with some minor adjustments: (i) remnantflag is replaced with nsflag (ii) code related to bh spin and pisn tracking is commented out
3. evolv1.f: tweaks to check with radial conversion.


