This directory contains python scripts used to solve Landau's dielectric funtion,
eqn 8.4.8 in Krall and Trivelpiece's Principles of Plasma Physics.

In `dispsolve.py`, a function, `D_h5data()`, is provided to compute this function for 
arbitrary distributions functions defined by an array of values that one might obtain from
an Osiris simulation, and which are stored inside an H5Data object from the library
PyVisOS, which can be found here https://github.com/UCLA-Plasma-Simulation-Group/pyVisOS.git.

It also contains a function, `get_D()`, which is used to compute the value of the
dielectric in a window of complex omega for given k0 and distribution function, which was
practical for the applications I was interested in.

See the `examples/` to get a sense for the usage.
