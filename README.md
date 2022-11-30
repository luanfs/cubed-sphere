# cubed-sphere
Cubed-sphere grid tools for geophysical fluid dynamics modelling.
Fully Fortran written, with outputs using Python.

Luan Santos
(luan.santos@usp.br)

- Work under constant development!

1) Use the Makefile to compile (just type 'make')

2) Run using "./main". These will call the necessary routines 
    for grid generation or modelling. Edit main.f90 
    for your purpuses. 

3) Mesh parameters must be set in par/mesh.par or other par/*.par files
   Other parameters must be set in par/ directory

4) Choose the simulation to be run in mesh.par (1, 2...)

5) Output is written in data/
 
6) Use Python scripts from pyscripts/ to plot output

This code is based on  <a href=https://github.com/pedrospeixoto/iModel/>iModel</a>.

