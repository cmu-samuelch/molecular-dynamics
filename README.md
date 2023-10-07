# MOLECULAR DYNAMICS
Author: "Samuel Chen"


# GLOBAL ASSUMPTIONS
- ε, σ, and m are all 1.
hopefully, this will hold true and I won't have to go back and insert these
everywhere.

# TODOs
PS #3:
- DONE: randomly initialize particle velocities with zero total momentum
- DONE: implement continuous force/energy with cutoff of 2.5 (dimless)
- DONE: calculates instantaneous temperature, pressure
    - done?: applies periodic boundary conditions and the nearest-image convention
    - done?: create side length as a variable set in the code

PS #4:
- ???

LONG-TERM:
- IN PROGRESS: improve variable names


# NOTES

- use N-by-3 arrays to store positions, velocities, forces
- please do not create N-by-N arrays (don't store all atomic separations/pair forces)

- please comment code :)

- background on mol sim file formats - will be useful later on
- for simple applications: XYZ is very simple to work with
- openBabel, avogadro can convert trajfile formats