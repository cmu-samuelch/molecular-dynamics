# MOLECULAR DYNAMICS
Author: "Samuel Chen"


# USAGE
- change only the values in `main` in {mc,md}.jl and in `setup` in Parameters.jl
to run different simulations.
    - to disable PBCs, set `L = 0`. this will ensure the code never alters 
    particle positions.
    - to disable cutoff radii, set `cut📏` to -1.
- run all of `VelocityVerlet.jl` and all of `ReadWrite.jl` before running `main.jl`;
this ensures that the functions from those two files are loaded in
- 

# GLOBAL ASSUMPTIONS
- ε, σ, and m are all 1.
hopefully, this will hold true and I won't have to go back and insert these
everywhere.

# TODOs
Psets: DONE
Final Project:
- tbd

LONG-TERM:
- IN PROGRESS: improve variable names
- IN PROGRESS: improve the style and interface, make the interface more opaque
and easier to continue running from the top level without having to change my
interfaces each problem set
    - DONE: implement a struct to pass along simulation parameters so i dont have to
    do like five runs just to debug my interface