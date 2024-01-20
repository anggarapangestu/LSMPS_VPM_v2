# VPM_Boundary_Layer
Fluid particle base simulation for Boundary Layer on a flat plate.
The code is based on the previous work by Ical and Adhika.

## Program Scenario
Boundary layer is a steady flow phenomena as long the flat plate not moving (static).
This simulation using explicit approach in which simulating the boundary layer from it's formation.
Initially the domain is a fluid with uniform flow by property of ρ, μ, and Re.
The flat plate solid will be impulsively started in the domain then the boundary layer start to develop.
The wall boundary condition is simulated by Brinkman Penalization.

## Output
The output of this simulation is the steady state is 
- The Boundary Layer Formation
- Boundary Layer Velocity Profile

## Method Note
- Grid node for adaptation, neighbor search, redistribution is very robust. The grid node data grouping method also delivers a fast calculation. <...>
- The velocity calculation by FMM provide a fast calculation (relatively). We need to get the best tuned value for tree level and the taylor expansion order for efficient simulation yet keep the accuracy level. 
    - TREE LEVEL put higher to reduce the particle number on cell for a faster calculation, but it will shrink the cell size, does it still meet the FMM criteria<?>
    - A higher EXPANSION ORDER have higher accuracy but it puts a lot of computational burden on M2L (multupole to local) process, since it takes at most $27*(p^2)$ order of calculation per cell.
- The data file reader (for particle or geometry data) work very well<!>
- The heaviest computation is particle saving (since it consisted lot of conditional). Please refine the particle writting code further <...>


## Programming Note
This note is author note regarding the code creation, debugging process, and physical interpretation of the simulation result.
This note helps author to make target on this program.
- Particle distribution near body surface is very crucial. A regular single resolution near the body surface is prefered for more accurate simulation.
- A symmetric particle distribution have a great computational stability. It takes lot of iteration time to start flow separation. *We can induce a slightly shifted distribution to get a small imbalance to trigger early (in term of time) flow separation.
- LSMPS having a problem calculating at domain edge. It may introduce large error at edge. Make sure no property calculation at domain edge.
- The vortex at downstream flow lasts very long. Make sure it not move to edge (Make a larger domain, so that the vortex motion and evolution at downstream not meet the domain edge).
- Proceed to 3D (^_^)

## May be need to be done
- Check the effect of iterative calculation in brinkman penalization
- Find a better resolve for LSMPS calculation

## Must be done
- Find the best force calculation
- Check the LSMPS on surface calculation (why need smoothing?) -> Read reference on penalization (more spesific)
- Check the effect of not flushing the vorticity value if not active (check remeshing)
- Only redistribute the active region but returns not change in simulation time (Why I dont really know) (Let me check in the 2D simulation)

## Occured problem
- Random spot pattern of vorticity [CAUSE]-> Fail to do vorticity penalization (*spesific: chi list ID get a wrong index order), [SOLUTION]-> Forget to reset the vector after calculation (chi data is stuffed in the later iteration), just add code line to reset the chi vector container.
- Broken vorticity interpolation at boundary [CAUSE]-> vortex get near to the domain edge. [SOLUTION]-> temporary: enlarge the domain, permanent: do something with the LSMPS calculation.
- Broken FMM calculation (The velocity field not continue, making a block discontinuity pattern) [CAUSE]-> Wrong loop indexing after code refactor (the looped cell ID at M2L not concludes all cell) -> [SOLUTION] rearrange the index to follow the loop criteria
- Broken FMM calculation (The velocity field not continue, having some block with discontinuity value) [CAUSE]-> Double cell ID in list 1 and 3 for the cell outside the vorticity source region -> [SOLUTION] Set a flag to prevent duplicate
- [2] Broken FMM calculation (Quite big discontinuity) [CAUSE]-> Wandering of nearest cell list for internal list 1 and 3. After taking the parent recursively, the code let a bypass to get child thus the cell can be separated apart from the evaluated cell (not near anymore) [SOLUTION]-> Put additional treatment, cannot find child for a cell that is the parent from candidate neighbor cell.
- Anticipated wrong starting cell ID at each level -> Add additional code for treatment near between level boundary.
- A mozaic pattern comes again into the velocity calculation [CAUSE]-> The double comparison problem in comparing the distance between cell in list evaluation (list_1 through list_4). The tiny deliberately boundary displacement adjustment using built in value of __DBL_EPSILON__ is BEYOND the comparison (turns out no effect). [SOLUTION]-> go back to the original user defined value of 10^-10 for EPSILON_TOLERANCE.

NOTE FOR MIGRATION to 3D
> Check the vortex strength ->
> vorticity become 3D

Data Collection
Simulation Timer
280,000 par   : 2.5 sec / iter
780,000 par   : 4.0 sec / iter     

FMM contain of 3 Error
> Direct from Multipole         -> Factor of 2.12
> Direct from M2M               -> Factor of 2.12
> Direct after translation      -> Factor of 1.82
Which gives error in order 3 for p = 10 (illust. 2^10 = 1000)

Also the error is combination of A -> sum of all vortex strength in the box [or] (box size * average vorticity)
> But: P > log (A/error) / log (factor), for A = 1, error = 1e-3, then P > ~10
> Along the way that average vorticity is below than 5 (O(10^1)) (But may reach 100 at the body surface, have a finest refinement)
> Whilst the box size at level 5 with domain size of 20 --> 0.625 (order of O(10-1)), the finest may come to (0.02 at level 10) -> Result more less tension on FMM expansion order


<!> Check in next code
std::numeric_limits<double>::epsilon()  #[DONE]

Determining the size of LSMPS bound

## LSMPS TEST
> Test for any local position calculate the LSMPS for variation: 
[1] Constant radius (Re*s) with different size
[2] Constant Re with different size
[3] Constant size with different Re
Create a chart between error.

What I want to do:
1. Do a simulation enhancement by adaptation:
   - See the error toward particle reducement
2. Built 3D simulation result.
   - FMM, Penalization, LSMPS
3. Check penalization, what is the problem with the force <?>
4. Check the LSMPS, what is the best practice for this simulation
5. 