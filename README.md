CIMPRESS TECH CHALLENGE SOLUTION README
ROB LORD 
oh.lord92 at gmail.com

1. SOLUTION MOTIVATION

The problem posed by the tech challenge is similar to the cutting stock problem, which is NP Complete.
As a result, a polynomial time solution would be worth a lot more than $10,000! Therefore, either some kind
of heuristic is required, or we parallelise the problem. After trying a number of heuristics with limited success, 
I opted for the latter by parallelising with OpenCL to allow powerful hardware like graphics cards to take a crack at
the problem.

2. OPENCL

OpenCL allows the solution to be heavily parallelised- whether on a graphics card, processor, supercomputer or FPGA.
It also means that the majority of the code is in an OpenCL kernel, with (in this case) python largely only calling
OpenCL library functions. As OpenCL libraries are available for most languages, this makes my solution very portable
to languages other than python (simply convert the API calls over). While OpenCl does pose a number of extra challenges
(e.g. no dynamic length arrays or lists, synchronisation between threads) I felt that the extra power that could be accessed
would be worth it.

3. USAGE
The algorithm is provided as a python class. Before it can be used to solve puzzles, init() must be called to read
and compile the OpenCL kernel. Sample python usage code:

solver = CLSolver()
solver.init("path/to/montecarlo.cl")
solution = solver.solve(puzzle,16384,35,128)

4. ALGORITHM

The algorithm I opted for takes inspiration from various optimisation techniques, with particular emphasis on Monte Carlo
and simulated annealing. The program operates in work groups of individual threads- each thread generates a random solution to
the problem, and the best solution for the group is taken as the starting point for the next round. Each individual thread then
randomly clears and refills part of the solution, and the best result for the group again serves as the starting point for the
next iteration. At the end of the execution, the group solutions are copied back from the graphics card, and the best overall
solution is selected.

Solutions are randomly generated as follows:
1. The maximum possible square size at each possible location is calculated and stored
2. Loop from the largest possible square overall to the smallest
3. Find all of the squares for the given size, and randomly select one
4. Either fill it with a square of that size, or randomly skip it (allowing it to be filled by smaller squares later)
5. Loop until all squares of the given size have been filled, and move on to the next smallest size
6. When we get to size 1, the puzzle will be full.

Solutions are randomly perturbed for the next iteration as follows:
1. loop through all squares in the solution
2. With a random probability, choose to ignore this square in the next iteration (keep the current solution in this area)

5. TUNING

The algorithm has quite a few hyperparameters that will need to be tuned depending on the hardware it is being run on.
The default parameters are the ones I used for the challenge with two AMD 280x graphics cards in crossfire. Your system
may require you to lower these parameters, or it may be able to handle higher ones (and so give faster/ better solutions).
These parameters also allow the user to tune the accuracy/time tradeoff. The parameters are currently tuned to give the best
solution without going over 10 seconds (on my system), but they could be tuned to give quick, relatively good solutions 
or spend a long time generating a very good solution.

The hyperparameters are as follows (note tuning should largely be directed towards the first 3 items unless you're running a
very short or very long solution):

-simulations: The number of random solutions that are generated at each stage
-iterations: the number of stages (the completely random first iteration + the subsequent refinement iterations)
-workGroupSize: the number of threads that the best solution is aggregated over (i.e. with this set to 128, at the end
		of the iteration all 128 threads in the group take the best solution over the whole group as the starting
		point for the next iteration).
-alpha (specified within the .cl file): Simulated annealing parameter, specifies how the temperature decays over time. The 
		greediness of the solution can be modified by increasing (less greedy) or decreasing (more greedy) the probability
		of skipping a large square and allowing it to be filled with small square. By default the solution is programmed
		to get more greedy over each iteration. alpha changes the rate that the greediness (temperature) increases at. 
-prob (within .cl in perturbPuzzle): specifies the probability of including a square from the best solution in the
		next round. Mentioned only for clarity as it should not need to be altered from 0.5.
