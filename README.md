# 4Quadrants
4 Quadrants of Dynamic Optimization

1: **Discrete-time**/**deterministic** solution using [SolveDSGE](https://github.com/RJDennis/SolveDSGE.jl).jl (@RJDennis). [Code link](https://github.com/azev77/4Quadrants/blob/main/DiscreteTime_Deterministic_SolveDSGE.jl). 
The corresponding simulations (using 6 different methods):

![image](https://user-images.githubusercontent.com/7883904/146826183-c2b4ddbf-eba4-4f06-a253-89a1fd4c9951.png)

2: **Discrete-time**/**deterministic** solution using [POMDPs](https://github.com/JuliaPOMDP/POMDPs.jl).jl (@zsunberg) 
-[Code link](https://github.com/azev77/4Quadrants/blob/main/DiscreteTime_Deterministic_POMDPs.jl). Note: this example is currently not working for me. See [discussion here](https://github.com/JuliaPOMDP/POMDPs.jl/discussions/351). 
-Please post correct code below. 

3: **Discrete-time**/**deterministic** solution using [JuMP](https://github.com/jump-dev/JuMP.jl).jl (@odow & co)
-Please post the code below. 

4: **Discrete Time**/**Stochastic**, Difference Equations [SolveDSGE](https://github.com/RJDennis/SolveDSGE.jl).jl (@RJDennis) 
-[Code link](https://github.com/azev77/4Quadrants/blob/main/DiscreteTime_Stochastic_SolveDSGE.jl).
The corresponding simulations (using perturbation methods):

![image](https://user-images.githubusercontent.com/7883904/146828428-48702b51-a0ac-4952-af3c-d6ea88b41292.png)

5: **Discrete Time**/**Stochastic**, Bellman Equations [POMDPs](https://github.com/JuliaPOMDP/POMDPs.jl).jl (@zsunberg)To Do. 
-Please post code below using e.g. a 2-state Markov Chain for z. 

6: **Continuous Time**/**Deterministic**, Bellman Equations (my HJB Solver, based on Ben Moll) 
-[Code link](https://github.com/azev77/4Quadrants/blob/main/ContinuousTime_Deterministic_Zevelev.jl).
The corresponding simulations:

![image](https://user-images.githubusercontent.com/7883904/146828590-7f2845fb-c916-48e1-a4a4-e8b2bbfc00de.png)

7: **Continuous Time**/**Deterministic**, Bellman Equations [EconPDEs](https://github.com/matthieugomez/EconPDEs.jl).jl (@matthieu)
To Do...

8: **Continuous Time**/**Deterministic**, Hamiltonian [InfiniteOpt](https://github.com/pulsipher/InfiniteOpt.jl).jl (@pulsipher)
-To Do. Only finite horizon currently.

9: **Continuous Time**/**Deterministic**, DE [SciML](https://github.com/SciML).jl (@ChrisRackauckas). 
-TO DO. Only finite horizon currently. 
If someone knows how to transform the problem to solve an infinite horizon system please post code below.

10: **Continuous Time**/**Stochastic**, Bellman Equations @matthieu (Custom HJB Solver & EconPDEs.jl). To Do.
