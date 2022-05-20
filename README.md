# 4Quadrants
4 Quadrants of Dynamic Optimization

$V(s_{t}) \text{ } \equiv \text{ } \underset{ c_{t} }{ \sup } \text{ } E_{t}\left[ \int_{s=t}^{s=T} e^{-\rho (s-t)} u(c_t, s_t) dt \right] $      
$ds_{t} = \mu(s_{t},c_{t})dt + \sigma(s_{t},c_{t})dZ_{t} $    
$c(s_{t}) \text{ } \equiv \text{ } \underset{ c_{t} }{ \arg\sup } \text{ } E_{t}\left[ \int_{s=t}^{s=T} e^{-\rho (s-t)} u(c_t, s_t) dt \right] $      


1: **Discrete-time**/**deterministic** solution using [SolveDSGE](https://github.com/RJDennis/SolveDSGE.jl).jl (@RJDennis). [Code link](https://github.com/azev77/4Quadrants/blob/main/DiscreteTime_Deterministic_SolveDSGE.jl). 
The corresponding simulations (using 6 different methods):

![image](https://user-images.githubusercontent.com/7883904/146826183-c2b4ddbf-eba4-4f06-a253-89a1fd4c9951.png)

2: **Discrete-time**/**deterministic** solution using [POMDPs](https://github.com/JuliaPOMDP/POMDPs.jl).jl (@zsunberg) 
-[Code link](https://github.com/azev77/4Quadrants/blob/main/DiscreteTime_Deterministic_POMDPs.jl). See [discussion here](https://github.com/JuliaPOMDP/POMDPs.jl/discussions/351). 
![image](https://user-images.githubusercontent.com/7883904/147296008-4afad2bb-a4ba-4934-8a23-1cf1d12f222d.png)
![image](https://user-images.githubusercontent.com/7883904/147296019-3b7d537e-2e75-4345-8bb9-4cb4e7a84df2.png)


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


To Do: Deep learning methods
https://discourse.julialang.org/t/solving-hjb-pde-using-deep-learning/60639/13


Transform finite horizon DE, into infinite Horizon DE:
https://discourse.julialang.org/t/solving-boundary-value-differential-equation-problems-for-economics/72871

TO Do: incorporate discrete time (deterministic/stochastic) w/ VFI 
https://github.com/shanemcmiken/PEinvestment
Except try to make generic
