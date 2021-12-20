using SolveDSGE, Plots;
if 1==1
    open("Det_INV.txt", "w") do io
        println(io,"# inv\n")
        println(io,"states:\nk, z \nend\n")
        println(io,"jumps:\ni \nend\n")
        println(io,"shocks:\nε \nend\n")
        println(io,"parameters:\nδ=0.10, ρ=0.15 \nend\n")
        println(io,"equations:")
        println(io,"k(+1) = (1-δ)*k + i")
        println(io,"i(+1) = ((ρ+δ- (0.5 + 0.5*z) )/(1-δ)) +((1+ρ)/(1-δ))*i")
        println(io,"z(+1) = 0.5 + 0.5*z + 1*ε")
        println(io,"end")
    end;
    
    ###########################################################################
    path = joinpath(@__DIR__,"Det_INV.txt")
    process_model(path)
    ngm = retrieve_processed_model(path)
    ###########################################################################
    # Steady State 
    ###########################################################################
    tol = 1e-8
    maxiters = 1000
    x = .25*ones(ngm.number_variables) #x = [0.25,0.25] #x = [0.25,0.25, .25] #x = [0.25,0.25, .25, .25] 
    ss = compute_steady_state(ngm, x, tol, maxiters)
    ss = ss.zero
    #ss ≈ [K_SS, I_SS]
    sss = ss[1:ngm.number_states] #ss for states only 
    ###########################################################################
    # Perturbation, Order = 1,2,3
    # soln_fo/soln_so/soln_to
    ###########################################################################
    if 1==1
        g0 = ss;          # point about which to perturb the model (generally the steady state) 
        cutoff = 1.0      # param that separates unstable from stable eigenvalues
        # order = "first" # order of the perturbation
        N   = PerturbationScheme(g0,cutoff,"first")   # Klein (2000)
        NN  = PerturbationScheme(g0,cutoff,"second")  # Gomme and Klein (2011)
        NNN = PerturbationScheme(g0,cutoff,"third")   # Binning (2013) w refinement from Levintal (2017)

        soln_fo = solve_model(ngm, N)
        soln_so = solve_model(ngm, NN)
        soln_to = solve_model(ngm, NNN)
    end

    # To Do: Stochastic Projection 
    # ProjectionSchemes: ChebyshevSchemes/SmolyakSchemes/PiecewiseLinearSchemes

    # Simulation 
    if 1 == 1
        T = 100 # num Sim pds 
        ω = 1.0 - 0.2; # distance from steady-state
        #ω = 1.0 - 0.75;

        sim_d1   = simulate(soln_fo, sss*ω,  T)
        sim_d2   = simulate(soln_so, sss*ω,  T)
        sim_d3   = simulate(soln_to, sss*ω,  T)
        # sim_d4   = simulate(soln_cha, sss*ω, T)
        # sim_d5   = simulate(soln_sa, sss*ω,  T)
        # sim_d6   = simulate(soln_pl, sss*ω,  T)
        # sim_d7   = simulate(soln_c1, sss*ω,  T)

        # Closed Form Sol 

        lv = reshape(ngm.variables,1,ngm.number_variables)
        p  = plot(legend=:bottomright);
        #p0 = plot!(1:T,cf,      title="Closed Form",    lab="");
        p0 = plot!()
        p1 = plot(p0, 1:T,sim_d1', title="Pert. 1st Order", lab=lv)
        p2 = plot(p0, 1:T,sim_d2', title="Pert. 2nd Order", lab=lv)
        p3 = plot(p0, 1:T,sim_d3', title="Pert. 3rd Order", lab=lv)
        # p4 = plot(p0, 1:T,sim_d4', title="Proj. Cheby 1",   lab=lv)
        # p5 = plot(p0, 1:T,sim_d5', title="Proj. Smol  1",   lab=lv)
        # p6 = plot(p0, 1:T,sim_d6', title="Proj. PL    1",   lab=lv)
        # p7 = plot(p0, 1:T,sim_d7', title="Proj. & Pert",    lab=lv)
        #p8 = plot!(1:T,sim_d8', title="Proj. & Pert",   lab=lv)
        plot(p1,p2,p3,size=1.25 .* (600, 400))
    end
    #
end    
