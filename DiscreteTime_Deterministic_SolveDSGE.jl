using SolveDSGE, Plots;
if 1==1
    open("Det_INV.txt", "w") do io
        println(io,"# inv\n")
        println(io,"states:\nk \nend\n")
        println(io,"jumps:\ni \nend\n")
        println(io,"shocks:\nend\n")
        println(io,"parameters:\nδ=0.10, ρ=0.15, z=0.27 \nend\n")
        println(io,"equations:")
        println(io,"k(+1) = (1-δ)*k + i")
        println(io,"i(+1) = ((ρ+δ-z)/(1-δ)) +((1+ρ)/(1-δ))*i")
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
    ss ≈ [K_SS, I_SS]
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

    # ProjectionSchemes: ChebyshevSchemes/SmolyakSchemes/PiecewiseLinearSchemes
    g0 = ss;                  #initial guess
    dom = [sss*2.0; sss*0.01] #domain for state vars. 2*n_s array.  #row 1 upper vals. row 2 lower vals
    tf = 1e-8;  #inner loop tol 
    tv = 1e-6;  #outer loop tol
    n_thr = 4; #num threads to use. 
    #nq = 9; #number quadrature points. Int. #Don't need if deterministic!!!!
    ###########################################################################
    # Projection, ChebyshevSchemeDet
    # soln_cha
    ###########################################################################
    if 1==1
        ng = chebyshev_nodes; #node generator: chebyshev_nodes/chebyshev_extrema
        #nn1 = [21]; #number of nodes for each state var #vec integers nn2 = [71];
        
        nn1 = 21 * ones(ngm.number_states) .|> Int64
        nn2 = 71 * ones(ngm.number_states) .|> Int64
        no1 = 6; #order of Cheby poly 

        P    = ChebyshevSchemeDet(g0,ng,nn1,no1,dom,tf,tv,maxiters)
        PP   = ChebyshevSchemeDet(g0,ng,nn1,4,  dom,tf,tv,maxiters)
        PPP  = ChebyshevSchemeDet(g0,ng,nn2,6,  dom,tf,tv,maxiters)
        
        #
        soln_cha = solve_model(ngm,P)
        soln_cha = solve_model(ngm,P,n_thr)
        soln_chb = solve_model(ngm,PP,n_thr)
        soln_chc = solve_model(ngm,PPP,n_thr)

    end 

    ###########################################################################
    # Projection, SmolyakSchemeDet            # soln_sa
    # Projection, PiecewiseLinearSchemeStoch  # soln_pl
    # Compose                                 # soln_c1
    ###########################################################################
    # Smolyak
    if 1==1 
        ng = chebyshev_gauss_lobatto; #node generator: chebyshev_gauss_lobatto/clenshaw_curtis_equidistant
        nl1 = 3; #number of layers to be used in the approximation
        no1 = 6; #order of Cheby poly 
        S = SmolyakSchemeDet(g0,ng,nl1,dom,tf,tv,maxiters)
        soln_sa = solve_model(ngm,S,n_thr)
    end 

    # PL 
    if 1==1
        nn1 = [21]; #number of nodes for each state var #vec integers 
        PL = PiecewiseLinearSchemeDet(g0,nn1,dom,tf,tv,maxiters)
        soln_pl = solve_model(ngm,PL,n_thr)
    end 
    # Compose 
    soln_c1 = solve_model(ngm,soln_to,PPP)
    soln_c2 = solve_model(ngm,soln_c1,PP)
    # Simulation 
    if 1 == 1
        T = 100 # num Sim pds 
        # T=10
        ω = 1.0 - 0.2; # distance from steady-state
        #ω = 1.0 - 0.75;

        sim_d1   = simulate(soln_fo, sss*ω,  T)
        sim_d2   = simulate(soln_so, sss*ω,  T)
        sim_d3   = simulate(soln_to, sss*ω,  T)
        sim_d4   = simulate(soln_cha, sss*ω, T)
        sim_d5   = simulate(soln_sa, sss*ω,  T)
        sim_d6   = simulate(soln_pl, sss*ω,  T)
        sim_d7   = simulate(soln_c1, sss*ω,  T)

        # Closed Form Sol 

        lv = reshape(ngm.variables,1,ngm.number_variables)
        p  = plot(legend=:bottomright);
        #p0 = plot!(1:T,cf,      title="Closed Form",    lab="");
        p0 = plot!([K_SS I_SS],  seriestype = :hline, lab="", color="grey", l=:dash)
        p1 = plot(p0, 1:T,sim_d1', title="Pert. 1st Order", lab=lv)
        p2 = plot(p0, 1:T,sim_d2', title="Pert. 2nd Order", lab=lv)
        p3 = plot(p0, 1:T,sim_d3', title="Pert. 3rd Order", lab=lv)
        p4 = plot(p0, 1:T,sim_d4', title="Proj. Cheby 1",   lab=lv)
        p5 = plot(p0, 1:T,sim_d5', title="Proj. Smol  1",   lab=lv)
        p6 = plot(p0, 1:T,sim_d6', title="Proj. PL    1",   lab=lv)
        p7 = plot(p0, 1:T,sim_d7', title="Proj. & Pert",    lab=lv)
        #p8 = plot!(1:T,sim_d8', title="Proj. & Pert",   lab=lv)
        plot(p1,p2,p3,p4,p5,p6, size=1.25 .* (600, 400))
    end
    #
end    
