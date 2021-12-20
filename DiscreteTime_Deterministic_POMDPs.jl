δ=0.10; ρ=0.15; # Param 
z= ρ + δ +.02     # Need: z > ρ + δ
I_SS = (z-ρ-δ)/(ρ+δ)
K_SS = I_SS/δ

if 1==1
    a_0(s;δ=δ) = δ*s
    min_s = eps(); max_s = 2.0*K_SS; n_s = 200;
    min_a = eps(); max_a = a_0(max_s); n_a = 400; # keeps a≥0 
    "states:"
    states = range(min_s, max_s, length=n_s)
    ceil_state(s) = states[searchsortedfirst(states, s)]
    "actions:" 
    valid_actions()  = range(min_a, max_a, length=n_a)        # possible actions any state.
    valid_actions(s) = filter(<(a0(s)), valid_actions())      # valid actions @ state=s 
    "Transition:" 
    μ(s,a;δ=δ)      = a - δ*s     
    "reward:" 
    r(s,a;z=z) = z*s -a -0.5*(a^2)
    "discount:"
    β = 1/(1+ρ);
    using QuickPOMDPs: QuickMDP           #QuickMDP()
    using POMDPModelTools: Deterministic
    m = QuickMDP(
        states     = states,
        actions    = valid_actions,
        transition = (s, a) -> Deterministic( ceil_state( μ(s,a) ) ),
        reward     = (s, a) -> r(s,a),
        discount   = β
    )
    # DiscreteValueIteration: Both work fast! 
    using DiscreteValueIteration
    s1 = DiscreteValueIteration.SparseValueIterationSolver()
    s2 = DiscreteValueIteration.ValueIterationSolver()
    @time sol1 = solve(s1, m) #
    @time sol2 = solve(s2, m) #
    # value(sol1, states[2]), action(sol1, states[2])
    # value(sol2, states[2]), action(sol2, states[2])
    value(sol1, states[end])

    using Plots
    # Value
    plot(legend=:bottomright, title="Value Functions");
    #plot!(states[2:end], i->A0 + A1 * log(i), lab="closed form") 
    plot!(states[2:end], i->value(sol1, i),   lab="sol1")
    plot!(states[2:end], i->value(sol2, i),   lab="sol2")
    plot!(states[2:end], i->value(sol4, i),   lab="sol4")
    plot!(states[2:end], i->value(sol3, i),   lab="sol3")

    # Simulation
    Tsim=150; s0=states[2]; 
    sim1 = []; push!(sim1, s0); 
    for tt in 1:Tsim
        tt==1 ? s = s0 : nothing 
        #a = valid_actions()[sol1.policy[searchsortedfirst(states, s)]]
        a = action(sol1, s)
        sp = μ(s,a) 
        sp = ceil_state(sp)
        push!(sim1, sp)
        s = sp
    end 

    sim4 = []; push!(sim4, s0); 
    for tt in 1:Tsim
        tt==1 ? s = s0 : nothing 
        a = action(sol4, s)
        sp = μ(s,a) 
        push!(sim4, sp)
        s = sp
    end 

    plot(legend=:bottomright, title="Simulation");
    #plot!(simcf,   lab="closed form")
    plot!(sim1,   lab="sol1")
    plot!(sim4,   lab="sol4")
end
