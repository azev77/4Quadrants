using LinearAlgebra, SparseArrays, Plots
if 1==1
    if 1==1
        ρ = 0.05; δ = 0.10; A = 0.20;
        r(s,a) = A*s -a -0.5*(a)^2.0
        μ(s,a) = a - δ*s
        dr(s,a) = -1 -a 
        FOC(s,Vs) = Vs - 1
        μ_inv(s,ṡ) = ṡ + δ*s 
        s_ss = (A-ρ-δ)/(δ*(ρ+δ))
        a_ss = δ*s_ss
        v_ss = r(s_ss,a_ss)/ρ        # ∫exp(-ρt)*r(s_ss,a_ss) dt = r(s_ss,a_ss)/ρ    
        #s_min = 0.00*s_ss  
        s_min = -.5  
        s_max = 2.000*s_ss
        # verify things are well defined @ corners 
        s_max, μ_inv(s_max,0), dr(s_max,μ_inv(s_max,0))
        s_min, μ_inv(s_min,0), dr(s_min,μ_inv(s_min,0))
        H = 10_000;
        s = collect(LinRange(s_min, s_max, H))
        ds = (s_max-s_min)/(H-1)
        dVf, dVb            = [zeros(H,1) for i in 1:2]
        dV_Upwind, a_Upwind = [zeros(H,1) for i in 1:2]
        dVf_end       = dr(s_max,μ_inv(s_max,0))  
        dVb_1         = dr(s_min,μ_inv(s_min,0))
        v0 = v_ss *ones(H)
        v0 = @. r(s, μ_inv(s,0))/ρ #initial guess for V
        v = v0
    end    
    Δ = 1_000; maxit = 10_000;ε = 10e-8; dist=[];

    # Generic HJB Solver
    # Parimonious: dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0
    for n=1:maxit
        V=v
        dV = (V[2:end] - V[1:end-1])/ds  

        # forward difference: if ṡ>0
        dVf = [dV; dVf_end]
        af  = FOC.(s,dVf)                 # choice w forward difference
        μ_f = μ.(s, af)                   # ṡ      w forward difference
        If  = μ_f .> 0

        # backward difference: if ṡ<0
        dVb = [dVb_1; dV]
        ab  = FOC.(s,dVb)                          # choice w backward difference
        μ_b = μ.(s, ab)
        Ib  = μ_b .< 0

        # neither difference: if ṡ=0
        a0  = μ_inv.(s,0)        # c if ṡ=0
        dV0 = dr.(s,a0)
        μ_0 = μ.(s, a0)          # μ_0 == zero(s)
        I0  = (1.0 .- If - Ib)   # choice betw forward & backward difference

        # I_concave = dVb .> dVf
        # scatter(I_concave) #1 everywhere EXCEPT @ last point H. 

        dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0   
        a_Upwind  = FOC.(s,dV_Upwind)    
        μ_Upwind  = μ.(s, a_Upwind)                                            
        u = r.(s,a_Upwind)

        # create the transition matrix AA
        X = -min.(μ_b,0)/ds
        Y = -max.(μ_f,0)/ds + min.(μ_b,0)/ds
        Z = max.(μ_f,0)/ds

        a1 = sparse(Diagonal((Y[:])))
        a2 = [zeros(1,H); sparse(Diagonal((X[2:H]))) zeros(H-1,1)]
        a3 = [zeros(H-1,1) sparse(Diagonal((Z[1:H-1]))); zeros(1,H)]
        AA = a1 + a2 + a3

        # Solve for new value 
        B = (ρ + 1/Δ)*sparse(I,H,H) - AA
        b = u + V./Δ
        V = B \ b

        # 8: stopping criteria 
        V_change = V-v
        v = V 

        push!(dist,findmax(abs.(V_change))[1])
        println(n, " ", dist[n])
        if dist[n] .< ε
            println("Value Function Converged Iteration=")
            println(n)
            break
        end
    end
    s_dot = μ.(s, a_Upwind)
    v_err = r.(s, a_Upwind) + dV_Upwind.*s_dot - ρ.*v # approx @ borrowing constraint
    plot(s, v, xlabel="s", ylabel="V(s)",legend=false, title="")
    plot!([s_ss],  seriestype = :vline, lab="", color="grey", l=:dash)
    plot!([v_ss],  seriestype = :hline, lab="", color="grey", l=:dash)

    # Simulate.
    using Interpolations
    â = LinearInterpolation(s, a_Upwind[:], extrapolation_bc = Line())
    Δt = 0.01; T = 150; time = 0.0:Δt:T
    s_sim, a_sim, ṡ_sim = [zeros(length(time),1) for i in 1:3]
    s_0 =0.5*s_ss 
    s_sim[1] = s_0 
    for i in 2:length(time)
        a_sim[i-1] = â(s_sim[i-1]) 
        ṡ_sim[i-1] = μ(s_sim[i-1], a_sim[i-1])
        s_sim[i] = s_sim[i-1] + Δt * ṡ_sim[i-1]
        # s_sim[i] = s_sim[i-1] + Δt * μ(s_sim[i-1], a_sim[i-1])
    end
    ix = 1:(length(time)-1)
    plot()
    plot!(time[ix], s_sim[ix], lab = "k")
    plot!(time[ix], a_sim[ix], lab = "i")
    plot!(time[ix], ṡ_sim[ix], lab = "k̇")
    plot!([s_min],  seriestype = :hline, lab="", color="grey", l=:dash)
    plot!([s_ss a_ss],  seriestype = :hline, lab="", color="grey", l=:dash)
end
