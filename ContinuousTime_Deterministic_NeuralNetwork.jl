# based on Gustavo Mellior's: https://github.com/GMellior/NGM_neuralnetwork
using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using ForwardDiff
import ModelingToolkit: Interval, infimum, supremum

@parameters x # x=K
@variables u(..), s(..) # u(x)=V(K) & s(x) = I(K) policy
Dx = Differential(x)  # Dx(u(x))  == V_K(K)
Vx(x) = Dx(u(x))
#
δ=0.10; ρ=0.05; z= ρ + δ +.10; Iss = (z-ρ-δ)/(ρ+δ); kss = Iss/δ;
eq = [s(x) ~ Vx(x)-1, ρ*u(x) ~ z*x - s(x) - 0.5*s(x)^2 + Vx(x)*(s(x) - δ*x)]
bcs  = [Vx(kss) ~ 1 + δ*kss]
domains = [x ∈ Interval(kss*0.1,kss*1.5)]
#
elems = size(eq)[1]
neurons = 20
chain             = [FastChain(FastDense(1,neurons,tanh)
                    ,FastDense(neurons,neurons,softplus)
                    ,FastDense(neurons,1)) for _ in 1:elems]
dx0               = kss*(1.5-0.1)/400
initθ             = map(c -> Float64.(c), DiffEqFlux.initial_params.(chain))
discretization    = PhysicsInformedNN(chain, GridTraining(dx0),init_params =initθ)
@named pde_system = PDESystem(eq,bcs,domains,[x],[u(x),s(x)])
prob = discretize(pde_system,discretization)
const losses = []
cb = function (p,l)
    push!(losses, l)
    if length(losses)%500==0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end
maxit1 = 1200; maxit2 = 4000;
res  = GalacticOptim.solve(prob, ADAM(0.06); atol=1e-7, cb = cb, maxiters=maxit1)
prob = remake(prob,u0=res.minimizer)
res  = GalacticOptim.solve(prob,BFGS(); cb = cb, maxiters=maxit2)
phi  = discretization.phi

using Plots
kgrid      = LinRange(kss*0.1,kss*1.5,1_000)
I          = Int.(collect(size(res.minimizer)./elems))
Vfunc(x)   = first(phi[1](x,res.minimizer[1:I[1]]))
Cfunc(x)   = first(phi[2](x,res.minimizer[I[1]+1:end]))

dVk = ForwardDiff.derivative.(Vfunc, Array(kgrid))
II  = dVk .- 1

# Plot V
plot(kgrid ,Vfunc.(kgrid),
    #line =:dashdot,
    linewidth=3,label="Neural net",
    legend=:topleft,xlabel = "k",ylabel="V(k)", 
    fg_legend = :false, xtick=0:2:10, ytick=-20:1:-8)
vline!([kss],lab="Kss")
hline!([(1/ρ)*(z*kss -Iss -.5*(Iss)^2 )],lab="Vss")

#display(plot!(kgrid,vfunc,label="Finite differences",legend=:topleft))

plot(kgrid ,Cfunc.(kgrid),
    line =:dashdot,linewidth=3,label="Neural net",
    legend=:topleft,xlabel = "k",ylabel="I(k)", 
    fg_legend = :false, xtick=0:2:10, ytick=0.5:0.5:2)
plot!(kgrid ,x -> Iss, lab="I_SS", ylims=(0,1))
#

# Compute HJB error
HJBerror = z*kgrid .- II .- .5*II.^2 + dVk .*(II-δ*kgrid) - ρ*Vfunc.(kgrid)
plot(kgrid,HJBerror,linewidth=3,ylims=(-1.2e-4,1e-4),
    xlabel = "k",ylabel="HJB error",label="Neural net", fg_legend = :false)
########
display(plot(losses,legend = false, yaxis=:log))






