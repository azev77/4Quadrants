if 1==1
    using EconPDEs, Plots
    stategrid = OrderedDict(:k => range(-7.0, 7.0, length = 1000))
    solend = OrderedDict(:V => ones(1000))
    function f(state::NamedTuple, sol::NamedTuple)
        #ρ = 0.05; δ = 0.10; z = 0.20;
        ρ = 0.15; δ = 0.10; z = 0.27;
        k = state.k
        V, Vk_up, Vk_down, Vkk = sol.V, sol.Vk_up, sol.Vk_down, sol.Vkk
        #
        Vk = Vk_up
        iter = 0
        @label start
        i = Vk - 1.0
        μk = i - δ*k
        if (iter == 0) & (μk <= 0)
            iter += 1
            Vk = Vk_down
            @goto start
        end 
        Vt = - (z*k - i -0.5*(i^2.0) + μk*Vk - ρ*V)
        (Vt = Vt,), (Vk=Vk,)
    end
    y, residual_norm =  pdesolve(f, stategrid, solend)
    plot(y[:V])
    #To Do: simulate
end 
