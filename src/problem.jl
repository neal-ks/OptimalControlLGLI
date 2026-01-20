function OPCFixedTimeProblem(cost_func, state_rate_func, boundary_init, boundary_final, state_guess, control_guess, params, disc_map)
    opt_func = OptimizationFunction(
        construct_cost_function(cost_func, disc_map), 
        ADTypes.AutoForwardDiff(),
        cons = construct_constraint_function(state_rate_func, disc_map; initial_bfn=boundary_init, final_bfn=boundary_final)
    )
    state0 = vec(state_guess)
    control0 = vec(control_guess)
    x0 = [state0; control0]
    lcons, ucons=get_constraint_bounds(disc_map)
    return OptimizationProblem(opt_func, x0, params, lcons=lcons, ucons=ucons)
end

function OPCMinimumTimeProblem(state_rate_func, boundary_init, boundary_final, x0, params, disc_map::LGLIGrid; time_guess=1.)
    x0 = [x0; time_guess]
    cost(opt_var, p) = opt_var[end]
    opt_func = OptimizationFunction(
        cost, 
        ADTypes.AutoForwardDiff(),
        cons = construct_free_time_constraint_function(
            state_rate_func, 
            disc_map; 
            initial_bfn=boundary_init, 
            final_bfn=boundary_final
        )
    )
    lcons, ucons = get_constraint_bounds(disc_map)
    return OptimizationProblem(opt_func, x0, params, lcons=lcons, ucons=ucons)
end