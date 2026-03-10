function construct_cost_function(cost_func, disc_map::LGLIGrid)
    return function (opt_var, params)
        states = get_states(opt_var, disc_map)
        controls = get_controls(opt_var, disc_map)
        return cost_func(states, controls, params)
    end
end

function OPCFixedTimeProblem(cost_func, state_rate_func, boundary_init, boundary_final, x0, params, disc_map::LGLIGrid)
    opt_func = OptimizationFunction(
        construct_cost_function(cost_func, disc_map), 
        ADTypes.AutoForwardDiff(),
        cons = construct_constraint_function(
            state_rate_func, 
            disc_map; 
            initial_bfn=boundary_init, 
            final_bfn=boundary_final
        )
    )
    lcons, ucons = get_constraint_bounds(disc_map)
    return OptimizationProblem(opt_func, x0, params, lcons=lcons, ucons=ucons)
end

function OPCMinimumTimeProblem(state_rate_func, boundary_init, boundary_final, x0, params, disc_map::LGLIGrid; control_lb=nothing, control_ub=nothing)
    x0 = [x0; disc_map.time.t_final]
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
    num_des_var = length(x0)
    lb = fill(-Inf, num_des_var)
    ub = fill(Inf, num_des_var)
    if !isnothing(control_lb)        
        lb[disc_map.control_idx] .= control_lb
    end
    if !isnothing(control_ub)
        ub[disc_map.control_idx] .= control_ub
    end
    return OptimizationProblem(opt_func, x0, params, lb=lb, ub=ub, lcons=lcons, ucons=ucons)
end