function integration_parameters(disc_map)
    A_reduced = lgl_integration_matrix(disc_map)
    A_sq = A_reduced[:, 2:end]
    A_inv = inv(A_sq)

    # Construct parameters for the dynamic constraints
    alpha = A_inv * A_reduced[:,1]
    E = A_inv*hcat(-ones(disc_map.order-1), I(disc_map.order-1))
    f_interp = hcat(alpha, I(disc_map.order-1))
    return f_interp, E
end

function construct_dynamic_constraints(_f, disc_map)
    f_interp, E = integration_parameters(disc_map)
    return function (states, controls, t_duration, params)
        state_rate = _f.(states, controls, params)
        nondim_state_rate = state_rate .* (0.5*t_duration)
        stacked_state_rate = stack(nondim_state_rate, dims=1)
        stacked_state = stack(states, dims=1)
        constraints = f_interp*stacked_state_rate - E*stacked_state
        return vec(constraints)
    end
end

function construct_constraint_function(_f, disc_map; initial_bfn, final_bfn)
    # Get the dynamic constraint function, uses states and controls
    dyn_con_func = construct_dynamic_constraints(_f, disc_map)

    return function constraint_func(res, opt_var, params)
        # Transform the 1D vector into a vector of state/control vectors at each discretisation point
        states = get_states(opt_var, disc_map)
        controls = get_controls(opt_var, disc_map)
        dyn_cons = dyn_con_func(states, controls, 2*disc_map.time.t_diff, params)        
        bcp_init = initial_bfn(states[1], params)
        bcp_final = final_bfn(states[end], params)
        res .= [dyn_cons; bcp_init; bcp_final]
    end
end

function construct_free_time_constraint_function(_f, disc_map; initial_bfn, final_bfn)
    # Get the dynamic constraint function, uses states and controls
    dyn_con_func = construct_dynamic_constraints(_f, disc_map)

    return function (res, opt_var, params)
        # Transform the 1D vector into a vector of state/control vectors at each discretisation point
        states = get_states(opt_var, disc_map)
        controls = get_controls(opt_var, disc_map)
        # last element of the design variable vector is the time taken, rescaling the derivatives as time taken
        # changes
        dyn_cons = dyn_con_func(states, controls, opt_var[end], params)        
        bcp_init = initial_bfn(states[1], params)
        bcp_final = final_bfn(states[end], params)
        res .= [dyn_cons; bcp_init; bcp_final]
    end
end

function get_constraint_bounds(disc_map)
    num_cons = disc_map.num_states*(disc_map.order + 1)
    bounds = zeros(num_cons)
    return bounds, bounds
end