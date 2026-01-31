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
    
    segment_con_func! = function (defects, state_rate_init, states, controls, t_duration, params)
        # Evaluate the state rates at points that aern't the first in the segment
        non_init_state = @view states[2:end]
        non_init_control = @view controls[2:end]
        state_rate = vec(_f.(non_init_state, non_init_control, params))

        # The state rate function returns a derivative wrt time, this needs to be converted to the derivatives wrt tau
        # which is the transformed time to between [-1, 1]
        pushfirst!(state_rate, state_rate_init)
        nondim_state_rate = state_rate .* (0.5*t_duration)

        # Constructing the dynamic constraint system
        stacked_state_rate = stack(nondim_state_rate, dims=1)
        stacked_state = stack(states, dims=1)
        constraints = f_interp*stacked_state_rate - E*stacked_state

        defects .= vec(constraints)
        state_rate_init .= nondim_state_rate[end]
    end
    t_duration = (disc_map.time.t_final - disc_map.time.t_initial)/disc_map.num_segments
    
    return function dynamic_constraints!(defects, states, controls, params)
        # calculate the initial state rate
        state_rate_init = vec(_f(states[1], controls[1], params) * 0.5 * t_duration)
        for i=1:disc_map.num_segments
            xu_idx = ((i-1)*(disc_map.order-1) + 1):(i*(disc_map.order-1) + 1)
            defect_idx = ((i-1)*(disc_map.order-1) + 1):(i*(disc_map.order-1))
            seg_states = @view states[xu_idx]
            seg_controls = @view controls[xu_idx]
            segment_defects = @view defects[defect_idx]
            segment_con_func!(segment_defects, state_rate_init, seg_states, seg_controls, t_duration, params) 
        end
    end
end

function construct_constraint_function(_f, disc_map; initial_bfn, final_bfn)
    # Get the dynamic constraint function, uses states and controls
    dyn_con_func! = construct_dynamic_constraints(_f, disc_map)
    dyn_con_idx = 1:(disc_map.num_states*(disc_map.order-1)*disc_map.num_segments)

    bfn_init_idx = dyn_con_idx[end] .+ (1:disc_map.num_states)
    bfn_final_idx = bfn_init_idx[end] .+ (1:disc_map.num_states)

    return function constraint_func!(res, opt_var, params)
        # Transform the 1D vector into a vector of state/control vectors at each discretisation point
        states = get_states(opt_var, disc_map)
        controls = get_controls(opt_var, disc_map)
        defects = @view res[dyn_con_idx]
        dyn_con_func!(defects, states, controls, params)        
        res[bfn_init_idx] .= initial_bfn(states[1], params)
        res[bfn_final_idx] .= final_bfn(states[end], params)
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
    num_cons = disc_map.num_states*((disc_map.order-1)*disc_map.num_segments + 2)
    bounds = zeros(num_cons)
    return bounds, bounds
end