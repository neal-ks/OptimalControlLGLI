@inline function __wrap_discretisation(des_var, num_params, num_nodes)
    wrapped_vector = [des_var[(1:num_params) .+ (i-1)*num_params] for i in 1:num_nodes]
    return wrapped_vector
end

@inline function get_states(opt_var, disc_map::LGLIGrid)
    num_state_var = disc_map.num_segments*(disc_map.order-1) + 1
    return __wrap_discretisation(opt_var[disc_map.state_idx], disc_map.num_states, num_state_var)
end

@inline function get_controls(opt_var, disc_map::LGLIGrid)
    num_state_var = disc_map.num_segments*(disc_map.order-1) + 1
    return __wrap_discretisation(opt_var[disc_map.control_idx], disc_map.num_controls, num_state_var)
end

@inline function set_controls!(opt_var, idx, control, disc_map)
    control_indicies = disc_map.control_idx[idx:disc_map.num_controls:end]
    opt_var[control_indicies] .= control
end

@inline function set_states!(opt_var, idx, state, disc_map)
    state_indicies = disc_map.state_idx[idx:disc_map.num_states:end]
    opt_var[state_indicies] .= state
end

@inline function zero_opt_var(grid::LGLIGrid)
    return zeros((grid.num_segments*(grid.order-1) + 1)*(grid.num_states + grid.num_controls))
end