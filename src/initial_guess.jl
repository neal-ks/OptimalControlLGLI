@inline function __wrap_discretisation(des_var, num_params, num_nodes)
    wrapped_vector = [des_var[(1:num_params) .+ (i-1)*num_params] for i in 1:num_nodes]
    return wrapped_vector
end

@inline function get_states(opt_var, disc_map::LGLIGrid)
    num_nodes = disc_map.num_segments*(disc_map.order-1) + 1
    return __wrap_discretisation(opt_var[disc_map.state_idx], disc_map.num_states, num_nodes)
end

@inline function get_states(opt_var, disc_map::LGLIGrid, idx::Int)
    idx = disc_map.state_idx[idx:disc_map.num_states:end]
    return opt_var[idx]
end

@inline function get_controls(opt_var, disc_map::LGLIGrid, idx::Int)
    idx = disc_map.control_idx[idx:disc_map.num_controls:end]
    return opt_var[idx]
end

@inline function get_controls(opt_var, disc_map::LGLIGrid)
    num_nodes = disc_map.num_segments*(disc_map.order-1) + 1
    return __wrap_discretisation(opt_var[disc_map.control_idx], disc_map.num_controls, num_nodes)
end

function set_state!(opt_var::Vector, grid::LGLIGrid, interpolation::AbstractInterpolation, idx::Int)
    states = @views get_states(opt_var, grid, idx)
    states .= interpolation(grid.time.(grid.collocation_nodes))
end

function set_control!(opt_var::Vector, grid::LGLIGrid, interpolation::AbstractInterpolation, idx::Int)
    states = @views get_controls(opt_var, grid, idx)
    states .= interpolation(grid.time.(grid.collocation_nodes))
end