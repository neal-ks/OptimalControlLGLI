struct GaussLobattoTime
    t_initial
    t_final
    t_diff
    t_mean
    function GaussLobattoTime(t_initial, t_final, num_segments=1)
        time_diff = t_final - t_initial
        time_mean = 0.5*(t_initial + t_final)
        return new(t_initial, t_final, time_diff, time_mean)
    end 
end

function (glt::GaussLobattoTime)(tau)
    return 0.5*glt.t_diff*tau + glt.t_mean
end

struct GaussLobatto
    num_nodes::Int
    nodes
    weights
    function GaussLobatto(num_nodes::Int)
        nodes, weights = gausslobatto(num_nodes)
        return new(
            num_nodes, 
            nodes, 
            weights
        )
    end
end

"""
    Calculates the location of the nodes of a multi-segment LGL-I given discretisation parameters
"""
function time_grid(nodes, num_segments::Int, num_disc::Int, quadrature_time::GaussLobattoTime)
    nodes_reduced = @view nodes[1:end-1]
    time_type = typeof(quadrature_time.t_diff)
    t_grid = Vector{time_type}(undef, num_disc)
    t_diff = 2/num_segments
    t_diff_half = 0.5*t_diff
    idx = 1
    
    for i in 1:num_segments
        t_seg_init = (i-1)*t_diff - 1
        for node in nodes_reduced
            t_grid[idx] = t_seg_init + t_diff_half*(node+1)
            idx += 1
        end
    end
    t_grid[end] = 1

    return t_grid
end

struct LGLIGrid
    num_states::Int
    num_controls::Int
    num_segments::Int
    order::Int
    state_idx
    control_idx
    time::GaussLobattoTime
    collocation_nodes
    quadrature::GaussLobatto
    function LGLIGrid(num_states::Int, num_controls::Int, time_initial::Float64, time_final::Float64, num_segments::Int=1, order::Int=3)
        order > 1 || error("Order on the LGL-I transciption must be greater than 1")
        num_segments > 0 || error("There are $num_segments in the transciption. Transcription must have at least 1 segment.")
        num_disc_points = (order-1)*num_segments + 1
        state_idx = 1:(num_states*num_disc_points)
        control_idx = state_idx[end] .+ (1:num_controls*num_disc_points)
        time = GaussLobattoTime(time_initial, time_final)
        quadrature = GaussLobatto(order)
        collocation_nodes = time_grid(quadrature.nodes, num_segments, num_disc_points, time)
        
        return new(
            num_states, 
            num_controls,
            num_segments, 
            order, 
            state_idx, 
            control_idx,
            time,
            collocation_nodes,
            quadrature
        )
    end
end