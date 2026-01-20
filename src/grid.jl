struct GaussLobattoTime
    t_initial
    t_final
    t_diff
    t_mean
    function GaussLobattoTime(t_initial, t_final, num_segments=1)
        time_diff = 0.5*(t_final - t_initial)
        time_mean = 0.5*(t_initial + t_final)
        return new(t_initial, t_final, time_diff, time_mean)
    end 
end

function (glt::GaussLobattoTime)(tau)
    return glt.t_diff*tau + glt.t_mean
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

struct LGLIGrid
    num_states::Int
    num_controls::Int
    num_segments::Int
    order::Int
    state_idx
    control_idx
    time::GaussLobattoTime
    quadrature::GaussLobatto
    function LGLIGrid(num_states::Int, num_controls::Int, time_initial, time_final, num_segments::Int=1, order::Int=3)
        order > 1 || error("Order on the LGL-I transciption must be greater than 1")
        num_segments > 0 || error("Transcription must have at least 1 segment. Check the order of the function arguments")
        state_idx = 1:num_states*order
        control_idx = state_idx[end] .+ (1:num_controls*order)
        quad = GaussLobatto(order)
        time = GaussLobattoTime(time_initial, time_final)
        return new(
            num_states, 
            num_controls,
            num_segments, 
            order, 
            state_idx, 
            control_idx,
            time,
            quad
        )
    end
end