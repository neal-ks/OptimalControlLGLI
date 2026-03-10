using FastGaussQuadrature: gausslobatto

@testset "Quadrature" begin
    num_states = 1 
    num_controls = 1
    num_segments = 1 
    order = 5
    time_initial = 8.0 
    time_final = 30.0
    function integrand(states, controls, params)
        tau, weights = gausslobatto(order)
        map_time = GaussLobattoTime(8, 30)
        time = map_time.(tau)
        terms = @. 2000*log(14e4/(14e4-2100*time)) - 9.8*time
        return terms
    end 
    grid = LGLIGrid(num_states, num_controls, time_initial, time_final, num_segments, order)
    int_func = IntegrationFunction(integrand, grid)
    integral = int_func(0, 0, 0)
    @test integral ≈ 11061.34 atol=1e-2
end