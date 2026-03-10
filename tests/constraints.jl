using LinearAlgebra

# num_states::Int, num_controls::Int, time_initial, time_final, num_segments::Int=1, order::Int=3
@testset "Integration parameters" begin
    num_states = 1
    num_controls = 1 
    num_segments = 1
    order = 3 
    time_initial = 0.0 
    time_final = 1.0
    grid = LGLIGrid(num_states, num_controls, time_initial, time_final, num_segments, order)
    int_mtx = lgl_integration_matrix(grid)
    f_interp, E = integration_parameters(grid)
    @test E[:, 2:end]*int_mtx[:, 2:end] ≈ I(grid.order-1) rtol=1e-9
end