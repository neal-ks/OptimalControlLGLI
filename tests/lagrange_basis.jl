using FastGaussQuadrature: gausslobatto

@testset "Integration Matrix Test" begin
    order = 5
    disc_map = LGLIGrid(1, 1, 0.0, 1.0, 1, order)
    integration_mtx = lgl_integration_matrix(disc_map)
    nodes, weights = gausslobatto(order)
    @test integration_mtx[end, :] ≈ weights rtol=1e-6
end