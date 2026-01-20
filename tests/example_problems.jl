@testset "Example Problems" begin
    @testset "Brachistochrone" begin
        using OptimizationBase, OptimizationIpopt

        function brachistone(states, controls, params)
            g = 9.805
            theta = controls[1]
            x, y, v = states
            return [v*sin(theta), -v*cos(theta), g*cos(theta)]
        end

        function init_bcfn(state, params)
            residual = [
                state[1],
                state[2] - 10,
                state[3]
            ]
            return residual
        end

        function final_bcfn(state, params)
            residual = [
                state[1] - 10,
                state[2] - 5,
                0
            ]
            return residual
        end

        num_states = 3
        num_controls = 1
        num_segments = 1
        order = 30
        t_init = 0
        t_final = 10

        grid = LGLIGrid(num_states, num_controls, t_init, t_final, num_segments, order)

        x0 = zero_opt_var(grid)
        set_controls!(x0, 1, range(0.5*pi/180, 99.5*pi/180, grid.order), grid)
        set_states!(x0, 1, range(0.0, 10.0, grid.order), grid)
        set_states!(x0, 2, range(10.0, 7.0, grid.order), grid)
        set_states!(x0, 3, range(0.0, 8.0, grid.order), grid)

        prob = OPCMinimumTimeProblem(brachistone, init_bcfn, final_bcfn, x0, 0, grid, time_guess=1.8)
        opt = IpoptOptimizer()
        sol = solve(prob, opt)
        @test sol.objective ≈ 1.801754704304293 atol=1e-2
        @test sol.retcode == SciMLBase.ReturnCode.Success
    end
end