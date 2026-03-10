@testset "Example Problems" begin
    @testset "Brachistochrone" begin
        # Example from the dymos documentation
        # https://openmdao.github.io/dymos/examples/brachistochrone/brachistochrone.html
        using OptimizationBase, OptimizationIpopt, Interpolations

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
        t_init = 0.0
        t_final = 10.0

        grid = LGLIGrid(num_states, num_controls, t_init, t_final, num_segments, order)

        x0 = zero_opt_var(grid)
        x_interp = linear_interpolation([0, 10], [0.0, 10.0])
        y_interp = linear_interpolation([0, 10], [10.0, 5.0])
        v_interp = linear_interpolation([0, 10], [0, 9.9])
        theta_interp = linear_interpolation([0, 10], [5, 100.5]*pi/180)
        set_state!(x0, grid, x_interp, 1)
        set_state!(x0, grid, y_interp, 2)
        set_state!(x0, grid, v_interp, 3)
        set_control!(x0, grid, theta_interp, 1)

        prob = OPCMinimumTimeProblem(brachistone, init_bcfn, final_bcfn, x0, 0, grid)
        opt = IpoptOptimizer()
        sol = solve(prob, opt)
        @test sol.objective ≈ 1.801754704304293 atol=1e-2
        @test sol.retcode == SciMLBase.ReturnCode.Success
    end

    @testset "1D Example" begin
        # Example from: 
        # Integral form of Legendre-Gauss-Lobatto collocation for optimal control,
        # Gabriela Abadia-Doyle, William W. Hager, Anil V. Rao
        # https://www.sciencedirect.com/science/article/pii/S001600322500612X?ref=pdf_download&fr=RR-2&rr=9c3c78a56db4942d

        using OptimizationBase, OptimizationIpopt
        
        # Define the dynamics and objective to maximise the final state
        ode_func(x, u, p) = @. 2.5*(-x+x*u-u*u)
        init_boundary(x, p) = x .- 1
        final_boundary(x, p) = 0
        cost(x, u, p) = -x[end][1]

        num_states = 1
        num_controls = 1
        num_segments = 1
        order = 30
        t_init = 0.0
        t_final = 2.0
        grid = LGLIGrid(num_states, num_controls, t_init, t_final, num_segments, order)
        x0 = zero_opt_var(grid)
        
        # True analytic solution
        alpha_verify(t) = 1 + 3*exp(2.5*t)
        x_verify(t) = 4/alpha_verify(t)
        u_verify(t) = 0.5*x_verify(t)
        nodes, weights = gausslobatto(order)
        t_nodes = grid.time.(nodes)

        prob = OPCFixedTimeProblem(cost, ode_func, init_boundary, final_boundary, x0, 0, grid)
        opt = IpoptOptimizer()
        sol = solve(prob, opt)
        state_sol = vcat(get_states(sol.u, grid)...)
        control_sol = vcat(get_controls(sol.u, grid)...)
        state_true = x_verify.(t_nodes)
        control_true = u_verify.(t_nodes)
        @test state_true ≈ state_sol rtol=1e-6
        @test control_true ≈ control_sol rtol=1e-6
    end
end