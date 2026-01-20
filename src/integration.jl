function IntegrationFunction(f_integrand, grid::LGLIGrid)
    return function (states, controls, params)
        integral = grid.time.t_diff*sum(grid.quadrature.weights.*f_integrand(states, controls, params))
        return integral
    end
end