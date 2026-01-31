struct LGLBasis
    nodes
    basis_node
    function LGLBasis(nodes, basis_idx::Int)
        node_mask = eachindex(nodes) .!= basis_idx
        return new(nodes[node_mask], nodes[basis_idx])
    end
end

function (basis::LGLBasis)(querry_node, params)
    return prod(@. (querry_node - basis.nodes) / (basis.basis_node - basis.nodes))
end

function lgl_integration_matrix(disc_map::LGLIGrid)
    integration_mtx = zeros(disc_map.order-1, disc_map.order)
    nodes = disc_map.quadrature.nodes
    for j in eachindex(nodes)
        lagrange_basis = LGLBasis(nodes, j)
        for i in eachindex(nodes)[2:end]
            integral_prob = IntegralProblem(lagrange_basis, [-1, nodes[i]])
            sol = solve(integral_prob, QuadGKJL())
            integration_mtx[i-1, j] = sol.u
        end
    end
    return integration_mtx
end
