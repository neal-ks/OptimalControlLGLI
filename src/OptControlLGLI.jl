module OptControlLGLI

using LinearAlgebra
using Optimization
using FastGaussQuadrature: gausslobatto
using Integrals
using ADTypes
using ForwardDiff
using Interpolations: AbstractInterpolation

export LGLIGrid
export lgl_integration_matrix
export GaussLobattoTime, GaussLobatto, IntegrationFunction
export construct_dynamic_constraints, construct_constraint_function, construct_free_time_constraint_function, 
integration_parameters, construct_cost_function
export OPCFixedTimeProblem, OPCMinimumTimeProblem
export get_states, get_controls, set_control!, set_state!, zero_opt_var, time_grid

include("grid.jl")
include("lagrange_basis.jl")
include("integration.jl")
include("constraints.jl")
include("problem.jl")
include("utils.jl")
end