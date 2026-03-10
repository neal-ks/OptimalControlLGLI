# OptControlLGLI

Package to solve optimum control problems using optimum control problems using LGL-Integral method (https://www.sciencedirect.com/science/article/pii/S001600322500612X?ref=pdf_download&fr=RR-2&rr=9bfe6d95ab2dcd95).

This document currently assumes knowledge of optimal control and general optimisation principles.

The following resources are recommended as approachable learning resources for engineering/STEM students, which are also freely available on the internet.

Optimal Control:
 - Dymos Documentation (https://openmdao.github.io/dymos/index.html)
 - Optimal Control (CMU 16-745) (https://www.youtube.com/watch?v=6rUdAOCNXAU)
 - A Primer on Pontryagin's Principle in Optimal Control, Ross M

Direct Collocation focused resources:
 - An Introduction to Trajectory Optimization: How to Do Your Own Direct Collocation (https://epubs.siam.org/doi/epdf/10.1137/16M1062569)
 - Survey of Direct Transcription for Low-Thrust Space Trajectory Optimization with Applications (https://epubs.siam.org/doi/epdf/10.1137/16M1062569)

General Optimisation:
 - Engineering Design Optimisation (https://mdobook.github.io/)

## Usage

### Brachistochrone Problem
A simple example taken from the dymos documentation is the brachistochrone problem (https://openmdao.github.io/dymos/examples/brachistochrone/brachistochrone.html). The objective is to minimise the time taken for a ball to travel between two points along a curve ignoring friction or other forces. In lieu of proper documentation here is a worked example.

#### Required Packages
This package is an extension of the Optimization.jl interface. For this example the base package and the IPOPT interface are used. Plots is used later on to show the results.

```julia
using OptimizationBase, OptimizationIpopt, OptControlLGLI, Plots
```

#### Dynamics
The ODE is described by the following ODE equations:
$\dot{x} = gsin(t)$
$\dot{y} = -vcos(t)$
$\dot{v} = gcos(t)$
This can be implemented in a similar manner to the DifferentialEquations.jl package, which the interface and design of this package is inspired by/shamelessly copied from (https://docs.sciml.ai/DiffEqDocs/stable/).

```julia
function brachistone(states, controls, params)
        g = 9.805
        theta = controls[1]
        x, y, v = states
        return [v*sin(theta), -v*cos(theta), g*cos(theta)]
end
```

#### Boundary Constraints

The user should define boundary functions at the initial and final time, these take the state and parameter values and calculate a series of boundary constraints as a vector with a length equal to the number of state variables. If a variable is free then return 0 for that element of the vector.

An example of an initial boundary function.

```julia
    function init_bcfn(state, params)
        residual = [
            state[1],
            state[2] - 10,
            state[3]
        ]
        return residual
    end
```

An example of a final boundary function.

```julia
    function final_bcfn(state, params)
        residual = [
            state[1] - 10,
            state[2] - 5,
            0
        ]
        return residual
    end
```

The user should then define a grid struct. This defines a map from the user defined function which work off state and control variables to the underlying optimisation problem which uses a 1D vector of design variables and constraits. For the LGL-I method it uses a constructor with the following signiture.

```julia
LGLIGrid(num_states::Int, num_controls::Int, time_initial, time_final, num_segments::Int=1, order::Int=3)
```

Where:
num_states - The number of state variables
num_controls - The number of control variables
time_initial - The initial time of the trajectory
time_final - The final time of the trajectory
num_segments - Placeholder, multiple segments are currently implmented
order - The order of the interpolating polynomial (if this doesn't make sense read this https://openmdao.github.io/dymos/getting_started/collocation.html)

The user should then define an initial guess for the states and controls. The zero_opt_var function creates an zero vector with the correct number of design variables defined by the discretistion.

```julia
x0 = zero_opt_var(grid)
```

The user can then set guesses for the control and state variables using the set_state! and set_control! methods. These take the following arguments:
opt_var - The initial guess of the design variables which is modified in place.
grid - A grid struct which defines the discretisaiton.
interp - An Abstract Inteprolation type as defined in Interpoaltions.jl
idx - The index of the state/control variables.

```julia
x_interp = linear_interpolation([0, 10], [0.0, 10.0])
y_interp = linear_interpolation([0, 10], [10.0, 5.0])
v_interp = linear_interpolation([0, 10], [0, 9.9])
theta_interp = linear_interpolation([0, 10], [5, 100.5]*pi/180)
set_state!(x0, grid, x_interp, 1)
set_state!(x0, grid, y_interp, 2)
set_state!(x0, grid, v_interp, 3)
set_control!(x0, grid, theta_interp, 1)
```

#### But what about the cost function!
But what about the actual objective function. In this case the OPCMinimumTimeProblem implements the cost function automatically. Using the already defined dynamics, boundary function and discretisation map a OPCMinimumTimeProblem can be defined as follows:

### Solving the Brachistochrone Problem

```julia
params = 0
prob = OPCMinimumTimeProblem(brachistone, init_bcfn, final_bcfn, x0, params, grid, time_guess=1.8)
```
No parameters are required in this case, but a guess for the end time is provided. The OPCMinimumTimeProblem takes the user differential equation functions and boundary functions and cosntructs a nonlinear optimisation problem, using the LGL-I method.

Then all the user has to do is to provide a optimiser to the problem, many are available but the author has had good results with IPOPT and solve the optimisation problem.

```julia
opt = IpoptOptimizer()
sol = solve(prob, opt)
```

### Plotting and verifying the results

Now that the optimisation has finished how can the results be verified. A good start is to check the optimisers return code. For all optimisers using the Optimization.jl interface display sol.retcode.

If the optimiser gives a sucessful return code this is a good start, however doesn't garauntee success. Particularly with optimal control problems it is recommended to plot the state and control history of the optimal solution. To help with this the package defines the utility functions get_states and get_controls. These take the design variable vector (sol.u for the optimised design variables) and extract out the state and control vectors at each state discretistion node.

```julia
states = get_states(sol.u, grid)
controls = get_controls(sol.u, grid)
```

Here is an example plot which plots the path of the ball. Visually inspect the results to see if they are what you expect. If not or the return code is failing

### Debugging steps
- Try providing a good intial guess
- Check the defined functions for ODEs and boundary constraints for bugs
- Try different auto-differentiation packages to see if there are different results (make this easy for the user in future currently requries source code modification, refactoring the OPCProblem interface would make this possible)
- Look at the optimiser settings, in particular design variable scaling (bad scaling can cause numerical issues, future feature should include co-state estimation/auto-scaling)
- Re-formulate the problem, optimisers are very good are exploiting loopholes in problem formulation

## TODO List
- ~~Minimum viable example Brachistrone or 1D analytic example~~
- A better name, a lot of the obvious ones are taken
- Add more example problems
- ~~Make defining initial guess easier~~
- Make code type stable
- Better documentation
- ~~Multiple segments~~
- Easy scaling/auto-scaling
- Plotting utilities to more easily plot solutions
- Extensibility to add other constraints and paramter optimisation.
- More flexible boundary functions
- Debug plotting animation
- Lagrange term objective function utility, partially implemented but not tested (partially implemented already see intergration.jl)
- Path constraints
- Add more collocation options (Hermite Simpson, LGL, LGR etc.) with a refactor into a core/sub-libraries package.