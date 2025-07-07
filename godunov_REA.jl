
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the functions to calculate Godunov's 1D REA algorithm given an initial profile
and other parameters described throught the "parameters" dictionary. Additionally, contains
the function for solving 2D advection using 1D REA alongside Dimensional-Splitting.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""


##################################
#@       PACKAGES & CALLS        
##################################

using OrderedCollections

include("slopes_and_limiters.jl")
include("numerical_integration.jl")



##################################
#$   PREPARING THE SIMULATION
##################################

"""
Reads the simulation parameters at SIMULATION1D_PARAMETERS, extracting the necessary
parameters as well as creating the numerical mesh used for the resolution.

1. [arg] simulation_params: OrderedDict storing all the simulation parameters.

[return] Temporal/spatial step, adv speed, CFL num, slope and limiter, t points and spatial mesh.
"""
function simulation_parameters1D(simulation_params:: OrderedDict) :: NTuple{9, Union{Real, Vector{<:Real}}}
    ti, tf, dt, xi, xf, u, cfl, lim, slp = [ simulation_params[param] for param in keys(simulation_params) ]
    dx = u * dt / cfl
    nt = trunc(Int, (tf - ti) / dt + 1)
    x_mesh = collect(xi: dx: xf)
    nx = length(x_mesh)
    return dt, dx, u, cfl, lim, slp, nt, nx, x_mesh
end


"""
Reads the simulation parameters at SIMULATION2D_PARAMETERS, extracting the necessary
parameters as well as creating the numerical mesh used for the resolution.

1. [arg] simulation_params: OrderedDict storing all the simulation parameters.

[return] Temporal/spatial step, adv speeds, CFL num, slope and limiter, t points and spatial mesh.
"""
function simulation_parameters2D(simulation_params:: OrderedDict) :: NTuple{13, Union{Real, Vector{<:Real}}}
    ti, tf, dt, xi, xf, yi, yf, u, v, cfl, lim, slp = 
    [ simulation_params[param] for param in keys(simulation_params) ]
    dx, dy = u * dt / cfl, v * dt / cfl
    nt = trunc(Int, (tf - ti) / dt + 1)
    x_mesh = collect(xi: dx: xf)
    y_mesh = collect(yi: dy: yf)
    nx, ny = length(x_mesh), length(y_mesh)
    return dt, dx, dy, u, v, cfl, lim, slp, nt, nx, ny, x_mesh, y_mesh
end



##################################
#^   COMPUTE THE REA ALGORITHM
##################################

"""
Calculates the 1D cell average of a given function with initial profile.

1. [arg] dx: Discretization of the spatial mesh.
2. [arg] func: Function of which calculate the cell average.
3. [arg] x_mesh: 1D spatial mesh over which the solution is propagated.

[return] The 1D cell average over the numerical mesh. Periodic boundary condition used.
"""
function _profile1D_cell_average(dx:: Real, func:: Function, x_mesh:: Vector{<:Real}) :: Vector{<:Real}
    integration = [ _integration_Quad1D(x - dx / 2, x + dx / 2, func) for x in x_mesh[2:end] ] ./ dx
    pushfirst!(integration, integration[end])  # Periodic boundary condition
    return integration
end


"""
Calculates the 2D cell average of a given function with initial profile.

1. [arg] dx: Discretization of the spatial x mesh.
2. [arg] dy: Discretization of the spatial y mesh.
3. [arg] func: Function of which calculate the cell average.
4. [arg] x_mesh: Spatial x mesh over which the solution is propagated.
5. [arg] y_mesh: Spatial y mesh over which the solution is propagated.

[return] The 2D cell average over the numerical mesh. Periodic boundary condition used.
"""
function _profile2D_cell_average(dx:: Real, dy:: Real, func:: Function, x_mesh:: Vector{<:Real}, 
    y_mesh:: Vector{<:Real}) :: Matrix{<:Real}
    
    integration = [ _integration_Quad2D(x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2, func)
    for y in y_mesh, x in x_mesh ] ./ (dx * dy)
    integration[:, 1] .= integration[:, end]  # Periodic boundary condition at x-axis
    integration[1, :] .= integration[end, :]  # Periodic boundary condition at y-axis
    return integration
end


"""
Calculates the cell average for the next time step using the general expression
of Godunov's REA algorithm for constant advection speed and linear piecewise 
reconstruction, given the slopes.

1. [arg] adv_speed: Eigenvalue for m = 1 A matrix of advection PDE.
2. [arg] dx: Discretization of the spatial mesh.
3. [arg] dt: Discretization of the temporal mesh.
4. [arg] CFL_number: Courant number (<= 1 for stability for three points schemes)
5. [arg] cell_averages: Q value of the actual time step.
6. [arg] slopes: Computed slopes of the desired method.

[return] The cell averages alongside a vector of the next time step.
"""
function _REA_algorithm(adv_speed:: Real, dx:: Real, dt:: Real, CFL_number:: Real, 
    cell_averages:: Vector{<:Real}, slopes:: Union{Real, Vector{<:Real}}) :: Vector{<:Real}
    
    new_cell_averages = [ (cell_averages[i] - CFL_number * (cell_averages[i] - cell_averages[i - 1]) 
    -0.5 * CFL_number * (dx - adv_speed * dt) * (slopes[i] - slopes[i - 1])) for i in 2 : length(cell_averages) ]
    pushfirst!(new_cell_averages, new_cell_averages[end])  # Periodic boundary condition
    return new_cell_averages
end



##################################
#!        GODUNOV REA
##################################

"""
Computes Godunov's REA algorithm to solve the 1D advection equation according to the
given parameters and inital profile.

1. [arg] simulation_params: Dictionary which contains the desired parameters for the resolution.
2. [arg] init_profile: Configured 1D profile which is going to be advected.

[return] The x and y spatial meshes, and the numerical solution as a Matrix were each row correspond
to a time step, while the columns to the numerical solution alongside the x-mesh.
"""
function godunov_REA1D(simulation_params:: OrderedDict, 
    init_profile:: Function) :: Tuple{Vector{<:Real}, Matrix{<:Real}}

    dt, dx, u, cfl, lim, slp, nt, nx, x_mesh = simulation_parameters1D(simulation_params)  # Load the simulation params
    numerical_solution = zeros(nt, nx)
    cell_average = _profile1D_cell_average(dx, init_profile, x_mesh)            # Compute the cell averages of the initial f(x) profile
    slope = (lim == false) & (slp == false) ? 0 : slp_lim_selection(slp, lim)   # Selects the slope or limiter to use

    if slope == 0  # Upwind scheme
        for n in 1 : nt - 1
            println("Computing the 1D advection: t = $(round(n * dt, digits = 2)) s")
            numerical_solution[n, :] = cell_average
            cell_average = _REA_algorithm(u, dx, dt, cfl, cell_average, zeros(nx))
        end
    else
        for n in 1 : nt - 1  # Other schemes which use slopes with or without TVD limiters
            println("Computing the 1D advection: t = $(round(n * dt, digits = 2)) s")
            numerical_solution[n, :] = cell_average
            slopes = slope(cell_average, dx)
            cell_average = _REA_algorithm(u, dx, dt, cfl, cell_average, slopes)
        end
    end
    return x_mesh, numerical_solution
end


"""
Computes Godunov's REA algorithm to solve the 2D advection equation using dimensional splitting,
according to the given parameters and inital profile. Additionaly, computes the analitical solution.

1. [arg] simulation_params: Dictionary which contains the desired parameters for the resolution.
2. [arg] init_profile: Configured 2D profile which is going to be advected.

[return] The x and y spatial meshes, and the numerical solution as an Array were each vertical level 
correspond to a time step, while the rows and columns to the numerical solution over the spatial mesh.
"""
function godunov_REA2D(simulation_params:: OrderedDict, 
    init_profile:: Function) :: Tuple{Vector{<:Real}, Vector{<:Real}, Array{<:Real}}

    dt, dx, dy, u, v, cfl, lim, slp, nt, nx, ny, x_mesh, y_mesh = simulation_parameters2D(simulation_params)  # Load the simulations params
    cell_averages = _profile2D_cell_average(dx, dy, init_profile, x_mesh, y_mesh)  # Compute the cell averages of initial profile f(x,y)
    slope = (lim == false) & (slp == false) ? 0 : slp_lim_selection(slp, lim)      # Selects the slope or limiter to use
    numerical_solution = zeros(ny, nx, nt - 1)

    if slope == 0  # Upwind scheme + Dimensional splitting
        for n in 1 : nt - 1
            println("Computing the 2D advection: t = $(round(n * dt, digits = 2)) s")
            numerical_solution[:, :, n] = cell_averages
            for j in 1 : ny  # Advection in x-direction
                cell_averages[j, :] = _REA_algorithm(u, dx, dt, cfl, cell_averages[j, :], zeros(nx))
            end
            for i in 1 : nx  # Advection in y-direction
                cell_averages[:, i] = _REA_algorithm(v, dy, dt, cfl, cell_averages[:, i], zeros(ny))
            end
        end
    else           # Other schemes + Dimensional splitting
        for n in 1 : nt - 1
            println("Computing the 2D advection: t = $(round(n * dt, digits = 2)) s")
            numerical_solution[:, :, n] = cell_averages
            for j in 1 : ny  # Advection in x-direction
                slopes = slope(cell_averages[j, :], dx)
                cell_averages[j, :] = _REA_algorithm(u, dx, dt, cfl, cell_averages[j, :], slopes)
            end
            for i in 1 : nx  # Advection in y-direction
                slopes = slope(cell_averages[:, i], dy)
                cell_averages[:, i] = _REA_algorithm(v, dy, dt, cfl, cell_averages[:, i], slopes)
            end
        end
    end
    return x_mesh, y_mesh, numerical_solution
end



# END OF SCRIPT ##################