
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the function with the general solution of the 1D and 2D advection equation.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""


using OrderedCollections



##################################
#@     ANALYTICAL SOLUTIONS
##################################

"""
Computes the analytical solution of the 1D advection equation given the parameters of 
the advection and numerical mesh.

1. [arg] x_mesh:
2. [arg] simulation_params:
3. [arg] init_profile:

[return] A Matrix with the analytical solution of the 1D Advection equation with periodic BC.
"""
function advection1D_solution(x_mesh:: Vector{<:Real}, simulation_params:: OrderedDict, 
    init_profile:: Function) :: Matrix{<:Real} 

    ti, tf, dt, u = [ simulation_params[param] for param in ["t_init", "t_final", "t_step", "adv_speed"] ]
    nt = trunc(Int64, (tf - ti) / dt + 1)
    nx = trunc(Int64, length(x_mesh))
    profile = init_profile(x_mesh)
    analytical_solution = zeros(nt, nx)

    for n in 1 : nt
        analytical_solution[n, :] = profile
        profile = init_profile(x_mesh[2:end] .- u * n * dt)  #! Fix periodic boundary condition.
        pushfirst!(profile, profile[end])   #! Boundary condition does not advect!
    end
    return analytical_solution
end



#! The 2D analytical solution was going to be done at the beggining, but was later discarded



# END OF SCRIPT ##################