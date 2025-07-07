
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the 2D initial profiles which will be used on the 2D resolution of the
advection equation PDE using Godunov's REA algorithm + dimensional splitting.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""

using Plots
using OrderedCollections



##################################
#@         2D PROFILES            
##################################

"""
A single 2D gaussian profile.

1. [arg] q0: Centered position in x and y of the profile.
2. [arg] x: Individual spatial point or spatial mesh of x.
3. [arg] y: Individual spatial point or spatial mesh of y.
4. [arg] sigmas: Standard desviation at x and y of the 2D gaussian.
5. [arg] height: Amplitude of the gaussian profile.

[return] A 2D Gaussian profile centered at (x0, y0), with a given std and amplitude.
"""
function _gaussian_profile2D(q0:: Vector{<:Real}, x:: Union{Real, Vector{<:Real}}, y:: Union{Real, Vector{<:Real}}, 
    sigmas:: Vector{<:Real}, height:: Real) :: Union{Real, Matrix{<:Real}}

    if isa(x, Vector) || isa(y, Vector)
        x, y = [ x for x in x, y in y ], [ y for x in x, y in y ]
    end
    return @. height * exp( -(x - q0[1])^2 / (2 *sigmas[1]) -(y - q0[2])^2 / (2 *sigmas[2]) )
end


"""
A single 2D rectangular profile

1. [arg] q0: Centered position in x and y of the profile.
2. [arg] x: Individual spatial point or spatial mesh of x.
3. [arg] y: Individual spatial point or spatial mesh of y.
4. [arg] widths: Widths of each side of the rectangle respect the center.
5. [arg] height: Amplitude of the rectangular profile.

[return] A 2D rectangular profile centered at (x0, y0), with given widths and amplitude.
"""
function _rectangular_profile2D(q0:: Vector{<:Real}, x:: Union{Real, Vector{<:Real}}, y:: Union{Real, Vector{<:Real}}, 
    widths:: Vector{<:Real}, height:: Real) :: Union{Real, Matrix{<:Real}}
    
    if isa(x, Vector) || isa(y, Vector)
        x, y = [ x for x in x, y in y ], [ y for x in x, y in y ]
    end
    return @. height * ((q0[1] - widths[1] / 2 <= x <= q0[1] + widths[1] / 2) * 
    (q0[2] - widths[2] / 2 <= y <= q0[2] + widths[2] / 2))
end



##################################
#$     CONFIGURATE PROFILES
##################################

"""
Vector with all possible 2D profiles.

    (1) _gaussian_profile2D(q0, x, y, sigmas, height)
    (2) _rectangular_profile2D(q0, x, y, widths, height)
"""
PROFILES2D = [
    _gaussian_profile2D,       # Simple 2D gaussian profile
    _rectangular_profile2D     # Simple 2D rectangular profile
]


"""
Configures the chosen 2D profile using the written parameters at PROFILE2D_PARAMETERS.
The use of this function has to follow an structure such us "" fun(x,y) = profile2D_configured(n,x,y)"".

1. [arg] x: Individual spatial point or spatial x-mesh.
2. [arg] y: Individual spatial point or spatial y-mesh.
3. [arg] profiles_params: OrderedDict with all the parameters.

[return] The profile avaluated at the given parameters and as a function of the x and y variables.
"""
function profile2D_configured(x:: Union{Real, Vector{<:Real}}, y:: Union{Real, Vector{<:Real}}, 
    profile_params:: OrderedDict) :: Union{Function, Real, Matrix{<:Real}}

    num_profile, q0, std_width, ampl_height = [ profile_params[param] for param in keys(profile_params) ]
    return PROFILES2D[num_profile](q0, x, y, std_width, ampl_height)
end



# END OF SCRIPT ##################