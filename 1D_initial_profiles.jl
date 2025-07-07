
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the 1D initial profiles which will be used on the 1D resolution of the
advection equation PDE using Godunov's REA method.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""

using OrderedCollections



##################################
#@         1D PROFILES            
##################################

"""
A single 1D gaussian profile.

1. [arg] x0: Centered position of the profile.
2. [arg] x: Individual spatial point or spatial mesh.
3. [arg] sigma: Standard desviation of the gaussian.
4. [arg] height: Amplitude of the gaussian profile.

[return] A 1D Gaussian profile centered at x0, with a given std and amplitude.
"""
function _gaussian_profile1D(x0:: Real, x:: Union{Real, Vector{<:Real}}, sigma:: Real, 
    height:: Real) :: Union{Real, Vector{<:Real}}

    return @. height * exp(-(x - x0)^2 / (2 * sigma))
end


"""
A 1D profile originated from superposition of gaussians.

1. [arg] x0s: Center position of all gaussians.
2. [arg] x: Individual spatial point or spatial mesh.
3. [arg] sigmas: Standard desviation of all gaussians.
4. [arg] heights: Amplitude of all gaussians.

[return] A 1D superposition of gaussian profiles.
"""
function _gaussian_superposition_profile1D(x0s:: Vector{<:Real}, x:: Union{Real, Vector{<:Real}}, 
    sigmas:: Vector{<:Real}, heights:: Vector{<:Real}) :: Union{Real, Vector{<:Real}}

    gaussians = [ _gaussian_profile1D(x0, x, sigma, height) for (x0, sigma, height) 
                in zip(x0s, sigmas, heights) ]
    return reduce(+, gaussians)
end


"""
A simple 1D wave packet profile.

1. [arg] x0: Centered position of the profile.
2. [arg] x: Individual spatial point or spatial mesh.
3. [arg] sigma: Standard desviation of the wave packet.
4. [arg] height: Max amplitude of the wave packet.
5. [arg] frequency: Frequency of the oscillations. 10 by default

[return] A 1D wave packet profile centered at x0, with a given std, amplitude and frequency.
"""
function _wave_packet1D(x0:: Real, x:: Union{Real, Vector{<:Real}}, sigma:: Real, 
    height:: Real, frequency:: Real = 10) :: Union{Real, Vector{<:Real}}
    
    return @. height * exp(-(x - x0)^2 / (2 * sigma)) * cos(frequency * x)
end


"""
A simple 1D rectangular profile.

1. [arg] x0: Centered position of the profile.
2. [arg] x: Individual spatial point of spatial mesh.
3. [arg] width: Width of the rectangular profile.
4. [arg] height: Height of the rectangular profile.

[return] Rectangular profile centered at x0 with a certain height and width.
"""
function _rectangular_profile1D(x0::Real, x::Union{Real, Vector{<:Real}}, width::Real, 
    height::Real) :: Union{Real, Vector{<:Real}}

    return height * ((x .>= x0 - width / 2) .& (x .<= x0 + width / 2))
end


"""
A 1D profile originated from superposition of rectangles.

1. [arg] x0s: Center position of all rectangles.
2. [arg] x: Individual spatial point or spatial mesh.
3. [arg] widths: Widths of all rectangles.
4. [arg] heights: Heights of all rectangles.

[return] A 1D superposition of rectangular profiles.
"""
function _rectangular_superposition_profile1D(x0s:: Vector{<:Real}, x:: Union{Real, Vector{<:Real}}, 
    widths:: Vector{<:Real}, heights:: Vector{<:Real}) :: Union{Real, Vector{<:Real}}
    
    rectangles = [ _rectangular_profile1D(x0, x, width, height) for (x0, width, height) 
                in zip(x0s, widths, heights) ]
    
    return reduce(+, rectangles)
end


"""
A 1D profile which mixes a gaussian and rectangular function.

1. [arg] x0s: Center position of gaussian and rectangle.
2. [arg] x: Individual spatial point or spatial mesh.
3. [arg] widths: Widths of gaussian and rectangle.
4. [arg] heights: Heights of guassian and rectangle.

[return] A 1D profile with a gaussian and a rectangle.
"""
function _gaussian_rectangular_profile1D(x0s:: Vector{<:Real}, x:: Union{Real, Vector{<:Real}}, 
    widths:: Vector{<:Real}, heights:: Vector{<:Real}) :: Union{Real, Vector{<:Real}}
    
    return _gaussian_profile1D(x0s[1], x, widths[1], heights[1]) .+ 
            _rectangular_profile1D(x0s[2], x, widths[2], heights[2])
end



##################################
#$     CONFIGURATE PROFILES
##################################

# Outputs of the script which should be included at the script for numerical resolution

"""
Vector with all possible 1D profiles.

    (1) _gaussian_profile1D(x0, x, sigma, height)
    (2) _rectangular_profile1D(x0, x, width, height)
    (3) _wave_packet1D(x0, x, sigma, height, frequency)
    (4) _gaussian_superposition_profile1D(Dict{Positions, sigmas, heights}, x)
    (5) _rectangular_superposition_profile1D(Dict{Position, widths, heights}, x)
    (6) _gaussian_rectangular_profile1D(Dict{Position, widths, heights},x)

For the 3rd, 4th and 5th profiles, the parameters must be vectors.
"""
PROFILES1D = [
    _gaussian_profile1D,                   # Simple gaussian profile
    _rectangular_profile1D,                # Simple rectangular profile
    _gaussian_superposition_profile1D,     # Superposition of gaussians profile
    _rectangular_superposition_profile1D,  # Superposition of rectangles profile
    _gaussian_rectangular_profile1D,       # Rectangle followed by gaussian profile
    _wave_packet1D,                        # Wave packet profile
]


"""
Configures the chosen 1D profile using the written parameters at PROFILE1D_PARAMETERS.
The use of this function has to follow an structure such us "" fun(x) = profile1D_configured(n,x)"".

1. [arg] x: Individual spatial point or spatial mesh.
2. [arg] profiles_params: OrderedDict with all the parameters. PROFILE1D_PARAMETERS by default.

[return] The profile avaluated at the given parameters and as a function of the x variable.
"""
function profile1D_configured(x:: Union{Real, Vector{<:Real}}, profile_params:: OrderedDict) :: Union{Function, Real, Vector{<:Real}}
    num_profile, x0, std_width, ampl_height, freq = [ profile_params[param] for param in keys(profile_params) ]
    if num_profile in 1:5
        return PROFILES1D[num_profile](x0, x, std_width, ampl_height)
    else num_profile == 6
        return PROFILES1D[num_profile](x0, x, std_width, ampl_height, freq)
    end
end



# END OF SCRIPT ##################