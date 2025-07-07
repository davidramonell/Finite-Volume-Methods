
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the functions to calculate different possible slopes to improve the
precission order of Godunov's REA algorithm for hyperbolic PDEs, given the cell
average over the domain and the spatial discretization (constant).

Additionally contains the functions to calculate several types of slope limitators 
which verify TDV condition for the second order approximation.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""


##################################
#@            SLOPES
##################################

"""
Calculates the Upwind (backward differences) slopes for Godunov's REA algorithm. Uses periodic
boundary conditions for the first slope.

1. [arg] cell_averages: Vector containing the numerical cell average over all spatial mesh.
2. [arg] dx: Discretization of the spatial mesh.

[return] The value of the suitable slopes according to Upwind construction.
"""
function _slope_upwind(cell_averages:: Vector{<:Real}, dx:: Real) :: Vector{<:Real}
    slopes = [ ((cell_averages[i] - cell_averages[i - 1]) / dx) for i in 2 : length(cell_averages) ]
    pushfirst!(slopes, slopes[end])  # Periodic boundary condition for the first slope
    return slopes
end


"""
Calculates the Downwind (forward differences) slopes for Godunov's REA algorithm. Uses periodic
boundary conditions for the last slope.

1. [arg] cell_averages: Vector containing the numerical cell average over all spatial mesh.
2. [arg] dx: Discretization of the spatial mesh.

[return] The value of the suitable slopes according to Downwind construction.
"""
function _slope_downwind(cell_averages:: Vector{<:Real}, dx:: Real) :: Vector{<:Real}
    slopes = [ ((cell_averages[i + 1] - cell_averages[i]) / dx) for i in 1 : length(cell_averages) - 1 ]
    append!(slopes, slopes[1])  # Periodic boundary condition
    return slopes
end


"""
Calculates the upwind (backward differences) slopes for Godunov's REA algorithm. Uses periodic
boundary conditions for the first and last slope.

1. [arg] cell_averages: Vector containing the numerical cell average over all spatial mesh.
2. [arg] dx: Discretization of the spatial mesh.

[return] The value of the suitable slopes according to centered construction.
"""
function _slope_centered(cell_averages:: Vector{<:Real}, dx:: Real) :: Vector{<:Real}
    slopes = [ ((cell_averages[i + 1] - cell_averages[i - 1]) / (2 * dx)) for i in 2 : length(cell_averages) - 1 ]
    append!(slopes, cell_averages[end] - cell_averages[end - 1])  # Downwind slope for final grid point
    pushfirst!(slopes, slopes[end])  # Periodic boundary condition for the first slope
end



##################################
#$        MINMOD & MAXMOD
##################################

"""
Finds the minimum module of two given variables.

1. [arg] a: First variable.
2. [arg] b: Second variable.

[return] The minimum module of two variables.
"""
function _minmod(a:: Real, b::Real) :: Real
    return a * b <= 0 ? 0 : sign(a) * min(abs(a), abs(b))
end


"""
Finds the maximum module of two given variables.

1. [arg] a: First variable.
2. [arg] b: Second variable.

[return] The maximum module of two variables.
"""
function _maxmod(a:: Real, b::Real) :: Real
    return a * b <= 0 ? 0 : sign(a) * max(abs(a), abs(b))
end



##################################
#^           LIMITERS
##################################

"""
Applies the basic minmod slope limiter for Godunov's REA algorithm.

1. [arg] cell_averages: Vector containing the numerical cell average over all spatial mesh.
2. [arg] dx: Discretization of the spatial mesh.

[return] The modified slopes using minmod slope limiter.
"""
function _limiter_minmod(cell_averages:: Vector{<:Real}, dx:: Real) :: Vector{<:Real}
    return _minmod.(_slope_upwind(cell_averages, dx), _slope_downwind(cell_averages, dx))
end


"""
Applies the Superbee slope limiter for Godunov's REA algorithm.

1. [arg] cell_averages: Vector containing the numerical cell average over all spatial mesh.
2. [arg] dx: Discretization of the spatial mesh.

[return] The modified slopes using Superbee slope limiter.
"""
function _limiter_superbee(cell_averages:: Vector{<:Real}, dx:: Real) :: Vector{<:Real}
    slope_upwind = _slope_upwind(cell_averages, dx)
    slope_downwind = _slope_downwind(cell_averages, dx)
    slope1 = _minmod.(slope_downwind, 2 .* slope_upwind)
    slope2 = _minmod.(2 .* slope_downwind, slope_upwind)
    return _maxmod.(slope1, slope2)
end


"""
Applies the monitorized center slope limiter for Godunov's REA algorithm.

1. [arg] cell_averages: Vector containing the numerical cell average over all spatial mesh.
2. [arg] dx: Discretization of the spatial mesh.

[return] The modified slopes using MC slope limiter.
"""
function _limiter_mc(cell_averages:: Vector{<:Real}, dx:: Real) :: Vector{<:Real}
    first_minmod = _minmod.(2 .* _slope_upwind(cell_averages, dx), 2 .* _slope_downwind(cell_averages, dx))
    return _minmod.(_slope_centered(cell_averages, dx), first_minmod)
end



##################################
#!       SLOPES & LIMITERS 
##################################

# Outputs of the script which should be included at the script for numerical resolution

"""
Vector with all defined ways to compute the slopes.

    (1) _slope_upwind(cell_averages, dx)
    (2) _slope_downwind(cell_averages, dx)
    (3) _slope_centered(cell_averages, dx)
"""
SLOPES = [
    _slope_upwind,    # Upwind construction   -> Beam-Warming scheme
    _slope_downwind,  # Downwind construction -> Lax-Wendroff scheme
    _slope_centered   # Centered construction -> Fromm scheme
]


"""
Vector with all defined ways to apply limiters on the slopes.

    (1) _limiter_minmod(cell_averages, dx)
    (2) _limiter_superbee(cell_averages, dx)
    (3) _limiter_mc(cell_averages, dx)
"""
LIMITERS = [
    _limiter_minmod,    # Minmod limiter
    _limiter_superbee,  # Superbee limiter
    _limiter_mc         # Monitored center limiter
]


"""
Loads the chosen slope or limiter for computing Godunov's REA algorithm. If no
arguments are passed to the function, it returns 0 (Upwind scheme). Selecting
a limiter will override the chosen slope.

1. [arg] num_slope: Position at SLOPES of the desired slope. False by default.
2. [arg] num_limiter: Position at LIMITERS of the desired limiter. False by default.

[return] Zero (Upwind) by default. The chosen slope or limiter elsewise.
"""
function slp_lim_selection(num_slope:: Union{Int, Bool} = false, 
    num_limiter:: Union{Int, Bool} = false) :: Union{Real, Function}

    if (num_limiter == false) && (num_slope == false)
        return 0                      # Upwind scheme
    elseif (num_slope in 1:length(SLOPES)) && (num_limiter == false)
        return SLOPES[num_slope]      # Upwind, Downwind or centered
    else
        return LIMITERS[num_limiter]  # Minmod, Superbee, MC
    end
end



# END OF SCRIPT ==================