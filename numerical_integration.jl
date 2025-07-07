
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the functions used for numerically integrate a given 1D function. The
functions on the script will be used to compute the cell averages for the different
FVM methods from the initial profile.

[WARNING] Although the 1D and 2D trapezoid rule are programmed, only the 1D and 2D 
quadrature methods are used as they increased simulation performance.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date 26/01/2024
Enconding: UFT-8
"""


using MultiQuad



"""
Numerically calculates 1D integral using trapezoidal rule.

1. [arg] x_init: Left limit of the interval.
2. [arg] x_final: Right limit of the interval.
3. [arg] func: Function f(x) to numerically integrate.
4. [arg] points: Number of divisions of the interval. 50 by default.

[return] Function integrated within the interval with trapezoidal rule.
"""
function _trapz_integration1D(x_init:: Real, x_final:: Real, func:: Function, 
    points:: Real = 50) :: Real
    
    h = (x_final - x_init) / points
    integration = h * (func(x_init) + func(x_final)) / 2
    for nx = 1 : points - 1
        midpoint = (x_final - x_init) * nx / points + x_init
        integration += h * func(midpoint)
    end
    return integration
end


"""
Numerically calculates the 1D integral using quadrature methods.

1. [arg] x_init: Left limit of the interval.
2. [arg] x_final: Right limit of the interval.
3. [arg] func: Function f(x) to numerically integrate.

[return] Function f(x) integrated within the interval with quadrature method.
"""
function _integration_Quad1D(x_init:: Real, x_final:: Real, func:: Function) :: Real
    return quad(func, x_init, x_final)[1]  # The second value [2] corresponds to the num. error
end


"""
Numerically calculates 2D integral using trapezoidal rule.

1. [arg] x_init: Left limit of the x interval.
2. [arg] x_final: Right limit of the x interval.
3. [arg] y_init: Left limit of the y interval.
4. [arg] y_final: Right limit of the y interval.
5. [arg] func: Function f(x, y) to numerically integrate.
6. [arg] points: Number of intervals used for trapezoidal integration. 25 by default.

[return] Function integrated within the 2D interval with trapezoidal rule.
"""
function _trapz_integration2D(x_init::Real, x_final::Real, y_init::Real, y_final:: Real,
    func::Function, points::Real = 25) :: Real

    hx = (x_final - x_init) / points
    hy = (y_final - y_init) / points
    x_midpoints = LinRange(x_init + hx / 2, x_final - hx / 2, points - 1)
    y_midpoints = LinRange(y_init + hy / 2, y_final - hy / 2, points - 1)
    integration = 0.25 * hx * hy * (func(x_init, y_init) + func(x_final, y_init) 
            + func(x_init, y_final) + func(x_final, y_final))
    integration += hx * hy * sum(func(x, y) for x in x_midpoints, y in y_midpoints)
    return integration
end


"""
Numerically calculates the 2D integral using quadrature methods.

1. [arg] x_init: Left limit of the x interval.
2. [arg] x_final: Right limit of the x interval.
3. [arg] y_init: Left limit of the y interval.
4. [arg] y_final: Right limit of the y interval.
5. [arg] func: Function f(x, y) to numerically integrate.

[return] Function f(x, y) integrated within the interval with quadrature method.
"""
function _integration_Quad2D(x_init:: Real, x_final:: Real, y_init:: Real, y_final:: Real, func:: Function) :: Real
    return dblquad(func, x_init, x_final, x -> y_init, x -> y_final)[1]  # The second value [2] corresponds to the num. error
end



# END OF SCRIPT ##################