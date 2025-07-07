
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the functions for plotting, animating and saving the results of solving the 1D advection 
PDE using Godunov's REA algorithm.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""


##################################
#@       PACKAGES & TITLES
##################################

using Plots
using ImageMagick
using LaTeXStrings
using OrderedCollections


"""
Name of each possible scheme used at the project to 
compute the slopes for linear piecewise reconstruction.
"""
TITLE_SCHEMES = [
    L"\textbf{Beam-Warming~scheme~:~Upwind~slope}", 
    L"\textbf{Lax-Wendroff~scheme~:~Downwind~slope}",
    L"\textbf{Fromm~scheme~:~Centered~slope}"
]


"""
Name of each possible TVD method used at the project to
compute the limiters over the slopes.
"""
TITLE_LIMITERS = [
    L"\textbf{Minmod~limiter}", 
    L"\textbf{Superbee~limiter}", 
    L"\textbf{Monitored~center~limiter}"
]



##################################
#$       INDIVIDUAL PLOT
##################################

"""
Plots a specified time step of the numerical resolution of Advection 1D PDE.

1. [arg] x_mesh: Numerical mesh for the x variable.
2. [arg] simulation_params: Dictionary with the parameters for 1D simulation.
3. [arg] adv1D_numerical: Matrix with the numerical solution of 1D Advection.
4. [arg] adv1D_analytical: Matrix with the analytical solution of 1D Advection.
5. [arg] frame: Time step of solution.

[return] The numerical and analytical solution plotted for a chosen time step.
"""
function _plot_advection1D(x_mesh:: Vector{<:Real}, simulation_params:: OrderedDict, 
    adv1D_numerical:: Matrix{<:Real}, adv1D_analytical:: Matrix{<:Real}, frame:: Real)

    dt, limiter, slope = [ simulation_params[param] for param in ["t_step", "Limiter", "Slope"] ]
    time = round(frame * dt, digits = 2)

    if (slope == false) && (limiter == false)
        title_graph = L"\textbf{Upwind~scheme~:~No~slope~nor~limiter}"
    elseif (slope != false) && (limiter == false)
        title_graph = TITLE_SCHEMES[slope]
    else
        title_graph = TITLE_LIMITERS[limiter]
    end

    plot(
        x_mesh, 
        adv1D_analytical[frame, :], 
        label = L"Analitical", 
        color = :blue, 
        dpi = 100
    )
    plot!(
        x_mesh,
        adv1D_numerical[frame, :], 
        label = L"Numerical",
        color = :red, 
        dpi = 100
    )
    plot!(
        legendcolumns = 2, 
        legend = :outerbottom, 
        xlabel = L"x", 
        ylabel = L"q(x,t)", 
        ylim = (-6, 6),                          # Change the limits of plot according to the amplitudes.
        title = "$(title_graph)  t = $time s", 
        titlefont = font(11, "Computer Modern"), 
        tickfont = font(11, "Computer Modern"), 
        size = (700, 500), 
        dpi = 100
    )
end



##################################
#^    ANIMATED GIF & SAVING
##################################

# Output of the script which should be included at the main script for numerical resolution

"""
Plots the numerical resolution of the 1D Advection PDE as an animated GIF, then saves it
into the given folder with the desired name.

1. [arg] x_mesh: Numerical mesh for the x variable.
2. [arg] simulation_params: Dictionary with the parameters for 1D simulation.
3. [arg] adv1D_numerical: Matrix with the numerical solution of 1D Advection.
4. [arg] adv1D_analytical: Matrix with the analytical solution of 1D Advection.
5. [arg] folder_path: Path to the folder were results will be saved.
6. [arg] name: Desired name for the saved animation.

[return] Animated GIF displaying the results of the numerical resolution, saved at desired folder.
"""
function animated_advection1D(x_mesh:: Vector{<:Real}, simulation_params:: OrderedDict, 
    adv1D_numerical:: Matrix{<:Real}, adv1D_analytical:: Matrix{<:Real}, folder_path:: String, 
    name:: String = "animation")

    animation = @animate for frame in 1 : size(adv1D_numerical, 1)
        _plot_advection1D(x_mesh, simulation_params, adv1D_numerical, adv1D_analytical, frame)
    end every 5
    gif(animation, joinpath(folder_path, "$(name).gif"), fps = 20)
end



# END OF SCRIPT ##################