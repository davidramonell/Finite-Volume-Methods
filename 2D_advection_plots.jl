
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Contains the functions for plotting, animating and saving the results of solving the 2D advection 
PDE using Godunov's REA algorithm and dimensional splitting.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""


##################################
#@      PACKAGES & TITLES
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
#$       INDIVIDUAL HEATMAP
##################################

"""
Plots a specified time step of the numerical resolution of the Advection2D PDE.

1. [arg] x_mesh: Numerical mesh for the x variable.
2. [arg] y_mesh: Numerical mesh for the y variable.
3. [arg] simulation_params: Dictionary with the parameters for 1D simulation.
4. [arg] adv1D_numerical: Matrix with the numerical solution of 1D Advection.
5. [arg] frame: Time step of solution.

[return] The numerical and analytical solution plotted for a chosen time step.
"""
function _heatmap_advection2D(x_mesh:: Vector{<:Real}, y_mesh:: Vector{<:Real}, 
    simulation_params:: OrderedDict, adv2D_numerical:: Array{<:Real}, frame:: Real)
    
    dt, limiter, slope = [ simulation_params[param] for param in ["t_step", "Limiter", "Slope"] ]
    time = round(frame * dt, digits = 2)

    if (slope == false) && (limiter == false)
        title_graph = L"\textbf{Upwind~scheme~:~No~slope~nor~limiter}"
    elseif (slope != false) && (limiter == false)
        title_graph = TITLE_SCHEMES[slope]
    else
        title_graph = TITLE_LIMITERS[limiter]
    end

    heatmap(
        x_mesh, 
        y_mesh, 
        adv2D_numerical[:, :, frame], 
        color = :turbo, 
        dpi = 100
    )
    heatmap!(
        xlabel = L"x", 
        ylabel = L"y",
        # zlim = (0, 5), 
        title = "$(title_graph)  t = $time s", 
        tickfont = font(11, "Computer Modern"), 
        titlefont = font(11, "Computer Modern"),  
        size = (700, 500), 
        dpi = 100, 
        clims = (0, 4)    # Change the limits of colormap according to the amplitude of profile
    )
end



##################################
#^    ANIMATED GIF & SAVING
##################################

# Output of the script which should be included at the main script for numerical resolution

"""
Plots the numerical resolution of the 2D Advection PDE as an animated GIF, then saves it
into the given folder with the desired name.

1. [arg] x_mesh: Numerical mesh for the x variable.
2. [arg] y_mesh: Numerical mesh for the y variable.
3. [arg] simulation_params: Dictionary with the parameters for 2D simulation.
4. [arg] adv2D_numerical: Matrix with the numerical solution of 2D Advection.
5. [arg] folder_path: Path to the folder were results will be saved.
6. [arg] name: Desired name for the saved animation.

[return] Animated GIF displaying the results of the numerical resolution, saved at desired folder.
"""
function animated_advection2D(x_mesh:: Vector{<:Real}, y_mesh:: Vector{<:Real}, 
    simulation_params:: OrderedDict, adv2D_numerical:: Array{<:Real}, folder_path:: String, 
    name:: String = "animation")

    animation = @animate for frame in 1 : size(adv2D_numerical, 3)
        _heatmap_advection2D(x_mesh, y_mesh, simulation_params, adv2D_numerical, frame)
    end every 5
    gif(animation, joinpath(folder_path, "$(name).gif"), fps = 20)
end



# END OF SCRIPT ##################