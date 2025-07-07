
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Main script to simulate the 2D advection equation using Godunov's REA method with 
dimensional splitting. To configure the simulation, the only parts which must be changed
are the parameters inside the ordered dictionaries. Then simply run the script.

[WARNING] The first average calculation over each cell for the given profile f(x,y) may take
from a few seconds to a minute.

Author: David Ramonell Pujol
Julia Version: 1.9.3
Date: 26/01/2024
Enconding: UFT-8
"""


##################################
#@  PACKAGES, CALLS & VARIABLES
##################################

using OrderedCollections

include("godunov_REA.jl")
include("2D_advection_plots.jl")
include("2D_initial_profiles.jl")
include("analytical_solutions.jl")

FOLDER_PATH = "Advection_Results2D"  # Folder were the results will be saved
SIMULATION_NAME = "animation"        # Desired GIF name. "animation" by default.



##################################
#$  CONFIGURE ADV2D SIMULATION
##################################

"""
Variables used to define the 2D initial profile. For individual profiles
now the given parameters must be vectors. p.e: [1, 3, ...].

    Profiles -----
    (1) _gaussian_profile2D
    (2) _rectangular_profile2D
"""
PROFILE2D_PARAMETERS = OrderedDict(
    "Profile" => 2,             # Chosen profile
    "q0" => [0, 0],             # Centered position of the profile
    "sigma/width" => [1, 1],    # Std/width of profile
    "amplitude/height" => 4     # Amplitude/height of the profile
)


"""
Variables used to define the the simulation. CFL must be under 1 for stability.
Use "false" at limiter and slope to specify if they'll be used or not. If both 
are false, the simulation will result in Upwind scheme. If a limiter is selected, 
it will override any parameter on "Slopes".

    Limiter -----
    (1) minmod
    (2) superbee
    (3) monitored centering
    Slopes -----
    (1) upwind    -> Beam-Warming
    (2) downwind  -> Lax-Wendroff
    (3) centered  -> Fromm
"""
SIMULATION2D_PARAMETERS = OrderedDict(
    "t_init" => 0,          # Starting time (s)
    "t_final" => 100,       # Final time (s)
    "t_step" => 0.1,        # Time step (s)
    "x_init" => -10,        # Starting x-point (m)
    "x_final" => 10,        # Final x-point (m)
    "y_init" => -10,        # Starting y-point (m)
    "y_final" => 10,        # Final y-point (m)
    "adv_speed_x" => 0.2,   # Speed of advection at x-axis (m/s)
    "adv_speed_y" => 0.5,   # Speed of advection at y-axis (m/s)
    "CFL_number" => 0.7,    # Courant number (a.u)
    "Limiter" => false,     # TVD limiter over slope
    "Slope" => false        # Slope to use
)



##################################
#!  RUN 2D ADVECTION SIMULATION 
##################################

"""
MAIN FUNCTION FOR 2D CASE

Numerically solves using Godunov's REA method and dimensional splitting, the 2D advection PDE 
according to the given initial profile and its parameters, as well as the simulation parameters. 
Then displays and saves the result as an animated GIF at the desired folder with the chosen name.
"""
function simulate_advection2D()

    profile(x, y) = profile2D_configured(               # Configure and load the 2D profile
        x, 
        y, 
        PROFILE2D_PARAMETERS
        ) 
    x_mesh, y_mesh, adv2D_numerical = godunov_REA2D(    # Numerically solves the adv2D equation
        SIMULATION2D_PARAMETERS, 
        profile
        )
    animated_advection2D(                               # Animated result of the simulation
        x_mesh, 
        y_mesh, 
        SIMULATION2D_PARAMETERS, 
        adv2D_numerical, 
        FOLDER_PATH, 
        SIMULATION_NAME
    )
end


simulate_advection2D()  # RUNS THE SIMULATION, GIF OF THE SOLUTION AND SAVES IT



# END OF SCRIPT ##################