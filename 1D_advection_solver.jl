
"""
Julia script for the Finite Volume Methods subject at UIB's FAMA master degree.

Main script to simulate the 1D advection equation using Godunov's REA algorithm.
To configure the simulation, the only parts which must be changed are the parameters
inside the ordered dictionaries. Then, simply run the script.

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
include("1D_advection_plots.jl")
include("1D_initial_profiles.jl")
include("analytical_solutions.jl")

FOLDER_PATH = "Advection_Results1D"  # Folder were the results will be saved
SIMULATION_NAME = "animation"        # Desired GIF name. "animation" by default



##################################
#$   CONFIGURE ADV1D SIMULATION
##################################

"""
Variables used to define the 1D initial profile. 
For superposition profiles use Vectors. p.e: [1, 4, ...]

    Profiles -----
    (1) _gaussian_profile1D
    (2) _rectangular_profile1D
    (3) _gaussian_superposition_profile1D
    (4) _rectangular_superposition_profile1D
    (5) _gaussian_rectangular_profile1D
    (6) _wave_packet1D
"""
PROFILE1D_PARAMETERS = OrderedDict(
    "Profile" => 6,           # Chosen profile
    "x0" => 15,               # Centered position of the profile/s
    "sigma/width" => 5,       # Std/width of profile/s
    "amplitude/height" => 5,  # Amplitude/height of the profiles
    "frequency" => 5,         # Frequency for wavepacket profile
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
SIMULATION1D_PARAMETERS = OrderedDict(
    "t_init" => 0,         # Starting time (s)
    "t_final" => 100,      # Final time (s)
    "t_step" => 0.1,       # Time step (s)
    "x_init" => 0,         # Starting point (m)
    "x_final" => 50,       # Final point (m)
    "adv_speed" => 0.75,   # Speed of the advection (m/s)
    "CFL_number" => 0.7,   # Courant number (a.u)
    "Limiter" => false,    # TVD limiter over slope.
    "Slope" => false       # Slope to use
)



##################################
#!  RUN 1D ADVECTION SIMULATION 
##################################

"""
MAIN FUNCTION FOR 1D CASE

Numerically solves using Godunov's REA method, the 1D advection PDE according to the given 
initial profile and its parameters, as well as the simulation parameters. Then displays and
saves the result as an animated GIF at the desired folder with the chosen name.
"""
function simulate_advection1D()

    profile(x) = profile1D_configured(        # Configure and load the 1D profile
        x, 
        PROFILE1D_PARAMETERS
        ) 
    x_mesh, adv1D_numerical = godunov_REA1D(    # Numerically solves the adv1D equation
        SIMULATION1D_PARAMETERS, 
        profile
        )
    adv1D_analytical = advection1D_solution(  # Analytical solution of adv1D equation
        x_mesh, 
        SIMULATION1D_PARAMETERS, 
        profile
        )
    animated_advection1D(                     # Numerical vs analytical animated result
        x_mesh, 
        SIMULATION1D_PARAMETERS, 
        adv1D_numerical, 
        adv1D_analytical, 
        FOLDER_PATH, 
        SIMULATION_NAME
    )
end


simulate_advection1D()  # RUNS THE SIMULATION, GIF OF THE SOLUTION AND SAVES IT



# END OF SCRIPT ##################