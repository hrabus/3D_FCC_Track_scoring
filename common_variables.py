# python version of track structure analysis from Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)
# This module contains the common variables used in the differen codes

from numpy import zeros, ones, sqrt

#     ---- Debugging options -----
DEBUG = zeros(5, 'b')  # '[False for i in range(5)]
DEBUG[0] = True  # Invoke message on Version date from Subroutine CLUSTR
DEBUG[1] = False  # Main program
DEBUG[2] = False  # CLUSTR main sections
DEBUG[3] = False  # CLUSTR IDIR loop
DEBUG[4] = False  # CLUSTR IPOS loop details
SRFLAG = True

#  ----- Site size
diameter_of_site = 0.                               # Diameter of spherical site in nm

#  ----- Parameters for lattice orientation
direction_of_motion = ones(3) / sqrt(3.)
n_azimuthal_rotations_of_track = 1                  # Max. number of azimuthal ratoations of track (phi1)
n_directions = 1
Reciprocal_base = zeros((3, 3, n_directions), 'f')
Real_space_base = zeros((3, 3, n_directions), 'f')

#
minimum_number_of_ionizations = 0
minimum_energy_imparted = 0.
KMAX = 4   #
#    ----- Parameters for track and radial distance histograms
max_number_of_energy_transfer_points = 100000
xyz = zeros((max_number_of_energy_transfer_points, 3), 'f')
site_xyz = zeros((max_number_of_energy_transfer_points, 3), 'f')
edep = zeros(max_number_of_energy_transfer_points, 'f')
ionizations = ones(max_number_of_energy_transfer_points, dtype='i')
edep_flag = False
type_flag = False
tpos_flag = False

interactions = zeros(max_number_of_energy_transfer_points, dtype='i')

# number_of_sites = 0
