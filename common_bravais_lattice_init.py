# python version of track structure analysis from Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)

from math import sqrt, acos, cos, sin, pi, exp, log
from numpy import matmul, zeros, transpose
from numpy.linalg import norm

# import common variables
import common_variables as cv

# ######################################################################################################################
Version_Date = '23-JUL-2022'  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ##############^^^^^^^^^^^#############################################################################################
# Version history:
# 23-JUL-2022: Accomplished transformation from Fortran code BLINIT as published in Ngcezu S and Rabus H, Radiat Environ
#              Biophys 60, 559-578 (2021)


#    ----- Parameters for lattice orientation
# When the ion trajectory passes the origin of the lattice, three Euler rotations are sufficient for defining
# relative position and alignment of track and lattice. The first (azimuthal) Euler rotation accounts for the
# rotational symmetry of the track structure around the primary particle trajectory, i.e. the positions of interactions
# are rotated around the straight line representing the trajectory. In principle, this rotation can be omitted as the
# rotational symmetry of track structure is already accounted for by scoring a large number of primary particle tracks.
# It may, however, be beneficial to include this azimuthal rotation to improve scoring statistics.
# The second and third rotation serve to establish the orientation in space of the primary particle trajectory with
# respect to the chosen coordinate frame. In the simulation setup, the trajectory of the primary particle is a natural
# choice for one of the cartesian coordinate axes, and we shall assume that it would be the z axis of this cartesian
# coordinate system associated with the particle track. In order to assess the targets receiving interactions, the
# coordinates (a,b,c) with respect to the lattice must be determined and then rounded to the nearest integer.
# This can be achieved by applying the composition of the three Euler rotations to each of the points of ionizing
# interactions in the track and then taking the scalar products with the adjoint vectors that represent the basis
# of the reciprocal space.
# It is computationally much more efficient to apply the inverse transformation once to the base and adjoint vectors
# of the lattice for a certain orientation of the primary particle trajectory and the lattice.
# To assess the dependence of the result on the relative orientation of track and lattice two different approaches can
# be followed. In the first approach, the values of \theta, \phi1 and \phi2 are randomly sampled from uniform
# distributions in the domains [0,\pi/2], [0,2\pi] and [-pi/3,\pi/3], respectively. (This spherical triangle is the 
# fundamental domain of the octahedral symmetry group.)
# In the second approach (which is implemented in module common_bravais_lattice_init.py, deterministic sampling is 
# applied such that the considered orientations of the primary particle trajectory is approximately uniformly 
# distributed in the spherical triangle defined by aforementioned domains for \theta and \phi2. 
# This is achieved by subdividing the spherical triangle into 4**MDIV equal triangles, of which some are pointng up
# (polar angles theta_plus) and some are pointing down (polar angles theta_minus).


# Code starts here -----------------------------------------------------------------------------------------------------
def Inverse_Euler_Matrix(phi1, theta, phi2):
    """
    #     ***************************************************************
    #     calculates combination of three inverse euler rotation matrices
    #     ***************************************************************
    :param phi1:  First azimuthal rotation angle (radians)
    :param theta: Polar rotation angle (radians)
    :param phi2:  Second azimuthal rotation angle (radians)
    :return: Euler matrix for inverse rotation
    """
    #
    rot1 = [[cos(phi1), sin(phi1), 0.],
            [-sin(phi1), cos(phi1), 0.],
            [0., 0., 1.]]
    # for the second rotation, we must take into account that y,z,x is righthanded, so that the sign of theta must
    # be reversed
    rot2 = [[cos(theta), 0, -sin(theta)],
            [0., 1., 0.],
            [sin(theta), 0., cos(theta)]]
    rot3 = [[cos(phi2), sin(phi2), 0.],
            [-sin(phi2), cos(phi2), 0.],
            [0., 0., 1.]]
    return matmul(rot3, matmul(rot2, rot1))


def FCC_basis(nearest_neighbor_distance, direction_of_motion=(1., 1., 1.)):
    """
    # Calculates the array of the three basis vectors of a face-centered cubic (FCC)
    # Bravais lattice from the given lattice constant and then transforms it such that
    # the direction of motion is along the z-axis of the cartesian coordinate system"""
   
    # In the next 2. lines the nearest neighbor distance is divided by sqrt(2.) and then the basis is given ín the
    # cartesian basis of the simple cubic lattice
    dcos45 = nearest_neighbor_distance/sqrt(2.)   # == lattice_constant / 2.
    basis = [[dcos45, dcos45, 0.],
             [0., dcos45, dcos45],
             [dcos45, 0., dcos45]]
    direction_unit_vector = matmul(basis, list(direction_of_motion))
    direction_unit_vector /= norm(direction_unit_vector)
    cos_theta = direction_unit_vector[2]
    theta = acos(cos_theta)
    if abs(cos_theta) == 1:
        phi1 = 0.
    else:
        sin_theta = sin(theta)
        cos_phi = direction_unit_vector[0] / sin_theta
        sin_phi = direction_unit_vector[1] / sin_theta
        phi1 = acos(cos_phi)
        if sin_phi < 0.:
            phi1 = 2. * pi - phi1
    return matmul(Inverse_Euler_Matrix(phi1, theta, 0.), basis)


def Reciprocal_bases(nearest_neighbor_distance, log_base2_of_n_directions):
    # Set dimensions and define arrays
    n_theta_plus = 2 ** log_base2_of_n_directions  # number of polar angles of upward pointing triangles
    cv.n_directions = int(cv.n_azimuthal_rotations_of_track * n_theta_plus**2)   # number of orientations
    cv.Reciprocal_base = zeros((3, 3, cv.n_directions), 'f')
    cv.Real_space_base = zeros((3, 3, cv.n_directions), 'f')
    #
    dphi1 = pi / cv.n_azimuthal_rotations_of_track   # step for rotation of track

    DcosTheta = 1. / n_theta_plus                 # step size in cos theta
    cos_theta_plus = 1. - 2. / 3. * DcosTheta     # cos theta at center of gravity of uppermost upward triangle 
    cos_theta_minus = 1. - 4. / 3. * DcosTheta    # cos theta at center of gravity of uppermost downward triangle 

    phi1 = 0.
    IDIR = 0
    for IAZ in range(cv.n_azimuthal_rotations_of_track):
        for K in range(n_theta_plus):         # Centers of upward pointing triangles
            dphi2 = 2./3.*pi/K                # step of trajectory orientation azimuth
            phi2 = (1./K-1.)*pi/3.            # first azimuth value
            theta = acos(cos_theta_plus)
            for J in range(K+1):
                if cv.DEBUG[1]:
                    print('BLINIT vor CALL EULERM')

                direction_of_motion = matmul(Inverse_Euler_Matrix(phi1, theta, phi2), [0., 0., 1.])
                cv.Reciprocal_base[:, :, IDIR], cv.Real_space_base[:, :, IDIR] = \
                    Reciprocal_basis(nearest_neighbor_distance, direction_of_motion)

                IDIR += 1                     # increment orientation
                phi2 += dphi2  # increment trajectory azimuth

            if cv.DEBUG[1]:
                print('BLINIT nach CALL EULERM')

            cos_theta_plus -= DcosTheta          # increment cosine of polar angle

        if cv.DEBUG[1]:
            print('BLINIT 3.2 IAZ'+str(IAZ))

        for K in range(n_theta_plus-1):           # Centers of downward pointing triangles
            dphi2 = 2./3. * pi / K                  # step of trajectory orientationazimuth
            phi2 = (1. / K - 1.) * pi / 3.          # first azimuth value
            theta = acos(cos_theta_minus)
            for J in range(K+1):
                direction_of_motion = matmul(Inverse_Euler_Matrix(phi1, theta,  phi2), [0., 0., 1.])
                cv.Reciprocal_base[:, :, IDIR], cv.Real_space_base[:, :, IDIR] = \
                    Reciprocal_basis(nearest_neighbor_distance, direction_of_motion)
                IDIR += 1                            # increment orientation
                phi2 += dphi2                        # increment trajectory azimuth
            cos_theta_minus -= cos_theta_minus   # increment cosine of polar angle
        phi1 = phi1 + dphi1                    # increment azimuth for track rotation
    return cv.Reciprocal_base, cv.Real_space_base


def Reciprocal_basis(nearest_neighbor_distance, direction_of_motion=(1., 1., 1.)):
    # Calculates the array of the three basis vectors of the reciprocal space for the
    # face-centered cubic (FCC) Bravais lattice determined in function fcc_basis
    from numpy.linalg import inv
    if len(direction_of_motion) == 3:   # simple case: unique axis
        # In the next line we transpose the basis since it will be used with row vectors
        bravais_base = transpose(FCC_basis(nearest_neighbor_distance, direction_of_motion))
        return inv(bravais_base), bravais_base
    elif len(direction_of_motion) == 1:   # sampling the fundamental domain
        log_base2_of_n_directions = direction_of_motion
        return Reciprocal_bases(nearest_neighbor_distance, log_base2_of_n_directions)
    else:
        print('number of elements in 2nd input variable must be 3 or 1 but is  '+str(len(direction_of_motion)))
        breakpoint()


def Reciprocal_basis_init(spherical_site_diameter, direction_of_motion=(1., 1., 1.)):
    nearest_neighbor_distance = spherical_site_diameter * sqrt(2.) * exp(log(pi / 12.) / 3.)  # lattice constant in nm
    cv.diameter_of_site = nearest_neighbor_distance
    cv.direction_of_motion = direction_of_motion
    return Reciprocal_basis(nearest_neighbor_distance, direction_of_motion)
