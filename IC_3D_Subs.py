# python version of track structure analysis from Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)
#
# Open points:
# - check whether filtering the relevant data can be done more efficient

from numpy import zeros, round, matmul, argsort, argwhere, int_
# import numpy

# import common variables
import common_variables as cv

# ##################################################################################################
Version_Date = '03-AUG-2022'  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ##############^^^^^^^^^^^#############################################################################################
# Version history:
# 03-AUG-2022: Pythonized the code by using the structure (list) "track_infomation"
# 01-AUG-2022: Moved filtering of data to a separate module to be called by the calling module
# 30-JUL-2022: Moved loop over directions to the calling module
# 26-JUL-2022: Accomplished transformation from Fortran code IC_3D_Subs as published in
#              Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)


# Code starts here -----------------------------------------------------------------------------------------------------
def CLUSTR(track_information, scoring_information=(0, 0.), sort_output=True):
    #     ***************************************************************
    #     Determine targets with ionization clusters
    #     ***************************************************************

    # >>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #
    number_of_energy_transfer_points = track_information[0]
    if number_of_energy_transfer_points == 0:
        print('User, we have a problem: Number of energy transfer points is 0')
        breakpoint()

    interactions = zeros(number_of_energy_transfer_points, 'i')

    # ionizations = track_information[1][0:number_of_energy_transfer_points, 0]
    # xyz = track_information[1][0:number_of_energy_transfer_points, 1: 4]
    # if len(track_information) == 4:   # energy information
    #     edep = track_information[3][0:number_of_energy_transfer_points]
    # else:
    #     edep = array([])
    site_info = track_information[1][0: number_of_energy_transfer_points]

    r2_sphere = 0.
    if type(scoring_information).__name__ == 'int':
        scoring_information = (scoring_information, 0.)
    if type(scoring_information).__name__ == 'float':
        scoring_information = (0, scoring_information)
    if ('list', 'tuple').count(type(scoring_information).__name__) == 0:
        print('User we have a problem: scoring information is not comprehensible')
        breakpoint()
    if type(scoring_information[0]).__name__ == 'int':
        IDIR = scoring_information[0]
        reciprocal_base = cv.Reciprocal_base[:, :, IDIR]
        real_space_base = cv.Real_space_base[:, :, IDIR]
        r2_sphere = scoring_information[1]
    else:
        reciprocal_base = scoring_information[0]
        real_space_base = scoring_information[1]
        if len(scoring_information) == 3:
            r2_sphere = scoring_information[2]

    sphere_flag = r2_sphere > 0.
    #
    #  <<< End declarations <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if cv.SRFLAG:
        print('Subroutine CLUSTR Version: '+Version_Date)
        cv.SRFLAG = False

    if cv.DEBUG[1]:
        print('CLUSTR vor Main Loop NDIR = '+str(cv.n_directions))
        print('CLUSTR vor Main Loop number_of_energy_transfer_points = '+str(number_of_energy_transfer_points))

    #  Find target volumes for all energy transfer points in track
    ITARGT = int_(round(matmul(site_info[1:3], reciprocal_base)))
    # breakpoint()
    if cv.DEBUG[2]:
        print('CLUSTR vor sort volume indices'+str(number_of_energy_transfer_points))

    if sphere_flag:
        centers = matmul(ITARGT, real_space_base)
        indices = argwhere((site_info[1:3]-centers)**2 < r2_sphere)
        number_of_energy_transfer_points = len(indices)
        if number_of_energy_transfer_points > 0:
            site_info = site_info[indices, :]
            ITARGT = ITARGT[indices]

    if number_of_energy_transfer_points > 0:
        #  Sort target volume indices ascending
        indices = argsort(ITARGT[:, 0].ravel())                      # Now 1st index is ascending
        indices = indices[argsort(ITARGT[indices, 1].ravel())]       # now 2nd index is ascending
        indices = indices[argsort(ITARGT[indices, 2].ravel())]       # now 3rd index is ascending
        ITARGT = ITARGT[indices, :]           # sort indices
        site_info = site_info[indices, :]  # sort site info

        if cv.DEBUG[2]:
            print('CLUSTR vor find unique volumes'+str(number_of_energy_transfer_points))
        #  Find unique target volumes and score ionizations
        number_of_sites = 0
        interactions[number_of_sites] = 1
        for i in range(1, number_of_energy_transfer_points):
            if (ITARGT[i, 0] != ITARGT[i-1, 0] or ITARGT[i, 1] != ITARGT[i-1, 1] or
                    ITARGT[i, 2] != ITARGT[i-1, 2]):
                # site is only accepted if cluster size is large enough
                number_of_sites = number_of_sites + 1
                interactions[number_of_sites] = 1
                site_info[number_of_sites, :] = site_info[i, :]
            else:
                interactions[number_of_sites] += 1
                site_info[number_of_sites, :] += site_info[i, :]

        # Normalize positions
        number_of_sites += 1
        for i in range(1, 4):
            if cv.tpos_flag:
                site_info[0:number_of_sites, i] = matmul(ITARGT[0:number_of_sites, :], real_space_base)
            elif cv.type_flag:
                site_info[0:number_of_sites, i] /= interactions[0:number_of_sites].ravel()
            else:
                site_info[0:number_of_sites, i] /= site_info[0:number_of_sites, 0].ravel()

        if sort_output:
            # Sort position ascending along track
            indices = argsort(site_info[0:number_of_sites, 3].ravel())  # sequence such that z is ascending
            site_info[0:number_of_sites, :] = site_info[indices, :]  # sort positions

        if cv.DEBUG[2]:
            print('CLUSTR vor EXIT')

        return number_of_sites, site_info[0:number_of_sites, :]
    else:
        print('No energy transfer points in region of interest')
        return 0, None, None
# ______________________________________________________________________________________________________________________
