# python version of track structure analysis from Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)
#
# Open points:
# - check whether filtering the relevant data can be done more efficient

from numpy import zeros, round, matmul, argsort, argwhere, int_
# import numpy

# import common variables
import common_variables as cv

# ##################################################################################################
Version_Date = '30-JUL-2022'  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ##############^^^^^^^^^^^#############################################################################################
# Version history:
# 01-AUG-2022: Moved filtering of data to the calling module
# 30-JUL-2022: Moved loop over directions to the calling module
# 26-JUL-2022: Accomplished transformation from Fortran code IC_3D_Subs as published in
#              Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)


# Code starts here -----------------------------------------------------------------------------------------------------
def filter_track_data(track_information, filter_information: tuple):
    #     ***************************************************************
    #     reduce the track information to those values atisfying the filter condirtoion
    #     first condition is on geometry, second on number of ionizatzions 
    #     (or whether an interaction is an ionization) and the energy deposit (or 
    #     energy imparted in site)
    #
    #     ***************************************************************

    # >>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #
    number_of_energy_transfer_points = track_information[0]
    xyz = track_information[1][0:number_of_energy_transfer_points, :]
    ionizations = track_information[2][0:number_of_energy_transfer_points]
    edep_flag = (len(track_information) == 4)
    if edep_flag:   # energy information
        edep = track_information[3][0:number_of_energy_transfer_points]
    
    # Check wehther filter information conatins several criteria
    if ['tuple','list'].count(type(filter_information).__name__) == 0:
        # only a signle value specified
        if type(filter_information) == 'int':  # only minimum number of ionizations specified
            minimum_number_of_ionizations = filter_information
            indices = argwhere(ionizations >= minimum_number_of_ionizations)
            number_of_energy_transfer_points = len(indices)
            if number_of_energy_transfer_points == 1:
                if indices[0] == -1:
                    number_of_energy_transfer_points = 0
            if number_of_energy_transfer_points > 0:
                xyz[0:number_of_energy_transfer_points, :] = xyz[indices, :]
                ionizations[0:number_of_energy_transfer_points] = ionizations[indices]
                if edep_flag:
                    edep[0:number_of_energy_transfer_points] = edep[indices]
        elif type(filter_information) == 'float':  # only minimum number of energy specified:
            if edep_flag:
                minimum_energy = filter_information
                indices = argwhere(edep >= minimum_energy)
                number_of_energy_transfer_points = len(indices)
                if number_of_energy_transfer_points == 1:
                    if indices[0] == -1:
                        number_of_energy_transfer_points = 0
                if number_of_energy_transfer_points > 0:
                    xyz[0:number_of_energy_transfer_points, :] = xyz[indices, :]
                    ionizations[0:number_of_energy_transfer_points] = ionizations[indices]
                    edep[0:number_of_energy_transfer_points] = edep[indices]
            else:
                print('filter_track_data % energy threshold specified, but no energy data ==> do FILTER APPPLIED!')
        else:
            print('incomprehensible filter information '+str(filter_information)+' ==> no filter applied')
    else:
        limits_region_of_interest = scoring_information[0]
    if ['tuple','list'].count(type(limits_region_of_interest).__name__) == 0:
        dummy = 0 # to prevent crash    
    
    if scoring_information[1].__class__.__name__ == 'int':
        IDIR = scoring_information[1]
        reciprocal_base = cv.Reciprocal_base[:, :, IDIR]
        # reciprocal_base = cv.Reciprocal_base[:, :, IDIR]
        real_space_base = cv.Real_space_base[:, :, IDIR]
        if len(scoring_information) == 3:
            minimum_number_of_ionizations = scoring_information[2]
        else:
            minimum_number_of_ionizations = 0
    else:
        reciprocal_base = scoring_information[1]
        # reciprocal_base = scoring_information[1]
        real_space_base = scoring_information[2]
        if len(scoring_information) == 4:
            minimum_number_of_ionizations = scoring_information[2]
        else:
            minimum_number_of_ionizations = 0
    #
    #  <<< End declarations <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if cv.SRFLAG:
        print('Subroutine CLUSTR Version: '+Version_Date)
        cv.SRFLAG = False

    if cv.DEBUG[1]:
        print('CLUSTR vor Main Loop NDIR = '+str(cv.n_directions))
        print('CLUSTR vor Main Loop number_of_energy_transfer_points = '+str(number_of_energy_transfer_points))
         
    #  Filter the energy transfer points within the region of interest
    indices = argwhere((xyz[:, 2] >= limits_region_of_interest[0]) *
                       (xyz[:, 2] <= limits_region_of_interest[1])).ravel()
    number_of_energy_transfer_points = len(indices)
    if number_of_energy_transfer_points > 0:
        xyz = xyz[indices, :]
        ionizations = ionizations[indices]
        edep = edep[indices]

        #  Find target volumes for all energy transfer points in track
        ITARGT = int_(round(matmul(xyz, reciprocal_base)))
        # breakpoint()
        if cv.DEBUG[2]:
            print('CLUSTR vor sort volume indices'+str(number_of_energy_transfer_points))

        #  Sort target volume indices ascending
        indices = argsort(ITARGT[:, 0].ravel())                      # Now 1st index is ascending
        indices = indices[argsort(ITARGT[indices, 1].ravel())]       # now 2nd index is ascending
        indices = indices[argsort(ITARGT[indices, 2].ravel())]       # now 3rd index is ascending
        ITARGT = ITARGT[indices, :]           # sort indices
        xyz = xyz[indices, :]                 # sort positions
        ionizations = ionizations[indices]
        if edep_flag:
            edep = edep[indices]

        if cv.DEBUG[2]:
            print('CLUSTR vor find unique volumes'+str(number_of_energy_transfer_points))
        #  Find unique target volumes and score ionizations
        number_of_sites = 0
        interactions[number_of_sites] = 1
        for i in range(1, number_of_energy_transfer_points):
            if (ITARGT[i, 0] != ITARGT[i-1, 0] or ITARGT[i, 1] != ITARGT[i-1, 1] or
                    ITARGT[i, 2] != ITARGT[i-1, 2]):
                # site is only accepted if cluster size is large enough
                if ionizations[number_of_sites] >= minimum_number_of_ionizations:
                    number_of_sites = number_of_sites + 1
                interactions[number_of_sites] = 1
                ionizations[number_of_sites] = ionizations[i]
                for j in range(3):
                    ITARGT[number_of_sites, j] = ITARGT[i, j]
                    xyz[number_of_sites, j] = xyz[i, j]
                if edep_flag:
                    edep[number_of_sites] = edep[i]
            else:
                interactions[number_of_sites] += 1
                ionizations[number_of_sites] += ionizations[i]
                for j in range(3):
                    xyz[number_of_sites, j] += xyz[i, j]
                if edep_flag:
                    edep[number_of_sites] += edep[i]

        # Normalize positions
        number_of_sites += 1
        for i in range(3):
            if cv.tpos_flag:
                xyz[0:number_of_sites, i] = matmul(ITARGT[0:number_of_sites, :], real_space_base)
            elif cv.type_flag:
                xyz[0:number_of_sites, i] /= interactions[0:number_of_sites].ravel()
            else:
                xyz[0:number_of_sites, i] /= ionizations[0:number_of_sites].ravel()

        if sort_output:
            # Sort position ascending along track
            indices = argsort(xyz[0:number_of_sites, 2].ravel())  # sequence such that z is ascending
            xyz[0:number_of_sites] = xyz[indices, :]  # sort positions
            ionizations[0:number_of_sites] = ionizations[indices]
            if edep_flag:
                edep[0:number_of_sites] = edep[indices]

        if cv.DEBUG[2]:
            print('CLUSTR vor EXIT')

        if edep_flag:
            return number_of_sites, xyz[0:number_of_sites, :], ionizations[0:number_of_sites], edep[0:number_of_sites]
        else:
            return number_of_sites, xyz[0:number_of_sites, :], ionizations[0:number_of_sites]
    else:
        print('No energy transfer points in region of interest')
        return 0, None, None
# ______________________________________________________________________________________________________________________
