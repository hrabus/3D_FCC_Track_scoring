# python version of track structure analysis from Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)
# For manual input call without input arguments  !!!!!!!!!!!!!!

from io import open
# import math
from numpy import zeros, ones
from os.path import exists as file_exists
from re import findall

# import common variables
import common_variables as cv

# import associated python modules
from common_bravais_lattice_init import Reciprocal_basis_init, Version_Date as Version_Date_BLINIT
from common_tools import read_bool, read_float, read_int, read_str, Read_format_init, TSTAMP
from filter_track_data import filter_sites_by_interactions, filter_sites_in_region_of_interest
from IC_3D_Subs import CLUSTR, Version_Date as Version_Date_CLUSTR

# ##################################################################################################
Version_Date = '03-AUG-2022'  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ##############^^^^^^^^^^^#############################################################################################
# Version history:
# 03-AUG-2022: Pythonized the code by using the structure (list) "track_infomation"
# 01-AUG-2022: Started modification of the code structure to extend to the other codes used in the REBS paper by
#              creating the module TRACK_3D.
#              Moved the filtering from the clustering routine to a separate subroutine called in TRACK_3D.
# 31-JUL-2022: Interchanged the input sequence of program file and input version date in the command file
# 30-JUL-2022: Accomplished transformation from Fortran code IC_3D as published in
#              Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)


# Code starts here -----------------------------------------------------------------------------------------------------
def save_track_results(outp_file, track_id, track_results):
    # Save track results to file
    number_of_sites = track_results[0]
    xyz = track_results[1]
    ionizations = track_results[2]
    if cv.edep_flag:
        edep = track_results[3]
        for k in range(number_of_sites):
            outp_file.write('{0:1d} {1:6.3f} {2:6.3f} {3:6.3f} {4:4.2f} {5:1d} \n'.format(track_id,
                                                                                          xyz[k, 0], xyz[k, 1],
                                                                                          xyz[k, 2], edep[k],
                                                                                          ionizations[k]))
    else:
        for k in range(number_of_sites):
            if ionizations[k] > 1:
                outp_file.write('{0:1d} {1:6.3f} {2:6.3f} {3:6.3f} {4:1d} \n'.format(track_id,
                                                                                     xyz[k, 0], xyz[k, 1],
                                                                                     xyz[k, 2], ionizations[k]))
# ______________________________________________________________________________________________________________________
#    END SUBROUTINE # save_track_results


def TRACK_3D(program_name='IC_3D', command_file_name='', output_file_initials=''):
    #    Python version of FORTRAN PROGRAM IC_3D
    # >>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #    ----- Parameters: Version Date and number
    Version_Dates = [Version_Date, Version_Date_BLINIT, Version_Date_CLUSTR]
    #    Version_Input_Format: Version number for command files YYMMDD.HHMM
    #            (date and time  file structure was last changed)
    Version_Input_Format = 220831.18
    # *******************************************************************

    this_line = ''

    DSITE = cv.diameter_of_site
    limits_region_of_interest = [0., 1.]
    IDIR = 0
    CODE = '%%%%%%%'
    cmd_file = ''
    outp_file = ''
    #    ------
    #     CODE: Character string encoding the meaning of the entries in a 
    #           line of the input file as follows:
    #           'T'#   - number of the primary particle track
    #           'X', 'Y', 'Z' - x, y, and z coordinates of the transfer point 
    #           'E'#   - energy deposit (if applicable)
    #           'I'#   - ionization cluster size (if applicable)
    #           '%'#   - additional data that are not used 
    #           Example: The string 'TZ%%EXY ' indicates that there are 7
    #                  entries in each line, of which the first is the 
    #                  track number, the second the z coordinate, the 
    #                  fifth the energy, and the sixth and seventh the x
    #                  and y coordinates
    #     column_ID: Array of the columns indices of (1-3) x, y, z coordinates 
    #                 (4) energy deposit (if present), (5) track number, 
    #                 (6) number of ionizations in cluster (if present)
    #                 (7-8) are there for future use
    NHEADL = 0
    # column_ID = -1

    # <<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
    
    manual_input = command_file_name == ''
    print('Program IC_3D Version: '+Version_Dates[0])
    #
    if manual_input:    # Manual input
        print('Enter debug options ')
        cv.DEBUG[1] = bool(input('Main program'))
        cv.DEBUG[2] = bool(input('TARGTS main sections'))
        cv.DEBUG[3] = bool(input('TARGTS IDIR loop'))
        cv.DEBUG[4] = bool(input('TARGTS IPOS loop details'))
        
        print()
        DSITE = float(input('Enter spherical site diameter in nm'))
        limits_region_of_interest[0] = float(input('z position of begin of region of interest'))
        limits_region_of_interest[1] = float(input('z position of end of region of interest'))

        print()
        CODE = input('Enter 8 character code for data file structure where' +
                     ' "T" indicates the track ID, ' +
                     ' "X", "Y" and "Z" the respective coordinates, ' +
                     ' "E" the energy deposit (if present) and "&" any other data')
        print()
        NHEADL = int(input('Number of header lines'))
        print()
        NFILES = int(input('Number of files to process'))
    elif file_exists(command_file_name):
        # Read intended program name (only for command files)
        command_file_name = read_str(cmd_file)
        if command_file_name[0:5] != 'IC_3D':
            print('Command file appears not to be for IC_3D but for ' + command_file_name + ' = > STOP.')
            exit()

        # Command file is there, so let's go ahead  ...

        # Check command file version
        cmd_file = open(command_file_name, "r")
        VCHECK = read_float(cmd_file)
        if VCHECK < Version_Input_Format:
            print('Command file structure '+str(VCHECK)+' older than ' +
                  'current version '+str(Version_Input_Format)+' ==> STOP.')
            exit()
        
        # Read debug options
        for i in range(1, 5):
            cv.DEBUG[i] = read_bool(cmd_file)

        # Read geometry parameters (site diameter and region of interest boundaries)
        DSITE = read_float(cmd_file)                          # spherical site diameter in nm'
        limits_region_of_interest[0] = read_float(cmd_file)   #
        limits_region_of_interest[1] = read_float(cmd_file)
        
        #   Get code for input file structure
        CODE = read_str(cmd_file)
        
        #   Number of header lines
        NHEADL = read_int(cmd_file)
        #   Number of files to process
        NFILES = read_int(cmd_file)
    else:
        print('Command file '+command_file_name+' does not exist. ==> STOP.')
        NFILES = 0

    if NFILES > 0:
        # Save value of site size
        cv.site_diameter = DSITE
        # interpret the code for the input file structure
        column_ID = Read_format_init(CODE)  
        
        # Init lattice
        cv.Reciprocal_base[:, :, IDIR], cv.Real_space_base[:, :, IDIR] = Reciprocal_basis_init(DSITE, (0, 0, 1))

        # Init arrays
        cv.xyz = zeros((cv.max_number_of_energy_transfer_points, 3), 'f')
        cv.ionizations = ones(cv.max_number_of_energy_transfer_points, 'i')
        #  Check whether energy is there
        cv.edep_flag = (column_ID[4] != -1)    # data column with energy deposit present
        if cv.edep_flag:
            cv.edep = zeros(cv.max_number_of_energy_transfer_points, 'f')
            site_info = [ones(cv.max_number_of_energy_transfer_points, 'i'),
                         ones(cv.max_number_of_energy_transfer_points, 'f'), 
                         ones(cv.max_number_of_energy_transfer_points, 'f'), 
                         ones(cv.max_number_of_energy_transfer_points, 'f'), 
                         ones(cv.max_number_of_energy_transfer_points, 'f')]
        else:
            # cv.edep = array([])
            site_info = [ones(cv.max_number_of_energy_transfer_points, 'i'),
                         ones(cv.max_number_of_energy_transfer_points, 'f'),
                         ones(cv.max_number_of_energy_transfer_points, 'f'),
                         ones(cv.max_number_of_energy_transfer_points, 'f')]
        cv.type_flag = (column_ID[5] != -1)    # data column with interaction type present
    #    
        for ifile in range(NFILES):    
            if manual_input:
                input_file_name = input('Enter file name '+str(ifile+1))
            else:
                input_file_name = read_str(cmd_file)
    
            if cv.DEBUG[1]:
                print('Begin of Block 4')
            #  Init counters#
            number_of_events = 1
    
            #  Open input data file
            input_file = open(input_file_name, "r")
    
            # Open output data file (if requested by calling routine)
            if len(output_file_initials) > 0:
                #   Include site diameter info in output file name
                if DSITE < 10.:
                    output_file_initials = output_file_initials+"{:.1f}nm_".format(DSITE)
                elif DSITE < 100.:
                    output_file_initials = output_file_initials+"{:.0f}nm_".format(DSITE)
                else:
                    output_file_initials = output_file_initials+"{:3d}nm_".format(DSITE)
    
                #        NFORMT = 0
                #        if NFORMT == 1:
                #            output_file_name = output_file_initials+input_file_name[6:80]
                #        else:
                output_file_name = output_file_initials+input_file_name
                outp_file = open(output_file_name, "w")
                outp_file.write('Output from TRACK_3D - Option '+program_name+' Version ' +
                                Version_Dates[0] + ' BLINIT:' + Version_Dates[1] + ' CLUSTR:' + Version_Dates[2]
                                + ' on ' + TSTAMP())
                outp_file.write('Ionization cluster positions' +
                                ' in ion tracks for DSITE = ' + str(DSITE) + ' nm \n')
                outp_file.write('Processed input data file: ' + input_file_name + '\n')
                if cv.edep_flag:
                    outp_file.write('Track x/nm y/nm z/nm Edep/eV ICS \n')
                else:
                    outp_file.write('Track x/nm y/nm z/nm ICS \n')
    
            #  read header lines from input file
            for i in range(NHEADL+1):
                this_line = input_file.readline()  # read the next data line
            read_line = findall(r"[-+]?\d*\.\d+|\d+", this_line)   # get the different values in the line
    
            #  get the first values and init arrays
            ITA = int(read_line[column_ID[0]])
            IT = -1
            # cv.xyz = numpy.zeros((1, 3))
            index_of_transfer_point = 0
    
            for i in range(1, len(site_info)):
                # cv.xyz[index_of_transfer_point, i] = float(read_line[column_ID[i]])
                site_info[i][index_of_transfer_point] = float(read_line[column_ID[i]])
            # if column_ID[4] >= 0:
                # cv.edep = numpy.zeros(1)
            #    cv.edep[index_of_transfer_point] = float(read_line[column_ID[4]])
            # cv.ionizations = numpy.zeros(1, 'i')
            if cv.type_flag:  # interaction type is specified
                # cv.ionizations[index_of_transfer_point] = (int(read_line[column_ID[5]]) % 10 == 3)  # Check
                site_info[0][index_of_transfer_point] = (int(read_line[column_ID[5]]) % 10 == 3)  # Check
            # else:
            #    cv.ionizations[index_of_transfer_point] = 1
    
            eof = False
            while not eof:
                the_line = input_file.readline()
                eof = not the_line
                if not eof:
                    read_line = findall(r"[-+]?\d*\.\d+|\d+", the_line)
                    IT = int(read_line[column_ID[3]])
                    index_of_transfer_point += 1
                    for i in range(1, len(site_info)):
                        site_info[i][index_of_transfer_point] = float(read_line[column_ID[i]])
                        # cv.xyz[index_of_transfer_point, i] = float(read_line[column_ID[i]])
                    # xyz = [read_line[column_ID[0]], read_line[column_ID[1]], read_line[column_ID[2]]]
                    # cv.xyz = numpy.vstack((cv.xyz, numpy.float_(xyz)))
                    # if cv.edep_flag:
                        #  cv.edep[index_of_transfer_point] = float(read_line[column_ID[4]])
                        # cv.edep = numpy.append((cv.edep, float(read_line[column_ID[4]])))
                    if cv.type_flag:
                        site_info[0][index_of_transfer_point] = (int(read_line[column_ID[5]]) % 10 == 3)  # Check
                        # cv.ionizations[index_of_transfer_point] = (int(read_line[column_ID[5]]) % 10 == 3)  # CHECK
                    #    cv.ionizations = numpy.append((cv.ionizations, int(read_line[column_ID[5]])) == 3)  # CHECK
                    # else:
                    #    cv.ionizations = numpy.append((cv.ionizations, 1))
                else:
                    IT += 1
                # COMMENT
                if IT != ITA:  # new track
                    number_of_sites = index_of_transfer_point  # last index belongs to next track
                    if cv.DEBUG[1]:
                        print('Before CLUSTR')
                    if number_of_events % 100 == 0:
                        print(number_of_events)
    
                    # call to SUBROUTINE CLUSTR which detects clusters
    
                    if number_of_sites > 0:  # this is needed for the hendling of the last track
                        track_results = [number_of_sites, site_info]
                        if cv.minimum_number_of_ionizations > 0 or cv.minimum_energy_imparted > 0.:
                            filter_sites_by_interactions(track_results, [cv.minimum_number_of_ionizations,
                                                                         cv.minimum_energy_imparted])
                        if track_results[0] > 0:
                            track_results = CLUSTR(track_results, IDIR)
    
                        if track_results[0] > 0:
                            filter_sites_in_region_of_interest(track_results, limits_region_of_interest)
                        
                        if track_results[0] == 0: 
                            number_of_events -= 1  # there were no results for this track
                        elif len(output_file_initials) > 0:
                            save_track_results(outp_file, ITA, track_results)
    
                    if cv.DEBUG[1]:
                        print('After CLUSTR')
                    if not eof:
                        cv.xyz[0, :] = cv.xyz[index_of_transfer_point, :]  # 1st transfer point in track
                        if cv.edep_flag:
                            cv.edep[0] = cv.edep[index_of_transfer_point]
                        cv.ionizations[0] = cv.ionizations[index_of_transfer_point]
                        ITA = IT                           # Remember new track number
                        number_of_events += 1
                        index_of_transfer_point = 0        # There was already one energy transfer point
            outp_file.close()
            input_file.close()

    if not manual_input:
        cmd_file.close()
# ______________________________________________________________________________________________________________________
#    END SUBROUTINE # TRACK_3D


def IC_3D(command_file_name: str = ''):
    TRACK_3D('IC_3D', command_file_name, 'IC_')
# ______________________________________________________________________________________________________________________
#    END SUBROUTINE # IC_3D


IC_3D('3.0nm_IC.inp')
