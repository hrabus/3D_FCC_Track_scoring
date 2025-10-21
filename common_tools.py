# python version of track structure analysis from Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)

from numpy import zeros
from re import findall

# ##################################################################################################
Version_Date = '30-JUL-2022'  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ##############^^^^^^^^^^^#############################################################################################
# Version history:
# 30-JUL-2022: Accomplished transformation from Fortran codes RFINIT (here Read_format_init) and TSTAMP as published in
#              Ngcezu S and Rabus H, Radiat Environ Biophys 60, 559-578 (2021)
#              The other functions have been added are for easier coding.


# Code starts here -----------------------------------------------------------------------------------------------------
def read_bool(file):
    return file.readline().split('#')[0].strip(' ') == 'True'


def read_int(file):
    text = file.readline()
    return int(findall(r"\d", text)[0])


def read_float(file):
    text = file.readline()
    return float(findall(r"[-+]?\d*\.\d+|\d+", text)[0])


def read_str(file):
    return file.readline().split('#')[0].strip(' ')


def Read_format_init(fcode: str):
    # Sets the information for the columns of the input file
    #
    # ------
    #  fcode: String encoding the meaning of the entries in a line of the input file as follows:
    #         'T' - number of the primary particle track
    #         'X', 'Y', 'Z' - x, y, and z coordinates of the transfer point
    #         'E' - energy deposit( if applicable)
    #         'I' - ionization cluster size( if applicable)
    #         '%' - additional data that are not used
    #         Example: The string 'TZ%%EXY ' indicates that there are 7
    #                  entries in each line, of which the first is the
    #                  track number, the second the z coordinate, the
    #                  fifth the energy, and the sixth and seventh the x
    #                  and y coordinates
    # index_data_columns: Array of the columns indices of
    #                     (0..2) x,y,z coordinates
    #                     (4) energy deposit (if present), (5) track number,
    #                     (6) number of ionizations in cluster (if present)
    #                     (7-8) are there for future use
    #
    index_data_columns = zeros(8, 'i') - 1
    number_data_columns = len(fcode)
    for j in range(number_data_columns):
        if (fcode[j:j+1] == 'T') and index_data_columns[0] == -1:
            index_data_columns[0] = j
        elif (fcode[j:j+1] == 'X') and (index_data_columns[1] == -1):
            index_data_columns[1] = j
        elif (fcode[j:j+1] == 'Y') and index_data_columns[2] == -1:
            index_data_columns[2] = j
        elif (fcode[j:j+1] == 'Z') and index_data_columns[3] == -1:
            index_data_columns[3] = j
        elif (fcode[j:j+1] == 'E') and index_data_columns[4] == -1:
            index_data_columns[4] = j
        elif (fcode[j:j+1] == 'I') and index_data_columns[5] == -1:
            index_data_columns[5] = j
        elif fcode[j:j+1] != ' ':
            number_data_columns = j

    if (index_data_columns[0] * index_data_columns[1] *
            index_data_columns[2] * index_data_columns[4]) == -1:
        if index_data_columns[0] == -1:
            print('No data column for track number')
        if index_data_columns[1] == -1:
            print('No data column for X')
        if index_data_columns[2] == -1:
            print('No data column for Y')
        if index_data_columns[3] == -1:
            print('No data column for Z')
        print('error')
        breakpoint()

    index_data_columns = tuple(index_data_columns[0: number_data_columns])
    return index_data_columns   # , number_data_columns


def TSTAMP():
    import time
    return time.asctime() + '  UTC ' + str(time.timezone) + ' s'
