      PROGRAM ROI_3D
# >>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
#     - - - - - Parameters: Version Date and number
      CHARACTER VDATE * 11
      REAL * 8 VINPUT    #  Version number for command files YYMMDD.HHMM 
#                      [date and time  file structure was last changed]
#     #################################################################
      PARAMETER[VDATE = '05 - JUN - 2021', VINPUT = 210605.16d0]   #  ############
#     #################################################################
#     05 - JUN - 2021 HR:
#       - Added timestap to output files 
#       - Reworked identification of input file structure
#     23 - MAY - 2021 HR: 
#       - Modified call to subroutines to get their version date
#     17 - APR - 2021 HR:
#       - Increased MXTRG2 owing to data set with larger values
#     13 - APR - 2021 HR:
#       - Fixed soft bug: Results were overwritten due to same filename
#       - Added reduction of output to non - zero values
#     02 - APR - 2021 HR:
#       - Added output of correlation matrix
#       - Detected and fixed bug in section with convolution
#     28 - MAR - 2021 HR:
#       - Changed normalization of traversing tracks
#     20 - MAR - 2021 HR:
#       - Added input of options file name
#       - Added discrimination between site size and lattice constant
#     14 - MAR - 2021 HR:
#       - Fixed bug with treatment of last track
#       - Added common block CLSTRS for info on cluster positions  
#     16 - OCT - 2020 HR:
#       - Added readin of geometry parameters via options file
#       - Added readin of multiple file names to process with options
#       - Added processing of multiple ROIs in track data set
#     15 - MAR - 2020 HR:
#       - Added VDATE and modified COMMON BLOCK VERBOSE
#       - Added variable NXYPOS to fix bug in TARG3D
#     - - - - - Parameters: General purpose constants
      INTEGER * 4 LUN
      REAL * 8 EPS, ONE, ZERO
      PARAMETER[LUN = 11, EPS = 1.0e - 8, ONE = 1.0, ZERO = 0.0]
#     - - - - - Parameters for lattice orientation
      INTEGER * 4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER[MAZMTH in range(MDIV = 0, MTHETA = 2 *  * MDIV,  
     &           MDIR = MAZMTH * [MTHETA * [MTHETA + 1]] / 2]
#     - - - - - Parameters for track and radial distance histograms
      INTEGER * 4 KMAX, MIONIZ, MXRPOS, MXTARG, MXTRG2, MXYPOS, MXZROI
      PARAMETER[KMAX = 9, MIONIZ = 100000, MXRPOS = 101, MXTARG = 1025, 
     &          MXTRG2 = 513, MXYPOS = 4 * MXRPOS * [MXRPOS - 1] + 1, MXZROI = 100]
#     - - - - - Input parameters for geometry
      REAL * 8 DBEAM    #  Diameter of beam in nm
      REAL * 8 DLATC    #  Cell lattice constant 
      REAL * 8 DROI     #  Diameter of region of interest [cell nucleus size]
      REAL * 8 DSITE    #  Diameter of spherical target 
      REAL * 8 DZROI    #  Increment in position of region of interest 
#     - - - - - Functions
      CHARACTER TSTAMP * 24
#     - - - - - Local scalars
      CHARACTER * 11 VDATES[2], SRDATE
      CHARACTER FILENM * 80, FILOUT * 85, HEADER * 80, PARAM * 2, PREFIX * 10
      INTEGER * 4 I, IFILE, IT, ITA, J, K, KMX, L, NAZMTH, NDIV, NFILES, 
     &          NHEADL, NZROI, NPHI, NTRACS
      INTEGER * 4 NIONIZ, NRPOS, NXYPOS
      LOGICAL ASKINP
      REAL * 8 COSPHI, DELTAR, DUMMY, FNORM, PIBY4, RROI, 
     &       SINPHI, SPOSIN, VCHECK
#     - - - - - Local arrays
      INTEGER * 4 NTARG[2], NTARGK[KMAX]
      REAL * 8 CORRSE[MXTARG, MXTRG2], CONVOL[MXTARG], FREQBV[MXTARG, MXTRG2]
      REAL * 8 FREQSE[MXTARG, KMAX], FREQTE[MXTARG, KMAX], ZROIC[MXZROI]
      REAL * 8 RLINE[8]
#     - - - - - Global variables
      INTEGER * 4 NDIR
      REAL * 8 ADJNTB[3, 3, MDIR]
      COMMON / LATTICE / ADJNTB, NDIR 
#     - - - - - - 
      INTEGER * 4 IRAD[MXYPOS]
      INTEGER * 4 NRROI2   #  Integer of Square of ROI radius 
      REAL * 8 XROI[MXYPOS], YROI[MXYPOS], ZROI[3]
      COMMON / ROICTR /  XROI, YROI, ZROI, IRAD, NRROI2
#     - - - - - - 
      REAL * 8 XYZ[MIONIZ, 3]
      COMMON / TRACKS / XYZ
#     - - - - - - 
      REAL * 8 CORR12[MXTARG, MXTRG2, MXRPOS]
      REAL * 8 RADIST[MXRPOS], FREQRD[MXTARG, MXRPOS, KMAX]
#  *       COMMON / HISTOG / RADIST, FREQRD
      COMMON / HISTOG / RADIST, FREQRD, CORR12
#            RADIST is the vector of radial distances 
#            FREQRD initially holds the sum, the sum of squares and the
#                   sum of variances per track over all tracks. 
#                   In the main program this is converted in the end to 
#                   mean and variance over all analyzed tracks plus the
#                   the average intra - track variance for to orientation
#     - - - - - - 
      LOGICAL DEBUG[4], SRFLAG
      COMMON / VERBOSE / DEBUG, SRFLAG
#     - - - - - - 
      INTEGER * 4 ICSIZE[MIONIZ], NSITES
      COMMON / CLSTRS / ICSIZE, NSITES
#     - - - - - - 
#      CODE: Character string encoding the meaning of the entries in a 
#            line of the input file as follows:
#            'T'     - number of the primary particle track
#            'X', 'Y', 'Z' - x, y, and z coordinates of the transfer point 
#            'E'     - energy deposit [if applicable]
#            'I'     - ionization cluster size [if applicable]
#            '%'     - additional data that are not used 
#            Example: The string 'TZ%%EXY ' indicates that there are 7
#                     entries in each line, of which the first is the 
#                     track number, the second the z coordinate, the 
#                     fifth the energy, and the sixth and seventh the x
#                     and y coordinates
#      IDT: Array of the columns indices of [1 - 3] x, y, z coordinates 
#                    [4] energy deposit [if present], [5] track number, 
#                    [6] number of ionizations in cluster [if present]
#                    [7 - 8] are there for future use
      CHARACTER CODE * 8
      INTEGER * 2 IDT[8], NDT
      COMMON / RFORMT / IDT, NDT, CODE
# <<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      print('Program ROI_3D Version: ' + VDATE
      SRFLAG = True   #  Invoke message on Version date from Subroutine TARG3D 

      if 1 == 1:    #  DUMMY block for readability: Read input 1111111
        DEBUG[1] = False   #  Main program
        DEBUG[2] = False   #  TARG3D main sections
        DEBUG[3] = False   #  TARG3D IDIR loop
        DEBUG[4] = False   #  TARG3D IPOS loop details
        
        PIBY4 = ATAN[ONE]
        PREFIX = '3D__0.0nm_'
        for K in range(KMAX
          NTARGK[K] = 0
        

#       Check command file version
        print('Enter 0 for manual input or version [YYMMDD.HHMM] ' + 
     &          'of command file structure ' 
        = input(VCHECK
        ASKINP = [VCHECK == ZERO] 
        if  not ASKINP: 
          if VCHECK < VINPUT: 
            print('Command file structure ', VCHECK, ' older than ' + 
     &              'current version ', VINPUT, ' = > STOP.'
            STOP
                 
          = input(FILENM
          if FILENM[1:6]   # = 'ROI_3D': 
            print('Command file appear not to be for ROI_3D but ' + 
     &              'for ' + FILENM + ' = > STOP.'
            STOP
              
        
        
#       Read debug options
        if ASKINP: print('Enter debug options file name or - for none' 
        = input(FILENM        
        open(LUN, FILE = FILENM, STATUS = 'OLD', ERR = 10]
        for I in range(4
          = input(DEBUG[I]
           #  I in range(4
        .close()
  10    CONTINUE

#       Read geometry parameters 
        if ASKINP: print('Enter site diameter in nm'
        = input(DSITE
        DLATC = DSITE * SQRT[2.] * EXP[LOG[PIBY4 / 3.] / 3.] 
        if DSITE < 10.0:  
          .write(PREFIX[5:7], '[f3.1]'] DSITE
        else::        
          .write(PREFIX[4:7], '[f4.1]'] DSITE
        
        
        if ASKINP: print('Enter ROI diameter in nm'
        = input(DROI
        RROI = DROI / 2.
        NRROI2 = int(RROI * RROI / DLATC / DLATC]   #  Note: This is correct as DLATC is the unit of length
        DELTAR = DROI / 40.
        
        if ASKINP: print('Enter beam diameter in nm <= ', 5. * DROI
        = input(DBEAM
        NRPOS = 1 + round(DBEAM / 2. / DELTAR]
        if NRPOS > MXRPOS: 
          NRPOS = MXRPOS
          print('Beam diameter ', DBEAM, ' nm is too large. <<<<<<<<<<<'
          DBEAM = 2. * DELTAR * float(MXRPOS - 1]
          print('>>> maximum possible value ', DBEAM, ' is used.'
        
        
        if ASKINP: print('Enter maximum # of targets in histogram'
        = input(NTARG[1], NTARG[2]
        if NTARG[1] > MXTARG.OR.NTARG[1] < 1] NTARG[1] = MXTARG
        if NTARG[2] > MXTRG2.OR.NTARG[2] == 1] NTARG[2] = MXTRG2

        if ASKINP: print('Enter maximum ionization cluster ' + 
     &                     'complexity [KMAX]'
        = input(KMX
        if KMX > KMAX.OR.KMX < 2] KMX = KMAX

        if ASKINP: print('Enter # of regions of interest along track'
        = input(NZROI
        IF [NZROI == 1: 
          if ASKINP: print('Enter z position of region of interest'
          = input(ZROIC[1]
          DZROI = ZERO
        else::
          if ASKINP: print('Enter position of first region of ' + 
     &            'interest [ROI] and increment in ROI position'
          = input(ZROIC[1], DZROI
          for I = 2, NZROI
            ZROIC[I] = ZROIC[I - 1] + DZROI
          
        

        if ASKINP: print('Enter 8 character code for data file ' + 
     &       'structure where ''T'' indicates the track ID, ' + 
     &       ' ''X'', ''Y'' and ''Z'' the respective coordinates, ' + 
     &       '''E'' the energy deposit [if present] and ''&'' any ' + 
     &       'other data'
        = input(CODE
        CALL RFINIT[]

        if ASKINP: print('Number of header lines'
        = input(NHEADL

        if ASKINP: print('Number of files to process'
        = input(NFILES
        
                  #  DUMMY block for readability: Read input 1111111

      if 2 == 2:    #  DUMMY block for readability: Initialize 2222222
 *         if DEBUG[1]] print('Hier v'
        NAZMTH = MAZMTH
        NDIV = MDIV
        if DEBUG[1]] print('Hier vor CALL BLINIT'
        CALL BLINIT[NAZMTH, NDIV, DLATC, SRDATE]   #  Init reziprocal lattice
        VDATES[1] = SRDATE
        if DEBUG[1]] print('Hier nach CALL BLINIT', NRPOS, NTARG
        
        for I in range(NRPOS
#         Define radial offsets of track w.r.t. ROI center 
          RADIST[I] = float(I - 1] * DELTAR
#         x&y positions of track w.r.t. ROI center [piecake method]
          if I == 1: 
            NXYPOS = 1
            NPHI = 0
            IRAD[NXYPOS] = 1
            XROI[NXYPOS] = ZERO
            YROI[NXYPOS] = ZERO
          else::
            NPHI = NPHI + 8
            NXYPOS = NXYPOS + 1
            IRAD[NXYPOS] = I
            XROI[NXYPOS] = RADIST[I]
            YROI[NXYPOS] = ZERO
            COSPHI = COS[PIBY4 / float(I - 1]]
            SINPHI = SIN[PIBY4 / float(I - 1]]
            for J = 2, NPHI
              NXYPOS = NXYPOS + 1
              IRAD[NXYPOS] = I
              XROI[NXYPOS] = XROI[NXYPOS - 1] * COSPHI - YROI[NXYPOS - 1] * SINPHI
              YROI[NXYPOS] = XROI[NXYPOS - 1] * SINPHI + YROI[NXYPOS - 1] * COSPHI
            
          
           #  for I in range(NRPOS
                  #  DUMMY block for readability: Initialize 2222222

      for IFILE in range(NFILES      
        if ASKINP: print('Enter file name ', IFILE
        = input(FILENM
        if  not ASKINP: print('Processing file ', FILENM
        
        if 4 == 4:    #  DUMMY block for readability: Main Loop  444444
          if DEBUG[1]] print('Hier beginnt Block 4'
          for J in range(NRPOS
#           Init counters       
            for I in range(NTARG[1]
              for K in range(KMAX
                FREQRD[I, J, K] = ZERO
              
              
#           Begin 02 - APR - 2021 >>>>>>>>>>>>>>>>>     
            for I in range(NTARG[1]
              for L in range(NTARG[2]
                CORR12[I, L, J] = ZERO
              
              
#           END 02 - APR - 2021 <<<<<<<<<<<<<<<<<<  
             #  for J in range(NRPOS
          NIONIZ = 1
          NTRACS = 1
          
          open(LUN, FILE = FILENM, STATUS = 'OLD']
          for I in range(NHEADL
            = input(HEADER
          
   90     = input([RLINE[I], I in range(NDT]
          if IDT[6] > 0: 
            if round(RLINE[IDT[6]]] < 2] GOTO 90   #  No ionization cluster
           
          ITA = round(RLINE[IDT[5]]]
          for I in range(3
            XYZ[NIONIZ, I] = RLINE[IDT[I]]
          
          
          NIONIZ = NIONIZ + 1
          if DEBUG[1]] print('Hier 4'
  100     CONTINUE   #  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          = input(LUN, * , END = 200, ERR = 200] [RLINE[I], I in range(NDT]
          if IDT[6] > 0: 
            if round(RLINE[IDT[6]]] < 2] GOTO 100   #  No ionization cluster
           
          IT = round(RLINE[IDT[5]]]
          for I in range(3
            XYZ[NIONIZ, I] = RLINE[IDT[I]]
          

          if IT == ITA: 
            NIONIZ = NIONIZ + 1
          else::   #  Entry from a new track was read
            NIONIZ = NIONIZ - 1   #  Decrease to number of ionizations in track
            if DEBUG[1]] print('Hier vor CALL TARG3D'
            if MOD[NTRACS, 10] == 0] print(NTRACS
            for J in range(NZROI 
              ZROI[2] = ZROIC[J]
              ZROI[1] = ZROI[2] - DROI / 2. - DLATC / 2.
              ZROI[3] = ZROI[2] + DROI / 2. + DLATC / 2.
              CALL TARG3D[NIONIZ, NRPOS, NTARG, NXYPOS, SRDATE]
               #  J in range(NZROI

            if DEBUG[1]] print('Hier nach CALL TARG3D'
            for J in range(3
              XYZ[1, J] = XYZ[NIONIZ + 1, J]   #  First ionization in new track
            
            ITA = IT   #  Remember new track number
            NTRACS = NTRACS + 1
            NIONIZ = 2   #  There was already one ionization
          
          GOTO 100   #  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  200     CONTINUE
            for J in range(NZROI 
              ZROI[2] = ZROIC[J] 
              ZROI[1] = ZROI[2] - DROI / 2. - DLATC / 2.
              ZROI[3] = ZROI[2] + DROI / 2. + DLATC / 2.
              CALL TARG3D[NIONIZ, NRPOS, NTARG, NXYPOS, SRDATE] 
               #  J in range(NZROI 
            VDATES[2] = SRDATE
          .close()
                    #  DUMMY block for readability: Main Loop  4444444

        if 6 == 6:    #  DUMMY block for readability: Normalize  666666
#         #Normalization
          FNORM = ONE / float(NTRACS]
          FNORM = ONE / float(NTRACS] / float(NZROI] 
 *           SUMCNT = ZERO
          for J in range(NRPOS
            for I in range(NTARG[1]
              for K in range(KMAX
                FREQRD[I, J, K] = FNORM * FREQRD[I, J, K]
                 #  K in range(KMAX
#             Start 02 - APR - 2021 new >>>>>>>>>>>>>>>>>>>>
              for L in range(NTARG[2]
                CORR12[I, L, J] = FNORM * CORR12[I, L, J]
                 #  NTARG[2]
#             End 02 - APR - 2021 new <<<<<<<<<<<<<<<<<<<<<<<<<
               #  I in range(NTARG[1]
             #  J in range(NRPOS
          for K = KMX - 1, 2, - 1
            for J in range(NRPOS
              for I in range(NTARG[1]
                CONVOL[I] = ZERO
                 #  I in range(NTARG[1]
              for I in range(NTARG[1]
                for L = 0, I - 1
                  CONVOL[I] = CONVOL[I] + FREQRD[I - L, J, K] * FREQRD[L + 1, J, K + 1]
                   #  L = 0, I - 1
                 #  I in range(NTARG[1]
              for I in range(NTARG[1]
                FREQRD[I, J, K] = CONVOL[I]
                 #  I in range(NTARG[1]
               #  J in range(NRPOS
             #  K = KMX, 1, - 1
                    #  DUMMY block for readability: Normalize  666666

        if 8 == 8:    #  DUMMY block for readability: Prepare Output  8888888
#         Calculate single - event distributions
          for K in range(KMX
            NTARGK[K] = 1
            for I in range(NTARG[1]
              FREQSE[I, K] = FREQRD[I, 1, K]
              SPOSIN = ONE
              DUMMY = ZERO
              for J = 2, NRPOS
                DUMMY = DUMMY + 8.
                FREQSE[I, K] = FREQSE[I, K] + DUMMY * FREQRD[I, J, K]
                SPOSIN = SPOSIN + DUMMY         
                IF [RADIST[J] <= RROI: 
                  FREQTE[I, K] = FREQSE[I, K] / SPOSIN
                
                 
              FREQSE[I, K] = FREQSE[I, K] / SPOSIN
              if FREQSE[I, K] >= EPS.OR.FREQTE[I, K] >= EPS] NTARGK[K] = I            
            
          
          
          if NTARG[2] > 0:  
            if NTARGK[2] > MXTRG2: 
              print('Major problem: 2nd array dimension MXTRG2 = ', 
     &                MXTRG2, ' < max. number of 2 + clusters NTARGK[2] = ', 
     &                NTARGK[2]
              STOP
            
            for I in range(NTARGK[1]
              for L in range(NTARGK[2]
                FREQBV[I, L] = CORR12[I, L, 1]  
                SPOSIN = ONE
                DUMMY = ZERO
                for J = 2, NRPOS
                  DUMMY = DUMMY + 8.
                  FREQBV[I, L] = FREQBV[I, L] + DUMMY * CORR12[I, L, J]
                  SPOSIN = SPOSIN + DUMMY         
                   
                FREQBV[I, L] = FREQBV[I, L] / SPOSIN
                if FREQSE[I, 1] * FREQSE[L, 2] > ZERO: 
                  CORRSE[I, L] = FREQBV[I, L] / [FREQSE[I, 1] * FREQSE[L, 2]]
                else::
                  CORRSE[I, L] = FREQBV[I, L]
                                
              
            
                      
                    #  DUMMY block for readability: Prepare Output 8888888

       if 9 == 9:    #  DUMMY block for readability: OUTPUT  9999999999
#         #Write results to output file 
          .write(PREFIX[2:2], '[a]'] 'D'
          FILOUT = PREFIX + FILENM
          print('Write output to ' + FILOUT
          open(LUN, FILE = FILOUT, STATUS = 'UNKNOWN']
          .write('  *  *  *  Output from PROGRAM ROI_3D Version ' + 
     &              VDATE + ' BLINIT:' + VDATES[1] + ' TARG3D:' + VDATES[2]
     &                   + ' on ' + TSTAMP[]
          .write(' Filename: ' + FILOUT
          .write(LUN, '[a10, 10i6]'] ' NTARG[K] = ', [NTARGK[I], I in range(KMAX]
          .write(LUN, '[6[a, f8.3]]'] ' DLATC = ', DLATC, ' nm    DROI = ', DROI, 
     &                             ' nm   DSITE = ', DSITE, ' nm   DBEAM = ', 
     &                             DBEAM, ' nm'
          PARAM = 'P '
          for K in range(KMX
            if K > 1] PARAM = 'F '
            .write(PARAM[2:2], '[I1]'] K
            .write(LUN, '[1X, 2a8, 103a15]'] 'Para - ', '#Sites', 'Average', 
     &                          'Average', ['Distance / nm', J in range(NRPOS]
            .write(LUN, '[1X, 2a8, 2a15, 101f15.6]'] 'meter', ' / track', 
     &                         'total', 'inside', [RADIST[J], J in range(NRPOS]
            for I in range(NTARGK[K]
              .write(LUN, '[1X, a8, i8, 103f15.8]'] PARAM, I - 1, 
     &                                     FREQSE[I, K], FREQTE[I, K], 
     &                                    [FREQRD[I, J, K], J in range(NRPOS]
            
            .write('____________________________________________'
            .write(' *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * '
          
          .close()
          
          if NTARG[2] > 0:    #  begin 02 - APR - 2021 >>>>>>>>>>>>
            .write(PREFIX[2:2], '[a]'] 'B'
            FILOUT = PREFIX + FILENM
            print('Write output to ' + FILOUT   
            open(LUN, FILE = FILOUT, STATUS = 'UNKNOWN']
            .write(' *  *  *  Output from PROGRAM ROI_3D Version ' + 
     &              VDATE + ' BLINIT:' + VDATES[1] + ' TARG3D:' + VDATES[2]
     &                   + ' on ' + TSTAMP[]
            .write('Filename: ' + FILOUT
            .write(LUN, '[6[a, f8.3]]'] ' DLATC = ', DLATC, ' nm    DROI = ', 
     &            DROI, ' nm   DSITE = ', DSITE, ' nm   DBEAM = ', DBEAM, ' nm'
            .write('Correlations P1 and F2'
            .write(NTARGK[1], NTARGK[2]
            for I in range(NTARGK[1]
              .write(LUN, '[1X, 1000e15.8]'] [FREQBV[I, L], L in range(NTARGK[2]]
            
            .close()

            .write(PREFIX[2:2], '[a]'] 'C'
            FILOUT = PREFIX + FILENM
            print('Write output to ' + FILOUT   
            open(LUN, FILE = FILOUT, STATUS = 'UNKNOWN']
            .write(' *  *  *  Output from PROGRAM ROI_3D Version ' + 
     &              VDATE + ' BLINIT:' + VDATES[1] + ' TARG3D:' + VDATES[2]
     &                   + ' on ' + TSTAMP[]
            .write('Filename: ' + FILOUT
            .write(LUN, '[6[a, f8.3]]'] ' DLATC = ', DLATC, ' nm    DROI = ', 
     &            DROI, ' nm   DSITE = ', DSITE, ' nm   DBEAM = ', DBEAM, ' nm'
            .write('Correlations P1 and F2'
            .write(NTARGK[1], NTARGK[2]
            for I in range(NTARGK[1]
              .write(LUN, '[1X, 1000e15.8]'] [CORRSE[I, L], L in range(NTARGK[2]]
            
            .close()
             #  [NTARG[2] > 0]
                    #  DUMMY block for readability: OUTPUT  999999999

         #  IFILE in range(NFILES
      
      END PROGRAM   #  ROI_3D
# ______________________________________________________________________      

      INCLUDE 'BLINIT.f'    #  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'ROI_3D_Subs.f'    #  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'RFINIT.f'
      INCLUDE 'TSTAMP.f'