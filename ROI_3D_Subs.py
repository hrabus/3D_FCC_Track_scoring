#  Wahrscheinlich geht es doch nicht, auf IC_3D_Subs zurückzugreiefen
#

import common_variables as cv
from IC_3D_Subs import CLUSTR

Version_Date = '´02-AUG-2022'


#   NOTE: Should use call to CLUSTR to simplify
#  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
def TARG3D(number_of_transfer_points, region_of_interest_infos): 
#      *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
#     Determine targets with ionizations
#      *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
# >>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      
    limits_region_of_interest = [region_of_interest_infos[0][0],region_of_interest_infos[0][2]]
    center_region_of_interest = region_of_interest_infos[0][1]
    XROI = region_of_interest_infos[1]
    YROI = region_of_interest_infos[2]
    NXYPOS = len(XROI)
    ICROI = [0, 1, 2]
#ITARGT[MIONIZ, 3]
      # ICSITE[KMAX]
      # CORREL[MXTARG, MXTRG2, MXRPOS]
      # FTGICS[MXTARG, MXRPOS, KMAX]
#            ICSITE is the histogram of absolute frequencies of  
#                   targets with 1, 2, ... >= KMAX ionizations 
#                   occuring for a particular orientation of the  
#                   lattice w.r.t. to the track
#            FTGICS initially is the sum of the absolute frequencies 
#                   of targets holding a certain ICS for an orientation
#                   in the particular track. 
#                   In the end it is normalized and trasferred to global 
#                   counters              
#     - - - - - Global variables
#     - - - - - - 
      # CORR12[MXTARG, MXTRG2, MXRPOS]
      # RADIST[MXRPOS], FREQRD[MXTARG, MXRPOS, KMAX]
#  *       COMMON / HISTOG / RADIST, FREQRD
#     COMMON / HISTOG / RADIST, FREQRD, CORR12
#            RADIST is the vector of radial distances 
#            FREQRD initially holds the sum, the sum of squares and the
#                   sum of variances per track over all tracks. 
#                   In the main program this is converted in the end to 
#                   mean and variance over all analyzed tracks plus the
#                   the average intra - track variance for to orientation
## <<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      
    if cv.SRFLAG: 
        print('Subroutine TARG3D Version: '+Version_Date)
    cv.SRFLAG = False
    
    if cv.DEBUG[2]:
        print('TARG3D Input', number_of_transfer_points, NRPOS, NTARG, NXYPOS)
    
    #     Init local histograms
    # FTGICS[I, J, K] = ZERO
    # CORREL[I, K, J] = ZERO
    
    #     End init local histograms
    
    if cv.DEBUG[2]:
        print('TARG3D vor Main Loop cv.n_directions = '+str(cv.n_directions))
        print('TARG3D vor Main Loop number_of_transfer_points = '+str(number_of_transfer_points))
    
    for IDIR in range(cv.n_directions):  #  Loop over all orientations
        if cv.DEBUG[3]:
             print('TARG3D Begin Loop IDIR'+str(IDIR))
            
        # call to SUBROUTINE CLUSTR which detects clusters
        
        if cv.edep_flag:
         track_results = CLUSTR([number_of_transfer_points, cv.xyz, cv.ionizations, cv.edep],
                                limits_region_of_interest)
        else:
         track_results = CLUSTR([number_of_transfer_points, cv.xyz, cv.ionizations],
                                limits_region_of_interest)
        
        number_of_sites = track_results[0]
        xyz = track_results[1]
        ionizations = track_results[2]
        edep_flag = len(track_results) >= 3
        if edep_flag:
            edep = track_results[3]
        
        NIONIS = 0
        
        if cv.DEBUG[3]:
            print('TARG3D Begin Loop IPOS '+str(NXYPOS))

        for IPOS in range(NXYPOS):  #  Loop over all track positions
              #   Calculate cell indices of ROI center
              for J in range(3):
                    SPROD = ADJNTB[1, J, IDIR] * XROI[IPOS] + ADJNTB[2, J, IDIR] * YROI[IPOS]
                             #+ ADJNTB[3, J, IDIR] * center_region_of_interest)
                    ICROI[J] = round(SPROD)

              if cv.DEBUG[4]:
                  print('TARG3D vor Score hit targets {0:d} {1:d} {2:5.2f} '.format(ICROI, IPOS, XROI[IPOS]))
#         # Score hit targets
          for I in range(KMAX   #  Zero local counter
            ICSITE[I] = 0
             #  I in range(KMAX   #  Zero local counter
          
          for I in range(NSITES   #  Count hit targets
            IDIST = 0
            for J in range(3
              for K = J, 3
                IDIST = IDIST + [ITARGT[I, J] - ICROI[J]]
     &                      * [ITARGT[I, K] - ICROI[K]]
              
            
            NCOUNT = 0
            if ABS[IDIST] <= NRROI2:    #  count if inside ROI
              NCOUNT = ICSIZE[I]   #   Ionization cluster size
              if NCOUNT > KMAX] NCOUNT = KMAX
              ICSITE[NCOUNT] = ICSITE[NCOUNT] + 1
 *               if cv.DEBUG[5]] print('Count hit targets', IDIST, NCOUNT, IPOS
            
            
             #  I in range(NSITES   #  Count hit targets
          
          if cv.DEBUG[3]] print('TARG3D vor add to sum arrays', cv.n_directions, IPOS
#       # Add this histogram to sum arrays
          IR = IRAD[IPOS]
          if cv.DEBUG[3]] print('TARG3D vor add to sum arrays', IR, ICSITE
          for K in range(KMAX   #  Update global counters 
            NCOUNT = ICSITE[K] + 1
            if NCOUNT > NTARG[1]] NCOUNT = NTARG[1]
            if NCOUNT > 0] FTGICS[NCOUNT, IR, K] = FTGICS[NCOUNT, IR, K] + 1.
          
#         Begin 02 - APR - 2021 >>>>>>>>>
          if NTARG[2] > 0: 
            ICSITE[1] = ICSITE[1] + 1 
            if ICSITE[1] > NTARG[1]] ICSITE[1] = NTARG[1]
            NCOUNT = 1
            for K = 2, KMAX
              NCOUNT = NCOUNT + ICSITE[K]
            
        if NCOUNT > NTARG[2]:
                NCOUNT = NTARG[2]
            CORREL[ICSITE[1], NCOUNT, IR] = CORREL[ICSITE[1], NCOUNT, IR] + 1.
          
#         End 02 - APR - 2021 <<<<<<<<<<
          
          if cv.DEBUG[3]] print('TARG3D nach sum arrays', cv.n_directions, IDIR
           #  IPOS in range(NXYPOS   #  Loop over all track positions
         #  IDIR in range(cv.n_directions   #  Loop over all orientations
         
    # Update global counters 
    for I in range(NRPOS):
        if I == 1: 
          WEIGHT = 1. / (cv.n_directions]
        else:
          WEIGHT = 1. / (8 * I - 8) / float(cv.n_directions)
    
        for J in range(NTARG[1]
            for K in range(KMAX):
                    FREQRD[J, I, K] = FREQRD[J, I, K] + WEIGHT * FTGICS[J, I, K]
            if NTARG[2] > 0: 
                for K in range(NTARG[2]
                CORR12[J, K, I] = CORR12[J, K, I] + WEIGHT * CORREL[J, K, I]        

      
      if cv.DEBUG[2]:
          print('TARG3D vor EXIT')
      
#      END SUBROUTINE TARG3D
# ______________________________________________________________________
