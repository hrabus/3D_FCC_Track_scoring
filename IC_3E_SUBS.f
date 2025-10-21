
C!**********************************************************************
      SUBROUTINE CLUSTR(NTRANS,ZROI,SRDATE) 
C!    ***************************************************************
C!    Determine targets with energy deposits
C!    ***************************************************************
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date
      CHARACTER VDATE*11
      PARAMETER(VDATE='04-JUL-2022')
C!    04-JUL-2022 HR: Added outout of version date to calling program
C!    20-MAY-2020 HR:  - Fixed severe bug with sorting the hit targets
C!    18-MAY-2020 HR:  - Adapted from IC_3D_Subs.f
C!    25-MAR-2020 HR:  - Removed unneeded variables
C!    20-MAR-2020 HR:  - Created from ROI_3D_Subs.f
C!    ----- Parameters: General purpose constants
      REAL*8 ONE, ZERO
      PARAMETER(ONE=1.0, ZERO=0.0)
C!    ----- Parameters for lattice orientation
      INTEGER*4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER(MAZMTH=1, MDIV=0, MTHETA=2**MDIV,  
     &           MDIR=MAZMTH*(MTHETA*(MTHETA+1))/2)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 KMAX, MTRANS, MXRPOS, MXYPOS, MXTARG
      PARAMETER(KMAX=4, MTRANS=100000, MXRPOS=41, MXTARG=1001,
     &          MXYPOS=4*MXRPOS*(MXRPOS-1)+1)
C!    ----- Local scalars
      INTEGER*4 I, IDIR, INTACT, IR, J, K
      INTEGER*4 IHOLD, NINTER
      REAL*8 SAVXYZ, SPROD
C!    ----- Local arrays
      INTEGER*4 ITARGT(MTRANS,3)
      INTEGER*4 ICSITE(KMAX)
C!           ICSITE is the histogram of absolute frequencies of  
C!                  targets with 1, 2, ... >=KMAX ionizations 
C!                  occuring for a particular orientation of the  
C!                  lattice w.r.t. to the track
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR  
C!    ------
      INTEGER*1 ITYPE(MTRANS)
      REAL*8 XYZE(MTRANS,4)
      COMMON /TRACKS/ XYZE, ITYPE
C!    ------
      REAL*8 RADIST(MXRPOS),FREQRD(MXTARG,MXRPOS,KMAX)
      COMMON /HISTOG/ RADIST, FREQRD
C!           RADIST is the vector of radial distances 
C!           FREQRD initially holds the sum, the sum of squares and the
C!                  sum of variances per track over all tracks. 
C!                  In the main program this is converted in the end to 
C!                  mean and variance over all analyzed tracks plus the
C!                  the average intra-track variance do to orientation
C!    ------
      LOGICAL DEBUG(4), SRFLAG
      COMMON /VERBOSE/ DEBUG, SRFLAG
C!    ------ Added 14-MAR-2021
      INTEGER*4 ICSIZE(MTRANS), NSITES
      COMMON /CLSTRS/ ICSIZE, NSITES
C!    ----- Scalars and arrays passed from the calling module
      INTEGER*4 NTRANS
      REAL*8  ZROI(2)
C!    ----- Scalars and arrays passed to the calling module
      CHARACTER*11 SRDATE
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      
      SRDATE=VDATE

      IF(SRFLAG) THEN
        PRINT*, 'Subroutine CLUSTR Version: ',VDATE
        SRFLAG=.FALSE.
      END IF
      IF(DEBUG(2)) PRINT*, 'CLUSTR Input',NTRANS
      
      DO K=1,KMAX ! Empty local histograms
        ICSITE(K)=0
      END DO ! I=1,KMAX ! Empty local histograms

      IF(DEBUG(2)) PRINT*, 'CLUSTR vor Main Loop NDIR=',NDIR
      IF(DEBUG(2)) PRINT*, 'CLUSTR vor Main Loop NTRANS=',NTRANS
      DO IDIR=1,NDIR ! Loop over all orientations
        IF(IDIR.GT.NDIR) GOTO 100
        IF(DEBUG(3)) PRINT*, 'CLUSTR Begin Loop IDIR',IDIR
C!      # Find target volume for all  ionizations in track      
        NINTER=0
        DO I=1,NTRANS 
          IF(XYZE(I,3).GE.ZROI(1).AND.XYZE(I,3).LE.ZROI(2)) THEN
            NINTER=NINTER+1
            DO J=1,3
              SPROD= ADJNTB(1,J,IDIR)*XYZE(I,1)
     &                +ADJNTB(2,J,IDIR)*XYZE(I,2)
     &                +ADJNTB(3,J,IDIR)*(XYZE(I,3))
              ITARGT(NINTER,J)=NINT(SPROD)
            END DO ! J=1,3
          END IF
        END DO ! I=1,NTRANS
      
        IF(DEBUG(3)) PRINT*, 'CLUSTR vor sort volume indices',NINTER
        IF(IDIR.GT.NDIR) STOP
C!      # Sort target volume indices ascending     
        DO I=1,NINTER 
          DO J=I+1,NINTER
            IF(    (ITARGT(I,1).GT.ITARGT(J,1))
     &         .OR.(ITARGT(I,1).EQ.ITARGT(J,1).AND.
     &              ITARGT(I,2).GT.ITARGT(J,2))
     &         .OR.(ITARGT(I,1).EQ.ITARGT(J,1).AND.
     &              ITARGT(I,2).EQ.ITARGT(J,2).AND.
     &              ITARGT(I,3).LE.ITARGT(J,3))) THEN
              DO IR=1,3 
                IHOLD=ITARGT(I,IR)
                ITARGT(I,IR)=ITARGT(J,IR)
                ITARGT(J,IR)=IHOLD
              END DO ! IR=1,3
              DO IR=1,4 
                SAVXYZ=XYZE(I,IR)     ! 14-MAR-2021
                XYZE(I,IR)=XYZE(J,IR)  ! 14-MAR-2021
                XYZE(J,IR)=SAVXYZ     ! 14-MAR-2021
              END DO ! IR=1,4
              IR=ITYPE(I)
              ITYPE(I)=ITYPE(J)
              ITYPE(J)=IR
            END IF
          END DO ! J=1,NINTER
        END DO ! I=1,NINTER
        
        IF(DEBUG(3)) PRINT*, 'CLUSTR vor find unique volumes',NINTER
C!      # Find unique target volumes and score ionizations   
        NSITES=1
        ICSIZE(NSITES)=ITYPE(1)
        INTACT=1
        
        DO I=2,NINTER 
          IF(    ITARGT(I,1).NE.ITARGT(I-1,1)
     &       .OR.ITARGT(I,2).NE.ITARGT(I-1,2)
     &       .OR.ITARGT(I,3).NE.ITARGT(I-1,3)) THEN
*           IF(ICSIZE(NSITES).GT.0) THEN
*             DO J=1,3 
*               XYZE(NSITES,J)=XYZE(NSITES,J)/REAL(ICSIZE(NSITES))
*             END DO
*           ELSE
*             DO J=1,3 
*               XYZE(NSITES,J)=XYZE(NSITES,J)/REAL(INTACT)
*             END DO
*           END IF
           DO J=1,3 
             XYZE(NSITES,J)=XYZE(NSITES,J)/XYZE(NSITES,4)
           END DO
            
            NSITES=NSITES+1
            ICSIZE(NSITES)=ITYPE(I)
            INTACT=1
            DO J=1,3
              ITARGT(NSITES,J)=ITARGT(I,J)
              XYZE(NSITES,J)=XYZE(I,J)
            END DO ! J=1,3
            J=4
            XYZE(NSITES,J)=XYZE(I,J)
          ELSE
            INTACT=INTACT+1
            ICSIZE(NSITES)=ICSIZE(NSITES)+ITYPE(I)
            DO J=1,3 
              XYZE(NSITES,J)=XYZE(NSITES,J)+XYZE(I,J)*XYZE(I,4)         ! Weight by energy deposition
            END DO
            DO J=4,4 
              XYZE(NSITES,J)=XYZE(NSITES,J)+XYZE(I,J)
            END DO
          END IF
        END DO ! I=2,NINTER

      END DO ! IDIR=1,NDIR ! Loop over all orientations
 100  CONTINUE

      IF(DEBUG(2)) PRINT*, 'CLUSTR vor EXIT'
      
      END SUBROUTINE CLUSTR
C!______________________________________________________________________
