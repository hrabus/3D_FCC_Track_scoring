
C!**********************************************************************
      SUBROUTINE TARG2D(NIONIZ,NRPOS,NTARG,NXYPOS,SRDATE) 
C!    ***************************************************************
C!    Determine targets with ionizations
C!    ***************************************************************
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date
      CHARACTER VDATE*11
      PARAMETER(VDATE='23-MAY-2021')
C!    23-MAY-2021 HR: Added outout of version date to calling program
C!    20-MAY-2020 HR: Fixed severe bug with sorting the hit targets
C!    02-APR-2021 HR: Created and adapted from ROI_3D_subs.f
C!    ----- Parameters: General purpose constants
      REAL*8 ONE, ZERO
      PARAMETER(ONE=1.0, ZERO=0.0)
C!    ----- Parameters for lattice orientation
      INTEGER*4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER(MAZMTH=1, MDIV=0, MTHETA=2**MDIV,  
     &           MDIR=MAZMTH*(MTHETA*(MTHETA+1))/2)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 KMAX, MIONIZ, MXRPOS, MXTARG, MXTRG2, MXYPOS
      PARAMETER(KMAX=9, MIONIZ=100000, MXRPOS=101, MXTARG=201,
     &          MXTRG2=49, MXYPOS=4*MXRPOS*(MXRPOS-1)+1)
C!    ----- Local scalars
      INTEGER*4 I, IDIR, IPOS, IR, J, K
      INTEGER*4 IHOLD, NCOUNT, NIONIS
      REAL*8 DISTSQ, WEIGHT, SAVXYZ, SPROD
C!    ----- Local arrays
      INTEGER*4 ICROI(3), ITARGT(MIONIZ,3)
      INTEGER*4 ICSITE(KMAX)
      REAL*8 CORREL(MXTARG,MXTRG2,MXRPOS)
      REAL*8 FTGICS(MXTARG,MXRPOS,KMAX)
C!           ICSITE is the histogram of absolute frequencies of  
C!                  targets with 1, 2, ... >=KMAX ionizations 
C!                  occuring for a particular orientation of the  
C!                  lattice w.r.t. to the track
C!           FTGICS initially is the sum of the absolute frequencies 
C!                  of targets holding a certain ICS for an orientation
C!                  in the particular track. 
C!                  In the end it is normalized and trasferred to global 
C!                  counters              
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR  
C!    ------
      INTEGER*4 IRAD(MXYPOS)
      REAL*8 RROI2 ! Integer of Square of ROI radius 
      REAL*8 XROI(MXYPOS), YROI(MXYPOS), ZROI(3)
      COMMON /ROICTR/  XROI, YROI, ZROI, RROI2, IRAD
C!    ------
      REAL*8 XYZ(MIONIZ,3)
      COMMON /TRACKS/ XYZ 
C!    ------
      REAL*8 CORR12(MXTARG,MXTRG2,MXRPOS)
      REAL*8 RADIST(MXRPOS),FREQRD(MXTARG,MXRPOS,KMAX)
C!*      COMMON /HISTOG/ RADIST, FREQRD
      COMMON /HISTOG/ RADIST, FREQRD, CORR12
C!           RADIST is the vector of radial distances 
C!           FREQRD initially holds the sum, the sum of squares and the
C!                  sum of variances per track over all tracks. 
C!                  In the main program this is converted in the end to 
C!                  mean and variance over all analyzed tracks plus the
C!                  the average intra-track variance do to orientation
C!    ------
      LOGICAL DEBUG(4), SRFLAG
      COMMON /VERBOSE/ DEBUG, SRFLAG
C!    ------ 
      INTEGER*4 ICSIZE(MIONIZ),NSITES
      COMMON /CLSTRS/ ICSIZE, NSITES
C!    ----- Scalars and arrays passed from the calling module
      INTEGER*4 NIONIZ, NRPOS, NTARG(2), NXYPOS
C!    ----- Scalars and arrays passed to the calling module
      CHARACTER*11 SRDATE
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      
      SRDATE=VDATE
      IF(SRFLAG) THEN
        PRINT*, 'Subroutine TARG2D Version: ',VDATE
        SRFLAG=.FALSE.
      END IF
      IF(DEBUG(2)) PRINT*, 'TARG2D Input',NIONIZ,NRPOS,NTARG,NXYPOS
      
C!    Init local histograms
      DO K=1,KMAX                   
        ICSITE(K)=0
        DO J=1,NRPOS
          DO I=1,NTARG(1)
            FTGICS(I,J,K)=ZERO
          END DO ! I=1,NTARG(1)
        END DO ! J=1,NRPOS
      END DO ! I=1,KMAX 
      
C!    Begin 02-APR-2021 >>>>>>>>>
      DO J=1,NRPOS
        DO I=1,NTARG(1)
          DO K=1,NTARG(2)
            CORREL(I,K,J)=ZERO
          END DO ! K=1,NTARG(2)
        END DO ! I=1,NTARG(1)
      END DO ! J=1,NRPOS
C!    End 02-APR-2021 <<<<<<<<<<
C!    End init local histograms
      IF(DEBUG(2)) PRINT*, 'TARG2D vor Main Loop NDIR=',NDIR
      IF(DEBUG(2)) PRINT*, 'TARG2D vor Main Loop NIONIZ=',NIONIZ
      DO IDIR=1,NDIR ! Loop over all orientations
        IF(IDIR.GT.NDIR) GOTO 100
        IF(DEBUG(3)) PRINT*, 'TARG2D Begin Loop IDIR',IDIR
C!      # Find target volume for all  ionizations in track      
        NIONIS=0
        DO I=1,NIONIZ 
          IF(XYZ(I,3).GE.ZROI(1).AND.XYZ(I,3).LE.ZROI(3)) THEN
            NIONIS=NIONIS+1
            DO J=1,3
              SPROD= ADJNTB(1,J,IDIR)*XYZ(I,1)
     &                +ADJNTB(2,J,IDIR)*XYZ(I,2)
     &                +ADJNTB(3,J,IDIR)*(XYZ(I,3)-ZROI(2))
              ITARGT(NIONIS,J)=NINT(SPROD)
            END DO ! J=1,3
          END IF
        END DO ! I=1,NIONIZ
      
        IF(DEBUG(3)) PRINT*, 'TARG2D vor sort volume indices',NIONIS
        IF(IDIR.GT.NDIR) STOP
C!      # Sort target volume indices ascending     
        DO I=1,NIONIS 
          DO J=I+1,NIONIS
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
                SAVXYZ=XYZ(I,IR)
                XYZ(I,IR)=XYZ(J,IR)
                XYZ(J,IR)=SAVXYZ
              END DO ! IR=1,3
            END IF
          END DO ! J=1,NIONIS
        END DO ! I=1,NIONIS
        
        IF(DEBUG(3)) PRINT*, 'TARG2D vor find unique volumes',NIONIS
C!      # Find unique target volumes and score ionizations   
        NSITES=1
        ICSIZE(NSITES)=1
        DO I=2,NIONIS 
          IF(    ITARGT(I,1).NE.ITARGT(I-1,1)
     &       .OR.ITARGT(I,2).NE.ITARGT(I-1,2)
     &       .OR.ITARGT(I,3).NE.ITARGT(I-1,3)) THEN
            DO J=1,3
              XYZ(NSITES,J)=XYZ(NSITES,J)/REAL(ICSIZE(NSITES))
            END DO
            NSITES=NSITES+1
            ICSIZE(NSITES)=1
            DO J=1,3
              ITARGT(NSITES,J)=ITARGT(I,J)
            END DO ! J=1,3
          ELSE
            ICSIZE(NSITES)=ICSIZE(NSITES)+1
            DO J=1,3
              XYZ(NSITES,J)=XYZ(NSITES,J)+XYZ(I,J)
            END DO
          END IF
        END DO ! I=2,NIONIS

        IF(DEBUG(3)) PRINT*, 'TARG2D Begin Loop IPOS',NXYPOS
        DO IPOS=1,NXYPOS ! Loop over all track positions 
C!        Calculate cell indices of ROI center
          DO J=1,3
            SPROD= ADJNTB(1,J,IDIR)*XROI(IPOS)
     &                +ADJNTB(2,J,IDIR)*YROI(IPOS)
*     &                +ADJNTB(3,J,IDIR)*ZROI(2)
            ICROI(J)=NINT(SPROD)
          END DO ! J=1,3
          
          IF(DEBUG(4)) PRINT*, 'TARG2D vor Score hit targets',ICROI,
     &                      IPOS, XROI(IPOS)
C!        # Score hit targets
          DO I=1,KMAX ! Zero local counter
            ICSITE(I)=0
          END DO ! I=1,KMAX ! Zero local counter
          
          DO I=1,NSITES ! Count hit targets
            DISTSQ=XYZ(I,1)*XYZ(I,1)+XYZ(I,2)*XYZ(I,2)
            NCOUNT=0
            IF(DISTSQ.LE.RROI2) THEN ! count if inside ROI
              NCOUNT=ICSIZE(I) !  Ionization cluster size
              IF(NCOUNT.GT.KMAX) NCOUNT=KMAX
              ICSITE(NCOUNT)=ICSITE(NCOUNT)+1
*              IF(DEBUG(5)) PRINT*, 'Count hit targets',IDIST,NCOUNT,IPOS
            END IF
            
          END DO ! I=1,NSITES ! Count hit targets
          
          IF(DEBUG(3)) PRINT*, 'TARG2D vor add to sum arrays',NDIR,IPOS
C!      # Add this histogram to sum arrays
          IR=IRAD(IPOS)
          IF(DEBUG(3)) PRINT*, 'TARG2D vor add to sum arrays',IR,ICSITE
          DO K=1,KMAX ! Update global counters 
            NCOUNT=ICSITE(K)+1
            IF(NCOUNT.GT.NTARG(1)) NCOUNT=NTARG(1)
            IF(NCOUNT.GT.0) FTGICS(NCOUNT,IR,K)=FTGICS(NCOUNT,IR,K)+ONE
          END DO
C!        Begin 02-APR-2021 >>>>>>>>>
          IF(NTARG(2).GT.0) THEN
            ICSITE(1)=ICSITE(1)+1 
            IF(ICSITE(1).GT.NTARG(1)) ICSITE(1)=NTARG(1)
            NCOUNT=1
            DO K=2,KMAX
              NCOUNT=NCOUNT+ICSITE(K)
            END DO
            IF(NCOUNT.GT.NTARG(2)) NCOUNT=NTARG(2)
            CORREL(ICSITE(1),NCOUNT,IR)=
     &                       CORREL(ICSITE(1),NCOUNT,IR)+ONE
          END IF
C!        End 02-APR-2021 <<<<<<<<<<
          IF(DEBUG(3)) PRINT*, 'TARG2D nach sum arrays',NDIR,IDIR
        END DO ! IPOS=1,NXYPOS ! Loop over all track positions
      END DO ! IDIR=1,NDIR ! Loop over all orientations
 100  CONTINUE

      
C!    # Update global counters 
      DO I=1,NRPOS 
        IF(I.EQ.1) THEN
          WEIGHT=ONE/REAL(NDIR)
        ELSE
          WEIGHT=ONE/REAL(8*I-8)/REAL(NDIR)
        END IF
        DO J=1,NTARG(1)
          DO K=1,KMAX
            FREQRD(J,I,K)=FREQRD(J,I,K)+WEIGHT*FTGICS(J,I,K)
          END DO  ! K=1,KMAX
C!        Begin 02-APR-2021 >>>>>>>>>
          IF(NTARG(2).GT.0) THEN
            DO K=1,NTARG(2)
              CORR12(J,K,I)=CORR12(J,K,I)+WEIGHT*CORREL(J,K,I)        
            END DO ! K=1,NTARG(2)
          END IF ! (NTARG(2).GT.0)
C!        End 02-APR-2021 <<<<<<<<<<
        END DO  ! J=1,NTARG(1)
      END DO ! I=1,NRPOS
      
      IF(DEBUG(2)) PRINT*, 'TARG2D vor EXIT'
      
      END SUBROUTINE TARG2D
C!______________________________________________________________________
