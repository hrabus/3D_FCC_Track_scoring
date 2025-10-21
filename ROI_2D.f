      PROGRAM ROI_2D
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='23-MAY-2021',VINPUT=210413.1643d0) ! ###########
C!    #################################################################
C!    23-MAY-2021 HR: 
C!      - Modified call to subroutines to get their version date
C!    02-APR-2021 HR:
C!      - Created from ROI_3D
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN
      REAL*8 ONE, ZERO
      PARAMETER(LUN=11, ONE=1.0, ZERO=0.0)
C!    ----- Parameters for lattice orientation
      INTEGER*4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER(MAZMTH=1, MDIV=0, MTHETA=2**MDIV,  
     &           MDIR=MAZMTH*(MTHETA*(MTHETA+1))/2)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 KMAX, MIONIZ, MXRPOS, MXTARG, MXTRG2, MXYPOS, MXZROI
      PARAMETER(KMAX=9, MIONIZ=100000, MXRPOS=101, MXTARG=201,
     &          MXTRG2=49, MXYPOS=4*MXRPOS*(MXRPOS-1)+1, MXZROI=100)
C!    ----- Input parameters for geometry
      REAL*8 DBEAM  ! Diameter of beam in nm
      REAL*8 DLATC  ! Cell lattice constant 
      REAL*8 DROI   ! Diameter of region of interest (cell nucleus size)
      REAL*8 DSITE  ! Diameter of spherical target 
      REAL*8 DZROI  ! Increment in position of region of interest 
C!    ----- Local scalars
      CHARACTER VDATES(2)*11, SRDATE
      CHARACTER FILENM*80, FILOUT*85, HEADER*80, PARAM*2, PREFIX*10
      INTEGER*4 I, IFILE, IT, ITA, J, K, L, NAZMTH, NDIV, NFILES, 
     &          NFORMT, NZROI, NPHI, NTRACS
      INTEGER*4 NIONIZ, NRPOS, NXYPOS
      REAL*8 COSPHI, DELTAR, DUMB1, DUMB2, DUMB3, FNORM, PIBY4, RROI,
     &       SINPHI, SPOSIN
C!    ----- Local arrays
      INTEGER*4 NTARG(2)
      REAL*8 CORRSE(MXTRG2)
      REAL*8 FREQSE(MXTARG), ZROIC(MXZROI)
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR 
C!    ------
      INTEGER*4 IRAD(MXYPOS)
      REAL*8 RROI2 ! Square of ROI radius 
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
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program ROI_2D Version: '//VDATE
      SRFLAG=.TRUE. ! Invoke message on Version date from Subroutine TARG2D 

      IF(1.EQ.1) THEN ! DUMMY block for readability: Read input 1111111
        DEBUG(1)=.FALSE. ! Main program
        DEBUG(2)=.FALSE. ! TARG2D main sections
        DEBUG(3)=.FALSE. ! TARG2D IDIR loop
        DEBUG(4)=.FALSE. ! TARG2D IPOS loop details
        
        PIBY4=ATAN(ONE)
        PREFIX='2D__0.0nm_'
        
        NRPOS=MXRPOS
*        PRINT*,'Enter NRPOS'
*        READ(*,*) NRPOS
        NTARG(1)=MXTARG
        FILENM='SPRE_test.dat'
        FILOUT='TestRes.dat'

C!      Read debug options
        PRINT*, 'Enter debug options file name' 
        READ(*,*) FILENM        
        OPEN(LUN,FILE=FILENM,STATUS='OLD',ERR=10)
        DO I=1,4
          READ(LUN,*) DEBUG(I)
        END DO ! I=1,4
        CLOSE(LUN)
  10    CONTINUE
  
C!      Read geometry parameters 
        PRINT*, 'Enter site diameter in nm'
        READ(*,*) DSITE
        DLATC=DSITE*SQRT(2.)*EXP(LOG(PIBY4/3.)/3.) 
        IF(DSITE.LT.10.0) THEN 
          WRITE(PREFIX(5:7),'(f3.1)') DSITE
        ELSE        
          WRITE(PREFIX(4:7),'(f4.1)') DSITE
        END IF
        
        PRINT*, 'Enter ROI diameter in nm'
        READ(*,*) DROI
        RROI=DROI/2.
        RROI2=RROI*RROI
        DELTAR=DROI/40.
              
        PRINT*, 'Enter beam diameter in nm <= ',5.*DROI
        READ(*,*) DBEAM
        NRPOS=1+NINT(DBEAM/2./DELTAR)
        IF(NRPOS.GT.MXRPOS) THEN
          NRPOS=MXRPOS
          PRINT*, 'Beam diameter ',DBEAM,' nm is too large. <<<<<<<<<<<'
          DBEAM=2.*DELTAR*REAL(MXRPOS-1)
          PRINT*, '>>> maximum possible value ',DBEAM,' is used.'
        END IF

        PRINT*, 'Enter maximum number of targets in histogram'
        READ(*,*) NTARG(1), NTARG(2)
        IF(NTARG(1).GT.MXTARG.OR.NTARG(1).LT.1) NTARG(1)=MXTARG
        IF(NTARG(2).GT.MXTRG2.OR.NTARG(2).EQ.1) NTARG(2)=MXTRG2
        PRINT*, 'Number of files to process'
        READ(*,*) NFILES
        PRINT*, 'Enter flag for data file structure 1=SPRE  2=Heidi'
        READ(*,*) NFORMT
        PRINT*, 'Enter number of regions of interest along track'
        READ(*,*) NZROI
        IF (NZROI.EQ.1) THEN
          PRINT*, 'Enter z position and length of region of interest'
          READ(*,*) ZROIC(1), DZROI
        ELSE
          PRINT*, 'Enter z position of first region of interest and '//
     &            'increment in z position of regions of interest'
          READ(*,*) ZROIC(1), DZROI
          DO I=2,NZROI
            ZROIC(I)=ZROIC(I-1)+DZROI
          END DO
        END IF
        
      END IF          ! DUMMY block for readability: Read input 1111111

      IF(2.EQ.2) THEN ! DUMMY block for readability: Initialize 2222222
*        IF(DEBUG(1)) PRINT*, 'Hier v'
        NAZMTH=MAZMTH
        NDIV=MDIV
        IF(DEBUG(1)) PRINT*, 'Hier vor CALL BLINIT'
        CALL BLINIT(NAZMTH,NDIV,DLATC,SRDATE) ! Init reziprocal lattice
        VDATES(1)=SRDATE
        IF(DEBUG(1)) PRINT*, 'Hier nach CALL BLINIT', NRPOS, NTARG
        
        DO I=1,NRPOS
C!        Define radial offsets of track w.r.t. ROI center 
          RADIST(I)=REAL(I-1)*DELTAR
C!        x&y positions of track w.r.t. ROI center (piecake method)
          IF(I.EQ.1) THEN
            NXYPOS=1
            NPHI=0
            IRAD(NXYPOS)=1
            XROI(NXYPOS)=ZERO
            YROI(NXYPOS)=ZERO
          ELSE
            NPHI=NPHI+8
            NXYPOS=NXYPOS+1
            IRAD(NXYPOS)=I
            XROI(NXYPOS)=RADIST(I)
            YROI(NXYPOS)=ZERO
            COSPHI=COS(PIBY4/REAL(I-1))
            SINPHI=SIN(PIBY4/REAL(I-1))
            DO J=2, NPHI
              NXYPOS=NXYPOS+1
              IRAD(NXYPOS)=I
              XROI(NXYPOS)=XROI(NXYPOS-1)*COSPHI-YROI(NXYPOS-1)*SINPHI
              YROI(NXYPOS)=XROI(NXYPOS-1)*SINPHI+YROI(NXYPOS-1)*COSPHI
            END DO
          END IF
        END DO ! DO I=1,NRPOS
      END IF          ! DUMMY block for readability: Initialize 2222222

      DO IFILE=1,NFILES      
        PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM

        IF(4.EQ.4) THEN ! DUMMY block for readability: Main Loop  444444
          IF(DEBUG(1)) PRINT*, 'Hier beginnt Block 4'
          DO J=1,NRPOS
C!          Init counters       
            DO I=1,NTARG(1)
              DO K=1,KMAX
                FREQRD(I,J,K)=ZERO
              END DO
            END DO  
C!          Begin 02-APR-2021 >>>>>>>>>>>>>>>>>     
            DO I=1,NTARG(1)
              DO L=1,NTARG(2)
                CORR12(I,L,J)=ZERO
              END DO
            END DO  
C!          END 02-APR-2021 <<<<<<<<<<<<<<<<<<  
          END DO ! DO J=1,NRPOS
          NIONIZ=1
          NTRACS=1
          
          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          IF (NFORMT.EQ.1) THEN 
            DO I=1,4
              READ(LUN,*) HEADER
            END DO
            READ(LUN,*) ITA,XYZ(NIONIZ,3),DUMB1,DUMB2,DUMB3,
     &                  XYZ(NIONIZ,1),XYZ(NIONIZ,2)
          ELSE
            READ(LUN,*) HEADER
            READ(LUN,*) ITA,DUMB1,DUMB2,XYZ(NIONIZ,3),
     &                   XYZ(NIONIZ,1),XYZ(NIONIZ,2)
          END IF
          
          NIONIZ=NIONIZ+1
          IF(DEBUG(1)) PRINT*, 'Hier 4'
  100     CONTINUE ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          IF (NFORMT.EQ.1) THEN 
            READ(LUN,*,END=200,ERR=200) IT,XYZ(NIONIZ,3),DUMB1,DUMB2,
     &                           DUMB3,XYZ(NIONIZ,1),XYZ(NIONIZ,2)
          ELSE
            READ(LUN,*,END=200,ERR=200) IT,DUMB1,DUMB2,XYZ(NIONIZ,3),
     &                   XYZ(NIONIZ,1),XYZ(NIONIZ,2)
          END IF
          IF(IT.EQ.ITA) THEN
            NIONIZ=NIONIZ+1
          ELSE ! Entry from a new track was read
            NIONIZ=NIONIZ-1 ! Decrease to number of ionizations in track
            IF(DEBUG(1)) PRINT*, 'Hier vor CALL TARG2D'
            IF(MOD(NTRACS,10).EQ.0) PRINT*, NTRACS
            DO J=1,NZROI 
              ZROI(2)=ZROIC(J)
              ZROI(1)=ZROI(2)-DROI/2.-DLATC/2.
              ZROI(3)=ZROI(2)+DROI/2.+DLATC/2.
              CALL TARG2D(NIONIZ,NRPOS,NTARG,NXYPOS,SRDATE)
              
            END DO ! J=1,NZROI

            IF(DEBUG(1)) PRINT*, 'Hier nach CALL TARG2D'
            DO J=1,3
              XYZ(1,J)=XYZ(NIONIZ+1,J) ! First ionization in new track
            END DO
            ITA=IT ! Remember new track number
            NTRACS=NTRACS+1
            NIONIZ=2 ! There was already one ionization
          END IF
          GOTO 100 ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  200     CONTINUE
            DO J=1,NZROI 
              ZROI(2)=ZROIC(J) 
              ZROI(1)=ZROI(2)-DZROI/2.-DLATC/2.
              ZROI(3)=ZROI(2)+DZROI/2.+DLATC/2.
              CALL TARG2D(NIONIZ,NRPOS,NTARG,NXYPOS,SRDATE) 
            END DO ! J=1,NZROI
            VDATES(2)=SRDATE
          CLOSE(LUN)
        END IF          ! DUMMY block for readability: Main Loop  4444444

        IF(6.EQ.6) THEN ! DUMMY block for readability: Normalize  666666
C!        #Normalization
          FNORM=ONE/REAL(NTRACS)
          FNORM=ONE/REAL(NTRACS)/REAL(NZROI) 
*          SUMCNT=ZERO
          DO J=1,NRPOS
            DO I=1,NTARG(1)
              DO K=1,KMAX
                FREQRD(I,J,K)=FNORM*FREQRD(I,J,K)
              END DO ! K=1,KMAX
C!            Start 02-APR-2021 new >>>>>>>>>>>>>>>>>>>>
              DO L=1,NTARG(2)
                CORR12(I,L,J)=FNORM*CORR12(I,L,J)
              END DO ! NTARG(2)
C!            End 02-APR-2021 new <<<<<<<<<<<<<<<<<<<<<<<<<
            END DO ! I=1,NTARG(1)
          END DO ! J=1,NRPOS
          DO K=KMAX-1,2,-1
            DO J=1,NRPOS
              DO I=1,NTARG(1)
                FREQSE(I)=ZERO
              END DO ! I=1,NTARG(1)
              DO I=1,NTARG(1)
                DO L=0, I-1
                  FREQSE(I)=FREQSE(I)+FREQRD(I-L,J,K)*FREQRD(L+1,J,K+1)
                END DO ! L=0, I-1
              END DO ! I=1,NTARG(1)
              DO I=1,NTARG(1)
                FREQRD(I,J,K)=FREQSE(I)
              END DO ! I=1,NTARG(1)
            END DO ! J=1,NRPOS
          END DO ! K=KMAX-1,1,-1
        END IF          ! DUMMY block for readability: Normalize  666666

        IF(9.EQ.9) THEN ! DUMMY block for readability: OUTPUT  9999999
C!        #Write results to output file 
          FILOUT=PREFIX//FILENM
          PRINT*, 'Write output to '//FILOUT
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ROI_2D Version '//
     &              VDATE//' BLINIT:'//VDATES(1)//' TARG3D:'//VDATES(2)
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,'(6(a,f8.3))') 'DLATC=',DLATC,' nm    DROI=',DROI,
     &                             ' nm   DSITE=',DSITE,' nm   DBEAM=',
     &                             DBEAM,' nm'
          PARAM='P '
          DO K=1,KMAX
            IF(K.GT.1) PARAM='F '
            WRITE(PARAM(2:2),'(I1)') K
            WRITE(LUN,'(1X,2a8,103a15)') 'Para-', '#Sites', 'Average',
     &                          'Average',('Distance/nm',J=1,NRPOS)
            WRITE(LUN,'(1X,2a8,2a15,101f15.6)') 'meter', '/track',
     &                         'total','inside',(RADIST(J),J=1,NRPOS)
            DO I=1,NTARG(1)
              FREQSE(1)=FREQRD(I,1,K)
              SPOSIN=ONE
              DUMB3=ZERO
              DO J=2,NRPOS
                DUMB3=DUMB3+8.
                FREQSE(1)=FREQSE(1)+DUMB3*FREQRD(I,J,K)
                SPOSIN=SPOSIN+DUMB3         
                IF (RADIST(J).LE.RROI) THEN
                  FREQSE(2)=FREQSE(1)/SPOSIN
                END IF
              END DO   
              FREQSE(1)=FREQSE(1)/SPOSIN    
              WRITE(LUN,'(1X,a8,i8,103f15.6)') PARAM, I-1,FREQSE(1),
     &                         FREQSE(2),(FREQRD(I,J,K),J=1,NRPOS)
            END DO
            WRITE(LUN,*) '____________________________________________'
            WRITE(LUN,*) '********************************************'
          END DO
          CLOSE(LUN)
          
          IF(NTARG(2).GT.0) THEN ! begin 02-APR-2021 >>>>>>>>>>>>
            WRITE(PREFIX(2:2),'(a)') 'C'
            FILOUT=PREFIX//FILENM
            PRINT*, 'Write output to '//FILOUT   
            OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
            WRITE(LUN,*) '*** Output from PROGRAM ROI_2D Version '//
     &              VDATE//' BLINIT:'//VDATES(1)//' TARG3D:'//VDATES(2)
            WRITE(LUN,*) 'Filename: '//FILOUT
            WRITE(LUN,'(6(a,f8.3))') 'DLATC=',DLATC,' nm    DROI=',DROI,
     &                             ' nm   DSITE=',DSITE,' nm   DBEAM=',
     &                             DBEAM,' nm'
            WRITE(LUN,*) 'Correlations P1 and F2'
            WRITE(LUN,*) NTARG
            DO I=1,NTARG(1)
              DO L=1,NTARG(2)
                CORRSE(L)=CORR12(I,L,1)  
                SPOSIN=ONE
                DUMB3=ZERO
                DO J=2,NRPOS
                  DUMB3=DUMB3+8.
                  CORRSE(L)=CORRSE(L)+DUMB3*CORR12(I,L,J)
                  SPOSIN=SPOSIN+DUMB3         
                END DO   
                CORRSE(L)=CORRSE(L)/SPOSIN    
              END DO
              WRITE(LUN,'(1X,103e15.6)') (CORRSE(L),L=1,NTARG(2))
            END DO
            CLOSE(LUN)
          END IF ! (NTARG(2).GT.0)
        END IF          ! DUMMY block for readability: OUTPUT     9999999

      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! ROI_2D
C!______________________________________________________________________      

      INCLUDE 'BLINIT.f'       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'ROI_2D_Subs.f'  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
