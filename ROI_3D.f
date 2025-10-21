      PROGRAM ROI_3D
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='05-JUN-2021',VINPUT=210605.16d0) ! ############
C!    #################################################################
C!    05-JUN-2021 HR:
C!      - Added timestap to output files 
C!      - Reworked identification of input file structure
C!    23-MAY-2021 HR: 
C!      - Modified call to subroutines to get their version date
C!    17-APR-2021 HR:
C!      - Increased MXTRG2 owing to data set with larger values
C!    13-APR-2021 HR:
C!      - Fixed soft bug: Results were overwritten due to same filename
C!      - Added reduction of output to non-zero values
C!    02-APR-2021 HR:
C!      - Added output of correlation matrix
C!      - Detected and fixed bug in section with convolution
C!    28-MAR-2021 HR:
C!      - Changed normalization of traversing tracks
C!    20-MAR-2021 HR:
C!      - Added input of options file name
C!      - Added discrimination between site size and lattice constant
C!    14-MAR-2021 HR:
C!      - Fixed bug with treatment of last track
C!      - Added common block CLSTRS for info on cluster positions  
C!    16-OCT-2020 HR:
C!      - Added readin of geometry parameters via options file
C!      - Added readin of multiple file names to process with options
C!      - Added processing of multiple ROIs in track data set
C!    15-MAR-2020 HR:
C!      - Added VDATE and modified COMMON BLOCK VERBOSE
C!      - Added variable NXYPOS to fix bug in TARG3D
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN
      REAL*8 EPS, ONE, ZERO
      PARAMETER(LUN=11, EPS=1.0e-8, ONE=1.0, ZERO=0.0)
C!    ----- Parameters for lattice orientation
      INTEGER*4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER(MAZMTH=1, MDIV=0, MTHETA=2**MDIV,  
     &           MDIR=MAZMTH*(MTHETA*(MTHETA+1))/2)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 KMAX, MIONIZ, MXRPOS, MXTARG, MXTRG2, MXYPOS, MXZROI
      PARAMETER(KMAX=9, MIONIZ=100000, MXRPOS=101, MXTARG=1025,
     &          MXTRG2=513, MXYPOS=4*MXRPOS*(MXRPOS-1)+1, MXZROI=100)
C!    ----- Input parameters for geometry
      REAL*8 DBEAM  ! Diameter of beam in nm
      REAL*8 DLATC  ! Cell lattice constant 
      REAL*8 DROI   ! Diameter of region of interest (cell nucleus size)
      REAL*8 DSITE  ! Diameter of spherical target 
      REAL*8 DZROI  ! Increment in position of region of interest 
C!    ----- Functions
      CHARACTER TSTAMP*24
C!    ----- Local scalars
      CHARACTER*11 VDATES(2), SRDATE
      CHARACTER FILENM*80, FILOUT*85, HEADER*80, PARAM*2, PREFIX*10
      INTEGER*4 I, IFILE, IT, ITA, J, K, KMX, L, NAZMTH, NDIV, NFILES, 
     &          NHEADL, NZROI, NPHI, NTRACS
      INTEGER*4 NIONIZ, NRPOS, NXYPOS
      LOGICAL ASKINP
      REAL*8 COSPHI, DELTAR, DUMMY, FNORM, PIBY4, RROI,
     &       SINPHI, SPOSIN, VCHECK
C!    ----- Local arrays
      INTEGER*4 NTARG(2), NTARGK(KMAX)
      REAL*8 CORRSE(MXTARG,MXTRG2),CONVOL(MXTARG),FREQBV(MXTARG,MXTRG2)
      REAL*8 FREQSE(MXTARG,KMAX),FREQTE(MXTARG,KMAX),ZROIC(MXZROI)
      REAL*8 RLINE(8)
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR 
C!    ------
      INTEGER*4 IRAD(MXYPOS)
      INTEGER*4 NRROI2 ! Integer of Square of ROI radius 
      REAL*8 XROI(MXYPOS), YROI(MXYPOS), ZROI(3)
      COMMON /ROICTR/  XROI, YROI, ZROI, IRAD, NRROI2
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
C!    ------
C!     CODE: Character string encoding the meaning of the entries in a 
C!           line of the input file as follows:
C!           'T'     - number of the primary particle track
C!           'X','Y','Z' - x, y, and z coordinates of the transfer point 
C!           'E'     - energy deposit (if applicable)
C!           'I'     - ionization cluster size (if applicable)
C!           '%'     - additional data that are not used 
C!           Example: The string 'TZ%%EXY ' indicates that there are 7
C!                    entries in each line, of which the first is the 
C!                    track number, the second the z coordinate, the 
C!                    fifth the energy, and the sixth and seventh the x
C!                    and y coordinates
C!     IDT: Array of the columns indices of (1-3) x,y,z coordinates 
C!                   (4) energy deposit (if present), (5) track number, 
C!                   (6) number of ionizations in cluster (if present)
C!                   (7-8) are there for future use
      CHARACTER CODE*8
      INTEGER*2 IDT(8),NDT
      COMMON /RFORMT/IDT,NDT,CODE
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program ROI_3D Version: '//VDATE
      SRFLAG=.TRUE. ! Invoke message on Version date from Subroutine TARG3D 

      IF(1.EQ.1) THEN ! DUMMY block for readability: Read input 1111111
        DEBUG(1)=.FALSE. ! Main program
        DEBUG(2)=.FALSE. ! TARG3D main sections
        DEBUG(3)=.FALSE. ! TARG3D IDIR loop
        DEBUG(4)=.FALSE. ! TARG3D IPOS loop details
        
        PIBY4=ATAN(ONE)
        PREFIX='3D__0.0nm_'
        DO K=1, KMAX
          NTARGK(K)=0
        END DO

C!      Check command file version
        PRINT*, 'Enter 0 for manual input or version (YYMMDD.HHMM) '//
     &          'of command file structure ' 
        READ(*,*) VCHECK
        ASKINP=(VCHECK.EQ.ZERO) 
        IF(.NOT.ASKINP) THEN
          IF(VCHECK.LT.VINPUT) THEN
            PRINT*, 'Command file structure ',VCHECK,' older than '//
     &              'current version ',VINPUT,'=> STOP.'
            STOP
          END IF       
          READ(*,*) FILENM
          IF(FILENM(1:6).NE.'ROI_3D') THEN
            PRINT*, 'Command file appear not to be for ROI_3D but '//
     &              'for '//FILENM//'=> STOP.'
            STOP
          END IF    
        END IF
        
C!      Read debug options
        IF(ASKINP) PRINT*, 'Enter debug options file name or - for none' 
        READ(*,*) FILENM        
        OPEN(LUN,FILE=FILENM,STATUS='OLD',ERR=10)
        DO I=1,4
          READ(LUN,*) DEBUG(I)
        END DO ! I=1,4
        CLOSE(LUN)
  10    CONTINUE

C!      Read geometry parameters 
        IF(ASKINP) PRINT*, 'Enter site diameter in nm'
        READ(*,*) DSITE
        DLATC=DSITE*SQRT(2.)*EXP(LOG(PIBY4/3.)/3.) 
        IF(DSITE.LT.10.0) THEN 
          WRITE(PREFIX(5:7),'(f3.1)') DSITE
        ELSE        
          WRITE(PREFIX(4:7),'(f4.1)') DSITE
        END IF
        
        IF(ASKINP) PRINT*, 'Enter ROI diameter in nm'
        READ(*,*) DROI
        RROI=DROI/2.
        NRROI2=INT(RROI*RROI/DLATC/DLATC) ! Note: This is correct as DLATC is the unit of length
        DELTAR=DROI/40.
        
        IF(ASKINP) PRINT*, 'Enter beam diameter in nm <= ',5.*DROI
        READ(*,*) DBEAM
        NRPOS=1+NINT(DBEAM/2./DELTAR)
        IF(NRPOS.GT.MXRPOS) THEN
          NRPOS=MXRPOS
          PRINT*, 'Beam diameter ',DBEAM,' nm is too large. <<<<<<<<<<<'
          DBEAM=2.*DELTAR*REAL(MXRPOS-1)
          PRINT*, '>>> maximum possible value ',DBEAM,' is used.'
        END IF
        
        IF(ASKINP) PRINT*, 'Enter maximum # of targets in histogram'
        READ(*,*) NTARG(1), NTARG(2)
        IF(NTARG(1).GT.MXTARG.OR.NTARG(1).LT.1) NTARG(1)=MXTARG
        IF(NTARG(2).GT.MXTRG2.OR.NTARG(2).EQ.1) NTARG(2)=MXTRG2

        IF(ASKINP) PRINT*, 'Enter maximum ionization cluster '//
     &                     'complexity (KMAX)'
        READ(*,*) KMX
        IF(KMX.GT.KMAX.OR.KMX.LT.2) KMX=KMAX

        IF(ASKINP) PRINT*, 'Enter # of regions of interest along track'
        READ(*,*) NZROI
        IF (NZROI.EQ.1) THEN
          IF(ASKINP) PRINT*, 'Enter z position of region of interest'
          READ(*,*) ZROIC(1)
          DZROI=ZERO
        ELSE
          IF(ASKINP) PRINT*, 'Enter position of first region of '//
     &            'interest (ROI) and increment in ROI position'
          READ(*,*) ZROIC(1), DZROI
          DO I=2,NZROI
            ZROIC(I)=ZROIC(I-1)+DZROI
          END DO
        END IF

        IF(ASKINP) PRINT*, 'Enter 8 character code for data file '//
     &       'structure where ''T'' indicates the track ID,'//
     &       ' ''X'',''Y'' and ''Z'' the respective coordinates, '//
     &       '''E'' the energy deposit (if present) and ''&'' any '//
     &       'other data'
        READ(*,*) CODE
        CALL RFINIT()

        IF(ASKINP) PRINT*, 'Number of header lines'
        READ(*,*) NHEADL

        IF(ASKINP) PRINT*, 'Number of files to process'
        READ(*,*) NFILES
        
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
        IF(ASKINP) PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM
        IF(.NOT.ASKINP) PRINT*, 'Processing file ',FILENM
        
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
          DO I=1,NHEADL
            READ(LUN,*) HEADER
          END DO
   90     READ(LUN,*) (RLINE(I),I=1,NDT)
          IF(IDT(6).GT.0) THEN
            IF(NINT(RLINE(IDT(6))).LT.2) GOTO 90 ! No ionization cluster
          END IF 
          ITA=NINT(RLINE(IDT(5)))
          DO I=1,3
            XYZ(NIONIZ,I)=RLINE(IDT(I))
          END DO
          
          NIONIZ=NIONIZ+1
          IF(DEBUG(1)) PRINT*, 'Hier 4'
  100     CONTINUE ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          READ(LUN,*,END=200,ERR=200) (RLINE(I),I=1,NDT)
          IF(IDT(6).GT.0) THEN
            IF(NINT(RLINE(IDT(6))).LT.2) GOTO 100 ! No ionization cluster
          END IF 
          IT=NINT(RLINE(IDT(5)))
          DO I=1,3
            XYZ(NIONIZ,I)=RLINE(IDT(I))
          END DO

          IF(IT.EQ.ITA) THEN
            NIONIZ=NIONIZ+1
          ELSE ! Entry from a new track was read
            NIONIZ=NIONIZ-1 ! Decrease to number of ionizations in track
            IF(DEBUG(1)) PRINT*, 'Hier vor CALL TARG3D'
            IF(MOD(NTRACS,10).EQ.0) PRINT*, NTRACS
            DO J=1,NZROI 
              ZROI(2)=ZROIC(J)
              ZROI(1)=ZROI(2)-DROI/2.-DLATC/2.
              ZROI(3)=ZROI(2)+DROI/2.+DLATC/2.
              CALL TARG3D(NIONIZ,NRPOS,NTARG,NXYPOS,SRDATE)
            END DO ! J=1,NZROI

            IF(DEBUG(1)) PRINT*, 'Hier nach CALL TARG3D'
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
              ZROI(1)=ZROI(2)-DROI/2.-DLATC/2.
              ZROI(3)=ZROI(2)+DROI/2.+DLATC/2.
              CALL TARG3D(NIONIZ,NRPOS,NTARG,NXYPOS,SRDATE) 
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
          DO K=KMX-1,2,-1
            DO J=1,NRPOS
              DO I=1,NTARG(1)
                CONVOL(I)=ZERO
              END DO ! I=1,NTARG(1)
              DO I=1,NTARG(1)
                DO L=0, I-1
                  CONVOL(I)=CONVOL(I)+FREQRD(I-L,J,K)*FREQRD(L+1,J,K+1)
                END DO ! L=0, I-1
              END DO ! I=1,NTARG(1)
              DO I=1,NTARG(1)
                FREQRD(I,J,K)=CONVOL(I)
              END DO ! I=1,NTARG(1)
            END DO ! J=1,NRPOS
          END DO ! K=KMX,1,-1
        END IF          ! DUMMY block for readability: Normalize  666666

        IF(8.EQ.8) THEN ! DUMMY block for readability: Prepare Output  8888888
C!        Calculate single-event distributions
          DO K=1,KMX
            NTARGK(K)=1
            DO I=1,NTARG(1)
              FREQSE(I,K)=FREQRD(I,1,K)
              SPOSIN=ONE
              DUMMY=ZERO
              DO J=2,NRPOS
                DUMMY=DUMMY+8.
                FREQSE(I,K)=FREQSE(I,K)+DUMMY*FREQRD(I,J,K)
                SPOSIN=SPOSIN+DUMMY         
                IF (RADIST(J).LE.RROI) THEN
                  FREQTE(I,K)=FREQSE(I,K)/SPOSIN
                END IF
              END DO   
              FREQSE(I,K)=FREQSE(I,K)/SPOSIN
              IF(FREQSE(I,K).GE.EPS.OR.FREQTE(I,K).GE.EPS) NTARGK(K)=I            
            END DO
          END DO
          
          IF(NTARG(2).GT.0) THEN 
            IF(NTARGK(2).GT.MXTRG2) THEN
              PRINT*, 'Major problem: 2nd array dimension MXTRG2=',
     &                MXTRG2,' < max. number of 2+ clusters NTARGK(2)=',
     &                NTARGK(2)
              STOP
            END IF
            DO I=1,NTARGK(1)
              DO L=1,NTARGK(2)
                FREQBV(I,L)=CORR12(I,L,1)  
                SPOSIN=ONE
                DUMMY=ZERO
                DO J=2,NRPOS
                  DUMMY=DUMMY+8.
                  FREQBV(I,L)=FREQBV(I,L)+DUMMY*CORR12(I,L,J)
                  SPOSIN=SPOSIN+DUMMY         
                END DO   
                FREQBV(I,L)=FREQBV(I,L)/SPOSIN
                IF(FREQSE(I,1)*FREQSE(L,2).GT.ZERO) THEN
                  CORRSE(I,L)=FREQBV(I,L)/(FREQSE(I,1)*FREQSE(L,2))
                ELSE
                  CORRSE(I,L)=FREQBV(I,L)
                END IF                
              END DO
            END DO
          END IF            
        END IF          ! DUMMY block for readability: Prepare Output 8888888

       IF(9.EQ.9) THEN ! DUMMY block for readability: OUTPUT  9999999999
C!        #Write results to output file 
          WRITE(PREFIX(2:2),'(a)') 'D'
          FILOUT=PREFIX//FILENM
          PRINT*, 'Write output to '//FILOUT
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM ROI_3D Version '//
     &              VDATE//' BLINIT:'//VDATES(1)//' TARG3D:'//VDATES(2)
     &                   //' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,'(a10,10i6)') ' NTARG(K)=',(NTARGK(I),I=1,KMAX)
          WRITE(LUN,'(6(a,f8.3))') ' DLATC=',DLATC,' nm    DROI= ',DROI,
     &                             ' nm   DSITE=',DSITE,' nm   DBEAM= ',
     &                             DBEAM,' nm'
          PARAM='P '
          DO K=1,KMX
            IF(K.GT.1) PARAM='F '
            WRITE(PARAM(2:2),'(I1)') K
            WRITE(LUN,'(1X,2a8,103a15)') 'Para-', '#Sites', 'Average',
     &                          'Average',('Distance/nm',J=1,NRPOS)
            WRITE(LUN,'(1X,2a8,2a15,101f15.6)') 'meter', '/track',
     &                         'total','inside',(RADIST(J),J=1,NRPOS)
            DO I=1,NTARGK(K)
              WRITE(LUN,'(1X,a8,i8,103f15.8)') PARAM, I-1,
     &                                     FREQSE(I,K),FREQTE(I,K),
     &                                    (FREQRD(I,J,K),J=1,NRPOS)
            END DO
            WRITE(LUN,*) '____________________________________________'
            WRITE(LUN,*) '********************************************'
          END DO
          CLOSE(LUN)
          
          IF(NTARG(2).GT.0) THEN ! begin 02-APR-2021 >>>>>>>>>>>>
            WRITE(PREFIX(2:2),'(a)') 'B'
            FILOUT=PREFIX//FILENM
            PRINT*, 'Write output to '//FILOUT   
            OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
            WRITE(LUN,*) '*** Output from PROGRAM ROI_3D Version '//
     &              VDATE//' BLINIT:'//VDATES(1)//' TARG3D:'//VDATES(2)
     &                   //' on '//TSTAMP()
            WRITE(LUN,*) 'Filename: '//FILOUT
            WRITE(LUN,'(6(a,f8.3))') ' DLATC=',DLATC,' nm    DROI= ',
     &            DROI,' nm   DSITE=',DSITE,' nm   DBEAM= ',DBEAM,' nm'
            WRITE(LUN,*) 'Correlations P1 and F2'
            WRITE(LUN,*) NTARGK(1),NTARGK(2)
            DO I=1,NTARGK(1)
              WRITE(LUN,'(1X,1000e15.8)') (FREQBV(I,L),L=1,NTARGK(2))
            END DO
            CLOSE(LUN)

            WRITE(PREFIX(2:2),'(a)') 'C'
            FILOUT=PREFIX//FILENM
            PRINT*, 'Write output to '//FILOUT   
            OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
            WRITE(LUN,*) '*** Output from PROGRAM ROI_3D Version '//
     &              VDATE//' BLINIT:'//VDATES(1)//' TARG3D:'//VDATES(2)
     &                   //' on '//TSTAMP()
            WRITE(LUN,*) 'Filename: '//FILOUT
            WRITE(LUN,'(6(a,f8.3))') ' DLATC=',DLATC,' nm    DROI= ',
     &            DROI,' nm   DSITE=',DSITE,' nm   DBEAM= ',DBEAM,' nm'
            WRITE(LUN,*) 'Correlations P1 and F2'
            WRITE(LUN,*) NTARGK(1),NTARGK(2)
            DO I=1,NTARGK(1)
              WRITE(LUN,'(1X,1000e15.8)') (CORRSE(I,L),L=1,NTARGK(2))
            END DO
            CLOSE(LUN)
          END IF ! (NTARG(2).GT.0)
        END IF          ! DUMMY block for readability: OUTPUT  999999999

      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! ROI_3D
C!______________________________________________________________________      

      INCLUDE 'BLINIT.f'  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'ROI_3D_Subs.f'  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'RFINIT.f'
      INCLUDE 'TSTAMP.f'