      PROGRAM IC_3E
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='05-JUN-2021',VINPUT=210605.16d0) ! ############
C!    #################################################################
C!    04-JUL-2022 HR: 
C!      - Added timestap to output files 
C!      - Reworked identification of input file structure
C!      - Added version check for input files
C!      - Optimized output
C!    20-MAY-2021 HR:
C!      - Finalized adaption from IC_3D.f for Ziad's data
C!    28-MAR-2021 HR:
C!      - Cleaned up unused variables 
C!    20-MAR-2021 HR:
C!      - Created this clone of ROI_3D.f for ionization CLUSTR output
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN, LUNT
      REAL*8 ONE, ZERO
      PARAMETER(LUN=11, LUNT=12, ONE=1.0, ZERO=0.0)
C!    ----- Parameters for lattice orientation
      INTEGER*4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER(MAZMTH=1, MDIV=0, MTHETA=2**MDIV,  
     &           MDIR=MAZMTH*(MTHETA*(MTHETA+1))/2)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 MTRANS
      PARAMETER(MTRANS=100000)
C!    ----- Input parameters for geometry
      REAL*8 DLATC  ! Cell lattice constant 
      REAL*8 DSITE  ! Diameter of spherical target 
C!    ----- Functions
      CHARACTER TSTAMP*24
C!    ----- Local scalars
      CHARACTER FILENM*80, FILOUT*84, HEADER*80, PREFIX*9
      CHARACTER*11 VDATES(2), SRDATE
      INTEGER*4 I, IFILE, IT, ITA, J, K, NAZMTH, NDIV, NFILES, 
     &          NFORMT, NHEADL, NTRACS
      INTEGER*4 NTRANS
      LOGICAL ASKINP
      REAL*8 DUMB1, PIBY4, VCHECK
      REAL*8 ZROI(2)
C!    ----- Local arrays
      REAL*8 RLINE(8)
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR 
C!    ------
      INTEGER*4 ITYPE(MTRANS)
      REAL*8 XYZE(MTRANS,4)
      COMMON /TRACKS/ XYZE, ITYPE
C!    ------
      LOGICAL DEBUG(4), SRFLAG
      COMMON /VERBOSE/ DEBUG, SRFLAG
C!    ------ 
      INTEGER*4 ICSIZE(MTRANS),NSITES
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

      PRINT*, 'Program IC_3E Version: '//VDATE
      SRFLAG=.TRUE. ! Invoke message on Version date from Subroutine CLUSTR 

      IF(1.EQ.1) THEN ! DUMMY block for readability: Read input 1111111
        DEBUG(1)=.FALSE. ! Main program
        DEBUG(2)=.FALSE. ! CLUSTR main sections
        DEBUG(3)=.FALSE. ! CLUSTR IDIR loop
        DEBUG(4)=.FALSE. ! CLUSTR IPOS loop details
        
        PIBY4=ATAN(ONE)
        PREFIX='IC_0.0nm_'
        
C!      Check command file version
        PRINT*, 'Enter 0 for manual input or version (YYMMDD.HHMM) '//
     &          'of command file structure ' 
        READ(*,*) FILENM
        READ(FILENM(1:11),*) VCHECK
        ASKINP=(VCHECK.EQ.ZERO)
        IF(.NOT.ASKINP) THEN
          IF(VCHECK.LT.VINPUT) THEN
            PRINT*, 'Command file structure ',VCHECK,' older than '//
     &              'current version ',VINPUT,'=> STOP.'
            STOP
          END IF       
          READ(*,*) FILENM
          IF(FILENM(1:5).NE.'IC_3D') THEN
            PRINT*, 'Command file appears not to be for IC_3D but '//
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
        DLATC=DSITE*SQRT(2.)*EXP(LOG(PIBY4/3)/3.)
        IF(DSITE.LT.10.0d0) THEN
          WRITE(PREFIX(4:6),'(f3.1)') DSITE
        ELSE
          IF(DSITE.LT.100.0d0) THEN
            WRITE(PREFIX(4:6),'(f3.0)') DSITE
          ELSE
            WRITE(PREFIX(4:6),'(i3)') INT(DSITE)
          END IF
        END IF
        
        IF(ASKINP) PRINT*, 'z position of begin of region of interest'
        READ(*,*) ZROI(1)
        IF(ASKINP) PRINT*, 'z position of end of region of interest'
        READ(*,*) ZROI(2)

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
        NAZMTH=MAZMTH
        NDIV=MDIV
        IF(DEBUG(1)) PRINT*, 'Before CALL BLINIT'
        CALL BLINIT(NAZMTH,NDIV,DLATC,SRDATE) ! Init reziprocal lattice
        VDATES(1)=SRDATE
        IF(DEBUG(1)) PRINT*, 'After CALL BLINIT'
      END IF          ! DUMMY block for readability: Initialize 2222222
      
      DO IFILE=1,NFILES      
        PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM  

        IF(4.EQ.4) THEN ! DUMMY block for readability: Main Loop  444444
          IF(DEBUG(1)) PRINT*, 'Begin of Block 4'
C!        Init counters       
          NTRANS=1
          NTRACS=1
          
          IF (NFORMT.EQ.1) THEN 
            FILOUT=PREFIX//FILENM(6:80)
          ELSE
            FILOUT=PREFIX//FILENM
          END IF
          
          OPEN(LUNT,FILE=FILOUT,STATUS='UNKNOWN') 
          WRITE(LUNT,'(a)') 'Output from IC_3E.f - Version '//VDATE 
          WRITE(LUNT,'(a,f8.3,a)') 'Ionization CLUSTR positions in '//
     &            'ion tracks for DSITE=',DSITE,' nm '   
          WRITE(LUNT,'(a)') 'Processed input data file: '//FILENM     

          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          READ(LUN,*) ITA,(XYZE(NTRANS,J),J=1,4),ITYPE(NTRANS) 
          HEADER='  event    x/nm     y/nm    z/nm    E_imp/eV  n_ion'     
          WRITE(LUNT,'(a)') HEADER   
          
          NTRANS=NTRANS+1
          IF(DEBUG(1)) PRINT*, 'Hier 4'
  100     CONTINUE ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          READ(LUN,*,END=200,ERR=200) IT,(XYZE(NTRANS,J),J=1,4),
     &                                     ITYPE(NTRANS)         
          IF(IT.EQ.ITA) THEN
            NTRANS=NTRANS+1
          ELSE ! Entry from a new track was read
            NTRANS=NTRANS-1 ! Decrease to number of ionizations in track
            IF(DEBUG(1)) PRINT*, 'Hier vor CALL CLUSTR'
            IF(MOD(NTRACS,10).EQ.0) PRINT*, NTRACS
            DO I=1,NTRANS
              DO J=1,3
                XYZE(I,J)=XYZE(I,J)*1000.
              END DO
            END DO
            CALL CLUSTR(NTRANS,ZROI) 
            DO K=1,NSITES-1 
              DO I=K+1,NSITES
                IF(XYZE(K,3).GT.XYZE(I,3)) THEN
                  DO J=1, 4
                    DUMB1=XYZE(K,J)
                    XYZE(K,J)=XYZE(I,J)
                    XYZE(I,J)=DUMB1
                  END DO ! J=1, 4
                  J=ICSIZE(K)
                  ICSIZE(K)=ICSIZE(I)
                  ICSIZE(I)=J
                END IF
              END DO ! I=K+1,NSITES
            END DO !  K=1,NSITES-1             
            DO K=1,NSITES 
              WRITE(LUNT,'(i8,4f13.3,i4)') NTRACS,(XYZE(K,J),J=1,4),
     &                                        ICSIZE(K)
            END DO
            
            IF(DEBUG(1)) PRINT*, 'Hier nach CALL CLUSTR'
            DO J=1,4
              XYZE(1,J)=XYZE(NTRANS+1,J) ! First ionization in new track
            END DO
            ITA=IT ! Remember new track number
            NTRACS=NTRACS+1
            NTRANS=2 ! There was already one ionization
          END IF
          GOTO 100 ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  200     CONTINUE
            CALL CLUSTR(NTRANS,ZROI) 
            DO K=1,NSITES-1 
              DO I=K+1,NSITES
                IF(XYZE(K,3).GT.XYZE(I,3)) THEN
                  DO J=1, 4
                    DUMB1=XYZE(K,J)
                    XYZE(K,J)=XYZE(I,J)
                    XYZE(I,J)=DUMB1
                  END DO ! J=1, 3
                  J=ICSIZE(K)
                  ICSIZE(K)=ICSIZE(I)
                  ICSIZE(I)=J
                END IF
              END DO ! I=K+1,NSITES
            END DO !  K=1,NSITES-1             
  300       DO K=1,NSITES 
              WRITE(LUNT,'(i8,4f13.3,i4)') NTRACS,(XYZE(K,J),J=1,4),
     &                                        ICSIZE(K)
              !PRINT*, NTRACS,(XYZE(K,J),J=1,4),  ICSIZE(K)
            END DO ! K=1,NSITES
            CLOSE(LUNT)
          CLOSE(LUN)
        END IF          ! DUMMY block for readability: Main Loop  4444444
          
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! IC_3E
C!______________________________________________________________________      

      INCLUDE 'BLINIT.f'       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'IC_3D_SUBS.f'   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'RFINIT.f'
      INCLUDE 'TSTAMP.f'
