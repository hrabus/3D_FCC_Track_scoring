      PROGRAM IC_3D
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='13-JUL-2023',VINPUT=210605.16d0) ! ############
C!    #################################################################
C!    20-JUL-2023 HR: Created this spawn that outputs also single ions
C!    13-JUL-2023 HR: Fixed bug in WRFLT subroutine, which caused cases
C!                    in which the numbers in the output file were not
C!                    separated by a blank.
C!    04-JUL-2022 HR: Changed file name of IC_3D_SUBS.f to uppercase
C!                    for compatibility with LINUX
C!    05-JUN-2021 HR:
C!      - Added timestap to output files 
C!      - Reworked identification of input file structure
C!      - Added version check for input files
C!    03-JUN-2021 HR: 
C!      - Optimized output
C!    02-JUN-2021 HR: 
C!      - Fixed soft bug with output of header line
C!      - Fixed problem with formatted output
C!    28-MAY-2021 HR: 
C!      - Fixed bug with length of SRDATE variable 
C!    23-MAY-2021 HR: 
C!      - Modified call to subroutines to get their version date
C!    20-MAY-2021 HR:
C!      - Adapted to potential site sizes >= 10 nm
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
      INTEGER*4 MIONIZ
      PARAMETER(MIONIZ=100000)
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
      INTEGER*4 NIONIZ
      LOGICAL ASKINP
      REAL*8 DUMMY, PIBY4, VCHECK
      REAL*8 ZROI(2)
C!    ----- Local arrays
      REAL*8 RLINE(8)
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR 
C!    ------
      REAL*8 XYZ(MIONIZ,3)
      COMMON /TRACKS/ XYZ
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

      PRINT*, 'Program IC_3D Version: '//VDATE
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
          NIONIZ=1
          NTRACS=1
          
          IF (NFORMT.EQ.1) THEN 
            FILOUT=PREFIX//FILENM(6:80)
          ELSE
            FILOUT=PREFIX//FILENM
          END IF
          
          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          DO I=1,NHEADL
            READ(LUN,'(a)') HEADER
          END DO
          READ(LUN,*) (RLINE(I),I=1,NDT)
          ITA=NINT(RLINE(IDT(5)))
          DO I=1,3
            XYZ(NIONIZ,I)=RLINE(IDT(I))
          END DO
          
          NIONIZ=NIONIZ+1

  100     CONTINUE ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          READ(LUN,*,END=200,ERR=200) (RLINE(I),I=1,NDT)
          IT=NINT(RLINE(IDT(5)))
          DO I=1,3
            XYZ(NIONIZ,I)=RLINE(IDT(I))
          END DO
          
          IF(IT.EQ.ITA) THEN
            NIONIZ=NIONIZ+1
          ELSE ! Entry from a new track was read
            NIONIZ=NIONIZ-1 ! Decrease to number of ionizations in track
            IF(DEBUG(1)) PRINT*, 'Before CALL CLUSTR'
            IF(MOD(NTRACS,10).EQ.0) PRINT*, NTRACS
            CALL CLUSTR(NIONIZ,ZROI,SRDATE)
            VDATES(2)=SRDATE
            IF(NTRACS.EQ.1) THEN
              OPEN(LUNT,FILE=FILOUT,STATUS='UNKNOWN') 
              WRITE(LUNT,'(a)') 'Output from IC_3D.f - Version '//
     &             VDATE//' BLINIT:'//VDATES(1)//' CLUSTR:'//VDATES(2)
     &                   //' on '//TSTAMP()
              WRITE(LUNT,'(a,f8.3,a)') 'Ionization cluster positions'//
     &            ' in ion tracks for DSITE=',DSITE,' nm '   
              WRITE(LUNT,'(a)') 'Processed input data file: '//FILENM     
              WRITE(LUNT,'(a)') 'Track x/nm y/nm z/nm ICS'   
            END IF            
            DO K=1,NSITES-1 
              DO I=K+1,NSITES
                IF(XYZ(K,3).GT.XYZ(I,3)) THEN
                  DO J=1, 3
                    DUMMY=XYZ(K,J)
                    XYZ(K,J)=XYZ(I,J)
                    XYZ(I,J)=DUMMY
                  END DO ! J=1, 3
                END IF
              END DO ! I=K+1,NSITES
            END DO !  K=1,NSITES-1             
            DO K=1,NSITES 
C!20-JUL-2023              IF(ICSIZE(K).GT.1) THEN
              IF(ICSIZE(K).GE.1) THEN
                CALL WRLINE(LUNT,IT,XYZ(K,1),XYZ(K,2),XYZ(K,3),
     &                           ICSIZE(K))
              END IF
            END DO ! K=1,NSITES
            
            IF(DEBUG(1)) PRINT*, 'After CALL CLUSTR'
            DO J=1,3
              XYZ(1,J)=XYZ(NIONIZ+1,J) ! First ionization in new track
            END DO
            ITA=IT ! Remember new track number
            NTRACS=NTRACS+1
            NIONIZ=2 ! There was already one ionization
          END IF
          GOTO 100 ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  200     CONTINUE
            CALL CLUSTR(NIONIZ,ZROI,SRDATE)
            DO K=1,NSITES-1 
              DO I=K+1,NSITES
                IF(XYZ(K,3).GT.XYZ(I,3)) THEN
                  DO J=1, 3
                    DUMMY=XYZ(K,J)
                    XYZ(K,J)=XYZ(I,J)
                    XYZ(I,J)=DUMMY
                  END DO ! J=1, 3
                END IF
              END DO ! I=K+1,NSITES
            END DO !  K=1,NSITES-1             
            DO K=1,NSITES 
C!20-JUL-2023              IF(ICSIZE(K).GT.1) THEN
              IF(ICSIZE(K).GE.1) THEN
                CALL WRLINE(LUNT,IT,XYZ(K,1),XYZ(K,2),XYZ(K,3),
     &                           ICSIZE(K))
              END IF
            END DO ! K=1,NSITES
            CLOSE(LUNT)
          CLOSE(LUN)
        END IF          ! DUMMY block for readability: Main Loop  4444444
          
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! IC_3D
C!______________________________________________________________________      

      INCLUDE 'BLINIT.f'       ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'IC_3D_SUBS.f'   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'RFINIT.f'
      INCLUDE 'TSTAMP.f'
      
      SUBROUTINE WRLINE(LUN,IT,X,Y,Z,ICS)
      INTEGER*4 ICS,IT,LUN,NCH
      REAL*8 X,Y,Z
      CHARACTER FORMTS*80
      COMMON /FRMT/FORMTS,NCH
      FORMTS='('
      NCH=1
      CALL WRINT(IT)
      CALL WRFLT(X,3)
      CALL WRFLT(Y,3)
      CALL WRFLT(Z,3)
      CALL WRINT(ICS)
      WRITE(FORMTS(NCH:NCH),'(a1)') ')'
      WRITE(LUN,FORMTS) IT,X,Y,Z,ICS
      END SUBROUTINE WRLINE

C!______________________________________________________________________
      SUBROUTINE WRINT(I)
      INTEGER*4 I, LINTEG, NCH, NST
      CHARACTER FORMTS*80, WFORMT*10
      COMMON /FRMT/FORMTS,NCH
      IF(I.EQ.0) THEN
        LINTEG=2
      ELSE
        IF(I.GT.0) THEN
          LINTEG=2+INT(LOG10(REAL(I)))
        ELSE
          LINTEG=3+INT(LOG10(ABS(REAL(I))))
        END IF
      END IF
      NST=NCH+1
      NCH=NCH+3
      WFORMT='(a1,i1,a1)'
      IF(LINTEG.GE.10) THEN
        NCH=NCH+1
        WFORMT='(a1,i2,a1)'
      END IF
      WRITE(FORMTS(NST:NCH),WFORMT) 'i',LINTEG,','
      END
      
C!______________________________________________________________________
      SUBROUTINE WRFLT(DF,IDIG)
      INTEGER*4 IDIG,NCH,NST
      REAL*8 DF
      CHARACTER FORMTS*80, WFORMT*16
      COMMON /FRMT/FORMTS,NCH
      IF(DF.EQ.0.0) THEN
        LFLOAT=3
      ELSE
C!        LFLOAT=4+INT(LOG10(ABS(DF))) ! HR 13-JUL-2023 
        LFLOAT=4+INT(LOG10(ABS(DF)+1e-3)) ! HR 13-JUL-2023 
        IF(LFLOAT.LT.4) LFLOAT=4
        IF(DF.GT.0) LFLOAT=LFLOAT-1
      END IF
      LFLOAT=LFLOAT+IDIG
      NST=NCH+1
      NCH=NCH+5
      WFORMT='(a1,i1,a1,i1,a1)'
      IF(LFLOAT.GE.10) THEN
        WRITE(WFORMT(6:6),'(i1)') 2
        NCH=NCH+1
      END IF
      IF(IDIG.GT.10) THEN
        WRITE(WFORMT(12:12),'(i1)') 2
        NCH=NCH+1
      END IF
      WRITE(FORMTS(NST:NCH),WFORMT) 'f',LFLOAT,'.',IDIG,','   
      END