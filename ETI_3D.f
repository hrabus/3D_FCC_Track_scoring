      PROGRAM ETI_3D ! Energy Transfer by Ionizations 
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date
      CHARACTER VDATE*11
      PARAMETER(VDATE='23-MAY-2021') !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C!    23-MAY-2021 HR: 
C!      - Modified call to subroutines to get their version date
C!    28-MAR-2021 HR:
C!      - Cleaned up unused variables 
C!    20-MAR-2021 HR:
C!      - Created this clone of ROI_3D.f for ionization ETRFPT output
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
C!    ----- Local scalars
      CHARACTER*11 VDATES(2), SRDATE
      CHARACTER FILENM*80, FILOUT*84, HEADER*80, PREFIX*9
      INTEGER*4 I, IFILE, IT, ITA, J, K, NAZMTH, NDIV, NFILES, 
     &          NFORMT, NTRACS
      INTEGER*4 NIONIZ
      REAL*8 DUMB1, DUMB2, DUMB3, PIBY4
      REAL*8 ZROI(2)
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR 
C!    ------
      REAL*8 XYZE(MIONIZ,4)
      COMMON /TRACKS/ XYZE
C!    ------
      LOGICAL DEBUG(4), SRFLAG
      COMMON /VERBOSE/ DEBUG, SRFLAG
C!    ------ 
      INTEGER*4 ICSIZE(MIONIZ),NSITES
      COMMON /CLSTRS/ ICSIZE, NSITES
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program ETI_3D Version: '//VDATE
      SRFLAG=.TRUE. ! Invoke message on Version date from Subroutine ETRFPT 

      IF(1.EQ.1) THEN ! DUMMY block for readability: Read input 1111111
        DEBUG(1)=.FALSE. ! Main program
        DEBUG(2)=.FALSE. ! ETRFPT main sections
        DEBUG(3)=.FALSE. ! ETRFPT IDIR loop
        DEBUG(4)=.FALSE. ! ETRFPT IPOS loop details
        
        PIBY4=ATAN(ONE)
        PREFIX='IC_0.0nm_'
        
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

        PRINT*, 'Number of files to process'
        READ(*,*) NFILES
        PRINT*, 'Enter flag for data file structure 1=SPRE  2=Heidi'
        READ(*,*) NFORMT
        PRINT*, 'Enter z position of begin of region of interest'
        READ(*,*) ZROI(1)
        PRINT*, 'Enter z position of end of region of interest'
        READ(*,*) ZROI(2)
      END IF          ! DUMMY block for readability: Read input 1111111

      IF(2.EQ.2) THEN ! DUMMY block for readability: Initialize 2222222
        NAZMTH=MAZMTH
        NDIV=MDIV
        IF(DEBUG(1)) PRINT*, 'Hier vor CALL BLINIT'
        CALL BLINIT(NAZMTH,NDIV,DLATC,SRDATE) ! Init reziprocal lattice
        VDATES(1)=SRDATE
        IF(DEBUG(1)) PRINT*, 'Hier nach CALL BLINIT'
      END IF          ! DUMMY block for readability: Initialize 2222222
      
C!Begin new code on 16-OCT-2020  
      DO IFILE=1,NFILES      
        PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM  
C!End new code on 16-OCT-2020  

        IF(4.EQ.4) THEN ! DUMMY block for readability: Main Loop  444444
          IF(DEBUG(1)) PRINT*, 'Hier beginnt Block 4'
C!        Init counters       
          NIONIZ=1
          NTRACS=1
          
          IF (NFORMT.EQ.1) THEN 
            FILOUT=PREFIX//FILENM(6:80)
          ELSE
            FILOUT=PREFIX//FILENM
          END IF
          
          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          IF (NFORMT.EQ.1) THEN 
            DO I=1,4
              READ(LUN,'(a)') HEADER
            END DO
            READ(LUN,*) ITA,XYZE(NIONIZ,3),DUMB1,DUMB2,XYZE(NIONIZ,4),
     &                  XYZE(NIONIZ,1),XYZE(NIONIZ,2)
          ELSE
            READ(LUN,'(a)') HEADER
            READ(LUN,*) ITA,DUMB1,DUMB2,XYZE(NIONIZ,3),
     &                   XYZE(NIONIZ,1),XYZE(NIONIZ,2)
            XYZE(NIONIZ,4)=ZERO
          END IF
          WRITE(LUNT,'(a)') HEADER   
          
          NIONIZ=NIONIZ+1
          IF(DEBUG(1)) PRINT*, 'Hier 4'
  100     CONTINUE ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          IF (NFORMT.EQ.1) THEN 
            READ(LUN,*,END=200,ERR=200) IT,XYZE(NIONIZ,3),DUMB1,DUMB2,
     &                    XYZE(NIONIZ,4),XYZE(NIONIZ,1),XYZE(NIONIZ,2)
          ELSE
            READ(LUN,*,END=200,ERR=200) IT,DUMB1,DUMB2,XYZE(NIONIZ,3),
     &                   XYZE(NIONIZ,1),XYZE(NIONIZ,2)
            XYZE(NIONIZ,4)=ZERO
          END IF
          IF(IT.EQ.ITA) THEN
            NIONIZ=NIONIZ+1
          ELSE ! Entry from a new track was read
            NIONIZ=NIONIZ-1 ! Decrease to number of ionizations in track
            IF(DEBUG(1)) PRINT*, 'Hier vor CALL ETRFPT'
            IF(MOD(NTRACS,10).EQ.0) PRINT*, NTRACS
            CALL ETRFPT(NIONIZ,ZROI,SRDATE)
            VDATES(2)=SRDATE 
            IF(NTRACS.EQ.1) THEN
              OPEN(LUNT,FILE=FILOUT,STATUS='UNKNOWN') 
              WRITE(LUNT,'(a)') 'Output from ETI_3D.f - Version '//
     &             VDATE//' BLINIT:'//VDATES(1)//' CLUSTR:'//VDATES(2)
              WRITE(LUNT,'(a,f8.3,a)') 'Targets with energy deposits'//
     &            ' in ion tracks for DSITE=',DSITE,' nm '   
              WRITE(LUNT,'(a)') 'Processed input data file: '//FILENM     
            END IF            
            
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
              WRITE(LUNT,'(2(i8,3f10.3))') IT,XYZE(K,3),ZERO,XYZE(K,4),
     &                                    ICSIZE(K),XYZE(K,1),XYZE(K,2)
            END DO
            
            IF(DEBUG(1)) PRINT*, 'Hier nach CALL ETRFPT'
            DO J=1,4
              XYZE(1,J)=XYZE(NIONIZ+1,J) ! First ionization in new track
            END DO
            ITA=IT ! Remember new track number
            NTRACS=NTRACS+1
            NIONIZ=2 ! There was already one ionization
          END IF
          GOTO 100 ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
  200     CONTINUE
            CALL ETRFPT(NIONIZ,ZROI,SRDATE) 
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
            DO K=1,NSITES 
              WRITE(LUNT,'(2(i8,3f10.3))') IT,XYZE(K,3),ZERO,XYZE(K,4),
     &                                    ICSIZE(K),XYZE(K,1),XYZE(K,2)
            END DO ! K=1,NSITES
            CLOSE(LUNT)
          CLOSE(LUN)
        END IF          ! DUMMY block for readability: Main Loop  4444444
          
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! ETI_3D
C!______________________________________________________________________      

      INCLUDE 'BLINIT.f'  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INCLUDE 'ETI_3D_Subs.f'   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
