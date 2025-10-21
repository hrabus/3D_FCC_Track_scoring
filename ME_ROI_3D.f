      PROGRAM ME_ROI_3D
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='06-JUN-2021',VINPUT=210605.1700d0) ! ###########
C!    06-JUN-2021 HR:
C!      - Added timestamp in output file
C!    05-JUN-2021 HR:
C!      - Reorganized initialization of counters.
C!    14-APR-2021 HR:
C!      - Adapted to changes done in ROI_3D (leaner output).
C!    31-MAR-2021 HR:
C!      - Fixed bug in initialization of convolution  
C!    28-MAR-2021 HR:
C!      - Added one-in-all output for easier transfer to Origin etc.
C!    26-MAR-2021 HR:
C!      - Added option for calculating fluence from dose (based on 
C!        new input info of proton energy)
C!    15-MAR-2020 HR:
C!      - Added VDATE and modified COMMON BLOCK VERBOSE
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN
      REAL*8 ONE, ZERO, PI, WMIN
      PARAMETER(LUN=11, ONE=1.0, ZERO=0.0, PI=3.1415926, WMIN=1d-12)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 KMAX, MXRPOS, MXTARG
      PARAMETER(KMAX=9, MXRPOS=101, MXTARG=65537)
C!    ----- Parameters for fluence calculation from power law 
C!          regression to PSTAR data between 1 MeV and 100 MeV
      REAL*8 DPERFL  ! Conversion factor from mass stopping power in
                     ! MeV cm²/g to dose per fluence in Gy nm²
      REAL*8 STPWR1  ! Stopping Power at 1 MeV
      REAL*8 STPEXP  ! Exponent of power law for stopping power
      PARAMETER(DPERFL=1.6e4, STPEXP=-0.7899, STPWR1=276.973)
C!    ----- Functions
      CHARACTER TSTAMP*24
C!    ----- Local scalars
      CHARACTER FILENM*120, FILOUT*120, HEADER*1200, PARAM*2, 
     &          RINFO(5,2)*11
      INTEGER*4 I, IFILE, IJ, ITARG, J, K, KMX, L, NFILES, NRPOS
      INTEGER*4 NTARG
      LOGICAL ASKINP, DEBUG, HIFLNC
      REAL*8 DBEAM   ! Diameter of beam in nm
      REAL*8 DOSE    ! Absorbed dose in Gray
      REAL*8 DROI    ! Diameter of ROI in nm
      REAL*8 ENERGY  ! Proton energy in MeV
      REAL*8 EVENTS  ! number of events
      REAL*8 FLUENC  ! average number of tracks per ROI cross section
      REAL*8 STPWRE  ! Stopping Power at energy E
      REAL*8 DLFLNC, DLWGHT, SPOSIN, SUMPOS, VCHECK, WGHT, TGMEAN(KMAX)
C!    ----- Local arrays
      INTEGER*4 JMAX(2),NTARGK(KMAX),NTARGT(KMAX)
      REAL*8 FREQSE(MXTARG,KMAX),FREQTE(MXTARG,KMAX),
     &       FREQME(MXTARG,KMAX),FNE(MXTARG,2,KMAX),
     &       FREQST(MXRPOS), FREQTR(MXTARG,KMAX,5), 
     &       FREQAR(MXTARG,KMAX,5), RELFLU(5)
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program ME_ROI_3D Version: '//VDATE
      
      IF(1.EQ.1) THEN ! DUMMY block for readability: Read input 1111111
         
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
          IF(FILENM(1:9).NE.'ME_ROI_3D') THEN
            PRINT*, 'Command file appear not to be for ME_ROI_3D but '//
     &              'for '//FILENM//'=> STOP.'
            STOP
          END IF    
        END IF
        
        IF(ASKINP) PRINT*, 'Run in debug mode? (1/0)'
        READ(*,*) I
        DEBUG=(I.EQ.1)

        IF(ASKINP) PRINT*, 'Absorbed dose in Gy'
        READ(*,*) DOSE

        IF(ASKINP) PRINT*, 'Number of files to process'
        READ(*,*) NFILES
        
      END IF          ! DUMMY block for readability: Read input 1111111

      IF(2.EQ.2) THEN ! DUMMY block for readability: Initialize 2222222
        RINFO(1,1)='        R=0'
        RINFO(2,1)='R=0.75R_ROI'
        RINFO(3,1)='    R=R_ROI'
        RINFO(4,1)=' R=1.5R_ROI'
        RINFO(5,1)='   R=2R_ROI'
        RINFO(1,2)='    R<R_ROI'
        RINFO(2,2)='1<R/R_ROI<2'
        RINFO(3,2)='2<R/R_ROI<3'
        RINFO(4,2)='3<R/R_ROI<4'
        RINFO(5,2)='4<R/R_ROI<5'
        NTARG=MXTARG
      END IF          ! DUMMY block for readability: Initialize 2222222
      
      DO IFILE=1,NFILES      
        IF(ASKINP) PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM  
        IF(.NOT.ASKINP) PRINT*, 'Processing file '//FILENM
        IF(ASKINP) PRINT*, 'Enter related proton energy '
        READ(*,*) ENERGY
        STPWRE=STPWR1*EXP(STPEXP*LOG(ENERGY))
        
        IF(3.EQ.3) THEN ! DUMMY block for readability: Init counters 3333
          DO K=1,KMAX
            DO I=1,MXTARG
              FREQSE(I,K)=ZERO
              FREQTE(I,K)=ZERO
              FREQME(I,K)=ZERO
              DO J=1,2
                FNE(I,J,K)=ZERO
              END DO
              DO J=1,5
                FREQTR(I,K,J)=ZERO
                FREQAR(I,K,J)=ZERO
              END DO
            END DO !I=1,MXTARG
            TGMEAN(K)=ZERO            
          END DO ! DO K=1,KMAX
        END IF          ! DUMMY block for readability: Init counters 3333

        IF(4.EQ.4) THEN ! DUMMY block for readability: Get data  444444
          IF(DEBUG) PRINT*, 'Begin of Block 4'
          
          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          DO I=1,2
            READ(LUN,'(A80)') HEADER
          END DO
          READ(LUN,'(10X,10i6)') (NTARGK(I),I=1,KMAX)
          DO I=1,KMAX
            IF(NTARGK(I).GT.0) KMX=I
          END DO
          
          READ(LUN,'(A80)') HEADER
          READ(HEADER(28:36),*) DROI
          READ(HEADER(69:77),*) DBEAM
          NRPOS=1+NINT(DBEAM/DROI*20.)
          IF(NRPOS.GE.41) THEN
            JMAX(1)=5
          ELSE
            JMAX(1)=4
            IF(NRPOS.LT.16) JMAX(1)=1
            IF(NRPOS.LT.21) JMAX(1)=2
            IF(NRPOS.LT.31) JMAX(1)=3
          END IF
          IF(NRPOS.EQ.101) THEN
            JMAX(2)=5
          ELSE
            JMAX(2)=4
            IF(NRPOS.LT.41) JMAX(2)=1
            IF(NRPOS.LT.61) JMAX(2)=2
            IF(NRPOS.LT.81) JMAX(2)=3
          END IF
C!          FLUENC=DROI*DROI*PI/4.*DOSE/DPERFL/STPWRE
C!        Single-event distribution uses 20*RROI as max. impact parameter
          FLUENC=0.25*DBEAM*DBEAM*PI*DOSE/DPERFL/STPWRE
          RELFLU(1)=DROI*DROI/(DBEAM*DBEAM)
          SPOSIN=ONE
          DO I=2,5
            SPOSIN=SPOSIN+2.
            RELFLU(I)=RELFLU(1)*SPOSIN
          END DO
          PRINT*, DOSE,ENERGY,STPWRE,FLUENC,DROI,DBEAM,NRPOS
          
          DO K=1,KMX
            DO J=1,2
              READ(LUN,*) PARAM
            END DO
            DO J=1,NTARGK(K)
              READ(LUN,*) PARAM,ITARG,FREQSE(J,K), FREQTE(J,K), 
     &                    (FREQST(L),L=1,NRPOS)  
              FREQTR(J,K,1)=FREQST(1)
              IF(NRPOS.GE.16) FREQTR(J,K,2)=FREQST(16)
              IF(NRPOS.GE.21) FREQTR(J,K,3)=FREQST(21)
              IF(NRPOS.GE.31) FREQTR(J,K,4)=FREQST(31)
              IF(NRPOS.GE.41) FREQTR(J,K,5)=FREQST(41)
              L=1
              SPOSIN=ZERO
              SUMPOS=ONE
              FREQAR(J,K,L)=FREQST(1)
              DO IJ=20*L-18,20*L+1
                SPOSIN=SPOSIN+8.
                FREQAR(J,K,L)=FREQAR(J,K,L)+SPOSIN*FREQST(IJ)
                SUMPOS=SUMPOS+SPOSIN
              END DO
              FREQAR(J,K,L)=FREQAR(J,K,L)*RELFLU(L)/SUMPOS
              DO L=2,JMAX(2)
                FREQAR(J,K,L)=ZERO
                SUMPOS=ZERO
                DO IJ=20*L-18,20*L+1
                  SPOSIN=SPOSIN+8.
                  FREQAR(J,K,L)=FREQAR(J,K,L)+SPOSIN*FREQST(IJ)
                  SUMPOS=SUMPOS+SPOSIN
                END DO
                FREQAR(J,K,L)=FREQAR(J,K,L)*RELFLU(L)/SUMPOS
              END DO
            END DO
            DO J=1,2
              READ(LUN,*) PARAM
            END DO
          END DO
          CLOSE(LUN)
        END IF          ! DUMMY block for readability: Get data  4444444
          
        IF(5.EQ.5) THEN ! DUMMY block for readability: Main loop  5555555
          IF(DEBUG) PRINT*, 'Begin of Block 5'
          
C!        Init counters       
          DO K=1,KMX
            FREQME(1,K)=EXP(-FLUENC) ! Kronecker's delta for J=1 (0 targets)
          END DO
          
          EVENTS=ONE
          HIFLNC=(FLUENC.GT.23.4)
          IF(HIFLNC) THEN
            DLFLNC=DLOG(FLUENC)
            DLWGHT=-FLUENC+DLFLNC
            WGHT=EXP(DLWGHT)
          ELSE
             WGHT=EXP(-FLUENC)*FLUENC 
          END IF
         
          DO K=1,KMX
            DO J=1,NTARGK(K)
              FNE(J,2,K)=FREQSE(J,K)
              IF(DEBUG) TGMEAN(K)=TGMEAN(K)+FNE(J,2,K)*(J-1)
              FREQME(J,K)=FREQME(J,K)+FNE(J,2,K)*WGHT
            END DO
            NTARGT(K)=NTARGK(K)
          END DO ! DO K=1,KMX
          
          IF(DEBUG) THEN
            IF(WGHT.GT.WMIN) THEN 
              PRINT*, '#',EVENTS,WGHT,TGMEAN(1)/EVENTS,TGMEAN(2)/EVENTS
            ELSE
              PRINT*, '#',EVENTS,WGHT
            END IF
          END IF
          
  10      CONTINUE 
          EVENTS=EVENTS+ONE
          IF(HIFLNC) THEN
            DLWGHT=DLWGHT+DLFLNC-DLOG(EVENTS)
            WGHT=EXP(DLWGHT)
          ELSE
             WGHT=WGHT*FLUENC/EVENTS
          END IF
          
          DO K=1,KMX 
            TGMEAN(K)=ZERO
            DO J=1,NTARGT(K)
              FNE(J,1,K)=FNE(J,2,K)
              FNE(J,2,K)=ZERO
            END DO
            DO J=1,NTARGT(K)
              DO I=1,NTARGK(K)
                IJ=I+J-1
                IF (IJ.GT.MXTARG) IJ=MXTARG
                FNE(IJ,2,K)=FNE(IJ,2,K)+FNE(J,1,K)*FREQSE(I,K)
              END DO
            END DO
            NTARGT(K)=NTARGT(K)+NTARGK(K)
            IF(NTARGT(K).GT.MXTARG) NTARGT(K)=MXTARG
            DO J=1,NTARGT(K)
              FREQME(J,K)=FREQME(J,K)+FNE(J,2,K)*WGHT
              IF(DEBUG) TGMEAN(K)=TGMEAN(K)+FNE(J,2,K)*(J-1)
            END DO
          END DO ! DO K=1,KMX
          
          IF(DEBUG) THEN
            IF(WGHT.GT.WMIN) THEN 
              PRINT*, '#',EVENTS,WGHT,TGMEAN(1),TGMEAN(2)
            ELSE
              PRINT*, '#',EVENTS,WGHT
            END IF
          END IF

          IF(EVENTS.LT.FLUENC.OR.WGHT.GT.WMIN) GOTO 10  ! >>>>>>>>>>
C!        End of inner loop
          
          DO K=1,KMX      
            TGMEAN(K)=ZERO
            DO J=1,NTARGK(K)
              TGMEAN(K)=TGMEAN(K)+FREQSE(J,K)*(J-1)
            END DO
          END DO ! DO K=1,KMX
          PRINT*, 'SE distributions averages (#hit targets) K=1,',KMX
          PRINT*, (TGMEAN(K),K=1,KMX)
          PRINT*, (TGMEAN(K)/EVENTS,K=1,KMX)
          
          DO K=1,KMX      
            TGMEAN(K)=ZERO
            DO J=1,NTARGK(K)
              TGMEAN(K)=TGMEAN(K)+FREQTE(J,K)*(J-1)
            END DO
          END DO ! DO K=1,KMX
          PRINT*, 'TE distributions averages (#hit targets) K=1,',KMX
          PRINT*, (TGMEAN(K),K=1,KMX)
          PRINT*, (TGMEAN(K)/EVENTS,K=1,KMX)
          
          DO K=1,KMX      
            TGMEAN(K)=ZERO
            DO J=1,NTARGT(K)
              TGMEAN(K)=TGMEAN(K)+FREQME(J,K)*(J-1)
            END DO
          END DO ! DO K=1,KMX
          PRINT*, 'ME distributions averages (#hit targets) K=1,',KMX
          PRINT*, (TGMEAN(K),K=1,KMX)
          PRINT*, (TGMEAN(K)/EVENTS,K=1,KMX)
          
        END IF          ! DUMMY block for readability: Main Loop  555555
          
        IF(9.EQ.9) THEN ! DUMMY block for readability: OUTPUT  9999999
C!        #Write results to output file 
          FILOUT='SE_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ME_ROI_3D Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,'(1X,2a8,9(a14,i1))') '#Sites','P1',('F',J,J=2,KMX)
          DO I=1,NTARGK(1)
            SUMPOS=ZERO
            DO K=1,KMX
              IF(SUMPOS.LT.FREQSE(I,K)) SUMPOS=FREQSE(I,K)
            END DO
            IF(SUMPOS.GT.WMIN) NTARG=I
          END DO
          DO I=1,NTARG
            WRITE(LUN,'(1X,i8,10f15.8)') I-1,(FREQSE(I,K),K=1,KMX)
          END DO
          CLOSE(LUN)
          
          FILOUT='TE_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ME_ROI_3D Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,'(1X,2a8,9(a14,i1))') '#Sites','P1',('F',J,J=2,KMX)
          DO I=1,NTARGK(1)
            SUMPOS=ZERO
            DO K=1,KMX
              IF(SUMPOS.LT.FREQTE(I,K)) SUMPOS=FREQTE(I,K)
            END DO
            IF(SUMPOS.GT.WMIN) NTARG=I
          END DO
          DO I=1,NTARG
            WRITE(LUN,'(1X,i8,10f15.6)') I-1,(FREQTE(I,K),K=1,KMX)
          END DO
          CLOSE(LUN)
          
          FILOUT='ME_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ME_ROI_3D Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'DOSE= ',DOSE,' Gy'
          WRITE(LUN,'(1X,2a8,9(a14,i1))') '#Sites','P1',('F',J,J=2,KMX)
          DO I=1,NTARGT(1)
            SUMPOS=ZERO
            DO K=1,KMX
              IF(SUMPOS.LT.FREQME(I,K)) SUMPOS=FREQME(I,K)
            END DO
            IF(SUMPOS.GT.WMIN) NTARG=I
          END DO
          DO I=1,NTARG
            WRITE(LUN,'(1X,i8,10f15.6)') I-1,(FREQME(I,K),K=1,KMX)
          END DO
          CLOSE(LUN)
          
          FILOUT='ALL_'//FILENM  
          KMX=2          
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ME_ROI_3D Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'DOSE= ',DOSE,' Gy'
          
          WRITE(LUN,'(1X,a8,50(a14,i1))') '#Sites',
     &                   'P',1,('F',K,K=2,KMX),'P',1,('F',K,K=2,KMX),
     &                   'P',1,('F',K,K=2,KMX),'P',1,('F',K,K=2,KMX),
     &                   'P',1,('F',K,K=2,KMX),'P',1,('F',K,K=2,KMX),
     &                   'P',1,('F',K,K=2,KMX),'P',1,('F',K,K=2,KMX),
     &                   'P',1,('F',K,K=2,KMX),'P',1,('F',K,K=2,KMX),
     &                   'P',1,('F',K,K=2,KMX),'P',1,('F',K,K=2,KMX),
     &                   'P',1,('F',K,K=2,KMX)

          WRITE(LUN,'(1X,a8,50a15)') 'Data:',('ME',K=1,KMX),
     &                ('SE',K=1,KMX),('TE',K=1,KMX),
     &                ((RINFO(J,2),K=1,KMX),J=1,JMAX(2)),
     &                ((RINFO(J,1),K=1,KMX),J=1,JMAX(1))

          DO I=1,NTARGT(1)
            SUMPOS=ZERO
            DO K=1,KMX
              IF(SUMPOS.LT.FREQSE(I,K)) SUMPOS=FREQSE(I,K)
              IF(SUMPOS.LT.FREQTE(I,K)) SUMPOS=FREQTE(I,K)
              IF(SUMPOS.LT.FREQME(I,K)) SUMPOS=FREQME(I,K)
              DO J=1,JMAX(1)
                IF(SUMPOS.LT.FREQTR(I,K,J)) SUMPOS=FREQTR(I,K,J)
              END DO
              DO J=1,JMAX(2)
                IF(SUMPOS.LT.FREQAR(I,K,J)) SUMPOS=FREQAR(I,K,J)
              END DO
            END DO
            IF(SUMPOS.GT.WMIN) NTARG=I
          END DO
          DO I=1,NTARG
            WRITE(LUN,'(1X,i8,50e15.6)') I-1,(FREQME(I,K),K=1,KMX),
     &                (FREQSE(I,K),K=1,KMX),(FREQTE(I,K),K=1,KMX),
     &                ((FREQAR(I,K,J),K=1,KMX),J=1,JMAX(2))     ,
     &                ((FREQTR(I,K,J),K=1,KMX),J=1,JMAX(1))     
          END DO
          CLOSE(LUN)
          
        END IF          ! DUMMY block for readability: OUTPUT     9999999
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! ME_ROI_3D
C!______________________________________________________________________      

      INCLUDE 'TSTAMP.f'