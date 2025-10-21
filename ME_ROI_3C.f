      PROGRAM ME_ROI_3C
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='06-JUN-2021',VINPUT=210605.1600d0) ! ###########
C!    06-JUN-2021 HR:
C!      - Added timestamp in output file
C!    05-JUN-2021 HR:
C!      - Reorganized command input file.
C!    17-APR-2021 HR:
C!      - Cleaned up debug messages
C!    16-APR-2021 HR:
C!      - Adapted to changes done in ROI_3D (leaner output).
C!    13-APR-2021 HR:
C!      - Created
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN
      REAL*8 ONE, ZERO, PI, WMIN, WMIN2
      PARAMETER(LUN=11, ONE=1.0, ZERO=0.0, PI=3.1415926, WMIN=1d-10,
     &          WMIN2=1d-40)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 MXTIC1, MXTIC2, MXTRG1, MXTRG2
      PARAMETER(MXTIC1=1025, MXTIC2=513, MXTRG1=129, MXTRG2=65)
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
      CHARACTER FILENM*120, FILOUT*120, HEADER*120
      INTEGER*4 I, I1, I2, IFILE, IJ, J, K, L, NFILES
      LOGICAL ASKINP, DEBUG, HIFLNC
      REAL*8 DBEAM   ! Diameter of beam in nm
      REAL*8 DOSE    ! Absorbed dose in Gray
      REAL*8 DROI    ! Diameter of ROI in nm
      REAL*8 ENERGY  ! Proton energy in MeV
      REAL*8 EVENTS  ! number of events
      REAL*8 FLUENC  ! average number of tracks per ROI cross section
      REAL*8 PVALUE  ! Probability that cluster volume is target 
      REAL*8 QVALUE  ! Probability that cluster volume is not a target 
      REAL*8 STPWRE  ! Stopping Power at energy E
      REAL*8 DLFLNC, DLWGHT, VCHECK, WGHT
C!    ----- Local arrays
      INTEGER*4 NMAXIC(2), NMAXME(2), NMAXSE(2), 
     &          NMCMAX(MXTIC1), NTARGT(8)
      REAL*8 ADEN(2), ANUM(2), APOT(2), PTOK(2), WEIGH(2)
      REAL*8 CORRIC(MXTIC1,MXTIC2)
      REAL*8 PMARG(MXTIC1,8)
      REAL*8 CORRSE(MXTRG1,MXTRG2),CVARSE(MXTRG1,MXTRG2)
      REAL*8 CORRME(MXTRG1,MXTRG2),CVARME(MXTRG1,MXTRG2)
      REAL*8 CTEMP1(MXTRG1,MXTRG2),CTEMP2(MXTRG1,MXTRG2)
      REAL*8 SUMS(8),TGMEAN(8),FNE(MXTRG1,2,2)
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program ME_ROI_3C Version: '//VDATE
      
      IF(1.EQ.1) THEN ! DUMMY block for readability: Read input 1111111
C!    11111111111111111111111111111111111111111111111111111111111111111
C!    Check command file version
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
        IF(FILENM(1:9).NE.'ME_ROI_3C') THEN
          PRINT*, 'Command file appear not to be for ME_ROI_3C but '//
     &              'for '//FILENM//'=> STOP.'
          STOP
        END IF    
      END IF
      
      IF(ASKINP) PRINT*, 'Run in debug mode? (1/0)'
      READ(*,*) I
      DEBUG=(I.EQ.1)
      ASKINP=ASKINP.OR.DEBUG

      PRINT*, 'Relative target density (0<PVALUE<1)'
      READ(*,*) PVALUE
      QVALUE=ONE-PVALUE
      IF(ASKINP) PRINT*, 'Absorbed dose in Gy'
      READ(*,*) DOSE

      IF(ASKINP) PRINT*, 'Size of multi-event bivariate distribution'
      READ(*,*) NMAXME(1), NMAXME(2)
      IF(NMAXME(1).GT.MXTRG1.OR.NMAXME(1).LE.1) NMAXME(1)=MXTRG1
      IF(NMAXME(2).GT.MXTRG2.OR.NMAXME(2).LE.1) NMAXME(2)=MXTRG2

      IF(ASKINP) PRINT*, 'Number of files to process'
      READ(*,*) NFILES
C!    11111111111111111111111111111111111111111111111111111111111111111
      END IF          ! DUMMY block for readability: Read input 1111111

      IF(2.EQ.2) THEN ! DUMMY block for readability: Initialize 2222222
C!    22222222222222222222222222222222222222222222222222222222222222222
C!    22222222222222222222222222222222222222222222222222222222222222222
      END IF          ! DUMMY block for readability: Initialize 2222222
      
      DO IFILE=1,NFILES      
        IF(ASKINP) PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM  
        IF(ASKINP) PRINT*, 'Enter related proton energy '
        READ(*,*) ENERGY
        STPWRE=STPWR1*EXP(STPEXP*LOG(ENERGY))
        
        IF(3.EQ.3) THEN ! DUMMY block for readability: Get data  333333
C!      333333333333333333333333333333333333333333333333333333333333333
        IF(DEBUG) PRINT*, 'Begin of Block 3'
        
        OPEN(LUN,FILE=FILENM,STATUS='OLD')
        DO I=1,2
          READ(LUN,'(A80)') HEADER
        END DO
        
        READ(LUN,'(A80)') HEADER
        PRINT*, HEADER
        READ(HEADER(28:36),*) DROI
        READ(HEADER(69:77),*) DBEAM                                      
        FLUENC=0.25*DBEAM*DBEAM*PI*DOSE/DPERFL/STPWRE
        PRINT*, DOSE,ENERGY,STPWRE,FLUENC,DROI,DBEAM
        
        READ(LUN,'(A)') FILOUT
        READ(LUN,*) NMAXIC(1), NMAXIC(2)
        IF(NMAXIC(1).GT.MXTIC1) THEN
          PRINT*, 'Problem number of lines in '//FILENM//' is ',
     &              NMAXIC(1),' > MXTIC1 =',MXTIC1
          STOP
        END IF
        IF(NMAXIC(2).GT.MXTIC2) THEN
          PRINT*, 'Problem number of columns in '//FILENM//' is ',
     &              NMAXIC(2),' > MXTIC2 =',MXTIC2
          STOP
        END IF
        
        DO J=1,NMAXIC(1)
          READ(LUN,*) (CORRIC(J,L),L=1,NMAXIC(2))
        END DO
        
        CLOSE(LUN)
C!      333333333333333333333333333333333333333333333333333333333333333
        END IF          ! DUMMY block for readability: Get data 3333333
          
        IF(4.EQ.4) THEN ! DUMMY block for readability: zero arrays 4444
C!      444444444444444444444444444444444444444444444444444444444444444
        IF(DEBUG) PRINT*, 'Begin of Block 4 '
        DO I=1,NMAXIC(1)
          NMCMAX(I)=0
          DO J=1,NMAXIC(2)
            IF(CORRIC(I,J).GT.WMIN2) NMCMAX(I)=J
          END DO            
        END DO
        
        DO I=1,MXTIC1
          DO K=1,8
            PMARG(I,K)=ZERO
          END DO
        END DO
        
        DO I=1,8
          NTARGT(I)=0
        END DO
        
        DO I=1,MXTRG1
          DO J=1,MXTRG2
            CVARSE(I,J)=ZERO
            CORRSE(I,J)=ZERO
            CVARME(I,J)=ZERO
            CORRME(I,J)=ZERO
            CTEMP1(I,J)=ZERO
            CTEMP2(I,J)=ZERO
          END DO
          DO J=1, 2
            DO K=1,2
              FNE(I,J,K)=ZERO
            END DO
          END DO
        
        END DO ! DO I=1,NMAXIC(1)
C!      444444444444444444444444444444444444444444444444444444444444444
        END IF          ! DUMMY block for readability: zero arrays 4444

        IF(5.EQ.5) THEN ! DUMMY block for readability: Convolve 5555555
C!      555555555555555555555555555555555555555555555555555555555555555
          IF(DEBUG) PRINT*, 'Begin of Block 5 '
          APOT(1)=ZERO
          PTOK(1)=ONE
          I1=0
          DO I=1,NMAXIC(1)
            IF(I1.LT.NMAXME(1)) I1=I1+1
            APOT(2)=ZERO
            PTOK(2)=ONE
            I2=0
            DO J=1,NMAXIC(2)
              IF(I2.LT.NMAXME(2)) I2=I2+1
              ANUM(1)=APOT(1)
              ADEN(1)=ZERO
              WEIGH(1)=PTOK(1)
              DO K=I,NMAXIC(1)
                IF(WEIGH(1).GT.WMIN2) THEN
                  ANUM(2)=APOT(2)
                  ADEN(2)=ZERO
                  WEIGH(2)=PTOK(2)*WEIGH(1)
                  DO L=J,NMCMAX(K)
                    IF(WEIGH(2).GT.WMIN2) THEN
                      CVARSE(I1,I2)=CVARSE(I1,I2)+WEIGH(2)*CORRIC(K,L)
                      ADEN(2)=ADEN(2)+ONE
                      ANUM(2)=ANUM(2)+ONE
                      WEIGH(2)=WEIGH(2)*QVALUE*ANUM(2)/ADEN(2)
                    END IF
                  END DO ! L=J,NMAXIC(2)
                  ADEN(1)=ADEN(1)+ONE
                  ANUM(1)=ANUM(1)+ONE
                  WEIGH(1)=WEIGH(1)*QVALUE*ANUM(1)/ADEN(1)
                END IF
              END DO ! K=I,NMAXIC(1)
              APOT(2)=APOT(2)+ONE
              PTOK(2)=PTOK(2)*PVALUE
            END DO !J=1,NMAXIC(2)
            APOT(1)=APOT(1)+ONE
            PTOK(1)=PTOK(1)*PVALUE
          END DO !I=1,NMAXIC(1)
C!      555555555555555555555555555555555555555555555555555555555555555
        END IF          ! DUMMY block for readability: Convolve 555555
        
        IF(6.EQ.6) THEN ! DUMMY block for readability: Margins 6666666
C!      66666666666666666666666666666666666666666666666666666666666666
        IF(DEBUG) PRINT*, 'Begin of Block 6'
C!        marginal distributions of targets with one IC or multi ICs
        DO I=1,NMAXIC(1)
          DO J=1,NMAXIC(2)
            PMARG(I,5)=PMARG(I,5)+CORRIC(I,J)
            PMARG(J,6)=PMARG(J,6)+CORRIC(I,J)
          END DO !I=1,NMAXIC(1)
        END DO !J=1,NMAXIC(2)
        NTARGT(5)=NMAXIC(1)
        NTARGT(6)=NMAXIC(2)
C!        marginal distributions of targets with one DSB or multi DSBs
        DO I=1,NMAXME(1)
          DO J=1,NMAXME(2)
            PMARG(I,3)=PMARG(I,3)+CVARSE(I,J)
            PMARG(J,4)=PMARG(J,4)+CVARSE(I,J)
            IF(CVARSE(I,J).GT.WMIN2) THEN
              NMAXSE(1)=I
              NMAXSE(2)=J
            END IF
          END DO !I=1,NMAXSE(1)
        END DO !J=1,NMAXSE(2)
        NTARGT(3)=NMAXSE(1)
        NTARGT(4)=NMAXSE(1)
C!        "Relative correlation matrix" - ratio bivariate frequency to
C!        product of marginal frequencies
        DO I=1,NMAXSE(1)
          DO J=1,NMAXSE(2)
            IF(PMARG(I,3).GT.WMIN2.AND.PMARG(J,4).GT.WMIN2) THEN
              CORRSE(I,J)=CVARSE(I,J)/(PMARG(I,3)*PMARG(J,4))
            ELSE
              CORRSE(I,J)=CVARSE(I,J)
            END IF
          END DO !I=1,NMAXSE(1)
        END DO !J=1,NMAXSE(2)
C!      66666666666666666666666666666666666666666666666666666666666666
        END IF          ! DUMMY block for readability: Margins 6666666


        IF(7.EQ.7) THEN ! DUMMY block for readability: Main loop  7777777
C!      777777777777777777777777777777777777777777777777777777777777777
        IF(DEBUG) PRINT*, 'Begin of Block 7'
        NTARGT(1)=NMAXSE(1)
        NTARGT(2)=NMAXSE(2)
        
        PMARG(1,1)=EXP(-FLUENC) ! Kronecker's delta for J=1 (0 targets)
        PMARG(1,2)=EXP(-FLUENC) ! Kronecker's delta for J=1 (0 targets)
        CVARME(1,1)=EXP(-FLUENC)
        
        EVENTS=ONE
        HIFLNC=(FLUENC.GT.23.4)
        IF(HIFLNC) THEN
          DLFLNC=DLOG(FLUENC)
          DLWGHT=-FLUENC+DLFLNC
          WGHT=EXP(DLWGHT)
        ELSE
          WGHT=EXP(-FLUENC)*FLUENC 
        END IF
       
        DO I=1,NTARGT(1)
          DO J=1,NTARGT(2)
            CTEMP2(I,J)=CVARSE(I,J)
            CVARME(I,J)=CVARME(I,J)+CTEMP2(I,J)*WGHT
          END DO
        END DO  
        
        DO K=1,2
          TGMEAN(K)=ZERO
          DO I=1,NTARGT(K)
            FNE(I,2,K)=PMARG(I,K+2)
            PMARG(I,K)=PMARG(I,K)+FNE(I,2,K)*WGHT
            IF(DEBUG) TGMEAN(K)=TGMEAN(K)+FNE(I,2,K)*(I-1)
          END DO
        END DO
                
        IF(DEBUG) THEN
          IF(WGHT.GT.WMIN) THEN 
            PRINT*, '#',EVENTS,WGHT,TGMEAN(1)/EVENTS,TGMEAN(2)/EVENTS
          ELSE
            PRINT*, '#',EVENTS,WGHT
          END IF
        END IF
        
          
  10    CONTINUE 
        EVENTS=EVENTS+ONE
        IF(HIFLNC) THEN
          DLWGHT=DLWGHT+DLFLNC-DLOG(EVENTS)
          WGHT=EXP(DLWGHT)
        ELSE
           WGHT=WGHT*FLUENC/EVENTS
        END IF
        
        DO K=1,2
          TGMEAN(K)=ZERO
          DO I=1,NTARGT(K)
            FNE(I,1,K)=FNE(I,2,K)
            FNE(I,2,K)=ZERO
          END DO
          DO I=1,NTARGT(K)
            DO J=1,NMAXSE(K)
              IJ=I+J-1
              IF (IJ.GT.NMAXME(K)) IJ=NMAXME(K)
              FNE(IJ,2,K)=FNE(IJ,2,K)+FNE(I,1,K)*PMARG(J,K+2)
            END DO
          END DO
          NTARGT(K)=NTARGT(K)+NMAXSE(K)
          IF(NTARGT(K).GT.NMAXME(K)) NTARGT(K)=NMAXME(K)
          DO I=1,NTARGT(K)
            PMARG(I,K)=PMARG(I,K)+FNE(I,2,K)*WGHT
            IF(DEBUG) TGMEAN(K)=TGMEAN(K)+FNE(I,2,K)*(I-1)
          END DO
        END DO

        DO I=1,NTARGT(1)
          DO J=1,NTARGT(2)
            CTEMP1(I,J)=CTEMP2(I,J)
            CTEMP2(I,J)=ZERO
          END DO
        END DO
        
        DO I=1,NTARGT(1)
          DO J=1,NTARGT(2)
            DO K=1,I
              DO L=1,J
                CTEMP2(I,J)=CTEMP2(I,J)+
     &                        CTEMP1(K,L)*CVARSE(I+1-K,J+1-L)
              END DO
            END DO
            CVARME(I,J)=CVARME(I,J)+CTEMP2(I,J)*WGHT
          END DO
        END DO

        IF(DEBUG) THEN
          IF(WGHT.GT.WMIN) THEN 
            PRINT*, '#',EVENTS,WGHT,TGMEAN(1),TGMEAN(2)
          ELSE
            PRINT*, '#',EVENTS,WGHT
          END IF
        END IF

        IF(EVENTS.LT.FLUENC.OR.WGHT.GT.WMIN) GOTO 10
C!      777777777777777777777777777777777777777777777777777777777777777
        END IF          ! DUMMY block for readability: Main Loop  777777
        
        IF(8.EQ.8) THEN ! DUMMY block for readability: Margins 8888888
C!      88888888888888888888888888888888888888888888888888888888888888
        IF(DEBUG) PRINT*, 'Begin of Block 8'
        DO I=1,NMAXME(1)
          DO J=1,NMAXME(2)
            PMARG(I,7)=PMARG(I,7)+CVARME(I,J)
            PMARG(J,8)=PMARG(J,8)+CVARME(I,J)
            IF(CVARME(I,J).GT.WMIN2) THEN
              NTARGT(1)=I
              NTARGT(2)=J
            END IF
          END DO !J=1,NMAXME(2)
        END DO !I=1,NMAXME(1)
        
        DO I=1,NTARGT(1)
          DO J=1,NTARGT(2)
            IF(PMARG(I,7).GT.WMIN2.AND.PMARG(J,8).GT.WMIN2) THEN
              CORRME(I,J)=CVARME(I,J)/(PMARG(I,7)*PMARG(J,8))
            ELSE
              CORRME(I,J)=CVARME(I,J)
            END IF
          END DO !I=1,NTARGT(1)
        END DO !J=1,NTARGT(2)

        DO K=1,8      
          TGMEAN(K)=ZERO
          SUMS(K)=ZERO
          IF(K.GT.6) NTARGT(K)=NTARGT(K-6)
          DO J=1,NTARGT(K)
            TGMEAN(K)=TGMEAN(K)+PMARG(J,K)*(J-1)
            SUMS(K)=SUMS(K)+PMARG(J,K)
          END DO
        END DO ! DO K=1,8
        PRINT*, 'Distributions averages (#hit targets)'
        WRITE(*,'(10f10.4)') (SUMS(J),J=1,8)
        WRITE(*,'(10f10.4)') (TGMEAN(J),J=1,8)
        WRITE(*,*) TGMEAN(3)/TGMEAN(7),TGMEAN(4)/TGMEAN(8)
C!      88888888888888888888888888888888888888888888888888888888888888
        END IF          ! DUMMY block for readability: Margins 8888888

        IF(9.EQ.9) THEN ! DUMMY block for readability: OUTPUT  9999999
C!      99999999999999999999999999999999999999999999999999999999999999
C!        #Write results to output file 
          FILOUT='SEB_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM ME_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) ' Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXSE(1), NMAXSE(2)
          DO I=1,NMAXSE(1)
            WRITE(LUN,'(1000D15.6)') (CVARSE(I,J),J=1,NMAXSE(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='SEC_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM SE_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) ' Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXSE(1), NMAXSE(2)
          DO I=1,NMAXSE(1)
            WRITE(LUN,'(1000D15.6)') (CORRSE(I,J),J=1,NMAXSE(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='MEA_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM ME_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) ' Multi and single event distribution for '//
     &                 'PVALUE= ',PVALUE,'  DOSE= ',DOSE,' Gy'
          WRITE(LUN,*) ' ***  Compared with distributions for ICs ***'
          
          WRITE(LUN,'(A9,12a15)') '#CVs','P1','F2','P1','F2','P1','F2',
     &                            'P1','F2'

          WRITE(LUN,'(1X,a8,50a15)') 'Data:',('ME',K=1,2),
     &                ('SE',K=1,2),('IC',K=1,2),('ME',K=1,2)
 
          DO I=1,NMAXIC(1)
            WRITE(LUN,'(1X,i8,10F15.10)') I-1,(PMARG(I,J),J=1,8)
          END DO
          CLOSE(LUN)

          FILOUT='MEB_'//FILENM             
          PRINT*, ' Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM ME_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) ' Bivariate freq. P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXME(1), NMAXME(2)
          DO I=1,NMAXME(1)
            WRITE(LUN,'(1000D15.6)') (CVARME(I,J),J=1,NMAXME(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='MEC_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM ME_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) ' Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXME(1), NMAXME(2)
          DO I=1,NMAXME(1)
            WRITE(LUN,'(1000D15.6)') (CORRME(I,J),J=1,NMAXME(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='MED_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) ' *** Output from PROGRAM ME_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) ' Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) ' Multi and single event distribution for '//
     &                 'PVALUE= ',PVALUE,'  DOSE= ',DOSE,' Gy'
          
          WRITE(LUN,'(A9,12a15)') '#CVs','P1','F2','P1','F2'
          WRITE(LUN,'(1X,a8,50a15)') 'Data:',('ME',K=1,2),
     &                ('SE',K=1,2)
 
          DO I=1,NTARGT(1)
            WRITE(LUN,'(1X,i8,10F15.10)') I-1,(PMARG(I,J),J=1,4)
          END DO
          CLOSE(LUN)
C!      99999999999999999999999999999999999999999999999999999999999999
        END IF          ! DUMMY block for readability: OUTPUT     9999999
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! ME_ROI_3C
C!______________________________________________________________________      

      INCLUDE 'TSTAMP.f'