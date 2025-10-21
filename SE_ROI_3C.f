      PROGRAM SE_ROI_3C
C!>>> Declarations >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IMPLICIT NONE
C!    ----- Parameters: Version Date and number
      CHARACTER VDATE*11
      REAL*8 VINPUT  ! Version number for command files YYMMDD.HHMM 
C!                     (date and time  file structure was last changed)
C!    #################################################################
      PARAMETER(VDATE='15-APR-2021',VINPUT=210415.0000d0) ! ###########
C!    15-APR-2021 HR:
C!      - Adapted to changes done in ROI_3D (leaner output).
C!    13-APR-2021 HR:
C!      - Created
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN
      REAL*8 ONE, ZERO, PI, WMIN, WMIN2
      PARAMETER(LUN=11, ONE=1.0, ZERO=0.0, PI=3.1415926, WMIN=1d-10,
     &          WMIN2=1d-40)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 MXTRG2, MXTARG
      PARAMETER(MXTRG2=513, MXTARG=1025)
C!    ----- Local scalars
      CHARACTER FILENM*120, FILOUT*120, HEADER*120
      INTEGER*4 I, IFILE, IJ, J, K, L, NFILES
      LOGICAL ASKINP, DEBUG, HIFLNC
      REAL*8 EVENTS  ! number of events
      REAL*8 PVALUE  ! Probability that cluster volume is target 
      REAL*8 QVALUE  ! Probability that cluster volume is not a target 
      REAL*8 SPOSIN, SUMPOS, VCHECK, WGHT
C!    ----- Local arrays
      INTEGER*4 NMAX(2), NMAXSE(2), NMCMAX(MXTARG), NTARGT(4)
      REAL*8 ADEN(2), ANUM(2), APOT(2), PTOK(2), WEIGH(2)
      REAL*8 CORRIC(MXTARG,MXTRG2), CORRSE(MXTARG,MXTRG2)
      REAL*8 CVARSE(MXTARG,MXTRG2), PMARG(MXTARG,4)
      REAL*8 SUMS(4),TGMEAN(4),FNE(MXTARG,2,2)
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program SE_ROI_3C Version: '//VDATE
      
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
        END IF
        
        IF(ASKINP) PRINT*, 'Run in debug mode? (1/0)'
        READ(*,*) I
        DEBUG=(I.EQ.1)
        ASKINP=ASKINP.OR.DEBUG
 
        IF(ASKINP) PRINT*, 'Relative target density (0<PVALUE<1)'
        READ(*,*) PVALUE
        QVALUE=ONE-PVALUE
        
        IF(ASKINP) PRINT*, 'Number of files to process'
        READ(*,*) NFILES
      END IF          ! DUMMY block for readability: Read input 1111111

      IF(2.EQ.2) THEN ! DUMMY block for readability: Initialize 2222222
      END IF          ! DUMMY block for readability: Initialize 2222222
      
      DO IFILE=1,NFILES      
        IF(ASKINP) PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM  
        
        IF(3.EQ.3) THEN ! DUMMY block for readability: Get data  333333
          IF(DEBUG) PRINT*, 'Hier beginnt Block 3'
          
          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          DO I=1,3
            READ(LUN,'(A80)') HEADER
          END DO
          
          READ(LUN,'(A)') FILOUT
          READ(LUN,*) NMAX(1), NMAX(2)
          
          DO J=1,NMAX(1)
            READ(LUN,*) (CORRIC(J,L),L=1,NMAX(2))
          END DO
          
          CLOSE(LUN)
          
          DO J=1,NMAX(1)
            NMCMAX(J)=0
            DO L=1,NMAX(2)
              CVARSE(J,L)=ZERO
              CORRSE(J,L)=ZERO
              IF(CORRIC(J,L).GT.WMIN2) THEN
                NMCMAX(J)=L
              END IF
            END DO            
          END DO
          
        END IF          ! DUMMY block for readability: Get data 3333333
          
        IF(4.EQ.4) THEN ! DUMMY block for readability: Convolve 4444444
          IF(DEBUG) PRINT*, 'Hier beginnt Block 4 '
          APOT(1)=ZERO
          PTOK(1)=ONE
          DO I=1,NMAX(1)
            IF(DEBUG.AND.MOD(I,50).EQ.0) PRINT*, 'Line ',I,' of ',
     &                                            NMAX(1)
            APOT(2)=ZERO
            PTOK(2)=ONE
            DO J=1,NMAX(2)
              ANUM(1)=APOT(1)
              ADEN(1)=ZERO
              WEIGH(1)=PTOK(1)
              DO K=I,NMAX(1)
                IF(WEIGH(1).GT.WMIN2) THEN
                  ANUM(2)=APOT(2)
                  ADEN(2)=ZERO
                  WEIGH(2)=PTOK(2)*WEIGH(1)
                  DO L=J,NMCMAX(K)
                    IF(WEIGH(2).GT.WMIN2) THEN
                      CVARSE(I,J)=CVARSE(I,J)+WEIGH(2)*CORRIC(K,L)
                      ADEN(2)=ADEN(2)+ONE
                      ANUM(2)=ANUM(2)+ONE
                      WEIGH(2)=WEIGH(2)*QVALUE*ANUM(2)/ADEN(2)
                    END IF
                  END DO ! L=J,NMAX(2)
                  ADEN(1)=ADEN(1)+ONE
                  ANUM(1)=ANUM(1)+ONE
                  WEIGH(1)=WEIGH(1)*QVALUE*ANUM(1)/ADEN(1)
                END IF
              END DO ! K=I,NMAX(1)
              APOT(2)=APOT(2)+ONE
              PTOK(2)=PTOK(2)*PVALUE
            END DO !J=1,NMAX(2)
            APOT(1)=APOT(1)+ONE
            PTOK(1)=PTOK(1)*PVALUE
          END DO !I=1,NMAX(1)
        END IF          ! DUMMY block for readability: Convolve 444444
        
        IF(5.EQ.5) THEN ! DUMMY block for readability: Margins 5555555
          IF(DEBUG) PRINT*, 'Hier beginnt Block 5'
          DO I=1,NMAX(1)
            DO K=1,4
              PMARG(I,K)=ZERO
            END DO
          END DO
          DO I=1,NMAX(1)
            DO J=1,NMAX(2)
              PMARG(I,1)=PMARG(I,1)+CVARSE(I,J)
              PMARG(J,2)=PMARG(J,2)+CVARSE(I,J)
              PMARG(I,3)=PMARG(I,3)+CORRIC(I,J)
              PMARG(J,4)=PMARG(J,4)+CORRIC(I,J)
              IF(CVARSE(I,J).GT.WMIN) THEN
                NMAXSE(1)=I
                NMAXSE(2)=J
              END IF
              IF(PMARG(I,1).GT.WMIN2.AND.PMARG(J,2).GT.WMIN2) THEN
                CORRSE(I,J)=CVARSE(I,J)/(PMARG(I,1)*PMARG(J,2))
              ELSE
                CORRSE(I,J)=CVARSE(I,J)
              END IF
            END DO !J=1,NMAX(2)
          END DO !I=1,NMAX(1)
          
          DO K=1,4      
            TGMEAN(K)=ZERO
            SUMS(K)=ZERO
            DO J=1,NMAX(1)
              TGMEAN(K)=TGMEAN(K)+PMARG(J,K)*(J-1)
              SUMS(K)=SUMS(K)+PMARG(J,K)
            END DO
          END DO ! DO K=1,4
          PRINT*, 'Distributions averages (#hit targets)'
          WRITE(*,'(10f10.4)') (SUMS(J),J=1,4)
          WRITE(*,'(10f10.4)') (TGMEAN(J),J=1,4)
          WRITE(*,*) TGMEAN(1)/TGMEAN(3),TGMEAN(2)/TGMEAN(4)
          
        END IF          ! DUMMY block for readability: Margins 5555555

        IF(9.EQ.9) THEN ! DUMMY block for readability: OUTPUT  9999999
C!        #Write results to output file 
          FILOUT='SEB_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM SE_ROI_3C Version '
     &                 //VDATE
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXSE(1), NMAXSE(2)
          DO I=1,NMAXSE(1)
            WRITE(LUN,'(1000F15.10)') (CVARSE(I,J),J=1,NMAXSE(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='SEC_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM SE_ROI_3C Version '
     &                 //VDATE
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXSE(1), NMAXSE(2)
          DO I=1,NMAXSE(1)
            WRITE(LUN,'(1000F15.10)') (CORRSE(I,J),J=1,NMAXSE(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='SED_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM SE_ROI_3C Version '
     &                 //VDATE
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'Single event distribution for PVALUE= ',PVALUE
          WRITE(LUN,'(A9,12a15)') '#CVs','P1_SEPIC','F2_SEPIC',
     &                             'P1_SEIC','F2_SEIC'
         DO I=1,NMAX(1)
            WRITE(LUN,'(1X,i8,10F15.10)') I-1,(PMARG(I,J),J=1,4)
          END DO
          CLOSE(LUN)
          
        END IF          ! DUMMY block for readability: OUTPUT     9999999
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! SE_ROI_3C
C!______________________________________________________________________      

