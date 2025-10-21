      PROGRAM ME_ROI_3P
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
C!    02-MAY-2021 HR:
C!      - Introduced PREFIX to output file name for inclusion of dose.
C!    15-APR-2021 HR:
C!      - Adapted to changes done in ROI_3D (leaner output).
C!    13-APR-2021 HR:
C!      - Created from ME_ROI_3D
C!    ----- Purpose 
C!    Calculate multi-ebent distribution of targets with single and 
C!    and multiple ionization clusters after filtering with Bernoulli
C!    process.    
C!    ----- Parameters: General purpose constants
      INTEGER*4 LUN
      REAL*8 ONE, ZERO, PI, WMIN, WMIN2
      PARAMETER(LUN=11, ONE=1.0, ZERO=0.0, PI=3.1415926, WMIN=1d-10,
     &          WMIN2=1d-40)
C!    ----- Parameters for track and radial distance histograms
      INTEGER*4 MXTARG, MXTRG1, MXTRG2, MXNDOS, MXNCOL
      PARAMETER(MXTARG=65636, MXTRG1=1025, MXTRG2=513, MXNDOS=8, 
     &          MXNCOL=2*MXNDOS+4)
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
      CHARACTER FILENM*120, FILOUT*120, HEADER*120, PREFIX*4
      INTEGER*4 I, IDOSE, IFILE, IJ, IMESSG(2), J, K, L, NCOLS, NDOSES, 
     &          NFILES
      LOGICAL ASKINP, DEBUG(2), HIFLNC
      REAL*8 DBEAM   ! Diameter of beam in nm
      REAL*8 DROI    ! Diameter of ROI in nm
      REAL*8 ENERGY  ! Proton energy in MeV
      REAL*8 EVENTS  ! number of events
      REAL*8 FLUENC  ! average number of tracks per ROI cross section
      REAL*8 PVALUE  ! Probability that cluster volume is target 
      REAL*8 QVALUE  ! Probability that cluster volume is not a target 
      REAL*8 STPWRE  ! Stopping Power at energy E
      REAL*8 DLFLNC, DLWGHT, VCHECK, WGHT
C!    ----- Local arrays
      CHARACTER CDOSES(MXNDOS)*7
      INTEGER*4 ICOL(2), ICOLIC(2), ICOLSE(2)
      INTEGER*4 NMAX(2), NMAXME(2), NMAXSE(2), NMCMAX(MXTRG1), NTARGT(8)
      REAL*8 ADEN(2), ANUM(2), APOT(2), PTOK(2), WEIGH(2)
      REAL*8 DOSE(MXNDOS)    ! Absorbed dose in Gray
      REAL*8 BVARIC(MXTRG1,MXTRG2)
      REAL*8 PMARG(MXTARG,MXNCOL)
      REAL*8 CORRSE(MXTRG1,MXTRG2),BVARSE(MXTRG1,MXTRG2)
      REAL*8 SUMS(MXNCOL),TGMEAN(MXNCOL),FNE(MXTARG,2,2)
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

      PRINT*, 'Program ME_ROI_3P Version: '//VDATE
      
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
          IF(FILENM(1:9).NE.'ME_ROI_3P') THEN
            PRINT*, 'Command file appear not to be for ME_ROI_3P but '//
     &              'for '//FILENM//'=> STOP.'
            STOP
          END IF    
        END IF
        
        IF(ASKINP) THEN
          PRINT*, 'Run in debug mode? (1/0) '
          READ(*,*) IMESSG(1)
          PRINT*, 'Interval between printing intermediate results '//
     &            '(0=suppress this output) '
          READ(*,*) IMESSG(2)
        ELSE
          READ(*,*) (IMESSG(I),I=1,2)
        END IF
        DEBUG(1)=(IMESSG(1).EQ.1)
        DEBUG(2)=(IMESSG(2).GE.1)
        ASKINP=ASKINP.OR.DEBUG(1)
 
        PRINT*, 'Relative target density (0<PVALUE<1)'
        READ(*,*) PVALUE
        QVALUE=ONE-PVALUE
        IF(ASKINP) PRINT*, 'Number of dose values (<=',MXNDOS,')'
        READ(*,*) NDOSES
        IF(NDOSES.GT.MXNDOS) NDOSES=MXNDOS
        IF(ASKINP) PRINT*, 'Absorbed doses in Gy'
        READ(*,*) (DOSE(I),I=1,NDOSES)
        
        IF(ASKINP) PRINT*, 'Number of files to process'
        READ(*,*) NFILES
      END IF          ! DUMMY block for readability: Read input 1111111

      IF(2.EQ.2) THEN ! DUMMY block for readability: Initialize 2222222
        DO I=1,2
          ICOLSE(I)=2*NDOSES+I
          ICOLIC(I)=ICOLSE(I)+2
        END DO
        NMAXME(1)=MXTARG
        NMAXME(2)=MXTARG
      END IF          ! DUMMY block for readability: Initialize 2222222
      
      DO IFILE=1,NFILES      
        IF(ASKINP) PRINT*, 'Enter file name ',IFILE
        READ(*,*) FILENM  
        IF(ASKINP) PRINT*, 'Enter related proton energy '
        READ(*,*) ENERGY
        STPWRE=STPWR1*EXP(STPEXP*LOG(ENERGY))
        
        IF(3.EQ.3) THEN ! DUMMY block for readability: Init counters 3333
          DO K=1,2
            DO I=1,MXTARG
              DO J=1,2
                FNE(I,J,K)=ZERO
              END DO
            END DO !I=1,MXTARG
            TGMEAN(K)=ZERO            
          END DO ! DO K=1,2
        END IF          ! DUMMY block for readability: Init counters 3333

        IF(4.EQ.4) THEN ! DUMMY block for readability: Get data  444444
          IF(DEBUG(1)) PRINT*, 'Begin of Block 4'
          
          OPEN(LUN,FILE=FILENM,STATUS='OLD')
          DO I=1,2
            READ(LUN,'(A80)') HEADER
          END DO
          
          READ(LUN,'(A80)') HEADER
          PRINT*, FILENM
          PRINT*, HEADER
          READ(HEADER(28:36),*) DROI
          READ(HEADER(69:77),*) DBEAM                                      
          PRINT*, (DOSE(J),J=1,NDOSES),ENERGY,STPWRE,DROI,DBEAM
          
          READ(LUN,'(A)') FILOUT
          READ(LUN,*) NMAX(1), NMAX(2)
          
          DO J=1,NMAX(1)
            READ(LUN,*) (BVARIC(J,L),L=1,NMAX(2))
          END DO
          
          CLOSE(LUN)
          
          DO J=1,NMAX(1)
            NMCMAX(J)=0
          END DO
          
          DO J=1,NMAX(1)
            DO L=1,NMAX(2)
              BVARSE(J,L)=ZERO
              CORRSE(J,L)=ZERO
              IF(BVARIC(J,L).GT.WMIN2) THEN
                NMCMAX(J)=L
              END IF
            END DO            
          END DO
          
        END IF          ! DUMMY block for readability: Get data 4444444
          
        IF(5.EQ.5) THEN ! DUMMY block for readability: Convolve 5555555
          IF(DEBUG(1)) PRINT*, 'Begin of Block 5 '
          APOT(1)=ZERO
          PTOK(1)=ONE
          DO I=1,NMAX(1)
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
                      BVARSE(I,J)=BVARSE(I,J)+WEIGH(2)*BVARIC(K,L)
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
        END IF          ! DUMMY block for readability: Convolve 555555
        
        IF(5.EQ.5) THEN ! DUMMY block for readability: Margins 6666666
          IF(DEBUG(1)) PRINT*, 'Begin of Block 6'
          DO I=1,MXTARG
            DO K=1,MXNCOL
              PMARG(I,K)=ZERO
            END DO
          END DO 

          DO I=1,NMAX(1)
            DO J=1,NMAX(2)
              PMARG(I,ICOLIC(1))=PMARG(I,ICOLIC(1))+BVARIC(I,J)
              PMARG(J,ICOLIC(2))=PMARG(J,ICOLIC(2))+BVARIC(I,J)
              IF(BVARSE(I,J).GT.WMIN2) THEN
                NMAXSE(1)=I
                NMAXSE(2)=J
              END IF
            END DO !I=1,NMAX(1)
          END DO !J=1,NMAX(2)
          
          DO I=1,NMAXSE(1)
            DO J=1,NMAXSE(2)
              PMARG(I,ICOLSE(1))=PMARG(I,ICOLSE(1))+BVARSE(I,J)
              PMARG(J,ICOLSE(2))=PMARG(J,ICOLSE(2))+BVARSE(I,J)
            END DO !I=1,NMAXSE(1)
          END DO !J=1,NMAXSE(2)
          
          DO I=1,NMAXSE(1)
            DO J=1,NMAXSE(2)
              IF(PMARG(I,ICOLSE(1)).GT.WMIN2.AND.
     &           PMARG(J,ICOLSE(2)).GT.WMIN2) THEN
                CORRSE(I,J)=BVARSE(I,J)
     &                      /(PMARG(I,ICOLSE(1))*PMARG(J,ICOLSE(2)))
              ELSE
                CORRSE(I,J)=BVARSE(I,J)
              END IF
            END DO !I=1,NMAXSE(1)
          END DO !J=1,NMAXSE(2)
          
        END IF          ! DUMMY block for readability: Margins 6666666

        IF(7.EQ.7) THEN ! DUMMY block for readability: Main loop  777777
          DO IDOSE=1,NDOSES
            IF(DEBUG(1)) PRINT*, 'Begin of Block 6, dose =',DOSE(IDOSE)
          
C!        Init counters       
          
            FLUENC=0.25*DBEAM*DBEAM*PI*DOSE(IDOSE)/DPERFL/STPWRE
            DO K=1,2
              NTARGT(K)=NMAXSE(K)
              ICOL(K)=2*(IDOSE-1)+K
              PMARG(1,ICOL(K))=EXP(-FLUENC) ! Kronecker's delta for J=1 (0 targets)
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
           
             DO K=1,2
              TGMEAN(ICOL(K))=ZERO
              DO I=1,NTARGT(K)
                FNE(I,2,K)=PMARG(I,ICOLSE(K))
                PMARG(I,ICOL(K))=PMARG(I,ICOL(K))+FNE(I,2,K)*WGHT
                IF(DEBUG(2)) TGMEAN(ICOL(K))=TGMEAN(ICOL(K))
     &                                       +FNE(J,2,K)*(J-1)
              END DO
            END DO
          
            IF(DEBUG(1)) THEN
              IF(WGHT.GT.WMIN) THEN 
                PRINT*, '#',EVENTS,WGHT,TGMEAN(ICOL(1))/EVENTS,
     &                  TGMEAN(ICOL(2))/EVENTS
              ELSE
                PRINT*, '#',EVENTS,WGHT
              END IF
            END IF
          
  10        CONTINUE 
            EVENTS=EVENTS+ONE
            IF(HIFLNC) THEN
              DLWGHT=DLWGHT+DLFLNC-DLOG(EVENTS)
              WGHT=EXP(DLWGHT)
            ELSE
               WGHT=WGHT*FLUENC/EVENTS
            END IF
            
            DO K=1,2
              DO I=1,NTARGT(K)
                FNE(I,1,K)=FNE(I,2,K)
                FNE(I,2,K)=ZERO
              END DO
              DO I=1,NTARGT(K)
                DO J=1,NMAXSE(K)
                  IJ=I+J-1
                  IF (IJ.GT.NMAXME(K)) IJ=NMAXME(K)
                  FNE(IJ,2,K)=FNE(IJ,2,K)+FNE(I,1,K)*PMARG(J,ICOLSE(K))
                END DO
              END DO
              NTARGT(K)=NTARGT(K)+NMAXSE(K)
              IF(NTARGT(K).GT.NMAXME(K)) NTARGT(K)=NMAXME(K)
              DO I=1,NTARGT(K)
                PMARG(I,ICOL(K))=PMARG(I,ICOL(K))+FNE(I,2,K)*WGHT
              END DO
            END DO

            IF(DEBUG(2).AND.MOD(NINT(EVENTS),IMESSG(2)).EQ.0) THEN
              IF(WGHT.GT.WMIN) THEN 
                DO K=1,2
                  TGMEAN(ICOL(K))=ZERO
                  DO J=1,NTARGT(K)
                    TGMEAN(ICOL(K))=TGMEAN(ICOL(K))+FNE(J,2,K)*(J-1)
                  END DO
                END DO
                PRINT*, '#',EVENTS,WGHT,TGMEAN(ICOL(1)),TGMEAN(ICOL(2))
              ELSE
                PRINT*, '#',EVENTS,WGHT
              END IF
            END IF

            IF(EVENTS.LT.FLUENC.OR.WGHT.GT.WMIN) GOTO 10
            PRINT*, DOSE
          
          END DO ! IDOSE=1,NDOSES
        END IF          ! DUMMY block for readability: Main Loop  777777
        
        IF(8.EQ.8) THEN ! DUMMY block for readability: Margins 8888888
          IF(DEBUG(1)) PRINT*, 'Begin of Block 8'
          
          NCOLS=2*NDOSES+4
          NMAX(1)=0
          DO K=1,NCOLS    
            TGMEAN(K)=ZERO
            SUMS(K)=ZERO
            IF(K.GT.2) NTARGT(K)=NTARGT(K-2)
            DO J=1,NTARGT(K)
              TGMEAN(K)=TGMEAN(K)+PMARG(J,K)*(J-1)
              SUMS(K)=SUMS(K)+PMARG(J,K)
              IF(PMARG(J,K).GT.WMIN.AND.J.GT.NMAX(1)) NMAX(1)=J
            END DO
          END DO ! DO K=1,NCOLS
          PRINT*, 'Distributions averages (#hit targets)'
          WRITE(*,'(20f10.4)') (SUMS(J),J=1,NCOLS)
          WRITE(*,'(20f10.4)') (TGMEAN(J),J=1,NCOLS)
          WRITE(*,*) TGMEAN(ICOLSE(1))/TGMEAN(ICOLIC(1)),
     &               TGMEAN(ICOLSE(2))/TGMEAN(ICOLIC(2))
          DO K=1,NCOLS-4
            IF(MOD(K,2).EQ.1) THEN
              TGMEAN(K)=TGMEAN(K)/TGMEAN(ICOLSE(1))
            ELSE
              TGMEAN(K)=TGMEAN(K)/TGMEAN(ICOLSE(2))
            END IF
          END DO
          WRITE(*,'(20f10.4)') (TGMEAN(J),J=1,NCOLS-4)
          
        END IF          ! DUMMY block for readability: Margins 8888888

        IF(9.EQ.9) THEN ! DUMMY block for readability: OUTPUT  9999999
C!        #Write results to output file 
          FILOUT='SEB_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ME_ROI_3P Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXSE(1), NMAXSE(2)
          DO I=1,NMAXSE(1)
            WRITE(LUN,'(1000F15.10)') (BVARSE(I,J),J=1,NMAXSE(2))
          END DO
          CLOSE(LUN)
          
          FILOUT='SEC_'//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM SE_ROI_3C Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'Corre1ations P1 and F2 for PVALUE= ',PVALUE
          WRITE(LUN,*) NMAXSE(1), NMAXSE(2)
          DO I=1,NMAXSE(1)
            WRITE(LUN,'(1000F15.10)') (CORRSE(I,J),J=1,NMAXSE(2))
          END DO
          CLOSE(LUN)
          
          PREFIX='MEP_'
          FILOUT=PREFIX//FILENM             
          PRINT*, 'Write output to '//FILOUT   
          OPEN(LUN,FILE=FILOUT,STATUS='UNKNOWN')
          WRITE(LUN,*) '*** Output from PROGRAM ME_ROI_3P Version '
     &                 //VDATE//' on '//TSTAMP()
          WRITE(LUN,*) 'Filename: '//FILOUT
          WRITE(LUN,*) HEADER
          WRITE(LUN,*) 'Single event distribution for PVALUE= ',PVALUE
          
          CDOSES(1)='     P1'
          CDOSES(2)='    P2+'
          WRITE(LUN,'(A9,50a15)') '#',((CDOSES(K),K=1,2),J=1,NDOSES),
     &                            'P1','P2+','P1','P2+'
          PRINT*, DOSE
          DO I=1,NDOSES
            WRITE(CDOSES(I),'(f5.1,a2)') DOSE(I),'Gy'
          END DO
          
          WRITE(LUN,'(1X,a8,50a15)') 'targets',
     &                ((CDOSES(J),K=1,2),J=1,NDOSES),
     &                ('SE',K=1,2),('IC',K=1,2)
          DO I=1,NMAX(1)
            WRITE(LUN,'(1X,i8,50F15.10)') I-1,(PMARG(I,J),J=1,NCOLS)
          END DO
          CLOSE(LUN)
          
        END IF          ! DUMMY block for readability: OUTPUT     9999999
      END DO ! IFILE=1,NFILES
      
      END PROGRAM ! ME_ROI_3P
C!______________________________________________________________________      

      INCLUDE 'TSTAMP.f'