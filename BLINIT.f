C!**********************************************************************
      SUBROUTINE BLINIT(NAZMTH,NDIV,DLATC,SRDATE)
C!    ***************************************************************
C!    Initialization of the base vectors of the reziprocal lattices
C!    ***************************************************************
C!    This routine calculates the elements of the array ADJNTB in the
C!    common block LATTICE, i.e. the set of base vectors of the 
C!    reciprocal lattice of a cubic fcc Bravais lattice for different 
C!    orientations of the primary particle's trajectory with respect to 
C!    the Bravais lattice of target volumes.
C!    25.02.2020: The algorithm for deterministic sampling has been 
C!    implemented such that for NDIV>0 the spherical triangle defined by  
C!    the north pole and the azimuths +/- 2/3*pi on the equator is 
C!    divided into four congruent triangles, then these smaller 
C!    triangles are again divided into four congruent triangles etc.
C!    until finally we have 4^NDIV triangles. The centers of gravity of 
C!    these triangles define directions of the primary particle's 
C!    trajectory that are homogenously distributed over the section of 
C!    the unit sphere that is sufficient to consider owing to the 
C!    symmetry of the cubic fcc Bravais lattice of target volumes. 
C!    In addition, also a rotation of the track in itself is considered
C!    to improve statistics.
C!    Declarations ---------------------------------------------------
      IMPLICIT NONE
C!    ----- Parameters: Version Date <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      CHARACTER VDATE*11
      PARAMETER(VDATE='23-MAY-2021') !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C!    23-MAY-2021 HR: Added outout of version date to calling program
C!    23-OCT-2020 HR: Fixed bug in initialization of reciprocal basis
C!    15-MAR-2020 HR: Added VDATE and modified COMMON BLOCK VERBOSE
C!    ----- Parameters: General purpose constants
      REAL*8 ONE, TWO, THREE, FOUR, ZERO
      PARAMETER(ONE=1.0, TWO=2.0, THREE=3.0, FOUR=4.0, ZERO=0.0)
C!    ----- Parameters for lattice orientation
      INTEGER*4 MAZMTH, MDIV, MTHETA, MDIR
      PARAMETER(MAZMTH=1, MDIV=0, MTHETA=2**MDIV,  
     &           MDIR=MAZMTH*(MTHETA*(MTHETA+1))/2)
C!    ----- Local scalars
      REAL*8 PIBY2, PIBY3, ROOT2, ROOT3 ! Constants
      REAL*8 COSTHM, COSTHP, DCOSTH, DPHI2, DPHI1, PHI1, PHI2, THETA
      INTEGER*4 IAZ, IDIR, IR, IS, J, K, L, NTHETA
C!    ----- Local arrays
      REAL*8 ADJBAS(3,3), REULER(3,3)
C!    ----- Global variables
      INTEGER*4 NDIR
      REAL*8 ADJNTB(3,3,MDIR)
      COMMON /LATTICE/ ADJNTB, NDIR 
C!    ----- Scalars and arrays passed from the calling module
      INTEGER*4 NAZMTH, NDIV
      REAL*8 DLATC
C!    ----- Scalars and arrays passed to the calling module
      CHARACTER*11 SRDATE
C!    ------
      LOGICAL DEBUG(4), SRFLAG
      COMMON /VERBOSE/ DEBUG, SRFLAG
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      
      IF(1.EQ.1) THEN ! Dummy block: Checks for consistency 
        PRINT*, 'SUBROUTINE BLINIT Version: '//VDATE
        IF(DEBUG(1)) PRINT*, 'BLINIT 1'
        IF(NAZMTH.GT.MAZMTH) THEN
          PRINT*, 'ERROR in BLINIT: NAZMTH > MAZMTH ',NAZMTH,MAZMTH
          STOP
        END IF
        IF(NDIV.GT.MDIV) THEN
          PRINT*, 'ERROR in BLINIT: NDIV > MDIV ',NDIV,MDIV
          STOP
        END IF
        SRDATE=VDATE
      END IF          ! Dummy block: Checks for consistency         
      
      IF(2.EQ.2) THEN ! Dummy block for readability: INITIALIZATION
        IF(DEBUG(1)) PRINT*, 'BLINIT 2'
        ROOT2=SQRT(2.0)
        ROOT3=SQRT(3.0)
        ADJBAS(1,1)=ONE/DLATC/ROOT3
        ADJBAS(2,1)=-ONE/DLATC
        ADJBAS(3,1)=-ONE/DLATC/ROOT2/ROOT3
        ADJBAS(1,2)=ONE/DLATC/ROOT3
        ADJBAS(2,2)=ONE/DLATC
        ADJBAS(3,2)=-ONE/DLATC/ROOT2/ROOT3
        ADJBAS(1,3)=ZERO
        ADJBAS(2,3)=ZERO                   ! fixed Bug 1st index was 1 23-OCT-2020
        ADJBAS(3,3)=ONE/DLATC*ROOT3/ROOT2  ! fixed Bug 1st index was 1 23-OCT-2020
        
        PIBY2=ATAN(ONE)*TWO
        PIBY3=PIBY2*TWO/THREE
      END IF          ! Dummy block for readability: INITIALIZATION

      PRINT*, ADJBAS
      IF(3.EQ.3) THEN ! Dummy block for readability: Get track direction
                      ! and adjust lattice orientation accordingly
        IF(DEBUG(1)) PRINT*, 'BLINIT 3'
        NTHETA=2**NDIV                ! number of polar angles
        NDIR=NAZMTH*NTHETA*NTHETA     ! number of orientations
        DCOSTH=ONE/REAL(NTHETA)       ! step size in cos theta
        COSTHP=ONE-TWO/THREE*DCOSTH   ! cos theta plus (see paper)
        COSTHM=ONE-FOUR/THREE*DCOSTH  ! cos theta minus (see paper)
        DPHI1=FOUR*PIBY2/NAZMTH       ! step for rotation of track
        
        PHI1=ZERO
        IDIR=0
        DO IAZ=1,NAZMTH
          DO K=1,NTHETA ! Centers of upward pointing triangles
            DPHI2=TWO*PIBY3/REAL(K) ! step of trajectory azimuth 
            PHI2=(ONE/REAL(K)-ONE)*PIBY3 ! first azimuth value
            THETA=ACOS(COSTHP)
            DO J=1,K
              IDIR=IDIR+1
          IF(DEBUG(1)) PRINT*, 'BLINIT vor CALL EULERM'
              CALL EULERM(REULER,THETA,PHI1,PHI2)
          IF(DEBUG(1)) PRINT*, 'BLINIT nach CALL EULERM',REULER
              DO IR=1,3
                 DO IS=1,3
                   ADJNTB(IR,IS,IDIR)=ZERO
                   DO L=1,3
                     ADJNTB(IR,IS,IDIR)=ADJNTB(IR,IS,IDIR)
     &                                  +REULER(IR,L)*ADJBAS(L,IS)
                   END DO ! IS=1,3
                 END DO ! IR=1,3
              END DO
              PHI2=PHI2+DPHI2  ! increment trajectory azimuth
            END DO ! J=1,K
            COSTHP=COSTHP+DCOSTH  ! increment cosine of polar angle
          END DO ! K=1,NTHETA
          
        IF(DEBUG(1)) PRINT*, 'BLINIT 3.2 IAZ', IAZ
          DO K=1,NTHETA-1 ! Centers of downward pointing triangles
            DPHI2=TWO*PIBY3/REAL(K) ! step of trajectory azimuth 
            PHI2=(ONE/REAL(K)-ONE)*PIBY3 ! first azimuth value
            THETA=ACOS(COSTHM)
            DO J=1,K
              IDIR=IDIR+1
              CALL EULERM(REULER,THETA,PHI1,PHI2)
              DO IR=1,3
                 DO IS=1,3
                   ADJNTB(IR,IS,IDIR)=ZERO
                   DO L=1,3
                     ADJNTB(IR,IS,IDIR)=ADJNTB(IR,IS,IDIR)
     &                                 +REULER(IR,L)*ADJBAS(L,IS)
                   END DO ! IS=1,3
                 END DO ! IR=1,3
              END DO
              PHI2=PHI2+DPHI2  ! increment trajectory azimuth
            END DO ! J=1,K
            COSTHM=COSTHM+DCOSTH  ! increment cosine of polar angle
          END DO ! K=1,NTHETA
          PHI1=PHI1+DPHI1 ! increment azimuth for track rotation
        END DO ! IAZ=1,NAZMTH       
      END IF          ! Dummy block: Adjust lattice orientation
      
      END SUBROUTINE BLINIT
C!______________________________________________________________________

C!**********************************************************************
      SUBROUTINE EULERM(REULER,THETA,PHI1,PHI2)
C!    ***************************************************************
C!    Calculates combination of three inverse Euler rotation matrices      
C!    ***************************************************************
C!    Declarations ---------------------------------------------------
      IMPLICIT NONE
C!    ----- Scalars and arrays passed from the calling module
      REAL*8 PHI1,PHI2,THETA
      REAL*8 REULER(3,3)
C!    ----- Local scalars
      INTEGER*4 I,J
      REAL*8 RCOS, RSIN
C!    ----- Local arrays
      REAL*8 ROT1(3,3), ROT2(3,3), ROT3(3,3), RAUXIL(3,3)
C!<<<<End declarations<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      DO I=1,3
        DO J=1,3
          IF(I.EQ.J) THEN
            REULER(I,J)=1.0
            ROT1(I,J)=1.0
            ROT2(I,J)=1.0
            ROT3(I,J)=1.0
          ELSE
            REULER(I,J)=0.0
            ROT1(I,J)=0.0
            ROT2(I,J)=0.0
            ROT3(I,J)=0.0
          ENDIF
          RAUXIL(I,J)=0.0
        END DO
      END DO
      
      RCOS=COS(PHI2)
      RSIN=SIN(PHI2)
      ROT1(1,1)=RCOS
      ROT1(2,2)=RCOS
      ROT1(1,2)=RSIN
      ROT1(2,1)=-RSIN
      
      RCOS=COS(THETA)
      RSIN=SIN(THETA)
      ROT2(1,1)=RCOS
      ROT2(3,3)=RCOS
      ROT2(1,3)=RSIN
      ROT2(3,1)=-RSIN
      
      RCOS=COS(PHI1)
      RSIN=SIN(PHI1)
      ROT3(1,1)=RCOS
      ROT3(2,2)=RCOS
      ROT3(1,2)=RSIN
      ROT3(2,1)=-RSIN
      
      CALL MMULT(ROT2,ROT1,RAUXIL)
      CALL MMULT(ROT3,RAUXIL,REULER)
      END SUBROUTINE EULERM
C!______________________________________________________________________      
      
C!**********************************************************************
      SUBROUTINE MMULT(A,B,C)
C!    ***************************************************************
C!    Matrix multiplication C=A*B                                   
C!    ***************************************************************
      REAL*8 A(3,3), B(3,3), C(3,3)
      INTEGER*4 I, J, K
      DO I=1,3
        DO K=1,3
          C(I,K)=0.0
          DO J=1,3
            C(I,K)=C(I,K)+A(I,J)*B(J,K)
          END DO
        END DO
      END DO
      END SUBROUTINE MMULT
C!______________________________________________________________________

C!**********************************************************************
      SUBROUTINE MVMULT(A,B,C)
C!    ***************************************************************
C!    Multiplication of matrix A and vector B: C=A*B         
C!    ***************************************************************
      REAL*8 A(3,3), B(3), C(3)
      INTEGER*4 I, J
      DO I=1,3
        C(I)=0.0
        DO J=1,3
          C(I)=C(I)+A(I,J)*B(J)
        END DO
      END DO
      END SUBROUTINE MVMULT
C!______________________________________________________________________
