      SUBROUTINE RFINIT()
C!    Sets the information for the columns of the input file
C!      
C!    ------ Global variables
C!    ------
C!     CODE: Character string encoding the meaning of the entries in a 
C!           line of the input file as follows:
C!           'T'     - number of the primary particle track
C!           'X','Y','Z' - x, y, and z coordinates of the transfer point 
C!           'E'     - energy deposit (if applicable)
C!           'I'     - ionization cluster size (if applicable)
C!           'C'     - GEANT4 code for interaction type
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
C!    Local variables
      INTEGER*1 J
      
      DO J=1,8
        IDT(J)=0
      END DO
      
      DO J=1,8
        IF(CODE(J:J).EQ.'X'.AND.IDT(1).EQ.0) IDT(1)=J
        IF(CODE(J:J).EQ.'Y'.AND.IDT(2).EQ.0) IDT(2)=J
        IF(CODE(J:J).EQ.'Z'.AND.IDT(3).EQ.0) IDT(3)=J
        IF(CODE(J:J).EQ.'E'.AND.IDT(4).EQ.0) IDT(4)=J
        IF(CODE(J:J).EQ.'T'.AND.IDT(5).EQ.0) IDT(5)=J
        IF(CODE(J:J).EQ.'I'.AND.IDT(6).EQ.0) IDT(6)=J
        IF(CODE(J:J).EQ.'C'.AND.IDT(7).EQ.0) IDT(7)=J
        IF(CODE(J:J).NE.' ') NDT=J
      END DO
      
      IF(IDT(1).EQ.0) PRINT*, 'No data column for X'
      IF(IDT(2).EQ.0) PRINT*, 'No data column for Y'
      IF(IDT(3).EQ.0) PRINT*, 'No data column for Z'
      IF(IDT(5).EQ.0) PRINT*, 'No data column for track number'
      IF(IDT(1)*IDT(2)*IDT(3)*IDT(5).EQ.0) STOP
      
      END