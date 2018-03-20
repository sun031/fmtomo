!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This is a simple program designed to generate a list
! of receiver locations on a regular grid
! in spherical coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM rgen
IMPLICIT NONE
INTEGER :: i,j,nr,nrth,nrph
REAL :: thn,ths,phw,phe
REAL :: th,ph,r
CHARACTER (LEN=20) :: ofile
!
! nr = Number of receivers
! nrth,nrph = Number of receivers in latitude and longitude
! thn,ths = North and South limits of grid
! phw,phe = West and East limits of grid
! th,ph,r = Latitude, longitude and depth of receiver
! ofile = Output file
!
OPEN(UNIT=10,FILE="arraygen.in",STATUS='old')
READ(10,'(a20)')ofile
READ(10,*)thn,ths
READ(10,*)phw,phe
READ(10,*)r
READ(10,*)nrth,nrph
CLOSE(10)
nr=nrth*nrph
OPEN(UNIT=20,FILE=ofile,STATUS='unknown')
WRITE(20,*)nr
DO i=1,nrth
   IF(nrth.EQ.1)THEN
      th=(thn+ths)/2.0
   ELSE
      th=ths+(i-1)*(thn-ths)/REAL(nrth-1)
   ENDIF
   DO j=1,nrph
      IF(nrph.EQ.1)THEN
         ph=(phw+phe)/2.0
      ELSE
         ph=phw+(j-1)*(phe-phw)/REAL(nrph-1)
      ENDIF
      WRITE(20,*)th,ph,r
   ENDDO
ENDDO
CLOSE(20)
STOP
END PROGRAM rgen
