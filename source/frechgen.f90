!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program will generate the frechet.in file
! required by fm3d; this file determines which
! partial derivatives are computed. The choice of
! which Frechet derivative is computed is made in
! invert3d.in and frechgen.in.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM frechgen
IMPLICIT NONE
INTEGER :: i,j,idum,cfd,cvd,cid,csd,dvt,nvg,ns
INTEGER, DIMENSION (50) :: gridi,gridv,vtype
INTEGER, DIMENSION (:), ALLOCATABLE :: sid
CHARACTER (LEN=25) :: invfile,frfile,vgfile,srcfile,srcdfile
!
! invfile = Input file for program invert3d
! frfile = Input file for program fm3d specifying frechet derivatives
! vgfile = File containing reference velocity grid
! srcfile = File containing reference source locations
! srcdfile = File specifying sources for inversion
! cfd = Compute frechet derivatives (0=no,1=yes)
! cvd = Compute velocity derivatives (0=no, >0=subset, -1=all)
! cid = Compute interface derivatives (0=no, >0=subset, -1=all)
! csd = Compute source derivatives (0=no, >0=subset, -1=all)
! dvt = Velocity derivative type (1=first, 2=second, 3=both)
! gridi = Ids of interface grids for inversion
! gridv = Ids of velocity grids for inversion
! vtype = Velocity type (1 or 2) for inversion
! nvg = Number of velocity grids (or layers)
! ns = Number of sources
! sid = Source id for partial derivative
!
! Since we want to compute Frechet derivatives, we set cfd=1
!
cfd=1
!
! Read in all necessary information
!
OPEN(UNIT=10,FILE='frechgen.in',STATUS='old')
READ(10,*)invfile
READ(10,*)frfile
READ(10,*)vgfile
READ(10,*)srcfile
OPEN(UNIT=20,FILE=invfile,STATUS='old')
DO i=1,22
   READ(20,*)
ENDDO
READ(20,*)cvd
READ(20,*)cid
READ(20,*)csd
CLOSE(20)
READ(10,*)idum
IF(cvd.NE.0)cvd=idum
READ(10,*)dvt
IF(cvd.GT.0)THEN
   READ(10,*)gridv(1:cvd)
   READ(10,*)vtype(1:cvd)
ELSE
   READ(10,*)
   READ(10,*)
ENDIF
READ(10,*)idum
IF(cid.NE.0)cid=idum
IF(cid.GT.0)THEN
   READ(10,*)gridi(1:cid)
ELSE
   READ(10,*)
ENDIF
READ(10,*)idum
IF(csd.NE.0)csd=idum
READ(10,*)srcdfile
CLOSE(10)
!
! Now create the frechet derivative file.
!
OPEN(UNIT=30,FILE=frfile,STATUS='unknown')
WRITE(30,*)cfd
IF(cvd.EQ.-1.OR.cid.EQ.-1)THEN
   OPEN(UNIT=40,FILE=vgfile,STATUS='old')
   READ(40,*)nvg
   CLOSE(40)
ENDIF
!
! Start with velocities
!
IF(cvd.EQ.0)THEN
   WRITE(30,*)cvd
ELSE IF(cvd.GT.0)THEN
   WRITE(30,*)cvd
   WRITE(30,*)gridv(1:cvd)
   WRITE(30,*)vtype(1:cvd)
ELSE
   IF(dvt.EQ.3)THEN
      WRITE(30,*)2*nvg
      DO i=1,nvg
         DO j=1,2
            gridv(2*(i-1)+j)=i
            vtype(2*(i-1)+j)=j
         ENDDO
      ENDDO
      WRITE(30,*)gridv(1:2*nvg)
      WRITE(30,*)vtype(1:2*nvg)
   ELSE
      WRITE(30,*)nvg
      DO i=1,nvg
         gridv(i)=i
         vtype(i)=dvt
      ENDDO
      WRITE(30,*)gridv(1:nvg)
      WRITE(30,*)vtype(1:nvg)
   ENDIF
ENDIF
!
! Now deal with interfaces
!
IF(cid.EQ.0)THEN
   WRITE(30,*)cid
ELSE IF(cid.GT.0)THEN
   WRITE(30,*)cid
   WRITE(30,*)gridi(1:cid)
ELSE
   WRITE(30,*)nvg-1
   DO i=2,nvg
      gridi(i)=i
   ENDDO
   WRITE(30,*)gridi(2:nvg)
ENDIF
!
! Finally, deal with sources
!
IF(csd.EQ.0)THEN
   WRITE(30,*)csd
ELSE
   OPEN(UNIT=50,FILE=srcfile,STATUS='old')
   READ(50,*)ns
   CLOSE(50)
   ALLOCATE(sid(ns))
   IF(csd.GT.0)THEN
      OPEN(UNIT=60,FILE=srcdfile,STATUS='old')
      READ(60,*)ns
      DO i=1,ns
         READ(60,*)sid(i)
      ENDDO
      CLOSE(60)
      WRITE(30,*)ns
      WRITE(30,*)sid(1:ns)
      CLOSE(60)
   ELSE
      WRITE(30,*)ns
      DO i=1,ns
         sid(i)=i
      ENDDO
      WRITE(30,*)sid(1:ns)
   ENDIF
   DEALLOCATE(sid)
ENDIF
CLOSE(30)
END PROGRAM frechgen
