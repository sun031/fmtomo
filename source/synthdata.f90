!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program simply adds a synthetic picking error to
! a set of "observed" traveltimes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM stimes
IMPLICIT NONE
INTEGER i,j,i1,i2,i3,i4,nr,telp,ns,ntels,ios,agn,rseed,wwsd
INTEGER telegn
INTEGER, DIMENSION(:), ALLOCATABLE :: tsid
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER, PARAMETER :: maxiter=100000000
REAL(KIND=i10) :: tt,cv,rd1,ttpert,sdgn,sdgnl,sdgnt,vgn
REAL(KIND=i10), DIMENSION(100) :: paths,patht
REAL, EXTERNAL :: gasdev
CHARACTER(LEN=30) :: ofile,ifile,scfile,rtfile
!
! ifile = Input "observed" time file
! ofile = Output "observed" time file
! scfile = Source file
! rtfile = Reference time file for telseisms
! nr = Number of rays
! tt = Traveltime
! ttpert = Perturbation added to traveltimes
! cv = Data covariance
! telp = Teleseisms present (0=no, 1=yes)
! ns = Number of sources
! tsid = teleseismic source (0=no,1=yes)
! paths = path segment information
! path type information
! gasdev = function which returns Gaussian noise
! agn = Add Gaussian noise (0=no,1=yes)
! sdgn = Standard deviation of Gaussian noise
! sdgnl = Standard deviation of Gaussian noise (local)
! sdgnt = Standard deviation of Gaussian noise (teleseismic)
! rseed = Random seed for noise generation
! wwsd = Weight with standard deviation (1) or noise value (2)
! vgn = Value of Gaussian noise added
!
OPEN(UNIT=10,FILE='synthdata.in',STATUS='old')
READ(10,'(a26)')ifile
READ(10,*)cv
READ(10,'(a26)')ofile
READ(10,'(a26)')scfile
READ(10,'(a26)')rtfile
READ(10,*)ttpert
READ(10,*)agn
READ(10,*)sdgnl
READ(10,*)sdgnt
READ(10,*)wwsd
READ(10,*)rseed
CLOSE(10)
!
! Determine number of receivers
!
OPEN(UNIT=10,FILE=ifile,STATUS='old')
nr=0
DO i=1,maxiter
   READ(10,*,IOSTAT=ios)i1,i2,i3,i4,tt
   IF(ios.LT.0)EXIT
   nr=nr+1
ENDDO
CLOSE(10)
OPEN(UNIT=10,FILE=ifile,STATUS='old')
OPEN(UNIT=20,FILE=ofile,STATUS='unknown')
!
! Determine if teleseisms are present
!
telp=0
tsid=0
OPEN(UNIT=30,FILE=scfile,STATUS='old')
READ(30,*)ns
ALLOCATE(tsid(ns))
DO i=1,ns
   READ(30,*)i1
   IF(i1.EQ.1)THEN
      telp=1
      tsid(i)=1
      READ(30,*)
   ENDIF
   READ(30,*)
   READ(30,*)i2
   DO j=1,i2
      READ(30,*)i3
      READ(30,*)paths(1:2*i3)
      READ(30,*)patht(1:i3)
   ENDDO
ENDDO
CLOSE(30)
OPEN(UNIT=10,FILE=ifile,STATUS='old')
OPEN(UNIT=20,FILE=ofile,STATUS='unknown')
IF(telp.EQ.1)THEN
   OPEN(UNIT=40,FILE=rtfile,STATUS='old')
ENDIF
WRITE(20,*)nr
DO i=1,nr
   sdgn=sdgnl
   READ(10,*)i1,i2,i3,i4,tt
   IF(tt.LT.0.0)tt=-100.0
   IF(telp.EQ.1)THEN
      READ(40,*)i1,i2,i3,i4,rd1
      IF(tsid(i2).EQ.1)THEN
         tt=tt-rd1
         sdgn=sdgnt
      ENDIF
   ENDIF
   tt=tt+ttpert
   IF(agn.EQ.1)THEN
      vgn=gasdev(rseed)*sdgn
      tt=tt+vgn
      IF(wwsd.EQ.2)THEN
         cv=ABS(vgn)
         IF(ABS(cv).LT.0.05*sdgn)cv=0.05*sdgn
      ELSE
         cv=sdgn
         IF(telegn.EQ.1)cv=sdgn
      ENDIF
   ENDIF
   WRITE(20,*)i1,i2,i3,i4,tt,cv
ENDDO
CLOSE(10)
CLOSE(20)
IF(telp.EQ.1)THEN
   CLOSE(40)
ENDIF
DEALLOCATE(tsid)
END PROGRAM stimes

REAL FUNCTION gasdev(idum)
IMPLICIT NONE
INTEGER :: iset,i
INTEGER, PARAMETER :: imax=100000
INTEGER :: idum
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: fac,rsq,v1,v2
REAL(KIND=i10), SAVE :: gset
REAL, EXTERNAL :: ran1
iset=0
IF(iset.EQ.0)then
   DO i=1,imax
      v1=2.0*ran1(idum)-1.
      v2=2.0*ran1(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.LT.1.AND.rsq.NE.0.)EXIT
   ENDDO
   fac=sqrt(-2.0*LOG(rsq)/rsq)
   gset=v1*fac
   gasdev=v2*fac
   iset=1
ELSE
   gasdev=gset
   iset=0
ENDIF
END FUNCTION gasdev

REAL FUNCTION ran1(idum)
IMPLICIT NONE
INTEGER :: idum
INTEGER, PARAMETER :: ia=16807,im=2147483647,iq=127773
INTEGER, PARAMETER :: ir=2836,ntab=32,ndiv=1+(im-1)/ntab
INTEGER :: j,k
INTEGER, SAVE :: iy
INTEGER, DIMENSION (:), SAVE ::iv(ntab)
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10), PARAMETER :: eps=1.2e-7,rnmx=1.0-eps,am=1./im
iv=ntab*0
iy=0
IF(idum.LE.0.OR.iy.EQ.0)THEN
   DO j=ntab+8,1,-1
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      IF(idum.LT.0)idum=idum+im
      IF(j.LE.ntab)iv(j)=idum
   ENDDO
   iy=iv(1)
ENDIF
k=idum/iq
idum=ia*(idum-k*iq)-ir*k
IF(idum.LT.0)idum=idum+im
j=1+iy/ndiv
iy=iv(j)
iv(j)=idum
ran1=min(am*iy,rnmx)
END FUNCTION ran1
