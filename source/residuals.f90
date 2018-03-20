!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program calculates the RMS traveltime residual and
! variance of the current model. The 3-D fast marching program
! fm3d is used to predict traveltimes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM residual
INTEGER :: i,id1,id2,id3,id4,nt,isum,telp,ne,ntels,iswtel,id2o
INTEGER :: rmftr
INTEGER, DIMENSION(:), ALLOCATABLE :: tsid
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
REAL(KIND=i10) :: eot,trms,tvar,rd1,rd2,rd3,chis
REAL(KIND=i10), DIMENSION(100) :: paths,patht
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: mtmean
CHARACTER (LEN=30) :: mtimes,otimes,scfile,rtfile,invfile
!
! trms = RMS traveltime residual
! tvar = Traveltime residual variance
! chis = Chi squared value
! eot = Error in observed traveltime
! nt = Number of traveltimes
! telp = Teleseisms present (0=no, 1=yes)
! ns = Number of sources
! tsid = teleseismic source (0=no,1=yes)
! paths = path segment information
! path type information
! scfile = Source file
! rtfile = Reference time file for telseisms
! invfile = File containing inversion parameters
! iswtel = Switch for teleseisms
! mtmean = Model traveltime residual mean
! rmftr = Remove mean from traveltime residual (0=no, 1=yes)
!
! Read in input parameters
!
OPEN(UNIT=10,FILE='residuals.in',STATUS='old')
READ(10,'(a27)')otimes
READ(10,'(a27)')mtimes
READ(10,'(a26)')scfile
READ(10,'(a26)')rtfile
READ(10,'(a26)')invfile
CLOSE(10)
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
!
! Determine if means are to be removed from model residuals
!
IF(telp.EQ.1)THEN
   OPEN(UNIT=10,FILE=invfile,STATUS='old')
   DO i=1,21
      READ(10,*)
   ENDDO
   READ(10,*)rmftr
   CLOSE(10)
ENDIF
!
! Compute model mean for all teleseismic sources
!
IF(telp.EQ.1)THEN
   OPEN(UNIT=10,FILE=rtfile,STATUS='old')
   OPEN(UNIT=20,FILE=mtimes,STATUS='old')
   OPEN(UNIT=40,FILE=otimes,STATUS='old')
   READ(40,*)nt
   ALLOCATE(mtmean(nt))
   mtmean=0
   isum=0
   DO i=1,nt
      READ(10,*)id1,id2,id3,id4,rd3
      READ(20,*)id1,id2,id3,id4,rd2 
      IF(rd2.GE.0.0)THEN
         isum=isum+1
         rd2=rd2-rd3
         mtmean(id2)=mtmean(id2)+rd2
         IF(i.GT.1)THEN
            IF(id2o.NE.id2)THEN
               mtmean(id2o)=mtmean(id2o)/REAL(isum-1)
               isum=1
            ELSE IF(i.EQ.nt)THEN
               mtmean(id2o)=mtmean(id2o)/REAL(isum)
            ENDIF
         ENDIF
      ENDIF
      id2o=id2
   ENDDO
   CLOSE(10)
   CLOSE(20)
   CLOSE(40)
ENDIF
!
! Read in observed and model traveltimes
!
OPEN(UNIT=10,FILE=otimes,STATUS='old')
OPEN(UNIT=20,FILE=mtimes,STATUS='old')
IF(telp.EQ.1)THEN
   OPEN(UNIT=40,FILE=rtfile,STATUS='old')
ENDIF
READ(10,*)nt
tvar=0.0
isum=0
iswtel=0
chis=0.0
DO i=1,nt
   READ(10,*)id1,id2,id3,id4,rd1,eot
   READ(20,*)id1,id2,id3,id4,rd2
   IF(telp.EQ.1)THEN
      READ(40,*)id1,id2,id3,id4,rd3
      IF(tsid(id2).EQ.1)THEN
         IF(rd2.GT.0.0)iswtel=1
         IF(rmftr.EQ.1)THEN
            rd2=rd2-rd3-mtmean(id2)
         ELSE
            rd2=rd2-rd3
         ENDIF
      ENDIF
   ENDIF
   IF(rd1.LT.-50.0)iswtel=0
   IF(rd1.GT.0.0.AND.rd2.GT.0.0.OR.iswtel.EQ.1)THEN
      isum=isum+1
      tvar=tvar+(rd1-rd2)**2
      chis=chis+((rd1-rd2)/eot)**2
   ENDIF
   iswtel=0
ENDDO
CLOSE(10)
CLOSE(20)
!
! Compute variance and RMS residual
!
nt=isum
trms=1000.0*SQRT(tvar/REAL(nt))
tvar=tvar/REAL(nt-1)
chis=chis/REAL(nt)
WRITE(6,'(f10.2,2X,f10.5,2X,f10.5)')trms,tvar,chis
IF(telp.EQ.1)THEN
   CLOSE(40)
   DEALLOCATE(mtmean)
ENDIF
DEALLOCATE(tsid)
END program residual
