!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program will generate source and receiver input
! files in the correct format for the 3D fast marching
! code fm3d. Input source and receiver files in a
! much simpler format are required. Note that this
! conversion program does not exploit the full
! flexibility of path combinations allowed by fm3d.
! Frechet derivatives will also be calculated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM moddata
IMPLICIT NONE
INTEGER :: i,j,k,nisf,nps,ns,nsif
INTEGER :: lot,npr,nr,istep,checkstat,mxnps
INTEGER :: nvg,csdif,nsi
INTEGER, DIMENSION (50) :: npsg
INTEGER, DIMENSION (:), ALLOCATABLE :: numps,soep,npsp,sid
INTEGER, DIMENSION (100,50) :: paths, velt
REAL :: depth,lat,long,covdep,covlat,covlon
CHARACTER (LEN=10) :: phase
CHARACTER (LEN=30) :: ifiles,ifiler,ofiles,ofiler,ofilesd
CHARACTER (LEN=30) :: ifilevg
!
! nisf = Number of input source files
! ns = Number of sources
! nsif = Number of sources in input file
! nr = Number of receivers
! npr = Number of paths to receiver
! lot = Local (0) or teleseismic (1) source
! nps = Number of paths from source
! mxnps = Maximum nps value
! npsg = Number of path segments
! ifiles = Input source file
! ifiler = Input receiver file
! ofiles = Output source file
! ofiler = Output receiver file
! ofilesd = Output source derivative file
! npsg = Number of path segments
! paths = Path sequence
! phase = Teleseismic path
! velt = Velocity type associated with segment
! depth,lat,long = Location of source or receiver
! covdep,covlat,covlon = Error associated with source location
! numps = Number of paths from each source
! soep = Source of each path to receiver
! npsp = Number of path in source path file
! nvg = Number of velocity grids
! csdif = Compute source derivative for this input file (0=no,1=yes)
! ifilevg = Input velocity grid file
! sid = Source id for partial derivative
! nsi = Number of sources for inversion
!
OPEN(UNIT=10, FILE='moddata.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ofiles
READ(10,1)ofiler
READ(10,1)ofilesd
READ(10,1)ifiler
READ(10,*)nisf
1 FORMAT(a26)
!
! Create the source file. The total number
! of path types for each source will be recorded during
! this process. The first step is to determine the total number
! of sources
!
ns=0
DO i=1,nisf
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,1)ifiles
   READ(10,*)
   READ(10,*)
   READ(10,*)nps
   DO j=1,nps
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)
   ENDDO
   OPEN(UNIT=20,FILE=ifiles,STATUS='old')
   READ(20,*)nsif
   ns=ns+nsif
   CLOSE(20)
ENDDO
CLOSE(10)
ALLOCATE(numps(ns), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM moddata: REAL numps'
ENDIF
ALLOCATE(sid(ns), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM moddata: REAL sid'
ENDIF
OPEN(UNIT=10, FILE='moddata.in',STATUS='old')
DO i=1,10
   READ(10,*)
ENDDO
OPEN(UNIT=20,FILE=ofiles,STATUS='unknown')
WRITE(20,*)ns
!
! Now create the source file
!
istep=0
nsi=0
mxnps=0
DO i=1,nisf
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,1)ifiles
   READ(10,*)lot
   READ(10,*)csdif
   READ(10,*)nps
   IF(mxnps.LT.nps)mxnps=nps
!
!  Loop through all paths and segments and record
!  their values
!
   DO j=1,nps
      READ(10,*)npsg(j)
      READ(10,*)
      READ(10,*)paths(1:2*npsg(j),j)
      READ(10,*)velt(1:npsg(j),j)
   ENDDO
!
!  Now read in each source and write path information
!  to file.
!
   OPEN(UNIT=30,FILE=ifiles,STATUS='old')
   READ(30,*)nsif
   DO j=1,nsif
      IF(csdif.EQ.1)THEN
         READ(30,*)lat,long,depth,covlat,covlon,covdep
      ELSE
         READ(30,*)lat,long,depth
      ENDIF
      IF(lot.EQ.1)READ(30,2)phase
2     FORMAT(a10)
      WRITE(20,*)lot
      IF(lot.EQ.1)WRITE(20,2)phase
      IF(csdif.EQ.1)THEN
         WRITE(20,3)depth,lat,long,covdep,covlat,covlon
      ELSE
         WRITE(20,*)depth,lat,long
      ENDIF
3     FORMAT(6f12.5)
      WRITE(20,*)nps
      DO k=1,nps
         WRITE(20,*)npsg(k)
         WRITE(20,*)paths(1:2*npsg(k),k)
         WRITE(20,*)velt(1:npsg(k),k)
      ENDDO
!
!     Record the number of paths from this source
!
      istep=istep+1
      numps(istep)=nps
      IF(csdif.EQ.1)THEN
         nsi=nsi+1
         sid(nsi)=istep
      ENDIF
   ENDDO
   CLOSE(30)
ENDDO
CLOSE(20)
CLOSE(10)
!
! Write source Frechet derivative file
!
IF(nsi.GT.0)THEN
   OPEN(UNIT=40,FILE=ofilesd,STATUS='unknown')
   WRITE(40,*)nsi
   WRITE(40,*)sid(1:nsi)
   CLOSE(40)
ENDIF
!
! Create arrays with path information
!
ALLOCATE(soep(ns*mxnps), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM moddata: REAL soep'
ENDIF
ALLOCATE(npsp(ns*mxnps), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM moddata: REAL npsp'
ENDIF
npr=0
DO i=1,ns
   DO j=1,numps(i)
      npr=npr+1
      soep(npr)=i
      npsp(npr)=j
   ENDDO
ENDDO
!
! Now create receiver file
!
OPEN(UNIT=10,FILE=ifiler,STATUS='old')
READ(10,*)nr
OPEN(UNIT=20,FILE=ofiler,STATUS='unknown')
WRITE(20,*)nr
DO i=1,nr
   READ(10,*)lat,long,depth
   WRITE(20,*)depth,lat,long
   WRITE(20,*)npr
   WRITE(20,*)soep(1:npr)
   WRITE(20,*)npsp(1:npr)
ENDDO
CLOSE(10)
CLOSE(20)
DEALLOCATE(numps,soep,npsp,sid, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM moddata: Final deallocate'
ENDIF
STOP
END PROGRAM moddata
