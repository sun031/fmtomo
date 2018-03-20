!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program will generate source and receiver input files, 
! and optionally, arrival times in the correct format for the 
! 3D fast marching code fm3d. Input source and receiver files 
! in a much simpler format are required. Note that this
! conversion program does not exploit the full
! flexibility of path combinations allowed by fm3d.
! Frechet derivatives will also be calculated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM obsdata
IMPLICIT NONE
INTEGER :: i,j,jj,k,nisf,nps,ns,nsif,nrif,extt
INTEGER :: lot,npr,nr,istep,checkstat,mxnps
INTEGER :: nvg,csdif,nsi
INTEGER :: idm1,idm2,siter,riter,npath,nval,npathr,isw
INTEGER, DIMENSION (50) :: npsg,pathids,pathidr
INTEGER, DIMENSION (:), ALLOCATABLE :: sid
INTEGER, DIMENSION (100,50) :: paths, velt
REAL :: depth,lat,long,covdep,covlat,covlon
REAL :: rlat,rlon,rdep,tt,tterr
CHARACTER (LEN=10) :: phase
CHARACTER (LEN=30) :: ifiles,ifiler,ofiles,ofiler,ofilesd
CHARACTER (LEN=30) :: ifilevg,ifilerd,ofilet,cdum
!
! nisf = Number of input source files
! ns = Number of sources
! nr = Number of receivers (=number of traveltimes)
! nsif = Number of sources in input file
! nrif = Number of receivers in input file
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
! ofilet = Output traveltime file
! ifilerd = Input receiver file directory
! extt = Extract traveltimes? (0=no, 1=yes)
! npsg = Number of path segments
! paths = Path sequence
! phase = Teleseismic path
! velt = Velocity type associated with segment
! depth,lat,long = Location of source or receiver
! covdep,covlat,covlon = Error associated with source location
! nvg = Number of velocity grids
! csdif = Compute source derivative for this input file (0=no,1=yes)
! ifilevg = Input velocity grid file
! sid = Source id for partial derivative
! nsi = Number of sources for inversion
! siter = Source iteration number
! riter = Receiver iteration number
! rlat,rlon,rdep = Coordinate of receiver
! tt,tterr = Traveltime and traveltime error
! npath = Number of paths for source
! npathr = Number of paths for receiver
! nval = null (0) value for arrival time file
! pathidr = path id for receiver output file
! pathids = path id for source input file
!
OPEN(UNIT=10, FILE='obsdata.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ofiles
READ(10,1)ofiler
READ(10,1)ofilesd
READ(10,1)ifilerd
READ(10,*)nisf
READ(10,*)extt
READ(10,1)ofilet
1 FORMAT(a26)
!
! Determine the number of sources and receivers. Note that
! the number of receivers is equal to the number of
! traveltimes. 
!
ns=0
nr=0
DO i=1,nisf
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,1)ifiles
   READ(10,*)lot
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
   IF(lot.EQ.1)THEN
      DO j=1,nsif
         READ(20,*)idm1,idm2,cdum
         WRITE(ifiler,*)TRIM(ifilerd),TRIM(cdum)
         ifiler=ADJUSTL(ifiler)
         OPEN(UNIT=30,FILE=ifiler,STATUS='old')
         READ(30,*)nrif
         nr=nr+nrif
         CLOSE(30)
      ENDDO
   ELSE
      DO j=1,nsif
         READ(20,*)
         IF(lot.EQ.1)READ(20,*)
         READ(20,*)npath
         DO k=1,npath
            READ(20,*)idm1,idm2,cdum
            WRITE(ifiler,*)TRIM(ifilerd),TRIM(cdum)
            ifiler=ADJUSTL(ifiler)
            OPEN(UNIT=30,FILE=ifiler,STATUS='old')
            READ(30,*)nrif
            nr=nr+nrif
            CLOSE(30)
         ENDDO
      ENDDO
      CLOSE(20)
   ENDIF
ENDDO
CLOSE(10)
ALLOCATE(sid(ns), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM obsdata: REAL sid'
ENDIF
OPEN(UNIT=10, FILE='obsdata.in',STATUS='old')
DO i=1,12
   READ(10,*)
ENDDO
!
! Now create the source file, the receiver file, and if required,
! the traveltime file.
!
OPEN(UNIT=20,FILE=ofiles,STATUS='unknown')
WRITE(20,*)ns
OPEN(UNIT=60,FILE=ofiler,STATUS='unknown')
WRITE(60,*)nr
IF(extt.EQ.1)THEN
   OPEN(UNIT=70,FILE=ofilet,STATUS='unknown')
   WRITE(70,*)nr
ENDIF
istep=0
siter=0
riter=0
nval=0
nsi=0
mxnps=0
npathr=1
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
      siter=siter+1
      IF(lot.EQ.0)THEN
         IF(csdif.EQ.1)THEN
            READ(30,*)lat,long,depth,covlat,covlon,covdep
         ELSE
            READ(30,*)lat,long,depth
         ENDIF
         READ(30,*)npath
      ELSE
         npath=1
      ENDIF
!
!     Read in receiver and path information
!
      DO jj=1,npath
         READ(30,*)pathids(jj),pathidr(jj),cdum
         IF(pathidr(jj).GT.nps)THEN
            WRITE(6,*)'ERROR!!!!!!!!'
            WRITE(6,*)'Requested path',pathidr(jj)
            WRITE(6,*)'does not exist for source file',i
            WRITE(6,*)'TERMINATING PROGRAM!!!'
            STOP
         ENDIF
         WRITE(ifiler,*)TRIM(ifilerd),TRIM(cdum)
         ifiler=ADJUSTL(ifiler)
         OPEN(UNIT=80,FILE=ifiler,STATUS='old')
         READ(80,*)nrif
         IF(lot.EQ.1)THEN
            READ(80,*)lat,long,depth
            READ(80,*)phase
         ENDIF
         DO k=1,nrif
            riter=riter+1
            IF(extt.EQ.1)THEN
               READ(80,*)rlat,rlon,rdep,tt,tterr
            ELSE
               READ(80,*)rlat,rlon,rdep
            ENDIF
            WRITE(60,*)rdep,rlat,rlon
            WRITE(60,*)npathr
            WRITE(60,*)siter
            WRITE(60,*)pathidr(jj)
            IF(extt.EQ.1)THEN
               WRITE(70,*)riter,siter,pathidr(jj),nval,tt,tterr
            ENDIF
         ENDDO
         CLOSE(80)
      ENDDO
      WRITE(20,*)lot
      IF(lot.EQ.1)WRITE(20,2)phase
2     FORMAT(a10)
      IF(csdif.EQ.1)THEN
         WRITE(20,3)depth,lat,long,covdep,covlat,covlon
      ELSE
         WRITE(20,*)depth,lat,long
      ENDIF
3     FORMAT(6f12.5)
      idm1=0
      DO k=1,nps
         DO jj=1,npath
            IF(k.EQ.pathids(jj))THEN
               idm1=idm1+1
               EXIT
            ENDIF
         ENDDO
      ENDDO
      WRITE(20,*)idm1
      DO k=1,nps
         isw=0
         DO jj=1,npath
            IF(k.EQ.pathids(jj))isw=1
         ENDDO
         IF(isw.EQ.1)THEN
            WRITE(20,*)npsg(k)
            WRITE(20,*)paths(1:2*npsg(k),k)
            WRITE(20,*)velt(1:npsg(k),k)
         ENDIF
      ENDDO
!
!     Record the number of paths from this source
!
      istep=istep+1
      IF(csdif.EQ.1)THEN
         nsi=nsi+1
         sid(nsi)=istep
      ENDIF
   ENDDO
   CLOSE(30)
ENDDO
CLOSE(20)
CLOSE(10)
CLOSE(60)
IF(extt.EQ.1)CLOSE(70)
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
DEALLOCATE(sid, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM obsdata: Final deallocate'
ENDIF
STOP
END PROGRAM obsdata
