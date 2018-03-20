!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: MODULE
! CODE: FORTRAN 90
! This module declares variable for global use, that is, for
! USE in any subroutine or function or other module.
! Variables whose values are SAVEd can have their most
! recent values reused in any routine.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE globalp
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: nint,ninp
INTEGER, DIMENSION(500) :: nnr,nnt,nnp
INTEGER, DIMENSION(:,:), ALLOCATABLE :: lay
REAL(KIND=i10) :: earthr
REAL(KIND=i10), PARAMETER :: pi=3.141592653589793
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: vela,inta
REAL(KIND=i10), DIMENSION(:,:,:), ALLOCATABLE :: intn
REAL(KIND=i10), DIMENSION(:,:,:,:), ALLOCATABLE :: veln
! vela = diced velocity values
! inta = diced interface values
! veln = velocity grid values
! intn = interface grid values
! nnr,nnt,nnp = number of diced nodes in r,theta,phi
! lay = pointer to layer in which point resides
! nint,ninp = number of interface nodes in theta,phi
! earthr = earth radius
END MODULE globalp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to take depth slices through
! the FMM computed traveltime field and put them in a
! form suitable for contour plotting by GMT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM slice
USE globalp
IMPLICIT NONE
INTEGER :: ii,i,j,k,extrds,extrnss,extrews,extrgcs
INTEGER :: checkstat,isw,parv,ni,nip,pard,ns,nr
INTEGER :: idm1,idm2,idm3,idm4,nnx,nnz,prp
INTEGER :: rnode,intid,indt,indp,ios,rnodeth,rnodeph,rnoder
INTEGER :: ddt,ddp,dnsr,dnst,dewr,dewp,nrs,nrp
INTEGER :: mnr,mnt,mnp,vgid,ppors,nvgt,ips,psr,maxpth
INTEGER :: idmnr,idmnt,idmnp,penp,ngch,ngcd
INTEGER, PARAMETER :: maxr=500000000
REAL(KIND=i10) :: lft,rgt,btm,top,sldep,slns,slew
REAL(KIND=i10) :: rd1,rd2,rd3,rd4,rd5,rd6,rd7,rd8,u,v,w
REAL(KIND=i10) :: rdep,rlat,rlong
REAL(KIND=i10) :: deltas,deltaskm,azim
REAL(KIND=i10) :: deltanp,thgcp,phgcp,goti,gopi
REAL(KIND=i10) :: slgclat1,slgclon1,slgclat2,slgclon2
REAL(KIND=i10) :: slgclat1d,slgclon1d,slgclat2d,slgclon2d
REAL(KIND=i10) :: mgsr,mgst,mgsp,rgsr,rgst,rgsp,gsit,gsip
REAL(KIND=i10), DIMENSION(500) :: gor,got,gop,gsr,gst,gsp,paths
REAL(KIND=i10), DIMENSION(:), ALLOCATABLE :: srcid
REAL(KIND=i10), PARAMETER :: dtol=1.0e-6
REAL(KIND=i10), PARAMETER :: stol=1.0e-8
CHARACTER (LEN=30) :: ifilev,ifilevr,ifilei,ifileir,sep
CHARACTER (LEN=30) :: ofiledb,ofilensb,ofileewb
CHARACTER (LEN=30) :: ofiledv,ofilensv,ofileewv,ofilecont
CHARACTER (LEN=30) :: ofilensi,ofileewi,ofileint,ofileib
CHARACTER (LEN=30) :: ifileray,ofilerd,ofilerew,ofilerns
CHARACTER (LEN=30) :: ifilesrc,ifilercv,ofilesd,ofilesew
CHARACTER (LEN=30) :: ofilesns,ofilercd,ofilercew,ofilercns
CHARACTER (LEN=30) :: ofilergc,ofilesgc,ofilercgc,ofilegcb
CHARACTER (LEN=30) :: ofilegcv,ofilegci
!
! sldep = slice depth
! slns = NS slice
! slew = EW slice
! ifilev = input velocity grid file
! ifilevr = input reference velocity grid file
! ifilei = input interface grid file
! ifileir = input reference interface grid file
! ofiledb = bounds for output depth slice file
! ofilensb = bounds for output N-S slice file
! ofileewb = bounds for output E-W slice file
! ofilegcb = bounds for output great circle file
! ofiledv = output velocity file for depth slice
! ofilensv = output velocity file for N-S slice
! ofileewv = output velocity file for E-W slice
! ofilegcv = output velocity file for great circle slice
! ofilensi = output interface file for N-S slice
! ofileewi = output interface file for E-W slice
! ofilegci = output interface file for great circle slice
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! gsr,gst,gsp = Node spacing in r,theta,phi
! lft,rgt,btm,top = plotting bounds
! ddt,ddp = dicing of velocity grid for depth slice
! dnsr,dnst = dicing of velocity grid for N-S slice
! dewr,dewp = dicing of velocity grid for E-W slice
! sep = character marker for separating interfaces
! extrds = extract depth slice? (0=no,1=yes)
! extrnss = extract N-S slice? (0=no,1=yes)
! extrews = extract E-W slice? (0=no,1=yes)
! extrgcs = extract great circle slice? (0=no, 1=yes)
! u,v,w = bspline independent parameters
! nnx,nnz = dimensions of vela
! rnode = reference node for slice
! dtol = tolerance for grid spacing difference
! parv = plot absolute (0) or relative (1) velocity
! ni = number of interfaces
! nip = number of interfaces for plotting
! intid = interface id. for surface plotting
! indt,indp = number of nodes in theta,phi for surface
! ofileint = output file for surfaces
! ofileib = output file for interface plotting bounds
! ofilecont = output interface depth file
! pard = plot absolute (0) or relative (1) depth
! mnr,mnt,mnp = maximum nnr,nnt,nnp
! idmnr,idmnt,idmnp = Layer ids for mnr,mnt,mnp
! rgsr,rgst,rgsp = refined grid spacing for dicing
! mgsr,mgst,mgsp = grid spacing for maximum nnr,nnt,nnp
! gsit,gsip = interface grid node spacing
! vgid = velocity grid identical (0) or not (1)
! ppors = plot P(0) or S(1) velocities
! nvgt = number of velocity grid types
! ips = counter for P(0) or S(1) velocities
! prp = plot ray paths (0=no, 1=yes)
! ifileray = input raypath file
! ofilerd = output ray file in depth projection
! ofilerew = output ray file in EW projection
! ofilerns = output ray file in NS projection
! ofilergc = output ray file in great circle projection
! ios = flag for end of file
! nrs = number of ray sections
! nrp = number of points in path
! rdep,rlat,rlong = Coordinates of ray path points
! maxr = arbitrary maximum number of rays
! ns = Number of sources
! nr = Number of receivers
! ifilesrc = Input source file
! ifilercv = Input receiver file
! ofilesd = Output source file in depth
! ofilesew  = Output source file in EW
! ofilesns = Output source file in NS
! ofilesgc = Output source file in great circle
! ofilercd = Output receiver file in depth
! ofilercew = Output receiver file in EW
! ofilercns = Output receiver file in NS
! ofilercgc = Output receiver file in great circle
! psr = Plot sources and receivers (0=no,1=yes)
! paths = Path segment information (ignored)
! srcid = Source ids for each receiver
! maxpth = Maximum number of paths for a source
! penp = Plot every nth ray point
! slgclat1,slgclon1 = First point for great circle slice
! slgclat2,slgclon2 = Second point for great circle slice
! slgclat1,slgclon1 = First point for great circle slice (in degrees)
! slgclat2,slgclon2 = Second point for great circle slice (in degrees)
! ngch,ngcd = Number of sample points in horizontal and depth 
!             for great circle slice.
! deltas = Great circle angular distance span of great circle slice
! deltaskm = Same as deltas, but in km
! stol = Slice tolerance for great circle slice
! azim = Azimuth of great circle slice
! deltanp = Angular distance to point along great circle
! thgcp,phgcp = Lat and long of point along great circle
! rnodeth,rnodeph,rnoder = Reference grid nodes for great circle slice
! goti,gopi = Grid origin of interface (theta,phi)
!
sep='>'
OPEN(UNIT=10,FILE='gmtslice.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a27)')ifilev
READ(10,'(a27)')ifilevr
READ(10,'(a27)')ifilei
READ(10,'(a27)')ifileir
READ(10,*)earthr
READ(10,*)ppors
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)prp
READ(10,'(a27)')ifileray
READ(10,'(a27)')ofilerd
READ(10,'(a27)')ofilerew
READ(10,'(a27)')ofilerns
READ(10,'(a27)')ofilergc
READ(10,*)penp
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)psr
READ(10,'(a27)')ifilesrc
READ(10,'(a27)')ifilercv
READ(10,'(a27)')ofilesd
READ(10,'(a27)')ofilesew
READ(10,'(a27)')ofilesns
READ(10,'(a27)')ofilesgc
READ(10,'(a27)')ofilercd
READ(10,'(a27)')ofilercew
READ(10,'(a27)')ofilercns
READ(10,'(a27)')ofilercgc
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
!
! Read in slice parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)extrds
READ(10,*)sldep
READ(10,'(a27)')ofiledb
READ(10,'(a27)')ofilecont
READ(10,*)extrnss
READ(10,*)slns
READ(10,'(a27)')ofilensb
READ(10,*)extrews
READ(10,*)slew
READ(10,'(a27)')ofileewb
READ(10,*)extrgcs
READ(10,*)slgclat1d,slgclon1d
READ(10,*)slgclat2d,slgclon2d
READ(10,'(a27)')ofilegcb
!
! Start off by reading in complete velocity grid.
!
OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)ni,nvgt
IF(nvgt.LT.1.OR.nvgt.GT.2)THEN
   WRITE(6,*)'There can only be one or two'
   WRITE(6,*)'sets of velocity grids corresponding'
   WRITE(6,*)'to P and S wavespeeds!!!'
   WRITE(6,*)'TERMINATING PROGRAM!!!'
   STOP
ENDIF
IF(ppors.EQ.0)nvgt=1
ni=ni+1
vgid=0
DO ips=1,nvgt
   IF(ips.EQ.2)THEN
      DEALLOCATE(veln, STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
      ENDIF
   ENDIF
   DO ii=1,ni-1
      READ(20,*)nnr(ii),nnt(ii),nnp(ii)
      READ(20,*)gsr(ii),gst(ii),gsp(ii)
      READ(20,*)gor(ii),got(ii),gop(ii)
      IF(ii.GT.1)THEN
         IF(nnr(ii)-2.NE.nnr(ii-1))vgid=1
         IF(nnt(ii)-2.NE.nnt(ii-1))vgid=1
         IF(nnp(ii)-2.NE.nnp(ii-1))vgid=1
      ENDIF
      nnr(ii)=nnr(ii)-2
      nnt(ii)=nnt(ii)-2
      nnp(ii)=nnp(ii)-2
      gst(ii)=gst(ii)*180.0/pi
      gsp(ii)=gsp(ii)*180.0/pi
      got(ii)=got(ii)*180.0/pi
      gop(ii)=gop(ii)*180.0/pi
      gor(ii)=gor(ii)+nnr(ii)*gsr(ii)-earthr
      got(ii)=got(ii)+gst(ii)*nnt(ii)
      gop(ii)=gop(ii)+gsp(ii)
      IF(ii.EQ.1)THEN
         ALLOCATE(veln(0:nnr(ii)+1,0:nnt(ii)+1,0:nnp(ii)+1,ni-1), STAT=checkstat)
         IF(checkstat > 0)THEN
            WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL veln'
         ENDIF
      ENDIF
      DO i=0,nnr(ii)+1
         DO j=0,nnt(ii)+1
            DO k=0,nnp(ii)+1
               IF(vgid.EQ.0)THEN
                  READ(20,*)veln(nnr(ii)+1-i,nnt(ii)+1-j,k,ii)
               ELSE
                  READ(20,*)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
CLOSE(20)
!
! If node spacing on velocity grid is not identical, then we
! need to read in the velocity grid again for memory
! allocation purposes.
!
idmnr=1
idmnt=1
idmnp=1
mnr=nnr(1)
mnt=nnt(1)
mnp=nnp(1)
mgsr=gsr(1)
mgst=gst(1)
mgsp=gsp(1)
IF(vgid.NE.0)THEN
   DEALLOCATE(veln, STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln'
   ENDIF
   DO ii=2,ni-1
      IF(nnr(ii).GT.mnr)THEN
         idmnr=ii
         mnr=nnr(ii)
         mgsr=gsr(ii)
      ENDIF
      IF(nnt(ii).GT.mnt)THEN
         idmnt=ii
         mnt=nnt(ii)
         mgst=gst(ii)
      ENDIF
      IF(nnp(ii).GT.mnp)THEN
         idmnp=ii
         mnp=nnp(ii)
         mgsp=gsp(ii)
      ENDIF
   ENDDO
   ALLOCATE(veln(0:mnr+1,0:mnt+1,0:mnp+1,ni-1), STAT=checkstat)
   IF(checkstat > 0)THEN
      WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL veln'
   ENDIF
   OPEN(UNIT=20,FILE=ifilev,status='old')
   READ(20,*)ni
   ni=ni+1
   DO ips=1,nvgt
      DO ii=1,ni-1
         READ(20,*)nnr(ii),nnt(ii),nnp(ii)
         READ(20,*)gsr(ii),gst(ii),gsp(ii)
         READ(20,*)gor(ii),got(ii),gop(ii)
         nnr(ii)=nnr(ii)-2
         nnt(ii)=nnt(ii)-2
         nnp(ii)=nnp(ii)-2
         gst(ii)=gst(ii)*180.0/pi
         gsp(ii)=gsp(ii)*180.0/pi
         got(ii)=got(ii)*180.0/pi
         gop(ii)=gop(ii)*180.0/pi
         gor(ii)=gor(ii)+nnr(ii)*gsr(ii)-earthr
         got(ii)=got(ii)+gst(ii)*nnt(ii)
         gop(ii)=gop(ii)+gsp(ii)
         DO i=0,nnr(ii)+1
            DO j=0,nnt(ii)+1
               DO k=0,nnp(ii)+1
                  READ(20,*)veln(nnr(ii)+1-i,nnt(ii)+1-j,k,ii)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Now read in complete interface grid
!
OPEN(UNIT=20,FILE=ifilei,status='old')
READ(20,*)idm1
READ(20,*)nint,ninp
nint=nint-2
ninp=ninp-2
IF(idm1.NE.ni)THEN
   WRITE(6,*)'Velocity and interface grids not'
   WRITE(6,*)'consistent.'
   WRITE(6,*)'TERMINATING PROGRAM!!'
   STOP
ENDIF
READ(20,*)gsit,gsip
gsit=gsit*180.0/pi
gsip=gsip*180.0/pi
READ(20,*)goti,gopi
goti=goti*180.0/pi
gopi=gopi*180.0/pi
ALLOCATE(intn(0:nint+1,0:ninp+1,ni), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM tslice: REAL intn'
ENDIF
DO ii=1,ni
   DO j=0,nint+1
      DO k=0,ninp+1
         READ(20,*)intn(nint+1-j,k,ii)
         intn(nint+1-j,k,ii)=intn(nint+1-j,k,ii)-earthr
      ENDDO
   ENDDO
ENDDO
CLOSE(20)
!
! Now read in velocity and interface grid parameters.
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)parv
READ(10,*)ddt,ddp
READ(10,'(a22)')ofiledv
READ(10,*)dnsr,dnst
READ(10,'(a22)')ofilensv
READ(10,*)dewr,dewp
READ(10,'(a22)')ofileewv
READ(10,*)ngch,ngcd
READ(10,'(a22)')ofilegcv
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,'(a22)')ofilensi
READ(10,'(a22)')ofileewi
READ(10,'(a22)')ofilegci
!
! Calculate GMT bounds file for depth slice if required.
!
IF(extrds.EQ.1)THEN
   lft=gop(1)
   rgt=gop(1)+(nnp(1)-1)*gsp(1)
   btm=got(1)-(nnt(1)-1)*gst(1)
   top=got(1)
   OPEN(UNIT=50,FILE=ofiledb,STATUS='unknown')
   WRITE(50,*)lft
   WRITE(50,*)rgt
   WRITE(50,*)btm
   WRITE(50,*)top
   nnx=(mnt-1)*ddt+1
   nnz=(mnp-1)*ddp+1
   WRITE(50,*)nnz
   WRITE(50,*)nnx
   CLOSE(50)
!
!  Write out interface depth file
!
   OPEN(UNIT=50,FILE=ofilecont,STATUS='unknown')
   WRITE(50,*)sldep,'  C'
   CLOSE(50)
ENDIF
!
! Calculate GMT bounds file for N-S slice if required.
!
IF(extrnss.EQ.1)THEN
   lft=got(1)-(nnt(1)-1)*gst(1)
   rgt=got(1)
   btm=gor(1)-(nnr(1)-1)*gsr(1)
   top=gor(1)
   OPEN(UNIT=60,FILE=ofilensb,STATUS='unknown')
   WRITE(60,*)lft
   WRITE(60,*)rgt
   WRITE(60,*)btm
   WRITE(60,*)top
   nnx=(mnr-1)*dnsr+1
   nnz=(mnt-1)*dnst+1
   WRITE(60,*)nnz
   WRITE(60,*)nnx
   CLOSE(60)
ENDIF
!
! Calculate GMT bounds file for E-W slice if required.
!
IF(extrews.EQ.1)THEN
   lft=gop(1)
   rgt=gop(1)+(nnp(1)-1)*gsp(1)
   btm=gor(1)-(nnr(1)-1)*gsr(1)
   top=gor(1)
   OPEN(UNIT=70,FILE=ofileewb,STATUS='unknown')
   WRITE(70,*)lft
   WRITE(70,*)rgt
   WRITE(70,*)btm
   WRITE(70,*)top
   nnx=(mnr-1)*dewr+1
   nnz=(mnp-1)*dewp+1
   WRITE(70,*)nnz
   WRITE(70,*)nnx
   CLOSE(70)
ENDIF
!
! Calculate GMT bounds file for great circle slice if required.
!
IF(extrgcs.EQ.1)THEN
   slgclat1=slgclat1d*pi/180.0
   slgclat2=slgclat2d*pi/180.0
   slgclon1=slgclon1d*pi/180.0
   slgclon2=slgclon2d*pi/180.0
   IF(slgclon1d.EQ.slgclon2d)slgclon2=slgclon2+dtol
!
!  Compute angular distance between end points of slice
!
   deltas=SIN(slgclat2)*SIN(slgclat1)
   deltas=deltas+COS(slgclat2)*COS(slgclat1)*COS(slgclon2-slgclon1)
   deltas=ACOS(deltas)
   deltaskm=deltas*earthr
   lft=0.0
   rgt=deltaskm
   btm=gor(1)-(nnr(1)-1)*gsr(1)
   top=gor(1)
   OPEN(UNIT=70,FILE=ofilegcb,STATUS='unknown')
   WRITE(70,*)lft
   WRITE(70,*)rgt
   WRITE(70,*)btm
   WRITE(70,*)top
   WRITE(70,*)ngch
   WRITE(70,*)ngcd
   CLOSE(70)
ENDIF
!
! Read in reference velocity grid file if required
!
IF(parv.EQ.1)THEN
   isw=0
   OPEN(UNIT=20,FILE=ifilevr,status='old')
   READ(20,*)idm1
   idm1=idm1+1
   IF(idm1.NE.ni)THEN
      WRITE(6,*)'Actual and reference velocity grids'
      WRITE(6,*)'have different numbers of layers!'
      WRITE(6,*)'TERMINATING PROGRAM!!!'
      STOP
   ENDIF
   DO ips=1,nvgt
      DO ii=1,ni-1
         READ(20,*)idm1,idm2,idm3
         idm1=idm1-2
         idm2=idm2-2
         idm3=idm3-2
         IF(nvgt.EQ.1.OR.nvgt.EQ.2.AND.ips.EQ.2)THEN
            IF(idm1.NE.nnr(ii).OR.idm2.NE.nnt(ii).OR.idm3.NE.nnp(ii))isw=1
         ENDIF
         READ(20,*)rd1,rd2,rd3
         rd2=rd2*180.0/pi
         rd3=rd3*180.0/pi
         IF(ABS(rd1-gsr(ii)).GT.dtol)isw=1
         IF(ABS(rd2-gst(ii)).GT.dtol)isw=1
         IF(ABS(rd3-gsp(ii)).GT.dtol)isw=1
         READ(20,*)rd1,rd2,rd3
         IF(isw.EQ.1)THEN
            WRITE(6,*)'ERROR! Actual velocity grid and reference'
            WRITE(6,*)'velocity grid have different dimensions or'
            WRITE(6,*)'different numbers of grid points!'
            WRITE(6,*)'TERMINATING PROGRAM!!!'
            STOP
         ENDIF
         DO i=0,nnr(ii)+1
            DO j=0,nnt(ii)+1
               DO k=0,nnp(ii)+1
                  READ(20,*)rd1
                  IF(nvgt.EQ.1.OR.nvgt.EQ.2.AND.ips.EQ.2)THEN
                     rd1=veln(nnr(ii)+1-i,nnt(ii)+1-j,k,ii)-rd1
                     veln(nnr(ii)+1-i,nnt(ii)+1-j,k,ii)=rd1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   CLOSE(20)
ENDIF
!
! Extract depth slice if required
!
IF(extrds.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(sldep.GT.gor(1).OR.sldep.LT.gor(1)-gsr(1)*(nnr(1)-1))THEN
     WRITE(6,*)'Requested depth slice lies outside'
     WRITE(6,*)'Model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Allocate memory to velocity grid array
!
  nnx=(mnt-1)*ddt+1
  nnz=(mnp-1)*ddp+1
  rgst=gst(1)*(nnt(1)-1)/REAL(nnx-1)
  rgsp=gsp(1)*(nnp(1)-1)/REAL(nnz-1)
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Now allocate memory to pointer array which
! indicates which layer each grid point lies in.
!
  ALLOCATE(lay(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL lay'
  ENDIF
!
! Loop through all the layers and determine which layer each
! point lies in.
!
  lay=1
  IF(ni.GT.2)THEN
     ALLOCATE(inta(nnx,nnz), STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
     DO i=2,ni-1
!
!       Call bspline subroutine to compute depth to interface i.
!
        DO j=1,nnz
           idm2=INT((j-1)*rgsp/gsip)+1
           v=(j-1)*rgsp-(idm2-1)*gsip
           v=v/gsip
           IF(idm2.EQ.ninp)THEN
              idm2=idm2-1
              v=1.0
           ENDIF
           DO k=1,nnx
              idm1=INT((k-1)*rgst/gsit)+1
              u=(k-1)*rgst-(idm1-1)*gsit
              u=u/gsit
              IF(idm1.EQ.nint)THEN
                 idm1=idm1-1
                 u=1.0
              ENDIF
              CALL ibspline(idm1,idm2,u,v,k,j,i)
           ENDDO
        ENDDO
        DO j=1,nnz
           DO k=1,nnx
              IF(sldep.LT.inta(k,j))lay(k,j)=i
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE(inta, STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
  ENDIF
!
! Now call a bspline subroutine to compute depth slice
!
  DO j=1,nnz
     DO k=1,nnx
        idm3=lay(k,j)
        rnode=INT((gor(idm3)-sldep)/gsr(idm3))+1
        u=ABS(gor(idm3)-sldep-(rnode-1)*gsr(idm3))/gsr(idm3)
        IF(rnode.EQ.nnr(idm3))THEN
           rnode=rnode-1
           u=1.0
        ENDIF
        idm1=INT((k-1)*rgst/gst(idm3))+1
        v=(k-1)*rgst-(idm1-1)*gst(idm3)
        v=v/gst(idm3)
        IF(idm1.EQ.nnt(idm3))THEN
           idm1=idm1-1
           v=1.0
        ENDIF
        idm2=INT((j-1)*rgsp/gsp(idm3))+1
        w=(j-1)*rgsp-(idm2-1)*gsp(idm3)
        w=w/gsp(idm3)
        IF(idm2.EQ.nnp(idm3))THEN
           idm2=idm2-1
           w=1.0
        ENDIF
        CALL vbspline(rnode,idm1,idm2,v,w,u,k,j,idm3)
     ENDDO
  ENDDO
!
! Now write out depth slice to file.
!
  OPEN(UNIT=30,FILE=ofiledv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
         WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,lay, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,lay,vi'
  ENDIF
ENDIF
!
! Extract N-S slice if required
!
IF(extrnss.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(slns.LT.gop(1).OR.slns.GT.gop(1)+gsp(1)*(nnp(1)-1))THEN
     WRITE(6,*)'Requested N-S slice lies outside'
     WRITE(6,*)'Model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Allocate memory to velocity grid array
!
  nnx=(mnr-1)*dnsr+1
  nnz=(mnt-1)*dnst+1
  rgsr=gsr(1)*(nnr(1)-1)/REAL(nnx-1)
  rgst=gst(1)*(nnt(1)-1)/REAL(nnz-1)
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Now allocate memory to pointer array which
! indicates which layer each grid point lies in.
!
  ALLOCATE(lay(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL lay'
  ENDIF
!
! Loop through all the layers and determine which layer each
! point lies in.
!
  lay=1
  IF(ni.GT.2)THEN
     ALLOCATE(inta(1,nnz), STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
     DO i=2,ni-1
!
!       Call bspline subroutine to compute depth to interface i.
!
        rnode=INT((slns-gop(i))/gsip)+1
        u=ABS(slns-gop(i)-(rnode-1)*gsip)/gsip
        IF(rnode.EQ.ninp)THEN
           rnode=rnode-1
           u=1.0
        ENDIF
        DO k=1,nnz
           idm1=INT((k-1)*rgst/gsit)+1
           v=(k-1)*rgst-(idm1-1)*gsit
           v=v/gsit
           idm1=INT((k-1)*rgst/gsit)+1
           IF(idm1.EQ.nint)THEN
              idm1=idm1-1
              v=1.0
           ENDIF
           CALL ibspline(idm1,rnode,v,u,1,k,i)
        ENDDO
        DO j=1,nnz
           DO k=1,nnx
              rd1=gor(idmnr)-gsr(idmnr)*(k-1)/REAL(dnsr)
              IF(rd1.LT.inta(1,j))lay(k,j)=i
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE(inta, STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
  ENDIF
!
! Now call a bspline subroutine to compute N-S slice
!
  DO j=1,nnz
     DO k=1,nnx
        idm3=lay(k,j)
        rnode=INT((slns-gop(idm3))/gsp(idm3))+1
        u=ABS(slns-gop(idm3)-(rnode-1)*gsp(idm3))/gsp(idm3)
        IF(rnode.EQ.nnp(idm3))THEN
           rnode=rnode-1
           u=1.0
        ENDIF
        idm1=INT((k-1)*rgsr/gsr(idm3))+1
        v=(k-1)*rgsr-(idm1-1)*gsr(idm3)
        v=v/gsr(idm3)
        IF(idm1.EQ.nnr(idm3))THEN
           idm1=idm1-1
           v=1.0
        ENDIF
        idm2=INT((j-1)*rgst/gst(idm3))+1
        w=(j-1)*rgst-(idm2-1)*gst(idm3)
        w=w/gst(idm3)
        IF(idm2.EQ.nnt(idm3))THEN
           idm2=idm2-1
           w=1.0
        ENDIF
        CALL vbspline(idm1,idm2,rnode,w,u,v,k,j,idm3)
     ENDDO
  ENDDO
!
! Now write out N-S slice to file.
!
  OPEN(UNIT=30,FILE=ofilensv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
        WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
!
! Compute interface cross-sections and write to file
!
  OPEN(UNIT=30,FILE=ofilensi,STATUS='unknown')
  ALLOCATE(inta(1,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
  ENDIF
  rnode=INT((slns-gop(1))/gsip)+1
  u=ABS(slns-gop(1)-(rnode-1)*gsip)/gsip
  IF(rnode.EQ.ninp)THEN
     rnode=rnode-1
     u=1.0
  ENDIF
  DO i=1,ni
     DO j=1,nnz
        idm1=INT((j-1)*rgst/gsit)+1
        v=(j-1)*rgst-(idm1-1)*gsit
        v=v/gsit
        IF(idm1.EQ.nint)THEN
           idm1=idm1-1
           v=1.0
        ENDIF
        CALL ibspline(idm1,rnode,v,u,1,j,i)
        rd1=got(1)-rgst*(j-1)
        WRITE(30,*)rd1,inta(1,j)
     ENDDO
     WRITE(30,'(a1)')sep
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,lay,inta, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,lay,inta'
  ENDIF
ENDIF
!
! Extract E-W slice if required
!
IF(extrews.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(slew.GT.got(1).OR.slew.LT.got(1)-gst(1)*(nnt(1)-1))THEN
     WRITE(6,*)'Requested E-W slice lies outside'
     WRITE(6,*)'Model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Allocate memory to velocity grid array
!
  nnx=(mnr-1)*dewr+1
  nnz=(mnp-1)*dewp+1
  rgsr=gsr(1)*(nnr(1)-1)/REAL(nnx-1)
  rgsp=gsp(1)*(nnp(1)-1)/REAL(nnz-1)
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Now allocate memory to pointer array which
! indicates which layer each grid point lies in.
!
  ALLOCATE(lay(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL lay'
  ENDIF
!
! Loop through all the layers and determine which layer each
! point lies in.
!
  lay=1
  IF(ni.GT.2)THEN
     ALLOCATE(inta(1,nnz), STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
     DO i=2,ni-1
!
!       Call bspline subroutine to compute depth to interface i.
!
        rnode=INT((got(i)-slew)/gsit)+1
        u=ABS(got(i)-slew-(rnode-1)*gsit)/gsit
        IF(rnode.EQ.nint)THEN
           rnode=rnode-1
           u=1.0
        ENDIF
        DO k=1,nnz
           idm1=INT((k-1)*rgsp/gsip)+1
           v=(k-1)*rgsp-(idm1-1)*gsip
           v=v/gsip
           IF(idm1.EQ.ninp)THEN
              idm1=idm1-1
              v=1.0
           ENDIF
           CALL ibspline(rnode,idm1,u,v,1,k,i)
        ENDDO
        DO j=1,nnz
           DO k=1,nnx
              rd1=gor(idmnr)-gsr(idmnr)*(k-1)/REAL(dewr)
              IF(rd1.LT.inta(1,j))lay(k,j)=i
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE(inta, STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
  ENDIF
!
! Now call a bspline subroutine to compute E-W slice
!
  DO j=1,nnz
     DO k=1,nnx
        idm3=lay(k,j)
        rnode=INT((got(idm3)-slew)/gst(idm3))+1
        u=ABS(got(idm3)-slew-(rnode-1)*gst(idm3))/gst(idm3)
        IF(rnode.EQ.nnt(idm3))THEN
           rnode=rnode-1
           u=1.0
        ENDIF
        idm1=INT((k-1)*rgsr/gsr(idm3))+1
        v=(k-1)*rgsr-(idm1-1)*gsr(idm3)
        v=v/gsr(idm3)
        IF(idm1.EQ.nnr(idm3))THEN
           idm1=idm1-1
           v=1.0
        ENDIF
        idm2=INT((j-1)*rgsp/gsp(idm3))+1
        w=(j-1)*rgsp-(idm2-1)*gsp(idm3)
        w=w/gsp(idm3)
        IF(idm2.EQ.nnp(idm3))THEN
           idm2=idm2-1
           w=1.0
        ENDIF
        CALL vbspline(idm1,rnode,idm2,u,w,v,k,j,idm3)
     ENDDO
  ENDDO
!
! Now write out E-W slice to file.
!
  OPEN(UNIT=30,FILE=ofileewv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
        WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
!
! Compute interface cross-sections and write to file
!
  OPEN(UNIT=30,FILE=ofileewi,STATUS='unknown')
  ALLOCATE(inta(1,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
  ENDIF
  rnode=INT((got(1)-slew)/gsit)+1
  u=ABS(got(1)-slew-(rnode-1)*gsit)/gsit
  IF(rnode.EQ.nint)THEN
     rnode=rnode-1
     u=1.0
  ENDIF
  DO i=1,ni
     DO j=1,nnz
        idm1=INT((j-1)*rgsp/gsip)+1
        v=(j-1)*rgsp-(idm1-1)*gsip
        v=v/gsip
        IF(idm1.EQ.ninp)THEN
           idm1=idm1-1
           v=1.0
        ENDIF
        CALL ibspline(rnode,idm1,u,v,1,j,i)
        rd1=gop(1)+rgsp*(j-1)
        WRITE(30,*)rd1,inta(1,j)
     ENDDO
     WRITE(30,'(a1)')sep
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,lay,inta, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,lay,inta'
  ENDIF
ENDIF
!
! Extract great circle slice if required
!
IF(extrgcs.EQ.1)THEN
!
! Make sure slice lies within model
!
  IF(slgclat1d.GT.got(1)+stol.OR.slgclat1d.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
     WRITE(6,*)'Requested latitude of point 1 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ELSE IF(slgclat2d.GT.got(1)+stol.OR.slgclat2d.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
     WRITE(6,*)'Requested latitude of point 2 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ELSE IF(slgclon1d.LT.gop(1)-stol.OR.slgclon1d.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
     WRITE(6,*)'Requested longitude of point 1 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ELSE IF(slgclon2d.LT.gop(1)-stol.OR.slgclon2d.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
     WRITE(6,*)'Requested longitude of point 2 for great circle'
     WRITE(6,*)'slice lies outside model bounds!'
     WRITE(6,*)'TERMINATING PROGRAM!!!'
     STOP
  ENDIF
!
! Compute azimuth of great circle cross-section
!
  azim=COS(slgclat1)*SIN(slgclat2)-COS(slgclat2)*SIN(slgclat1)*COS(slgclon2-slgclon1)
  azim=ACOS(azim/SIN(deltas))
!
! Allocate memory to velocity grid array
!
  nnx=ngcd
  nnz=ngch
  ALLOCATE(vela(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vela'
  ENDIF
!
! Now allocate memory to pointer array which
! indicates which layer each grid point lies in.
!
  ALLOCATE(lay(nnx,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL lay'
  ENDIF
!
! Loop through all the layers and determine which layer each
! point lies in.
!
  lay=1
  IF(ni.GT.2)THEN
     ALLOCATE(inta(1,nnz), STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
     DO i=2,ni-1
!
!       Call bspline subroutine to compute depth to interface i. For each
!       point, we need to find the corresponding (u,v) values and
!       in which b-spline cell it lies in.
!
        DO j=1,nnz
           deltanp=deltas*(j-1)/REAL(nnz-1)
           thgcp=SIN(slgclat1)*COS(deltanp)+COS(slgclat1)*COS(azim)*SIN(deltanp)
           IF(thgcp.LT.-1.0)THEN
              thgcp=-1.0
           ELSE IF(thgcp.GT.1.0)THEN
              thgcp=1.0
           ENDIF
           thgcp=ASIN(thgcp)
           phgcp=COS(deltanp)-SIN(thgcp)*SIN(slgclat1)
           phgcp=phgcp/(COS(thgcp)*COS(slgclat1))
           IF(phgcp.LT.-1.0)THEN
              phgcp=-1.0
           ELSE IF(phgcp.GT.1.0)THEN
              phgcp=1.0
           ENDIF
           phgcp=ACOS(phgcp)
           IF(azim.GE.0.0.AND.azim.LE.pi)THEN
              phgcp=ABS(phgcp)
           ELSE
              phgcp=-ABS(phgcp)
           ENDIF
           IF(slgclon1.LT.slgclon2)THEN
              phgcp=slgclon1+phgcp
           ELSE
              phgcp=slgclon1-phgcp
           ENDIF
           thgcp=thgcp*180.0/pi
           phgcp=phgcp*180.0/pi
           !
           ! Make sure point lies within model
           !
           IF(thgcp.GT.got(1)+stol.OR.thgcp.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
              WRITE(6,*)'Great circle slice latitude'
              WRITE(6,*)'lies outside model bounds!'
              WRITE(6,*)'TERMINATING PROGRAM!!!'
              STOP
           ELSE IF(phgcp.LT.gop(1)-stol.OR.phgcp.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
              WRITE(6,*)'Great circle slice longitude'
              WRITE(6,*)'lies outside model bounds!'
              WRITE(6,*)'TERMINATING PROGRAM!!!'
              STOP
           ENDIF
!
!          Now identify the cell and corresponding (u,v) values
!
           rnodeth=INT((got(i)-thgcp)/gsit)+1
           u=ABS(got(i)-thgcp-(rnodeth-1)*gsit)/gsit
           IF(rnodeth.EQ.nint)THEN
              rnodeth=rnodeth-1
              u=1.0
           ENDIF
           rnodeph=INT((phgcp-gop(i))/gsip)+1
           v=ABS(phgcp-gop(i)-(rnodeph-1)*gsip)/gsip
           IF(rnodeph.EQ.ninp)THEN
              rnodeph=rnodeph-1
              v=1.0
           ENDIF
           CALL ibspline(rnodeth,rnodeph,u,v,1,j,i)
        ENDDO      
        DO j=1,nnz
           DO k=1,nnx
              rd1=gor(idmnr)-(k-1)*(nnr(i)-1)*gsr(i)/REAL(nnx-1)
              IF(rd1.LT.inta(1,j))lay(k,j)=i
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE(inta, STAT=checkstat)
     IF(checkstat > 0)THEN
        WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL inta'
     ENDIF
  ENDIF
!
! Now call a bspline subroutine to compute great circle slice
!
  DO j=1,nnz
     deltanp=deltas*(j-1)/REAL(nnz-1)
     thgcp=SIN(slgclat1)*COS(deltanp)+COS(slgclat1)*COS(azim)*SIN(deltanp)
     IF(thgcp.LT.-1.0)THEN
        thgcp=-1.0
     ELSE IF(thgcp.GT.1.0)THEN
        thgcp=1.0
     ENDIF
     thgcp=ASIN(thgcp)
     phgcp=COS(deltanp)-SIN(thgcp)*SIN(slgclat1)
     phgcp=phgcp/(COS(thgcp)*COS(slgclat1))
     IF(phgcp.LT.-1.0)THEN
        phgcp=-1.0
     ELSE IF(phgcp.GT.1.0)THEN
        phgcp=1.0
     ENDIF
     phgcp=ACOS(phgcp)
     IF(azim.GE.0.0.AND.azim.LE.pi)THEN
        phgcp=ABS(phgcp)
     ELSE
        phgcp=-ABS(phgcp)
     ENDIF
     IF(slgclon1.LT.slgclon2)THEN
        phgcp=slgclon1+phgcp
     ELSE
        phgcp=slgclon1-phgcp
     ENDIF
     thgcp=thgcp*180.0/pi
     phgcp=phgcp*180.0/pi
     !
     ! Make sure point lies within model
     !
     IF(thgcp.GT.got(1)+stol.OR.thgcp.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
        WRITE(6,*)'Great circle slice latitude'
        WRITE(6,*)'lies outside model bounds!'
        WRITE(6,*)'TERMINATING PROGRAM!!!'
        STOP
     ELSE IF(phgcp.LT.gop(1)-stol.OR.phgcp.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
        WRITE(6,*)'Great circle slice longitude'
        WRITE(6,*)'lies outside model bounds!'
        WRITE(6,*)'TERMINATING PROGRAM!!!'
        STOP
     ENDIF
     DO k=1,nnx
        idm3=lay(k,j)
        rnodeth=INT((got(idm3)-thgcp)/gst(idm3))+1
        u=ABS(got(idm3)-thgcp-(rnodeth-1)*gst(idm3))/gst(idm3)
        IF(rnodeth.EQ.nnt(idm3))THEN
           rnodeth=rnodeth-1
           u=1.0
        ENDIF
        rnoder=INT((k-1)*(nnr(idm3)-1)/REAL(nnx-1))+1
        v=(k-1)*(nnr(idm3)-1)/REAL(nnx-1)-(rnoder-1)
        IF(rnoder.EQ.nnr(idm3))THEN
           rnoder=rnoder-1
           v=1.0
        ENDIF
        rnodeph=INT((phgcp-gop(idm3))/gsp(idm3))+1
        w=ABS(phgcp-gop(idm3)-(rnodeph-1)*gsp(idm3))/gsp(idm3)
        IF(rnodeph.EQ.nnp(idm3))THEN
           rnodeph=rnodeph-1
           w=1.0
        ENDIF
        CALL vbspline(rnoder,rnodeth,rnodeph,u,w,v,k,j,idm3)
     ENDDO
  ENDDO
!
! Now write out great circle slice to file.
!
  OPEN(UNIT=30,FILE=ofilegcv,STATUS='unknown')
  DO i=1,nnz
     DO j=nnx,1,-1
        WRITE(30,*)vela(j,i)
     ENDDO
  ENDDO
  CLOSE(30)
!
! Compute interface cross-sections and write to file
!
  OPEN(UNIT=30,FILE=ofilegci,STATUS='unknown')
  ALLOCATE(inta(1,nnz), STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
  ENDIF
  DO i=1,ni
     goti=got(i)
     gopi=gop(i)
     IF(i.EQ.ni)THEN
        goti=got(i-1)
        gopi=gop(i-1)
     ENDIF
     DO j=1,nnz
        deltanp=deltas*(j-1)/REAL(nnz-1)
        thgcp=SIN(slgclat1)*COS(deltanp)+COS(slgclat1)*COS(azim)*SIN(deltanp)
        IF(thgcp.LT.-1.0)THEN
           thgcp=-1.0
        ELSE IF(thgcp.GT.1.0)THEN
           thgcp=1.0
        ENDIF
        thgcp=ASIN(thgcp)
        phgcp=COS(deltanp)-SIN(thgcp)*SIN(slgclat1)
        phgcp=phgcp/(COS(thgcp)*COS(slgclat1))
        IF(phgcp.LT.-1.0)THEN
           phgcp=-1.0
        ELSE IF(phgcp.GT.1.0)THEN
           phgcp=1.0
        ENDIF
        phgcp=ACOS(phgcp)
        IF(azim.GE.0.0.AND.azim.LE.pi)THEN
           phgcp=ABS(phgcp)
        ELSE
           phgcp=-ABS(phgcp)
        ENDIF
        IF(slgclon1.LT.slgclon2)THEN
           phgcp=slgclon1+phgcp
        ELSE
           phgcp=slgclon1-phgcp
        ENDIF
        thgcp=thgcp*180.0/pi
        phgcp=phgcp*180.0/pi
!
!       Now identify the cell and corresponding (u,v) values
!
        rnodeth=INT((goti-thgcp)/gsit)+1
        u=ABS(goti-thgcp-(rnodeth-1)*gsit)/gsit
        IF(rnodeth.EQ.nint)THEN
           rnodeth=rnodeth-1
           u=1.0
        ENDIF
        rnodeph=INT((phgcp-gopi)/gsip)+1
        v=ABS(phgcp-gopi-(rnodeph-1)*gsip)/gsip
        IF(rnodeph.EQ.ninp)THEN
           rnodeph=rnodeph-1
           v=1.0
        ENDIF
        CALL ibspline(rnodeth,rnodeph,u,v,1,j,i)
        rd1=(j-1)*deltaskm/(nnz-1)
        WRITE(30,*)rd1,inta(1,j)
     ENDDO
     WRITE(30,'(a1)')sep
  ENDDO
  CLOSE(30)
  DEALLOCATE(vela,lay,inta, STAT=checkstat)
  IF(checkstat > 0)THEN
     WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL vela,lay,inta'
  ENDIF
ENDIF
!
! Now compute interface depths for contour plots if required.
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)pard
READ(10,*)nip
READ(10,*)indt,indp
READ(10,'(a26)')ofileib
IF(nip.GE.1)THEN
!
! Compute relative depth if required
!
   IF(pard.EQ.1)THEN
      OPEN(UNIT=20,FILE=ifileir,status='old')
      READ(20,*)idm1
      READ(20,*)idm2,idm3
      idm2=idm2-2
      idm3=idm3-2
      IF(idm1.NE.ni.OR.idm2.NE.nint.OR.idm3.NE.ninp)THEN
         WRITE(6,*)'True and reference interface grids'
         WRITE(6,*)'not consistent.'
         WRITE(6,*)'TERMINATING PROGRAM!!'
         STOP
      ENDIF
      READ(20,*)rd1,rd2
      READ(20,*)rd1,rd2
      DO ii=1,ni
         DO j=0,nint+1
            DO k=0,ninp+1
               READ(20,*)rd1
               intn(nnt+1-j,k,ii)=intn(nint+1-j,k,ii)-(rd1-earthr)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(20)
   ENDIF
!
! Calculate GMT bounds file for surface
!
   lft=gop(1)
   rgt=gop(1)+(ninp-1)*gsip
   btm=got(1)-(nint-1)*gsit
   top=got(1)
   OPEN(UNIT=80,FILE=ofileib,STATUS='unknown')
   WRITE(80,*)lft
   WRITE(80,*)rgt
   WRITE(80,*)btm
   WRITE(80,*)top
   nnx=(nint-1)*indt+1
   nnz=(ninp-1)*indp+1
   WRITE(80,*)nnz
   WRITE(80,*)nnx
   CLOSE(80)
   DO i=1,nip
      READ(10,*)intid
      READ(10,'(a26)')ofileint
      nnx=(nint-1)*indt+1
      nnz=(ninp-1)*indp+1
      ALLOCATE(inta(nnx,nnz), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL inta'
      ENDIF
!
!     Call subroutine to compute interface depths and write to file
!
      CALL ibspline2d(intid,indt,indp)
      OPEN(UNIT=90,FILE=ofileint,STATUS='unknown')
      DO j=1,nnz
         DO k=nnx,1,-1
            WRITE(90,*)inta(k,j)
         ENDDO
      ENDDO
      CLOSE(90)
      DEALLOCATE(inta, STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL inta'
      ENDIF
   ENDDO
ENDIF
DEALLOCATE(veln,intn, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM slice: REAL veln,intn'
ENDIF
!
! Plot raypaths if required
!
IF(prp.EQ.1)THEN
   OPEN(UNIT=20,FILE=ifileray,STATUS='old',IOSTAT=ios)
   OPEN(UNIT=30,FILE=ofilerd,STATUS='unknown')
   OPEN(UNIT=40,FILE=ofilerew,STATUS='unknown')
   OPEN(UNIT=50,FILE=ofilerns,STATUS='unknown')
   IF(extrgcs.EQ.1)OPEN(UNIT=60,FILE=ofilergc,STATUS='unknown')
   DO i=1,maxr
      READ(20,*,IOSTAT=ios)idm1,idm2,idm3,idm4,nrs
      IF(ios.LT.0)EXIT
      IF(nrs.EQ.0)CYCLE
      DO j=1,nrs
         READ(20,*)nrp
         DO k=1,nrp
            READ(20,*)rdep,rlat,rlong
            rdep=rdep-earthr
            rlat=rlat*180./pi
            rlong=rlong*180./pi
            IF(MOD(k-1,penp).EQ.0.OR.k.EQ.nrp)THEN
               WRITE(30,*)rlong,rlat
               WRITE(40,*)rlong,rdep
               WRITE(50,*)rlat,rdep
            ENDIF
!
!           Projecting the ray paths onto the great circle
!           plane is more complicated.
!
            IF(extrgcs.EQ.1)THEN
!
!              First, define the great circle that is perpendicular to
!              the required great cricle slice
!  
               deltanp=deltas
               rlat=rlat*pi/180.0
               rlong=rlong*pi/180.0
               thgcp=SIN(rlat)*COS(deltanp)+COS(rlat)*COS(azim-pi/2)*SIN(deltanp)
               IF(thgcp.LT.-1.0)THEN
                  thgcp=-1.0
               ELSE IF(thgcp.GT.1.0)THEN
                  thgcp=1.0
               ENDIF
               thgcp=ASIN(thgcp)
               phgcp=COS(deltanp)-SIN(thgcp)*SIN(rlat)
               phgcp=phgcp/(COS(thgcp)*COS(rlat))
               IF(phgcp.LT.-1.0)THEN
                  phgcp=-1.0
               ELSE IF(phgcp.GT.1.0)THEN
                  phgcp=1.0
               ENDIF
               phgcp=ACOS(phgcp)
               IF(azim-pi/2.GE.0.0.AND.azim-pi/2.LE.pi)THEN
                  phgcp=ABS(phgcp)
               ELSE
                  phgcp=-ABS(phgcp)
               ENDIF
               IF(slgclon1.LT.slgclon2)THEN
                  phgcp=rlong+phgcp
               ELSE
                  phgcp=rlong-phgcp
               ENDIF
               rd1=rlat
               rd2=rlong
               rd3=thgcp
               rd4=phgcp
               rd5=slgclat1
               rd6=slgclon1
               rd7=slgclat2
               rd8=slgclon2
               CALL gcint(rd5,rd6,rd7,rd8,rd1,rd2,rd3,rd4)
!
!              Now determine which of the candidate points lies in the region
!
               IF(rd1.GT.got(1)+stol.OR.rd1.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
                  rd1=rd3
                  IF(rd1.GT.got(1)+stol.OR.rd1.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
                     CYCLE
                  ENDIF
               ENDIF
               IF(rd2.LT.gop(1)-stol.OR.rd2.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
                  rd2=rd4
                  IF(rd2.LT.gop(1)-stol.OR.rd2.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
                     CYCLE
                  ENDIF
               ENDIF
!
!              Finally, convert to distance
!
               rd1=rd1*pi/180.0
               rd2=rd2*pi/180.0
               deltanp=SIN(rd1)*SIN(slgclat1)
               deltanp=deltanp+COS(rd1)*COS(slgclat1)*COS(rd2-slgclon1)
               deltanp=ACOS(deltanp)
               rd1=deltanp*earthr
               IF(MOD(k-1,penp).EQ.0.OR.k.EQ.nrp)WRITE(60,*)rd1,rdep
            ENDIF
         ENDDO
      ENDDO
      WRITE(30,'(a1)')sep
      WRITE(40,'(a1)')sep
      WRITE(50,'(a1)')sep
      IF(extrgcs.EQ.1)WRITE(60,'(a1)')sep
   ENDDO
   CLOSE(20)
   CLOSE(30)
   CLOSE(40)
   CLOSE(50)
   IF(extrgcs.EQ.1)CLOSE(60)
ENDIF
!
! Plot sources and receivers if required (teleseismic sources
! not plotted).
!
IF(psr.EQ.1)THEN
!
! Start with sources
!
   maxpth=0
   OPEN(UNIT=20,FILE=ifilesrc,STATUS='old')
   READ(20,*)ns
   OPEN(UNIT=30,FILE=ofilesd,STATUS='unknown')
   OPEN(UNIT=40,FILE=ofilesew,STATUS='unknown')
   OPEN(UNIT=50,FILE=ofilesns,STATUS='unknown')
   IF(extrgcs.EQ.1)OPEN(UNIT=60,FILE=ofilesgc,STATUS='unknown')
   DO i=1,ns
      READ(20,*)idm1
      IF(idm1.EQ.1)READ(20,*)
      READ(20,*)rdep,rlat,rlong
      IF(idm1.EQ.0)THEN
         WRITE(30,*)rlong,rlat
         WRITE(40,*)rlong,-rdep
         WRITE(50,*)rlat,-rdep
      ENDIF
      READ(20,*)idm2
      IF(maxpth.LT.idm2)maxpth=idm2 
      DO j=1,idm2
         READ(20,*)idm3
         READ(20,*)paths(1:2*idm3)
         READ(20,*)paths(1:idm3)
      ENDDO
!
!     Projecting the sources onto the great circle
!     plane is more complicated.
!
      IF(extrgcs.EQ.1)THEN
!
!        First, define the great circle that is perpendicular to
!        the required great cricle slice
!
         deltanp=deltas
         rlat=rlat*pi/180.0
         rlong=rlong*pi/180.0
         thgcp=SIN(rlat)*COS(deltanp)+COS(rlat)*COS(azim-pi/2)*SIN(deltanp)
         IF(thgcp.LT.-1.0)THEN
            thgcp=-1.0
         ELSE IF(thgcp.GT.1.0)THEN
            thgcp=1.0
         ENDIF
         thgcp=ASIN(thgcp)
         phgcp=COS(deltanp)-SIN(thgcp)*SIN(rlat)
         phgcp=phgcp/(COS(thgcp)*COS(rlat))
         IF(phgcp.LT.-1.0)THEN
            phgcp=-1.0
         ELSE IF(phgcp.GT.1.0)THEN
            phgcp=1.0
         ENDIF
         phgcp=ACOS(phgcp)
         IF(azim-pi/2.GE.0.0.AND.azim-pi/2.LE.pi)THEN
            phgcp=ABS(phgcp)
         ELSE
            phgcp=-ABS(phgcp)
         ENDIF
         IF(slgclon1.LT.slgclon2)THEN
            phgcp=rlong+phgcp
         ELSE
            phgcp=rlong-phgcp
         ENDIF
         rd1=rlat
         rd2=rlong
         rd3=thgcp
         rd4=phgcp
         rd5=slgclat1
         rd6=slgclon1
         rd7=slgclat2
         rd8=slgclon2
         CALL gcint(rd5,rd6,rd7,rd8,rd1,rd2,rd3,rd4)
!
!        Now determine which of the candidate points lies in the region
!
         IF(rd1.GT.got(1)+stol.OR.rd1.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
            rd1=rd3
            IF(rd1.GT.got(1)+stol.OR.rd1.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
               CYCLE
            ENDIF
         ENDIF
         IF(rd2.LT.gop(1)-stol.OR.rd2.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
            rd2=rd4
            IF(rd2.LT.gop(1)-stol.OR.rd2.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
               CYCLE
            ENDIF
         ENDIF
!
!        Finally, convert to distance
!
         rd1=rd1*pi/180.0
         rd2=rd2*pi/180.0
         deltanp=SIN(rd1)*SIN(slgclat1)
         deltanp=deltanp+COS(rd1)*COS(slgclat1)*COS(rd2-slgclon1)
         deltanp=ACOS(deltanp)
         rd1=deltanp*earthr
         WRITE(60,*)rd1,-rdep
      ENDIF
   ENDDO
   CLOSE(20)
   CLOSE(30)
   CLOSE(40)
   CLOSE(50)
   IF(extrgcs.EQ.1)CLOSE(60)
!
!  Now do receivers
!
   OPEN(UNIT=20,FILE=ifilercv,STATUS='old')
   READ(20,*)nr
   ALLOCATE(srcid(ns*maxpth))
   OPEN(UNIT=30,FILE=ofilercd,STATUS='unknown')
   OPEN(UNIT=40,FILE=ofilercew,STATUS='unknown')
   OPEN(UNIT=50,FILE=ofilercns,STATUS='unknown')
   IF(extrgcs.EQ.1)OPEN(UNIT=60,FILE=ofilercgc,STATUS='unknown')
   DO i=1,nr
      READ(20,*)rdep,rlat,rlong
      WRITE(30,*)rlong,rlat
      WRITE(40,*)rlong,-rdep
      WRITE(50,*)rlat,-rdep
      READ(20,*)idm2
      READ(20,*)srcid(1:idm2)
      READ(20,*)srcid(1:idm2)
!
!     Projecting the receivers onto the great circle
!     plane is more complicated.
!
      IF(extrgcs.EQ.1)THEN
!
!        First, define the great circle that is perpendicular to
!        the required great cricle slice
!
         deltanp=deltas
         rlat=rlat*pi/180.0
         rlong=rlong*pi/180.0
         thgcp=SIN(rlat)*COS(deltanp)+COS(rlat)*COS(azim-pi/2)*SIN(deltanp)
         IF(thgcp.LT.-1.0)THEN
            thgcp=-1.0
         ELSE IF(thgcp.GT.1.0)THEN
            thgcp=1.0
         ENDIF
         thgcp=ASIN(thgcp)
         phgcp=COS(deltanp)-SIN(thgcp)*SIN(rlat)
         phgcp=phgcp/(COS(thgcp)*COS(rlat))
         IF(phgcp.LT.-1.0)THEN
            phgcp=-1.0
         ELSE IF(phgcp.GT.1.0)THEN
            phgcp=1.0
         ENDIF
         phgcp=ACOS(phgcp)
         IF(azim-pi/2.GE.0.0.AND.azim-pi/2.LE.pi)THEN
            phgcp=ABS(phgcp)
         ELSE
            phgcp=-ABS(phgcp)
         ENDIF
         IF(slgclon1.LT.slgclon2)THEN
            phgcp=rlong+phgcp
         ELSE
            phgcp=rlong-phgcp
         ENDIF
         rd1=rlat
         rd2=rlong
         rd3=thgcp
         rd4=phgcp
         rd5=slgclat1
         rd6=slgclon1
         rd7=slgclat2
         rd8=slgclon2
         CALL gcint(rd5,rd6,rd7,rd8,rd1,rd2,rd3,rd4)
!
!        Now determine which of the candidate points lies in the region
!
         IF(rd1.GT.got(1)+stol.OR.rd1.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
            rd1=rd3
            IF(rd1.GT.got(1)+stol.OR.rd1.LT.got(1)-gst(1)*(nnt(1)-1)-stol)THEN
               CYCLE
            ENDIF
         ENDIF
         IF(rd2.LT.gop(1)-stol.OR.rd2.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
            rd2=rd4
            IF(rd2.LT.gop(1)-stol.OR.rd2.GT.gop(1)+gsp(1)*(nnp(1)-1)+stol)THEN
               CYCLE
            ENDIF
         ENDIF
!
!        Finally, convert to distance
!
         rd1=rd1*pi/180.0
         rd2=rd2*pi/180.0
         deltanp=SIN(rd1)*SIN(slgclat1)
         deltanp=deltanp+COS(rd1)*COS(slgclat1)*COS(rd2-slgclon1)
         deltanp=ACOS(deltanp)
         rd1=deltanp*earthr
         WRITE(60,*)rd1,-rdep
      ENDIF
   ENDDO
   DEALLOCATE(srcid)
   CLOSE(20)
   CLOSE(30)
   CLOSE(40)
   CLOSE(50)
   IF(extrgcs.EQ.1)CLOSE(60)
ENDIF
CLOSE(10)
END PROGRAM slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine uses cubic B-spline interpolation return a
! value from a 2-D interface grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ibspline(gt,gp,u,v,stt,stp,iid)
USE globalp
IMPLICIT NONE
INTEGER :: i1,j1,stp,stt,iid,gt,gp
REAL(KIND=i10) :: u,v,sumi,sumj
REAL(KIND=i10), DIMENSION(4) :: ui,vi
!
! iid = Interface id.
! u,v = independent surface parameters
! ui,vi = bspline basis functions
! stt,stp = global grid coordinate of diced node
! sumi,sumj = Summation variables for spline
! gt,gp =  i,j coordinate of current node cell
!
! Compute the values of the basis functions
!
ui(1)=(1.0-u)**3/6.0
ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
ui(4)=u**3/6.0
vi(1)=(1.0-v)**3/6.0
vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
vi(4)=v**3/6.0
sumi=0.0
DO i1=1,4
   sumj=0.0
   DO j1=1,4
      sumj=sumj+ui(j1)*intn(gt-2+j1,gp-2+i1,iid)
   ENDDO
   sumi=sumi+vi(i1)*sumj
ENDDO
inta(stt,stp)=sumi
END SUBROUTINE ibspline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine uses cubic B-spline interpolation to take
! N-S or E-W slices through a surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ibspline1d(iid,ddt,rnode,w,isl)
USE globalp
IMPLICIT NONE
INTEGER :: i,m,i1,j1,stt,iid,ddt
INTEGER :: cont,rnode,isl,checkstat,nn
REAL(KIND=i10) :: u,w,sumi,sumj
REAL(KIND=i10), DIMENSION(:) :: wi(4)
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui
!
! iid = Interface id.
! ddt = dicing of interface grid
! u,w = independent surface parameters
! ui,wi = bspline basis functions
! stt = global grid coordinate of diced node
! cont,conp = dicing factor for node division
! sumi,sumj = Summation variables for spline
! isl = slice type (1=N-S, 2=E-W)
! nn = number of nodes
!
! Compute the values of the basis functions
!
ALLOCATE(ui(ddt+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
ENDIF
DO i=1,ddt+1
   u=ddt
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
wi(1)=(1.0-w)**3/6.0
wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
wi(4)=w**3/6.0
!
! Set parameters for slice type
!
IF(isl.EQ.1)THEN
  nn=nint
ELSE
  nn=ninp
ENDIF
WRITE(6,*)nint,ninp
DO i=1,nn-1
   cont=ddt
   IF(i==nn-1)cont=ddt+1
   DO m=1,cont
      WRITE(6,*)'hello'
      stt=ddt*(i-1)+m
      sumi=0.0
      DO i1=1,4
         sumj=0.0
         DO j1=1,4
            IF(isl.EQ.1)THEN
               sumj=sumj+ui(m,j1)*intn(i-2+j1,rnode-2+i1,iid)
            ELSE
               sumj=sumj+wi(j1)*intn(rnode-2+j1,i-2+i1,iid)
            ENDIF
         ENDDO
         IF(isl.EQ.1)THEN
            sumi=sumi+wi(i1)*sumj
         ELSE
            sumi=sumi+ui(m,i1)*sumj
         ENDIF
      ENDDO
      WRITE(6,*)stt,'hello1'
      inta(1,stt)=sumi
   ENDDO
ENDDO
DEALLOCATE(ui, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE ibspline1d: REAL ui'
ENDIF
END SUBROUTINE ibspline1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine uses cubic B-spline interpolation to
! dice up an interface grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ibspline2d(iid,ddt,ddp)
USE globalp
IMPLICIT NONE
INTEGER :: i,j,l,m,i1,j1,stp,stt,iid,ddt,ddp
INTEGER :: cont,conp,checkstat
REAL(KIND=i10) :: u,sumi,sumj
REAL(KIND=i10), DIMENSION(:,:), ALLOCATABLE :: ui,vi
!
! iid = Interface id.
! ddt,ddp = dicing of interface grid
! u = independent surface parameter
! ui,vi = bspline basis functions
! stt,stp = global grid coordinate of diced node
! cont,conp = dicing factor for node division
! sumi,sumj = Summation variables for spline
!
! Compute the values of the basis functions
!
ALLOCATE(ui(ddt+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL ui'
ENDIF
DO i=1,ddt+1
   u=ddt
   u=(i-1)/u
   ui(i,1)=(1.0-u)**3/6.0
   ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   ui(i,4)=u**3/6.0
ENDDO
ALLOCATE(vi(ddp+1,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM slice: REAL vi'
ENDIF
DO i=1,ddp+1
   u=ddp
   u=(i-1)/u
   vi(i,1)=(1.0-u)**3/6.0
   vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
   vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
   vi(i,4)=u**3/6.0
ENDDO
!
! Work out the value of u for the given depth
!
DO i=1,ninp-1
   conp=ddp
   IF(i==ninp-1)conp=ddp+1
   DO j=1,nint-1
      cont=ddt
      IF(j==nint-1)cont=ddt+1
      DO l=1,conp
         stp=ddp*(i-1)+l
         DO m=1,cont
            stt=ddt*(j-1)+m
            sumi=0.0
            DO i1=1,4
               sumj=0.0
               DO j1=1,4
                  sumj=sumj+ui(m,j1)*intn(j-2+j1,i-2+i1,iid)
               ENDDO
               sumi=sumi+vi(l,i1)*sumj
            ENDDO
            inta(stt,stp)=sumi
         ENDDO
      ENDDO
   ENDDO
ENDDO
DEALLOCATE(ui,vi, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: SUBROUTINE ibspline2d: REAL ui,vi'
ENDIF
END SUBROUTINE ibspline2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine uses cubic B-spline interpolation to sample
! a point wthiin a 3-D velocity grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vbspline(rnode,gt,gp,u,v,w,stt,stp,lid)
USE globalp
IMPLICIT NONE
INTEGER :: i1,j1,k1,stp,stt,gt,gp
INTEGER :: rnode,lid
REAL(KIND=i10), DIMENSION(4) :: ui,vi,wi
REAL(KIND=i10) :: u,v,w,sumi,sumj,sumk,rdm
!
! u,v,w = independent surface parameters
! ui,vi,wi = bspline basis functions
! stt,stp = global grid coordinate of diced node
! sumi,sumj,sumk = summation variables for spline
! lid= layer id.
! gt,gp = i,j coordinate of current cell
! rnode = Depth coordinate of current cell
!
! Compute the values of the basis functions
!
ui(1)=(1.0-u)**3/6.0
ui(2)=(4.0-6.0*u**2+3.0*u**3)/6.0
ui(3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
ui(4)=u**3/6.0
vi(1)=(1.0-v)**3/6.0
vi(2)=(4.0-6.0*v**2+3.0*v**3)/6.0
vi(3)=(1.0+3.0*v+3.0*v**2-3.0*v**3)/6.0
vi(4)=v**3/6.0
wi(1)=(1.0-w)**3/6.0
wi(2)=(4.0-6.0*w**2+3.0*w**3)/6.0
wi(3)=(1.0+3.0*w+3.0*w**2-3.0*w**3)/6.0
wi(4)=w**3/6.0
!
! Set parameters for slice type.
!
sumi=0.0
DO i1=1,4
   sumj=0.0
   DO j1=1,4
      sumk=0.0
      DO k1=1,4
         rdm=wi(k1)*veln(rnode-2+k1,gt-2+j1,gp-2+i1,lid)
         sumk=sumk+rdm
      ENDDO
      sumj=sumj+ui(j1)*sumk
   ENDDO
   sumi=sumi+vi(i1)*sumj
ENDDO
vela(stt,stp)=sumi
END SUBROUTINE vbspline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: SUBROUTINE
! CODE: FORTRAN 90
! This subroutine finds the intersection points between two
! great circle planes on the surface of the Earth.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gcint(lat0,lon0,lat1,lon1,lat2,lon2,lat3,lon3)
USE globalp
IMPLICIT NONE
REAL(KIND=i10) :: lat0,lon0,lat1,lon1,lat2,lon2,lat3,lon3
REAL(KIND=i10) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3
REAL(KIND=i10) :: a,b,c,d,e,f,g,h,k
!
! Convert to co-latitude
!
lat0=pi/2-lat0
lat1=pi/2-lat1
lat2=pi/2-lat2
lat3=pi/2-lat3
!
! Convert to Cartesian coordinates
!
x0=earthr*COS(lon0)*SIN(lat0)
y0=earthr*SIN(lon0)*SIN(lat0)
z0=earthr*COS(lat0)
x1=earthr*COS(lon1)*SIN(lat1)
y1=earthr*SIN(lon1)*SIN(lat1)
z1=earthr*COS(lat1)
x2=earthr*COS(lon2)*SIN(lat2)
y2=earthr*SIN(lon2)*SIN(lat2)
z2=earthr*COS(lat2)
x3=earthr*COS(lon3)*SIN(lat3)
y3=earthr*SIN(lon3)*SIN(lat3)
z3=earthr*COS(lat3)
!
! Find intersection points
!
a=(y0*z1-y1*z0)
b=-(x0*z1-x1*z0)
c=(x0*y1-x1*y0)
d=(y2*z3-y3*z2)
e=-(x2*z3-x3*z2)
f=(x2*y3-x3*y2)
h=(d*c-f*a)/(e*a-d*b)
g=(-b*h-c)/a
k=SQRT(earthr**2/(g**2+h**2+1))
!
! Convert back to longitude
!
lat2=ASIN(k/earthr)*180.0/pi
lon2=ATAN(h/g)*180.0/pi
lat3=-lat2
lon3=lon2+180.0
END SUBROUTINE gcint
