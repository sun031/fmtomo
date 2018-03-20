!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN PROGRAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TYPE: PROGRAM
! CODE: FORTRAN 90
! This program is designed to generate a layered velocity
! medium in 3-D spherical coordinates. For n layers, there will
! be n+1 interfaces. Each layer is described by a uniform
! velocity grid, and each interface is described by a uniform
! depth grid. Although the background velocity can only
! vary with depth, random noise and checkerboard patterns
! can be superimposed. Likewise, interfaces can only be planar,
! but random noise or checkerboard patterns may be
! superimposed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM gridder
IMPLICIT NONE
INTEGER, PARAMETER :: i10=SELECTED_REAL_KIND(10,100)
INTEGER :: ii,i,j,k,l,m,ars,umocvg,checkstat
INTEGER :: nvr,nvt,nvp,pors,aapmc,adch,nch,usp
INTEGER :: vusp,vusp1,vusp2,vuspo,vusp1o,vusp2o
INTEGER :: asp,nsp,ni,nvgt,rseed,ipsw,ips,idum
INTEGER :: npr,npt,npp,npgr,npcell,ueif,dovm
INTEGER, DIMENSION (100) :: ispr, ispt,ispp
REAL(KIND=i10) :: gor,got,gop,gsr,gst,gsp,rdum
REAL(KIND=i10) :: glr,glt,glp,veltol
REAL(KIND=i10) :: vtop,vbot,vgr,vel,rssf
REAL(KIND=i10) :: d1,d2,v1,v2,rdep,decm,mpc
REAL(KIND=i10) :: chvp,chvp1,chvp2,pnchd
REAL(KIND=i10) :: ggor,ggot,ggop,dep,earthr
REAL(KIND=i10) :: dnw,dne,dsw,xnw,xne,xsw,znw,zne,zsw
REAL(KIND=i10) :: aplc,bplc,cplc,xpl,zpl,grcf,rd1
REAL(KIND=i10) :: gorp,gotp,gopp,gsrp,gstp,gspp
REAL(KIND=i10), DIMENSION (:), ALLOCATABLE :: velm
REAL(KIND=i10), DIMENSION (:,:), ALLOCATABLE :: pdep
REAL(KIND=i10), DIMENSION (100) :: spa,spr,spt,spp
REAL(KIND=i10), PARAMETER :: pi=3.141592653589793
REAL, EXTERNAL :: gasdev
CHARACTER (LEN=20) :: ovfile,oifile,opfile,cdum
CHARACTER (LEN=20) :: extifile,extvfile
!
! ni = number of interfaces
! nvgt = number of velocity grid types
! nvr = number of vertices in radial direction
! nvt = number of vertices in theta (N-S) direction
! nvp = number of vertices in phi (E-W) direction
! gor = grid origin in radius
! got = grid origin in theta (N-S)
! gop = grid origin in phi (E-W)
! glr = grid limit in radius
! glt = grid limit in theta (N-S)
! glp = grid limit in phi (E-W)
! gsr = grid separation in radius
! gst = grid separation in theta (N-S)
! gsp = grid separation in phi (E-W)
! ovfile = Name of output velocity grid file
! oifile = Name of output interface grid file
! vtop,vbot = Velocity at top and bottom of grid
! vgr = Vertical velocity gradient
! vel = velocity of grid node
! dep = velocity node depth
! ars = Add random structure (0=no,1=yes)
! rssf = Random structure scaling factor
! umocvg = Use velocity model (0) or constant gradient (1)
! velm = Velocity discretized at depth from a 1-D model
! pors = P (0) or S (1) velocity model
! r1,r2 = radius at which v1,v2 occur for 1-D model
! v1,v2 = model velocities at r1,r2
! rdep = Required depth at which node requires model velocity
! aapmc = Add a priori model covariance (0=no, 1=yes)
! decm = Diagonal elements of covariance matrix
! adch = Add checkerboard (0=no,1=yes)
! mpc = Maximum perturbation of checkerboard
! nch = size of checkerboard cell
! chvp = Checkerboard velocity perturbation
! chvp1,chvp2 = Dummy checkerboard variables
! usp = use spacing for checkerboard (0=no,1=yes)
! vusp = Checkerboard spacing variable
! vusp1,vusp2 = Dummy spacing variables
! vuspo,vusp1o,vusp2o = Previous values of above
! asp = Apply spike (0=no,1=yes)
! nsp = Number of spikes
! spa = Amplitude of spikes
! spr,spr,spp = Coordinates of spikes
! ispr,ispt,ispp = Grid locations of spikes
! earthr = Earth radius in km
! ggor,ggot,ggop = Global gor,got,gop
! dnw,dne,dsw = Depth of interface at NW, NE and SW grid
! xnw,xne,xsw = Same as above but for longitude
! znw,zne,zsw = Same as above but for latitude
! aplc,bplc,cplc = Linear interpolation coefficients
! xpl,zpl = lat and long points on plane
! gasdev = Function for returning random noise
! rseed = Random seed
! pnchd = Pinchout distance for minimum layer thickness
! pdep = Depth to previous (overlying) interface
! ipsw = Interface pinchout switch (1=pinchouts created)
! ips = Counter for P(1) and S(2) velocity fields
! idum = Dummy variable
! opfile = Output file for propagation grid
! grcf = Grid cushion factor for propagation gird
! npr,npt,npp = Number of prop. grid points in rad, lat, long
! npgr = Refinement factor for source grid
! npcell = Number of cells in refined grid
! gorp,gotp,gopp = Grid origin of propagation grid (rad,lat,long)
! gsrp,gsrp,gspp = Grid separation of propagatipon grid (rad,lat,long)
! ueif = Use external interface file (0=no, 1=yes)
! extifile = External interface grid file
! extvfile = External velocity grid file
! dovm = Dimension of velocity model
! veltol = Velocity tolerance (minimum permitted velocity)
! 
OPEN(UNIT=10,FILE='grid3dg.in',STATUS='old')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)ni
ni=ni+1
READ(10,*)nvgt
READ(10,*)pnchd
IF(pnchd.LT.0.0)THEN
   WRITE(6,*)'Pinchout distance must be greater than'
   WRITE(6,*)'or equal to zero'
   WRITE(6,*)'Terminating program'
   STOP
ENDIF
READ(10,1)ovfile
READ(10,1)oifile
1 FORMAT(a20)
READ(10,*)rseed
READ(10,*)veltol
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)gor,glr
READ(10,*)got,glt
READ(10,*)gop,glp
READ(10,*)earthr
!
! First, set up the propagation grid
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)opfile
READ(10,*)npr,npt,npp
READ(10,*)npgr,npcell
READ(10,*)grcf
OPEN(UNIT=20,FILE=opfile,STATUS='unknown')
WRITE(20,*)npr,npt,npp
gorp=gor
gotp=glt+grcf
gopp=gop+grcf
gsrp=(gorp-glr)/REAL(npr-1)
gstp=(got-grcf-gotp)/REAL(npt-1)
gspp=(glp-grcf-gopp)/REAL(npp-1)
WRITE(20,*)gsrp,gstp,gspp
WRITE(20,*)gorp,gotp,gopp
WRITE(20,*)npgr,npcell
CLOSE(20)
got=got*pi/180.
gop=gop*pi/180.
glt=glt*pi/180.
glp=glp*pi/180.
!
! Now loop through all the layers
!
OPEN(UNIT=20,FILE=ovfile,STATUS='unknown')
WRITE(20,*)ni-1,nvgt
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
DO ips=1,nvgt
   IF(ips.EQ.2)THEN
      CLOSE(10)
      OPEN(UNIT=10,FILE='grid3dg.in',STATUS='old')
      DO i=1,29
         READ(10,*)
      ENDDO
   ENDIF
   DO ii=1,ni-1
      READ(10,*)
      READ(10,*)
      READ(10,*)
      IF(ips.EQ.1)THEN
         READ(10,*)nvr
         READ(10,*)nvt
         READ(10,*)nvp
      ELSE
         READ(10,*)idum,nvr
         READ(10,*)idum,nvt
         READ(10,*)idum,nvp
      ENDIF
      ALLOCATE(velm(nvr+2), STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with ALLOCATE: PROGRAM gridder: REAL velm'
      ENDIF
      gsr=(gor-glr)/REAL(nvr-1)
      gst=(got-glt)/REAL(nvt-1)
      gsp=(glp-gop)/REAL(nvp-1)
      ggor=gor-nvr*gsr+earthr
      ggot=got-gst*nvt
      ggop=gop-gsp
      WRITE(20,*)nvr+2,nvt+2,nvp+2
      WRITE(20,*)gsr,gst,gsp
      WRITE(20,*)ggor,ggot,ggop
      READ(10,*)umocvg
      READ(10,'(a1)')cdum
      IF(cdum.EQ.'P')THEN
        pors=0
      ELSE
        pors=1
      ENDIF
      READ(10,1)extvfile
      READ(10,*)dovm
      IF(ips.EQ.1)THEN
         READ(10,*)vtop
         READ(10,*)vbot
      ELSE
         READ(10,*)rdum,vtop
         READ(10,*)rdum,vbot
      ENDIF
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)ars
      READ(10,*)rssf
      READ(10,*)aapmc
      READ(10,*)decm
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)adch
      READ(10,*)mpc
      READ(10,*)nch
      READ(10,*)usp
      READ(10,*)
      READ(10,*)
      READ(10,*)
      READ(10,*)asp
      READ(10,*)nsp
      DO i=1,nsp
         READ(10,*)spa(i)
         READ(10,*)spr(i),spt(i),spp(i)
         spt(i)=spt(i)*pi/180.
         spp(i)=spp(i)*pi/180.
      ENDDO
!
!     If spikes are required, compute grid locations
!
      IF(asp.EQ.1)THEN
         DO i=1,nsp
            ispr(i)=(gor-spr(i))/gsr+1
            ispr(i)=nvr+1-ispr(i)
            ispt(i)=(got-spt(i))/gst+1
            ispt(i)=nvt+1-ispt(i)
            ispp(i)=(spp(i)-gop)/gsp+1
         ENDDO
      ENDIF
!
!     Determine velocity gradient (umocvg=1)
!     or velocity values from model (umocvg=0)
!
      IF(umocvg.EQ.1)THEN
         vgr=(vbot-vtop)/((nvr-1)*gsr)
      ELSE IF(dovm.EQ.1)THEN
!
!        Open file containing 1-D reference model and read in first
!        two values.
!
         OPEN(UNIT=40,FILE=extvfile,STATUS='old')
         IF(nvgt.EQ.1)THEN
            IF(pors.EQ.0)THEN
               READ(40,*)d1,v1
               READ(40,*)d2,v2
            ELSE
               READ(40,*)d1,rdum,v1
               READ(40,*)d2,rdum,v2
            ENDIF
         ELSE
            IF(ips.EQ.1)THEN
               READ(40,*)d1,v1
               READ(40,*)d2,v2
            ELSE
               READ(40,*)d1,rdum,v1
               READ(40,*)d2,rdum,v2
            ENDIF
         ENDIF
!
!        Make sure that the grid surface does not lie
!        above the 1-D model
!
         rdep=-(gor+gsr)
         IF(rdep.LT.d1)THEN
            WRITE(6,*)'Top of 3-D grid lies above radial extent of 1-D model!!'
            WRITE(6,*)'Terminating program'
            STOP
         ENDIF
!
!        Now calculate the velocity at each depth node
!
         DO i=1,nvr+2
            rdep=-(gor+gsr-(i-1)*gsr)
!
!           Ensure that rdep lies within the bounds of d1 and d2
!
            DO
               IF(rdep.GT.d2)THEN
                  d1=d2
                  v1=v2
                  IF(pors.EQ.0)THEN
                     READ(40,*)d2,v2
                  ELSE
                    READ(40,*)d2,rdum,v2
                  ENDIF
               ELSE
                  EXIT
               ENDIF
            ENDDO
!
!           Now calculate the velocity at the specified depth
!
            velm(i)=v1+(v2-v1)*(rdep-d1)/(d2-d1)
         ENDDO
         CLOSE(40)
      ELSE
!
!        Here, we open a 3-D grid file. This must have dimensions
!        equal to (nvr+2)x(nvt+2)x(nvp+2).
!
         OPEN(UNIT=40,FILE=extvfile,STATUS='old')
         IF(pors.EQ.1.AND.nvgt.EQ.1.OR.ips.EQ.2)THEN
            DO i=0,nvr+1
               DO j=0,nvt+1
                  DO k=0,nvp+1
                     READ(40,*)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
!     Superimpose checkerboard if required.
!
      IF(adch.EQ.1)THEN
         chvp1=mpc
         chvp2=mpc
         chvp=mpc
         vusp1=-1
         vusp2=-1
         vusp=-1
      ENDIF
      DO i=0,nvr+1	! radius
         IF(adch.EQ.1)THEN
            IF(MOD(i,nch).EQ.0)THEN
               chvp1=-chvp1
               IF(usp.EQ.1)THEN
                  IF(vusp1.EQ.0)THEN
                     IF(vusp1o.EQ.-1)THEN
                        vusp1=1
                     ELSE
                        vusp1=-1
                     ENDIF
                  ELSE
                     vusp1o=vusp1
                     vusp1=0
                  ENDIF
               ENDIF
            ENDIF
            chvp2=chvp1
            vusp2=1
            vusp2o=1
         ENDIF
         DO j=0,nvt+1	! latitude
            IF(adch.EQ.1)THEN
               IF(MOD(j,nch).EQ.0)THEN
                  chvp2=-chvp2
                  IF(usp.EQ.1)THEN
                     IF(vusp2.EQ.0)THEN
                        IF(vusp2o.EQ.-1)THEN
                           vusp2=1
                        ELSE
                           vusp2=-1
                        ENDIF
                     ELSE
                        vusp2o=vusp2
                        vusp2=0
                     ENDIF
                  ENDIF
               ENDIF
               chvp=chvp2
               vusp=1
               vuspo=1
            ENDIF
            DO k=0,nvp+1	! longitude
               IF(umocvg.EQ.1)THEN
                  vel=vtop+gsr*(nvr-i)*vgr
               ELSE IF(dovm.EQ.1)THEN
                  vel=velm(nvr+2-i)
               ELSE
                  READ(40,*)vel
               ENDIF
!
!              Add random structure if required.
!
               IF(ars.EQ.1)THEN
                  vel=vel+gasdev(rseed)*rssf
               ENDIF
!
!              Add checkerboard if required
!
               IF(adch.EQ.1)THEN
                  IF(MOD(k,nch).EQ.0)THEN
                     chvp=-chvp
                     IF(usp.EQ.1)THEN
                        IF(vusp.EQ.0)THEN
                           IF(vuspo.EQ.-1)THEN
                              vusp=1
                           ELSE
                              vusp=-1
                           ENDIF
                        ELSE
                           vuspo=vusp
                           vusp=0
                        ENDIF
                     ENDIF
                  ENDIF
                  vel=vel+vusp1*vusp2*vusp*chvp
               ENDIF
!
!              Apply spikes if required
!
               IF(asp.EQ.1)THEN
                  DO l=1,nsp
                     IF(k.EQ.ispp(l).OR.k.EQ.ispp(l)+1)THEN
                        IF(j.EQ.ispt(l).OR.j.EQ.ispt(l)+1)THEN
                           IF(i.EQ.ispr(l).OR.i.EQ.ispr(l)+1)THEN
                              vel=vel+spa(l)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
!
!              Write out specified covariance value
!              if required
!
               IF(vel.LT.veltol)vel=veltol
               IF(aapmc.EQ.1)THEN
                  WRITE(20,'(2f12.8)')vel,decm
               ELSE
                  WRITE(20,'(f12.8)')vel
               ENDIF
            ENDDO
            WRITE(20,'(1X)')
         ENDDO
         WRITE(20,'(1X)')
      ENDDO
      DEALLOCATE(velm, STAT=checkstat)
      IF(checkstat > 0)THEN
         WRITE(6,*)'Error with DEALLOCATE: PROGRAM gridder: velm'
      ENDIF
      IF(umocvg.NE.1.AND.dovm.NE.1)CLOSE(40)
   ENDDO
ENDDO
CLOSE(20)
!
! Now loop through all the interfaces
!
OPEN(UNIT=20,FILE=oifile,STATUS='unknown')
WRITE(20,*)ni
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)nvt
READ(10,*)nvp
gst=(got-glt)/REAL(nvt-1)
gsp=(glp-gop)/REAL(nvp-1)
WRITE(20,*)nvt+2,nvp+2
WRITE(20,*)gst,gsp
ggot=got-gst*nvt
ggop=gop-gsp
WRITE(20,*)ggot,ggop
ALLOCATE(pdep(0:nvp+1,0:nvt+1), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM gridder: REAL pdep'
ENDIF
ipsw=0
DO ii=1,ni
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,*)ueif
   READ(10,1)extifile
   READ(10,*)dnw
   READ(10,*)dne
   READ(10,*)dsw
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,*)ars
   READ(10,*)rssf
   READ(10,*)aapmc
   READ(10,*)decm
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,*)adch
   READ(10,*)mpc
   READ(10,*)nch
   READ(10,*)usp
   READ(10,*)
   READ(10,*)
   READ(10,*)
   READ(10,*)asp
   READ(10,*)nsp
   DO i=1,nsp
      READ(10,*)spa(i)
      READ(10,*)spt(i),spp(i)
      spt(i)=spt(i)*pi/180.0
      spp(i)=spp(i)*pi/180.0
   ENDDO
!
!  Determine equation of plane.
!
   xnw=1.0
   znw=1.0
   xne=REAL(nvp)
   zne=1.0
   xsw=1.0
   zsw=REAL(nvt)
   aplc=(zne-znw)*(dsw-dnw)-(dne-dnw)*(zsw-znw)
   bplc=(xne-xnw)*(dsw-dnw)-(dne-dnw)*(xsw-xnw)
   cplc=(xne-xnw)*(zsw-znw)-(zne-znw)*(xsw-xnw)
!
!  If spikes are required, compute grid locations
!
   IF(asp.EQ.1)THEN
      DO i=1,nsp
         ispt(i)=(got-spt(i))/gst+1
         ispt(i)=nvt+1-ispt(i)
         ispp(i)=(spp(i)-gop)/gsp+1
      ENDDO
   ENDIF
!
!  Superimpose checkerboard if required.
!
   IF(adch.EQ.1)THEN
      chvp1=mpc
      chvp2=mpc
      chvp=mpc
      vusp1=-1
      vusp2=-1
      vusp=-1
   ENDIF
   IF(ueif.EQ.1)THEN
      OPEN(UNIT=30,FILE=extifile,STATUS='old')
   ENDIF
   DO j=0,nvt+1
      IF(adch.EQ.1)THEN
         IF(MOD(j,nch).EQ.0)THEN
            chvp2=-chvp2
            IF(usp.EQ.1)THEN
               IF(vusp2.EQ.0)THEN
                  IF(vusp2o.EQ.-1)THEN
                     vusp2=1
                  ELSE
                     vusp2=-1
                  ENDIF
               ELSE
                  vusp2o=vusp2
                  vusp2=0
               ENDIF
            ENDIF
         ENDIF
         chvp=chvp2
         vusp=1
         vuspo=1
      ENDIF
      DO k=0,nvp+1
!
!        Determine depth to interface by linear
!        interpolation or read in from file.
!
         IF(ueif.EQ.0)THEN
            xpl=REAL(k)
            zpl=REAL(nvt+1-j)
            dep=(bplc*(zpl-znw)-aplc*(xpl-xnw))/cplc+dnw
         ELSE
            READ(30,*)dep
         ENDIF
!
!        Add random structure if required.
!
         IF(ars.EQ.1)THEN
            dep=dep+gasdev(rseed)*rssf
         ENDIF
!
!        Add checkerboard if required
!
         IF(adch.EQ.1)THEN
            IF(MOD(k,nch).EQ.0)THEN
               chvp=-chvp
               IF(usp.EQ.1)THEN
                  IF(vusp.EQ.0)THEN
                     IF(vuspo.EQ.-1)THEN
                        vusp=1
                     ELSE
                        vusp=-1
                     ENDIF
                  ELSE
                     vuspo=vusp
                     vusp=0
                  ENDIF
               ENDIF
            ENDIF
            dep=dep-vusp1*vusp2*vusp*chvp
         ENDIF
!
!        Apply spikes if required
!
         IF(asp.EQ.1)THEN
            DO l=1,nsp
               IF(k.EQ.ispp(l).OR.k.EQ.ispp(l)+1)THEN
                  IF(j.EQ.ispt(l).OR.j.EQ.ispt(l)+1)THEN
                        dep=dep+spa(l)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         dep=earthr+dep
!
!        Check that depth to interface node is not less
!        than interface above
!
         IF(ii.GT.1)THEN
            IF(dep.GE.pdep(k,j)-pnchd)THEN
               dep=pdep(k,j)-pnchd
               ipsw=1
            ENDIF
         ENDIF
!
!        Write out specified covariance value
!        if required
!
         IF(aapmc.EQ.1)THEN
            WRITE(20,'(2f12.5)')dep,decm
         ELSE
            WRITE(20,'(f12.5)')dep
         ENDIF
!
!        Record depth to all interface nodes
!
         IF(ii.NE.ni)THEN
            pdep(k,j)=dep
         ENDIF
      ENDDO
      WRITE(20,'(1X)')
   ENDDO
   WRITE(20,'(1X)')
   IF(ueif.EQ.1)THEN
      CLOSE(30)
   ENDIF
ENDDO
IF(ipsw.EQ.1)THEN
   WRITE(6,*)'Program has corrected for intersecting interfaces'
   WRITE(6,*)'by creating layer pinchouts where necessary.'
ENDIF
CLOSE(20)
CLOSE(10)
DEALLOCATE(pdep, STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with DEALLOCATE: PROGRAM gridder: pdep'
ENDIF
STOP
END PROGRAM gridder

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
