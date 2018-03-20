!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine  teleseismic_initialization(isec,source_id)

  use mod_3dfm
  implicit none
!
!     Determines ak135 traveltimes for a specified phase to
!     all bottom interface nodes  
!
!     Marthijn de Kool
!     Based on aktsurf by Nick Rawlinson 
!     Australian National University
!
!    
!     NOTE: Makes use of subroutines from the freeware
!           program ttimes.
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  type(Tintersection) :: isec
  integer             :: source_id

  integer             :: i,j,k,n,phid,iter
  integer,parameter   :: maxnp=60
  real                :: tt(maxnp),dtdd(maxnp),dtdh(maxnp),dddp(maxnp)
  real                :: evlat,evlon,evdep,ndlat,ndlon,rslat,usrc(2),nddep
  real                :: lon1,lon2,lon3,lat1,lat2,lat3,deltas1,deltas2,deltas3,delta
  real                :: deltan1,deltan2,deltan3,gt_norm,f1,f2,f3
  real                :: stt1,stt2,stt3,sdtdd1,sdtdh1,sdddp1,sdtdd2,sdtdh2,sdddp2,sdtdd3,sdtdh3,sdddp3
  real                :: ntt1,ntt2,ntt3,ndtdd1,ndtdh1,ndddp1,ndtdd2,ndtdh2,ndddp2,ndtdd3,ndtdh3,ndddp3
  real                :: deltas,cazim,bazim,edist,bazr,etcor,tph,azima
  real                :: odep,olat,olon,ddep,dlat,dlon,angi,kmpd
  real(kind=dp)       :: dodep,dolat,dolon,dddep,ddlat,ddlon,vels,deg_to_rad,deg_per_km
  real,parameter      :: pi=3.1415926535
  character(len=8)    :: phcd(maxnp),phlst(maxnp),speph,sname,loc_phase
  character(len=25)   :: modnam,sfile,efile,ofile
  logical             :: prnt(2)

!c     tt = traveltime of phase
!c     dtdd, dtdh,dddp = partial derivatives
!c     phcd = Phase id
!c     phlst = Phase list
!c     prnt = debugging print flag
!c     modnam = Name of velocity model
!c     evlat = event latitude
!c     evlon = event longitude
!c     evdep = event depth
!c     ndlat = node latitude
!c     ndlon = node longitude
!c     rslat = station co-latitude
!c     deltas = event-station angular distance
!c     cazim,bazim,azima = azimuth information
!c     m = number of phases found
!c     edist = event-station distance
!c     bazr = adjusted bazim
!c     etcor = elliptical correction for travel time
!c     phid = Phase id of specific phase
!c     tph = traveltime of specific phase
!c     maxnp = maximum number of ak135 phases
!c     speph = specific phase name for given event
!c     pi = pi
!c     er = Earth radius
!c     kmpd = Number of km per great circle degree


  if (.not.source(source_id)%is_teleseismic) stop 'teleseismic_initialization called with local source'

  if (.not. associated(isec%arrivaltime) allocate(isec%arrivaltime(isec%nnode))

  deg_to_rad = acos(-1.0_dp)/180.0_dp
  deg_per_km = earth_radius*deg_to_rad

!      Specify model type

  modnam='ak135'

!     Open travel time tables

  prnt(1)=.false.
  prnt(2)=.false.
  phlst(2)='  '
  call tabin(10,modnam)

!     get source location and phase type

  evlat = source(source_id)%lat/deg_to_rad
  evlon = source(source_id)%long/deg_to_rad
  evdep = earth_radius - source(source_id)%r
  speph = source(source_id)%teleseismic_phase

  loc_phase = speph

  phlst(1)=speph
  call brnset(1,phlst,prnt)
  rslat = (90.-evlat)*0.017453292
  call ellref(rslat)




!  Now loop through all the bottom interface points

  do n=1,isec%nnode

 ! convert position of the node to input taken by ttimes
     
     ndlat = isec%lat(n)/deg_to_rad
     ndlon = isec%long(n)/deg_to_rad
     nddep = earth_radius - isec%r(n)

!  This is a bit of a hack fix for the case when event and grid node longitude are equal. 
!  For some reason the ellipticity corrections are wrong for this case

     if(ndlon.eq.evlon)then
        ndlon=ndlon+isec%dlong0/(50.0*deg_to_rad)
     endif


! get angular distance between source and node
 
     call ydist(evlat,evlon,ndlat,ndlon,deltas,cazim,bazim,azima)


 ! initial guesses for position of points bracketting the solution

! the travel time and arrival time gradient from the source to the two initial guess points

     delta  = nddep*deg_per_km
     deltas1 = deltas

     phlst(1)=speph
     call brnset(1,phlst,prnt)
     call depset(evdep,usrc)

     deltas2 = deltas + delta    

     call trtm(deltas1,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
     call extract_phase
     stt1=tt(phid) ; sdtdd1=dtdd(phid) ; sdtdh1=dtdh(phid) ; sdddp1=dddp(phid)


     call trtm(deltas2,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
     call extract_phase
     stt2=tt(phid) ; sdtdd2=dtdd(phid) ; sdtdh2=dtdh(phid) ; sdddp2=dddp(phid)


 ! the travel time and arrival time gradient from the node to the two initial guess points

   
     phlst(1)=loc_phase
     call brnset(1,phlst,prnt)
     call depset(nddep,usrc)

     deltan1 = 0.0

     call trtm(deltan1,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
     call extract_phase
     ntt1=tt(phid) ; ndtdd1=dtdd(phid) ; ndtdh1=dtdh(phid) ; ndddp1=dddp(phid)

     gt_norm = ndtdh1  ! store the magnitude of the surface time gradient for later use

     deltan2=delta

     call trtm(deltan2,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
     call extract_phase
     ntt2=tt(phid) ; ndtdd2=dtdd(phid) ; ndtdh2=dtdh(phid) ; ndddp2=dddp(phid)

     f1 = ndtdh1 - sdtdh1
     f2 = ndtdh2 - sdtdh2

     if (f1 < 0.0) stop 'error 1 in tele_init'


 ! if the root is not bracketted in the first interval, try extending it

     iter = 1

     do while (f1*f2 > 0.0)

        iter = iter + 1
        if (iter > 3) stop 'tele_init: cannot find bracketted solution on surface'

        delta = delta*2.0
        deltas2 = deltas + delta  
        deltan2 = delta

        phlst(1)=speph
        call brnset(1,phlst,prnt)
        call depset(evdep,usrc)
        call trtm(deltas2,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
        call extract_phase
        stt2=tt(phid) ; sdtdd2=dtdd(phid) ; sdtdh2=dtdh(phid) ; sdddp2=dddp(phid)

        phlst(1)=loc_phase
        call brnset(1,phlst,prnt)
        call depset(nddep,usrc)
        call trtm(deltan2,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
        call extract_phase
        ntt2=tt(phid) ; ndtdd2=dtdd(phid) ; ndtdh2=dtdh(phid) ; ndddp2=dddp(phid)

        f2 = ndtdh2 - sdtdh2

     end do


!  start the root finding iteration 

     iter = 1
     delta_min=isec%dlong0/(50.0*deg_to_rad)
     f3 = gt_norm

     do while (abs(f3) > 1.e-5*gt_norm .and. delta > delta_min) 


        iter = iter + 1
        if (iter > 20) stop 'tele_init: cannot converge to solution'


        deltan3 = deltan1 +((f1/(f1-f2))*(deltan2-deltan1))
        deltas3 = deltas + deltan3  

        phlst(1)=speph
        call brnset(1,phlst,prnt)
        call depset(evdep,usrc)
        call trtm(deltas3,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
        call extract_phase
        stt3=tt(phid) ; sdtdd3=dtdd(phid) ; sdtdh3=dtdh(phid) ; sdddp3=dddp(phid)

        phlst(1)=loc_phase
        call brnset(1,phlst,prnt)
        call depset(nddep,usrc)
        call trtm(deltan3,maxnp,m,tt,dtdd,dtdh,dddp,phcd)
        call extract_phase
        ntt3=tt(phid) ; ndtdd3=dtdd(phid) ; ndtdh3=dtdh(phid) ; ndddp3=dddp(phid)

        f3 = ndtdh3 - sdtdh3

        if (f3 > 0.0) then
           f1 = f3 ; deltan1=deltan3 ; deltas1=deltas3 ; ntt1 = ntt3
        else
           f2 = f3 ; deltan2=deltan3 ; deltas2=deltas3 ; ntt2 = ntt3
        endif

        delta =abs(deltan2-deltan1)

     end do

!    Apply elliptical corrections
!

     bazr=bazim*0.017453292

     phlst(1)=speph
     call brnset(1,phlst,prnt)
     call depset(evdep,usrc)
     rslat = (90.-evlat)*0.017453292
     call ellref(rslat)
     edist=deltas3*0.017453292
     call ellcor(edist, bazr, evdep, speph, etcor)
     stt3=stt3+etcor

     phlst(1)=loc_phase
     call brnset(1,phlst,prnt)
     call depset(nddep,usrc)
     rslat = (90.-ndlat)*0.017453292
     call ellref(rslat)
     edist=deltan3*0.017453292
     call ellcor(edist, bazr, nddep, loc_phase, etcor)
     ntt3=ntt3+etcor

! arrival time at the bottom interface node is the difference in arrival times 

     isec%arrivaltime(n) = stt3 - ntt3

  enddo

  return

  contains

    subroutine extract_phase

!   Extract specified phase.

        k=1
        phid=0
        do while(k.le.m.and.phid.eq.0)
           if(phcd(k).eq.speph)then
              phid=k
           else
              k=k+1
           endif
        enddo

        if(phid.eq.0)then
           print *,'Requested phase ',speph, 'does not'
           print *,'exist for delta = ',deltas
           print *,'Source number = ',source%id
           print *,'Terminating program'
           stop
        endif

        return

    end subroutine extract_phase

end subroutine teleseismic_initialization
