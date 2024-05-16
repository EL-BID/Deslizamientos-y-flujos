MODULE SPH_FEM_topo_2019
 
! Changes in MODULE SPH_FEM_topo_2014 March 21 2014-----------------------------
 
!(a) New variables / uses ...........

!  ic_topo_props = 0 nothing	      
!                = 1 bottom friction (obs)
!  	             = 2 erosion         (obs)
!  	             = 3 both	     (obs)
!                = 10 New system. Keep compatibility as far as possible
!
!  ic_basal_friction = 1 if present, zero otherwise 
!  ic_basal_erosion  = 1 if present, zero otherwise
!  ic_basal_PwpBCs   = 1 if present, zero otherwise             
!
!  ntopo_zones   is the number of topo zones. Can be delimited in X, Y,..
!  ntopo_props   10 topo props stored in the new version
!  basal_type (ntopo_zones)
!                 = 1 friction  2 erosion 3 Basal Pwps  4 Obstacle field
!  basal_law (ntopo_zones) type of law
!              Erosion:  (1) Hungr
!                        (2) Tom
!                        (3) MP new
!              obstacle type (1) 1D (2) 2D
!
!  basal_props (ntopo_props,ntopo_zones) Assigned properties
!        (1)   Type of geomm... 1 Zmin Zmax
!                           ... 2 Xmin Xmax
!                           ... 3 Ymin Ymax
!                           ...	4 Xmin Ymin Xmax Ymax (vertex rect.)
!                           ...	5 same than 4 and Theta (incl Xmin Xmax axis)
!                           ... 101 read pts from file
!                           ... 102 uses finite elements to define zone
!                           ... 200 obstacle fiel, read in get_obstacles(ic_obs)
!        (2-6)  Geomm parameters (see above) 
!        (7-11) Material parameters, depend on phenomenon and law
!               erosion + Hungr: Er Max_erodible depth
!               friction       : delta basal
!               Pwp            : 

! (b) New subroutine Which_Topo_Zone
!
!       Given (xx,yy) of a point, 
!             with basal_type btype (1 friction, 2 erosion, 3 Pwp...)
!       Obtains the topo zone at which it belongs index_topo_zone
!              which is -999 if the points doesn't belong to any

!     --------------------------------------------------------------------------------
!
!   Added  code 10 and 11...11 includes some options
!          code 5  cone for Stromboli
!   Added terrain props 2 July 2007
!         (1) Zmin  (2) Zmax  (3) Tan delta (4) max erosion depth
!   Added  Subroutine Update_topol 19 July 2007
!           from topol(1,:) updates Zx Zy Zxx,..
!   Added  Subroutine Get_DTM_FE   uses dtm in fe format (Malpasset)
!          Sept 2008.

! Modifications to include erosion in the code

! drum added Ictop is 6 Data is read as:
!
!        read  (topdat_file,1) text
!        read  (topdat_file,*) x0, x1, R, d1
!
!.................................................
! 
!     topol ( 1,ipoig) = Z
!     topol ( 2,ipoig) = Zx
!     topol ( 3,ipoig) = Zy
!     topol ( 4,ipoig) = 1 + Zxx^2  
!     topol ( 5,ipoig) = Zx.Zy  
!     topol ( 6,ipoig) = 1 + Zyy^2
!     topol ( 7,ipoig) = Zxx / rt 
!     topol ( 8,ipoig) = Zxy / rt 
!     topol ( 9,ipoig) = Zyy / rt 
!     topol (10,ipoig) = rt = (1 + Zx^2 + Zyy^2)**0.5   Zxx/(1+Zx^2+Zy^2)**0.5
!     topol (11,ipoig) = Zxx
!     topol (12,ipoig) = Zxy
!     topol (13,ipoig) = Zyy
!     topol (14,ipoig) = eroded topo (updated)
!     topol (15,ipoig) = available erosion depth
!     topol (16,ipoig) = Geo zone to which ipoig belongs (MP April 2022)
!.................................................
! 

USE SPH_FEM_variable_types_2019
USE SPH_FEM_Driver_time_vars_2019    ! new CHGD MP March 28 2020

implicit none
private

integer(ink), public::  npoigx, npoigy, npoig, ntopo
real   (irk), public::  Xming, Xmaxg, Yming, Ymaxg, deltxg, deltyg, xshiftg, yshiftg
real   (irk), public::  Xming_org, Xmaxg_org, Yming_org, Ymaxg_org ! Saeidmt 12Nov2023

integer(ink), public::  ic_topo_props 
integer(ink), public::  ic_basal_friction, ic_basal_erosion, ic_basal_PwpBCs
integer(ink), public::  mass_zones           ! Saeidmt 12Nov2023
integer(ink), public::  ic_basal_flux        !  MP August 2021 CHGD 
integer(ink), public::  ntopo_props, ntopo_zones 
integer(ink), public, allocatable::  basal_type (:), basal_law (:)     
real   (irk), public, allocatable::  topo_props (:,:) ! obsolete, substituted by basal_props  
real   (irk), public, allocatable::  basal_props(:,:)  
real   (irk), public, allocatable::  auxd1(:,:)       ! Saeidmt 12Nov2023

integer(ink),public, allocatable:: poig_in_model (:)
real   (irk),public, allocatable:: coorg  (:,:) 
real   (irk),public, allocatable:: topol  (:,:), topol0(:,:) 

integer(ink)  Topchk_file, Topdat_file, Topgid_msh, Topgid_res, Topgid_msh3D, Topgid_res3D, len1 

integer(ink) ic_FDiff  ! mp 21st Jan 2019 set to 0 default  1 compute Zxx etc using FDiffs formulas

TYPE TZdata                                           ! MP CH 2019 06 02 June
   integer (ink)                        :: nplist     ! nr of topo nodes in list
   integer (ink)                        :: npmin      ! lower bound for searching
   integer (ink)                        :: npmax      ! upper bound 
   integer (ink)                        :: if_list    ! 0 , 1 not allocated  2 allocated   !  MP 2 march 2022
   integer (ink), dimension (:), pointer:: NodeList   ! list of nodes belonging to TopoZone       
END TYPE TZdata

TYPE (TZdata),  allocatable:: TopoZone_Info (:)

public:: Topo_Main, Get_Z_topo, Update_Topol, Which_Topo_Zone,Get_Z_plus_hs 
public:: Get_Topo_restart

CONTAINS

!-------------------------------------------------------------------

       SUBROUTINE Topo_main

!-------------------------------------------------------------------

 
!      ------  Gets Topographic data for SW 

integer(ink) ictop, ipoig, jpoig, idimn, inodg, jnodg
integer(ink) i, ip   , ipoigx, ipoigy, nptsz, iptsz
integer(ink) neleg, ieleg, ieleg2, ielgx, ielgy, itopo_zones

integer(ink) ic_obst 
integer(ink) ic_chk_dam   !  MP June 2023 for check dams

real   (irk) xg, yg, Zconst, T1, R, Theta, T3, xx, zz, ThetaRd, beta
real   (irk) x0,d0,x1,d1,x2,d2, xn0, xn1, x3, h1, h0, xn2, f0, f1
real   (irk) xeros_min, xeros_max, yeros_min, yeros_max, a0z, Lz, Zxmax

real   (irk), allocatable:: xptsZ   (:)     ! -- x's for Z segments
real   (irk), allocatable:: zptsZ   (:)     ! -- z's for Z segments
real   (irk), allocatable:: geome   (:,:)   ! -- derivatives of shape functions 7*neleg 

integer(ink), allocatable:: intmag  (:,:)   ! -- triangle mesh /1  (3*neleg)

character(16) TOPO_problem_name       ! name of problem being solved.
character(60) text                    ! general purpose char string

Topchk_file   = 61
topDat_file   = 71
Topgid_msh    = 81
Topgid_res    = 82
! Topgid_msh3D  = 83
! Topgid_res3D  = 84

ntopo       = 16     ! CHGD MP 13 April 2022 for Geo Zones used in trigger
ntopo_props = 4

!print *, ' ------------------------------------ ' 

!print *, ' '
!print *, '   TOPO module '
!print *, ' '
!print *,'Input the problem name (TOPO mesh) ?'
!print *, '  '

if (if_Read_TopoDat.eq.0) then 
    print *, ' '
    print *,'Input the problem name (TOPO mesh) ?'
    read  *, TOPO_problem_name
    if_Read_TopoDat = if_Read_TopoDat + 1
else
!    print *,'Already available (TOPO mesh) ?'
    RETURN
endif

print *, ' '
print *, '  '
print *, ' '

len1 = len_trim(TOPO_problem_name)

open (  Topdat_file,  file = TOPO_problem_name(1:len1)//'.top'         )
open (  Topchk_file,  file = TOPO_problem_name(1:len1)//'.top.chk'     )
open (   Topgid_msh,  file = TOPO_problem_name(1:len1)//'.top.post.msh')
open (   Topgid_res,  file = TOPO_problem_name(1:len1)//'.top.post.res')

! open (   Topgid_msh,  file = TOPO_problem_name(1:len1)//'.flavia.dat'  )
! open (   Topgid_res,  file = TOPO_problem_name(1:len1)//'.flavia.res'  )
! open ( Topgid_msh3D,  file = TOPO_problem_name(1:len1)//'3d.flavia.bon')
! open ( Topgid_res3D,  file = TOPO_problem_name(1:len1)//'3d.flavia.res')


!      ------  First thing is to know wheteher topo is read from a DTM or generated here. 
!              0 means DTM  npoigx, npoigy and a set of npoigx*npoigy points (x,y,z) are read.  

read (topDat_file,1) text
write(Topchk_file,1) text
read (topDat_file,*) ictop
write(Topchk_file,*) ictop

1 format(a60)

ic_Fdiff = 0             ! ** MP 21st Jan 2019  ic_fdiff gt 100 means compute Zxx using FDiffs
if (ictop.GE.100) then
    ic_Fdiff = 1
    ictop    = ictop - 100
endif

xshiftg = 0.0 ; yshiftg = 0.0

IF (ictop.eq.0) THEN         		!  ------  Read from file  2D meshes only!!
    read (topdat_file,1) text
    write(Topchk_file,1) text
    read (topdat_file,*) npoigx, npoigy
    write(Topchk_file,*) npoigx, npoigy
    read (topdat_file,1) text
    write(Topchk_file,1) text
    npoig    = npoigx*npoigy
    
    if (.NOT.allocated (coorg) ) allocate (coorg  (2,npoig))     ; coorg = 0.0
    if (.NOT.allocated (topol) ) allocate (topol  (ntopo,npoig)) ; topol = 0.0
    if (.NOT.allocated (poig_in_model) ) allocate (poig_in_model(npoig)) ; poig_in_model = 0

    coorg    = 0.0
    topol    = 0.0  
    poig_in_model = 1
    xming    =  1.e9
    yming    =  1.e9
    xmaxg    = -1.e9
    ymaxg    = -1.e9
    DO ip    =  1,npoig
       read (topdat_file,*) ipoig, (coorg(idimn,ipoig), idimn=1,2), topol(1,ipoig)
       write(Topchk_file,*) ipoig, (coorg(idimn,ipoig), idimn=1,2), topol(1,ipoig)
       if (coorg(1,ipoig).le.xming) xming = coorg(1,ipoig)
       if (coorg(1,ipoig).ge.xmaxg) xmaxg = coorg(1,ipoig)
       if (coorg(2,ipoig).le.yming) yming = coorg(2,ipoig)
       if (coorg(2,ipoig).ge.ymaxg) ymaxg = coorg(2,ipoig)
    ENDDO
    deltxg = (xmaxg-xming)/(npoigx-1)
    deltyg = (ymaxg-yming)/(npoigy-1)

ELSEIF (ictop.ge.10.and.ictop.le.16) THEN   !  12=10 and 13=11 plus reading a rect. and its H
                                            !  15 is to select a window (heavy dtms)
    call Get_DTM (ictop)

ELSEIF (ictop.eq.20.or.ictop.eq.21) THEN

    call Get_DTM_FEM (ictop)
    
ELSE            	!  ------  All other cases. We first generate X,Y 

    read (topdat_file,1) text
    write(Topchk_file,1) text
    read (topdat_file,*) Xming, Xmaxg, Yming, Ymaxg
    write(TopChk_file,*) Xming, Xmaxg, Yming, Ymaxg  
    read (topdat_file,1) text
    write(Topchk_file,1) text
    read (TopDat_file,*) npoigx, npoigy
    write(TopChk_file,*) npoigx, npoigy
    npoig  = npoigx*npoigy
       
    if (.NOT.allocated (coorg) ) allocate (coorg  (2,npoig))     ; coorg = 0.0
    if (.NOT.allocated (topol) ) allocate (topol  (ntopo,npoig)) ; topol = 0.0
    if (.NOT.allocated (poig_in_model) ) allocate (poig_in_model(npoig)) ; poig_in_model = 0.0

    coorg    = 0.0
    topol    = 0.0 
    poig_in_model = 1

    deltxg = (Xmaxg-Xming)/(npoigx-1)
    deltyg = (Ymaxg-Yming)/(npoigy-1)
    ipoig  = 0
    DO ipoigy = 1,npoigy 
       yg     = yming + (ipoigy-1)*deltyg
       Do ipoigx = 1,npoigx 
          ipoig  = ipoig + 1
          xg     = xming + (ipoigx-1)*deltxg
          coorg(1,ipoig) = xg
          coorg(2,ipoig) = yg
       enddo
    ENDDO

    if (ictop.eq.1) then  	! --- Bottom elev. is constant
        read (topdat_file,1) text
        write(Topchk_file,1) text
        read (topdat_file,*) Zconst
        write(Topchk_file,*) Zconst
        do  ipoig = 1,npoig
            topol(1,ipoig) = Zconst
            write(topchk_file,*) ipoig,(coorg(idimn,ipoig),idimn=1,2)
        enddo 

    elseif (ictop.eq.2) then 		!--- slope with cylindr. transition MP 04-2003  
	 
	read (topdat_file,1) text
	write(Topchk_file,1) text
	read (topdat_file,*) T1, R, Theta, T3 ! Long trm 1; rad, lope angle, long horiz
	write(Topchk_file,*) T1, R, Theta, T3
	thetard = theta*3.1416/180
	x0 = 0.0
	x1 = x0 + T1
	x2 = x1 + R*sin(thetard)
	x3 = x2 + T3
	H1 = R*(1.-cos(thetard))
	H0 = T1*tan(thetard)
	do ipoig = 1,npoig
	   xx  = coorg(1,ipoig)
	   if (xx.ge.x0.and.xx.lt.x1) then
	       zz   = H1+(x1-xx)*tan(thetard)
           elseif (xx.ge.x1.and.xx.lt.x2) then
              beta  = asin((x2-xx)/R)
              zz    = R*(1-cos(beta))
           elseif (xx.ge.x2.and.xx.le.x3) then 
              zz    = 0.0
           else
              write(*,*) 'cyl trans. ip',ipoig,'outofrange'
           endif
           topol(1,ipoig) = zz
           write(Topchk_file,*) ipoig,(coorg(idimn,ipoig),idimn=1,2),topol(1,ipoig)
	enddo 	
    elseif (ictop.eq.3) then 	  ! --- Bilinear law in X		

        read (topdat_file,1) text
        write(Topchk_file,1) text
        read (topdat_file,*) x0,d0,x1,d1,x2,d2
        write(Topchk_file,*) x0,d0,x1,d1,x2,d2
        do ipoig = 1,npoig
           xx = coorg(1,ipoig)
           if (xx.ge.x0.and.xx.le.x1) then
               xn0 = (xx-x1)/(x0-x1)
               xn1 = 1.-xn0
               zz  = xn0*d0+xn1*d1
           elseif (xx.ge.x1.and.xx.le.x2) then
               xn1 = (xx-x2)/(x1-x2)
               xn2 = 1.-xn1
               zz  = xn1*d1+xn2*d2
           endif
           topol(1,ipoig) = zz
        enddo 
        
    elseif (ictop.eq.4) then  ! --- Multilinear law in X
       
        read  (topdat_file,1) text
        write (Topchk_file,1) text
        read  (topdat_file,*) nptsZ
        write (Topchk_file,*) nptsz
        if (.NOT.allocated (xptsZ) ) allocate ( xptsZ(nptsZ) )
        if (.NOT.allocated (ZptsZ) ) allocate ( ZptsZ(nptsZ) )
        xptsZ = 0.0
        ZptsZ = 0.0
        read  (topdat_file,1) text
        write (Topchk_file,1) text
        read  (topdat_file,*) (xptsZ(i),i=1,nptsZ)
        write (Topchk_file,*) (xptsZ(i),i=1,nptsZ)
        read  (topdat_file,*) (ZptsZ(i),i=1,nptsZ)
        write (Topchk_file,*) (ZptsZ(i),i=1,nptsZ)

        do ipoig = 1,npoig
           xx = coorg(1,ipoig)
           do iptsZ = 1,nptsZ-1
              x0 = xptsZ(iptsZ)
              x1 = xptsZ(iptsZ+1)
              if (xx.ge.x0.and.xx.le.x1) then
                 f0 = (xx-x1)/(x0-x1)
                 f1 = 1.-f0
                 zz = f0*ZptsZ(iptsZ)+f1*ZptsZ(iptsZ+1)
                 topol(1,ipoig) = zz
                 exit
              endif
           enddo
        enddo 
        deallocate(xptsZ)
        deallocate(ZptsZ)
               
    elseif (ictop.eq.5) then  ! --- Cone: give: r at top, R at bottom, Height
    
        call get_cone
        
    elseif (ictop.eq.6) then  ! --- DRUM: give: Xleft Xright R dx  *** 20 oct 2013
    
        read  (topdat_file,1) text
        write (Topchk_file,1) text
        read  (topdat_file,*) x0, x1, R, d1
        write (Topchk_file,*) x0, x1, R, d1
        DO ipoig = 1, npoig
           xx    = coorg(1,ipoig)
           Theta = asin (xx/R)
           zz    = R*(1-cos(Theta))
           topol(1,ipoig) = zz
        enddo 
         
    elseif (ictop.eq.7) then  ! --- MP test trigger may 2022
    
        read  (topdat_file,1) text
        write (Topchk_file,1) text
        read  (topdat_file,*) Zxmax, Lz 
        write (Topchk_file,*) Zxmax, Lz 
        a0z   = Zxmax / (2.*Lz)        !a0z = Zxmax 
        ipoig = 0
        DO ipoigx = 1, npoigx
        DO ipoigy = 1, npoigy
           ipoig = ipoig + 1
           xx    = coorg(1,ipoig)
           zz    = a0z *xx*xx      !       zz    = a0z *xx 
           topol(1,ipoig) = zz
        enddo 
        enddo
         
    elseif (ictop.eq.-7) then  ! --- MP test trigger may 2022
    
        read  (topdat_file,1) text
        write (Topchk_file,1) text
        read  (topdat_file,*) a0z, Lz
        write (Topchk_file,*) a0z, Lz
        ipoig = 0
        DO ipoigx = 1, npoigx
        DO ipoigy = 1, npoigy
           ipoig = ipoig + 1
           xx    = coorg(1,ipoig)
           zz    = 0.5*a0z*(1. + cos (coorg(1,ipoig)*3.141592/Lz))
           topol(1,ipoig) = zz
        enddo 
        enddo
    else
        write(*,*) ' Error: topography not implemented'
        stop
    endif
    
ENDIF

!      ------  Read properties associated to terrain, such as basal friction and max erodible depth
!
!     ntopo_props = 4,6 or 7
!               1 ... Zmin
!               2 ... Zmax 
!               3 ... Delta (friction with bottom id diff from Fi)
!               4 ... Hungr's constant which modifies c4
!               5 ... max erosion depth
!               6 ... BC type
!               7 ... BC value

read  (topdat_file,1) text
write (Topchk_file,1) text
read  (topdat_file,*) ic_topo_props
write (Topchk_file,*) ic_topo_props

ic_chk_dam = 0                  ! default MP June 2023
if (ic_topo_props.eq.11) then   !  11  means check_dams present
   ic_chk_dam = 1
   ic_topo_props = 0
elseif (ic_topo_props.eq.12) then   !  12  means check_dams present
   ic_chk_dam = 1                   ! and other goodies
   ic_topo_props = 10
endif   
 
if (ic_topo_props.eq.1) then
    ntopo_props = 4		! keep compatibility with old versions
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    write (Topchk_file,*) ' take care, 4 constants'
    read  (topdat_file,*) ntopo_zones
    write (Topchk_file,*) ntopo_zones
    if (.NOT.allocated (topo_props) ) allocate (topo_props(ntopo_props,ntopo_zones))
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    write(*,*) ' ntopo_props   ntopo_zones '
    write(*,*) ntopo_props, ntopo_zones
    DO itopo_zones = 1, ntopo_zones
       read  (topdat_file,*) (topo_props(i,itopo_zones),i=1,ntopo_props)
       write (Topchk_file,*) (topo_props(i,itopo_zones),i=1,ntopo_props)  
    enddo
elseif (ic_topo_props.gt.1.and.ic_topo_props.lt.4) then
    ntopo_props = 6
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    read  (topdat_file,*) ntopo_zones
    write (Topchk_file,*) ntopo_zones
    write (Topchk_file,*) ' take care, read 6 constants and not 4 '
    if (.NOT.allocated (topo_props) ) allocate (topo_props(ntopo_props,ntopo_zones))
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    write(*,*) ' ntopo_props   ntopo_zones '
    write(*,*) ntopo_props, ntopo_zones
    DO itopo_zones = 1, ntopo_zones
       read  (topdat_file,*) (topo_props(i,itopo_zones),i=1,ntopo_props)
       write (Topchk_file,*) (topo_props(i,itopo_zones),i=1,ntopo_props)  
    enddo
elseif (ic_topo_props.eq.4) then 
    ntopo_props = 7
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    read  (topdat_file,*) ntopo_zones
    write (Topchk_file,*) ntopo_zones
    write (Topchk_file,*) ' take care, read 7 constants and not 4 nor 6'
    if (.NOT.allocated (topo_props) ) allocate (topo_props(ntopo_props,ntopo_zones))
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    write(*,*) ' ntopo_props   ntopo_zones '
    write(*,*) ntopo_props, ntopo_zones
    DO itopo_zones = 1, ntopo_zones
       read  (topdat_file,*) (topo_props(i,itopo_zones),i=1,ntopo_props)
       write (Topchk_file,*) (topo_props(i,itopo_zones),i=1,ntopo_props)  
    enddo
elseif (ic_topo_props.eq.10) then
    ntopo_props = 16
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    read  (topdat_file,*) ntopo_zones
    write (Topchk_file,*) ntopo_zones
    ic_basal_friction = 0
    ic_basal_erosion  = 0
    ic_basal_PwpBCs   = 0
    ic_basal_flux     = 0
    mass_zones = 0                ! Saeidmt 12Nov2023
    if (.NOT.allocated (TopoZone_Info) ) then
       allocate (TopoZone_Info (ntopo_zones))
       topoZone_info (:)%if_list = 1
       allocate (basal_type    (ntopo_zones))             ; basal_type  = 0
       allocate (basal_law     (ntopo_zones))             ; basal_law   = 0    
       allocate (basal_props (ntopo_props, ntopo_zones))  ; basal_props = 0
    endif
    DO itopo_zones = 1, ntopo_zones       
       read  (topdat_file,1) text
       write (Topchk_file,1) text
       read  (topdat_file,*)  basal_type(itopo_zones) 
       write (Topchk_file,*)  basal_type(itopo_zones) 
           i = basal_type(itopo_zones) + 0.01
           if (i.eq.1)  ic_basal_friction = 1 
           if (i.eq.2)  ic_basal_erosion  = 1 
           if (i.eq.3)  ic_basal_PwpBCs   = 1 
           if (i.eq.4)  then                     ! use 4 for basl grids desaturation
               ic_basal_flux           = 1       ! flag will be ic_basal_flux = 1
               ic_basal_PwpBCs         = 1       !  MP CHGD August 2021
               basal_type(itopo_zones) = 3
           endif
           if (i.eq.5)  mass_zones   = 1             ! Saeidmt 12Nov2023
           if (i.eq.55)  mass_zones   = 2             ! Saeidmt 12Nov2023
       read  (topdat_file,1) text
       write (Topchk_file,1) text 
       read  (topdat_file,*)  basal_law (itopo_zones) 
       write (Topchk_file,*)  basal_law (itopo_zones)      
       read  (topdat_file,1) text
       write (Topchk_file,1) text  
       read  (topdat_file,*) (basal_props(i,itopo_zones),i=1,6)
       write (Topchk_file,*) (basal_props(i,itopo_zones),i=1,6)     
       read  (topdat_file,1) text
       write (Topchk_file,1) text  
       read  (topdat_file,*) (basal_props(i,itopo_zones),i=7,ntopo_props)
       write (Topchk_file,*) (basal_props(i,itopo_zones),i=7,ntopo_props) 
       i = basal_props(1,itopo_zones) + 0.01
       if (i.eq.101) then				! changed Jan 2016 
          call Get_Eros_Pts ( TOPO_problem_name, itopo_zones )
       elseif (i.eq.102) then
          call Get_Structure_Pts ( TOPO_problem_name, itopo_zones )
       elseif (i.eq.200) then
          ic_obst = basal_law (itopo_zones)
          call get_obstacles (ic_obst, ictop)
       endif   
    enddo

elseif (ic_topo_props.eq.100) then  ! we will define Geo zones 
    call Get_Geo_Zones              ! to use in Trigger MP April 2022
elseif (ic_topo_props.eq.101.or.ic_topo_props.eq.-101) then  ! we will define Geo zones 
    call Get_Geo_Zones_NX            ! Non convex zones MP 22 June 2022  !   MP & Saeid April 2024
endif 
 
if (ic_chk_dam.eq.1) then      ! MP include check dams
    call get_check_dams
endif

!      ------   Now we will obtain the conectivities of the triangles /1 background mesh 

IF (ic_FDiff.eq.0) THEN
    neleg = 2*(npoigx-1)*(npoigy-1)
    if (.NOT.allocated (intmag) ) then 
       allocate ( intmag(3,neleg) )
       allocate ( geome (7,neleg) )
    endif
    geome = 0
    do ipoigy = 1,npoigy-1
       do ipoigx = 1, npoigx-1
          ieleg  =   (ipoigy-1)*(npoigx-1) + ipoigx 
          ieleg2 = 2*ieleg
          ipoig  = (ipoigy-1)*npoigx + ipoigx
          intmag(1,ieleg2-1) = ipoig   
          intmag(2,ieleg2-1) = ipoig  + 1
          intmag(3,ieleg2-1) = ipoig  + 1 + npoigx
          intmag(1,ieleg2)   = ipoig   
          intmag(2,ieleg2)   = ipoig  + 1 + npoigx
          intmag(3,ieleg2)   = ipoig      + npoigx
       enddo
    enddo

!      ------  Next, we obtain GEOME and TOPOL

    call GetGeom2D (neleg,intmag,geome)
    call Get_Topol (neleg,intmag,geome)
    call ouplot (TOPO_problem_name, neleg,intmag) !      ------  write results
    !call ouplot_deposit (TOPO_problem_name)      ! Saeidmt 12Nov2023
   
    deallocate (intmag)
    deallocate (geome)

ELSE

    call FDiff_derivs   !   compute fast Zxx etc by Finite diffs
    call ouplot_FDiffs (TOPO_problem_name)
    !call ouplot_deposit (TOPO_problem_name)      ! Saeidmt 12Nov2023

ENDIF

if (.NOT.allocated(topol0)) allocate ( topol0(ntopo,npoig) )
topol0 = topol                         !     ------  Store initial topo in topol0

write(Topchk_file,*) ' ***** shift is -----------' 
write(Topchk_file,*) xshiftg, yshiftg

close ( Topdat_file)
close ( Topchk_file)
! close (  Topgid_msh)
! close (  Topgid_res)
! close (Topgid_msh3D)
! close (Topgid_res3D)

END SUBROUTINE Topo_main


!--------------------------------------------------------

  subroutine get_check_dams
  
!--------------------------------------------------------  

!      reads info concerning chk dams  MP June 2023. See input reservoirs
!             
implicit none
character(60) text
integer(ink)  n_chk_dams, i_chk_dams
integer(ink)  inx, iny, inx0, iny0, inx1, iny1, ipoig, ipoigx , ipoigy
integer(ink)  ipoigx0, ipoigy0, ipoigx1, ipoigy1 , idumm 
real   (irk)  xmin_chk_dam, xmax_chk_dam, ymin_chk_dam, ymax_chk_dam
real   (irk)  Z_swl, x0_dam, y0_dam, xf_dam, yf_dam, tandelt
real   (irk)  xm, ym                                    ! center of chkdam, max slope at it
real   (irk)  RR, nx_dam, ny_dam, nabs, xx, yy, Zg      ! dist from (x,y) to dam along normal axis
real   (irk)  Zp, Zpx, Zpy, Zpgrad , deltZg             ! dz/dx, dZ/dy   

1 format(a60)

read  (topdat_file,1) text
write (Topchk_file,1) text 
read  (topdat_file,*) n_chk_dams
write (Topchk_file,*) n_chk_dams

DO i_chk_dams = 1, n_chk_dams

   read  (topdat_file,1) text
   write (Topchk_file,1) text 
   read  (topdat_file,*) xmin_chk_dam, xmax_chk_dam, ymin_chk_dam, ymax_chk_dam
   write (Topchk_file,*) xmin_chk_dam, xmax_chk_dam, ymin_chk_dam, ymax_chk_dam 

   read  (topdat_file,1) text
   write (Topchk_file,1) text 
   read  (topdat_file,*) Z_swl, x0_dam, y0_dam, xf_dam, yf_dam, tandelt
   write (Topchk_file,*) Z_swl, x0_dam, y0_dam, xf_dam, yf_dam, tandelt

   xm = ( xmin_chk_dam + xmax_chk_dam )/ 2    ! center of chkdam line
   ym = ( ymin_chk_dam + ymax_chk_dam )/ 2
   call Get_Z_Topo (idumm, xm, ym, Zp, Zpx, Zpy)
   nx_dam = -(yf_dam-y0_dam)                    ! this is the normal pointing inwards reservoir
   ny_dam =   xf_dam-x0_dam
   nabs   =  (nx_dam*nx_dam +  ny_dam*ny_dam)**0.5
   nx_dam = nx_dam /nabs
   ny_dam = ny_dam /nabs

   ipoigx0 = (xmin_chk_dam - xming)/deltxg + 1   !   ---  limits in geo grid for dam defined region
   ipoigy0 = (ymin_chk_dam - yming)/deltyg + 1  
   ipoigx1 = (xmax_chk_dam - xming)/deltxg + 1
   ipoigy1 = (ymax_chk_dam - yming)/deltyg + 1 
   inx0      = max(1,ipoigx0-1) 
   iny0      = max(1,ipoigy0-1)
   inx1      = min(npoigx, ipoigx1 + 1) 
   iny1      = min(npoigy, ipoigy1 + 1)                                                        
   
   DO inx    = inx0,inx1   ! topo points in the region containing the dam
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      xx     = coorg (1,ipoig)
      yy     = coorg (2,ipoig) 
      if(tandelt.LT.1.e-6) then
         RR     = (xx-x0_dam)*nx_dam + (yy-y0_dam)*ny_dam
         Zg     = topol(1,ipoig)
         if (RR.GT.0.AND.Zg.LT.Z_swl) then
            topol(1,ipoig) = Z_swl
         endif
      else
         RR     = (xx-xm)*nx_dam + (yy-ym)*ny_dam
         deltZg =  RR*tandelt
         Zg     = topol(1,ipoig)  
         if (RR.GT.0.AND.Zg.LT.(Z_swl + deltZg)) then
            topol(1,ipoig) = Z_swl + deltZg
         endif
      endif
   ENDDO
   ENDDO
    
ENDDO  ! chk dams

END SUBROUTINE get_check_dams

!--------------------------------------------------------

  subroutine Get_Structure_Pts ( TOPO_problem_name , itopo_zones)
  
!--------------------------------------------------------  


!      Reads a fem mesh describing a structures, then modify
!            topol(1,ipoig)

implicit none

integer(ink) TopEros_pts,  TopEros_chk 	!      input and chk files
integer(ink) npts_Eros  ,  nelem_Eros   !      nr pts in file, 
integer(ink) type_Eros                  !      type_Eros = +1 Z includes erodible depth
                                        !      type_Eros = -1 Z does NOT includes it : structure on top                                       
                                        !                = -2 structure intersects terrain Zs are total                                       
                                        !                = -3 basin
integer(ink) itopo_zones 
integer(ink) ip, ie, le, inode, lnode		  ! dummy nod, element, elem node
integer(ink) inx0, inx1, iny0, iny1, inx, iny ! limits in column and rows of grid for element
integer(ink) ipoig, jdimn, ipoigx, ipoigy	  ! point of topo grid and dummy index for dimension
integer(ink) node1, node2, node3		      ! nodes 1, 2 and 3 of element studied

integer(ink), allocatable:: tmp_list(:)          ! stores list of topo pts belonging to structure   ! MP 02 June  
integer(ink)                nplist, npmin, npmax ! nr of pts, min and max bounds
integer(ink)                iplist

real   (irk)  xmin_e, xmax_e, ymin_e, ymax_e	 ! each element's limits
real   (irk)  xx, yy				! x and y of a grid point
real   (irk)  x1, x2, x3, y1, y2, y3, z1, z2, z3! x,y and zs of triangle nodes
real   (irk)  x31, x32, y31, y32, x30, y30, detx! x1-x3, x2-x3,..  used to get shape functions
real   (irk)  sh1, sh2, sh3			! values of shape functions at xx and yy
real   (irk)  xg, yg				! temp values for x and y at topo mesh

real   (irk), allocatable:: coord_eros  (:,:)	! fem nodes coords (x,y,z)
integer(ink), allocatable:: intmat_eros (:,:)	! fem element nodes
real   (irk), allocatable:: aux2D       (:)     ! auxiliar, to store topo(1,:)

character(16) TOPO_problem_name
character(60) text                              ! general purpose char string
1 format(a60) 
logical if_yes                                  ! auxiliar MP feb 2022

TopEros_pts = 185
TopEros_chk = 186
len1 = len_trim(TOPO_problem_name)

open (  TopEros_pts,  file = TOPO_problem_name(1:len1)//'.eros.pts' )
open (  TopEros_chk,  file = TOPO_problem_name(1:len1)//'.eros.chk' )

read  (TopEros_pts,1) text
write (TopEros_chk,1) text       
read  (TopEros_pts,*) npts_Eros, nelem_Eros, type_Eros   
write (TopEros_chk,*) npts_Eros, nelem_Eros, type_Eros   

if (.NOT.allocated(coord_eros))  allocate ( coord_eros  (3,npts_Eros)); coord_eros = 0.0   !  X, Y and erodible depth (or Zstruct)
if (.NOT.allocated(intmat_eros)) allocate ( intmat_eros(3,nelem_eros)); intmat_eros = 0
if (.NOT.allocated(aux2D))       allocate (aux2D(npoig))              ; aux2D       = 0.0

topol (15,:) = 0.0

!      ------  read nodes: X, Y, erodible depth

read (TopEros_pts,1) text
write(TopEros_chk,1) text     

Do ip = 1, npts_Eros
   read (TopEros_pts,*) ( coord_eros(jdimn,ip), jdimn=1,3 )
   write(TopEros_chk,*) ( coord_eros(jdimn,ip), jdimn=1,3 )
ENDDO

!      ------  read elements

read (TopEros_pts,1) text
write(TopEros_chk,1) text
DO ie = 1, nelem_eros
   read (TopEros_pts,*) le, ( intmat_eros(inode,le), inode=1,3 )
   write(TopEros_chk,*) le, ( intmat_eros(inode,le), inode=1,3 )
ENDDO

!      ------  Loop in elements to find ipoigs belonging to them and interpolate

xmin_e = xmaxg; xmax_e = xming  !  we get min and max x and y of element
ymin_e = ymaxg; ymax_e = yming  !  initialize to max and mins in grid

if (.NOT.allocated(tmp_list)) allocate (tmp_list(2*npoig))                          ! MP 19 July updates to 2*npoig 02 June ***
tmp_list = 0 ; npmin = 1.e08  ; npmax = -npmin        ! MP 02 June ***
nplist = 0

DO ie = 1, nelem_eros    !  Loop on elements
   DO inode = 1, 3
      lnode = intmat_eros(inode,ie)
      if (coord_eros(1,lnode).lt.xmin_e) xmin_e = coord_eros(1,lnode) 
      if (coord_eros(1,lnode).gt.xmax_e) xmax_e = coord_eros(1,lnode)
      if (coord_eros(2,lnode).lt.ymin_e) ymin_e = coord_eros(2,lnode) 
      if (coord_eros(2,lnode).gt.ymax_e) ymax_e = coord_eros(2,lnode) 
   ENDDO
   
   node1 = intmat_eros(1,ie)   
   node2 = intmat_eros(2,ie)   
   node3 = intmat_eros(3,ie)
   x1    = coord_eros(1,node1) ; y1 = coord_eros(2,node1) ; z1 = coord_eros(3,node1) 
   x2    = coord_eros(1,node2) ; y2 = coord_eros(2,node2) ; z2 = coord_eros(3,node2)
   x3    = coord_eros(1,node3) ; y3 = coord_eros(2,node3) ; z3 = coord_eros(3,node3)
   x31   = x1 - x3 ; y31 = y1 - y3
   x32   = x2 - x3 ; y32 = y2 - y3 
   detx  = (x31*y32-y31*x32)
   
   inx0 = (xmin_e - xming)/deltxg +1
   iny0 = (ymin_e - yming)/deltyg +1
   if (inx0.lt.1) inx0 = 1
   if (iny0.lt.1) iny0 = 1

   inx1 = (xmax_e - xming)/deltxg + 2
   iny1 = (ymax_e - yming)/deltyg + 2
   if (inx1.gt.npoigx) inx1 = npoigx
   if (iny1.gt.npoigy) iny1 = npoigy
   
   DO inx   = inx0, inx1	!  We loop around points in the inscribing rectangle
   DO iny   = iny0, iny1
      ipoig = (iny-1)*npoigx + inx
      xx    = coorg (1,ipoig); yy = coorg(2,ipoig)
      x30   = xx - x3 ; y30 = yy - y3
      sh1   = (x30*y32-y30*x32)/detx
      sh2   = (y30*x31-x30*y31)/detx
      sh3   = 1.0 - sh1 - sh2
      if ( (sh1.ge.0.and.sh1.le.1).AND.(sh2.ge.0.and.sh2.le.1).AND.(sh3.ge.0.and.sh3.le.1) ) then
         topol (15, ipoig) = sh1*z1 + sh2*z2 + sh3*z3 ! interpolate Z using shape functions
         nplist = nplist + 1
         tmp_list (nplist) = ipoig
         if (ipoig.LT.npmin) npmin = ipoig
         if (ipoig.GT.npmax) npmax = ipoig
      endif
   ENDDO
   ENDDO 

ENDDO 

basal_props (2, itopo_zones) = xmin_e
basal_props (3, itopo_zones) = xmax_e
basal_props (4, itopo_zones) = ymin_e
basal_props (5, itopo_zones) = ymax_e

!  if (.NOT.allocated (topozone_info) ) allocate (topozone_info (itopo_zones)%NodeList(nplist)) 
! MP  o2 June

iplist = topozone_info (itopo_zones)%if_list   !  MP 3 March 2022
if (iplist.eq.1) then
    allocate (topozone_info (itopo_zones)%NodeList(nplist)) 
    topozone_info (itopo_zones)%if_list = 2
elseif(iplist.eq.2) then
    deallocate (topozone_info (itopo_zones)%NodeList)
    allocate (topozone_info (itopo_zones)%NodeList(nplist))
endif

do iplist = 1, nplist
   topozone_info (itopo_zones)%NodeList(iplist) = tmp_list(iplist)  !  MP 02 June
enddo
topozone_info (itopo_zones)%npmin  = npmin
topozone_info (itopo_zones)%npmax  = npmax
topozone_info (itopo_zones)%nplist = nplist

aux2D (:) = 0.0

if (type_Eros.eq.-1) then
    topol(1,:) = topol(1,:) + topol(15,:)
elseif (type_Eros.eq.-2) then
    aux2D(:) = topol(1,:)
    topol(1,:) = max ( topol(1,:) , topol(15,:) )
    topol(15,:) = topol(1,:) - aux2D(:) 
elseif (type_Eros.eq.-3) then
    topol(1,:) = topol(1,:) - topol(15,:)
endif

close ( TopEros_pts)
close ( TopEros_chk)

deallocate (coord_eros, intmat_eros, aux2D, tmp_list)


END SUBROUTINE Get_Structure_Pts



!--------------------------------------------------------

  subroutine Get_Eros_Pts ( TOPO_problem_name , itopo_zones)
  
!--------------------------------------------------------  


!      Reads erodible pts from file, store erodible depth
!            in topol(15,ipoig)

implicit none

integer(ink) TopEros_pts,  TopEros_chk 	!      ------  input and chk files
integer(ink) npts_Eros  ,  type_Eros    !      ------  nr pts in file, 
                                        !      type_Eros = +1 Z includes erodible depth
                                        !      type_Eros = -1 Z does NOT includes it
                                        
integer(ink)  ip_Src, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i   
integer(ink)  n1, n2, n3, n4, np2D, ieros, id, ndivSx, ndivSy, itopo_zones 
integer(ink)  ixmin_eros, ixmax_eros, iymin_eros, iymax_eros                                      

integer(ink), allocatable:: tmp_list(:)          ! stores list of topo pts belonging to structure   ! MP 02 June  
integer(ink)                nplist, npmin, npmax ! nr of pts, min and max bounds
integer(ink)                iplist

real   (irk) dx_s                            ! representative spacing of erodible pts 
real   (irk), allocatable:: x_Eros     (:,:) ! coordinates of erodible pts (ndimn,npts_Eros) 
real   (irk), allocatable:: depth_Eros (:)   ! depth of erodible pts (npts_Eros) 
real   (irk) xeros_min, xeros_max, yeros_min, yeros_max   ! erosion points limits

real   (irk), allocatable::  aux2D(:,:)
 
real   (irk)  dist2D, xx, yy ,xe, ye 
real   (irk)  distx2, disty2, dist, delt2D

character(16) TOPO_problem_name
character(60) text                    ! general purpose char string
1 format(a60) 

TopEros_pts = 185
TopEros_chk = 186
len1 = len_trim(TOPO_problem_name)

open (  TopEros_pts,  file = TOPO_problem_name(1:len1)//'.eros.pts' )
open (  TopEros_chk,  file = TOPO_problem_name(1:len1)//'.eros.chk' )

read (TopEros_pts,1) text
write(TopEros_chk,1) text       
 
read  (TopEros_pts,*) npts_Eros, type_Eros, dx_s 
write (TopEros_chk,*) npts_Eros, type_Eros, dx_s  

delt2D = max (deltxg, deltyg)                ! we use deltxg, deltyg
dist2D = (deltxg*deltxg + deltyg*deltyg)**0.5
np2D   = npoigx * npoigy                     ! is npoig

if (.NOT.allocated (x_Eros) ) then 
   allocate (x_Eros(2,npts_Eros))     ; x_Eros     = 0.0
   allocate (depth_Eros  (npts_Eros)) ; depth_Eros = 0.0
   allocate (aux2D(4,np2D))           ; aux2D      = 0.0
   allocate (tmp_list(npoig))                          ! MP 02 June ***
endif
tmp_list = 0 ; npmin = 1.e08  ; npmax = -npmin      ! MP 02 June ***
nplist = 0

topol (15,:) = 0.0

read (TopEros_pts,1) text
write(TopEros_chk,1) text 

xeros_min = 1.e10 ; xeros_max = -1.e10; yeros_min = 1.e10 ; yeros_max = -1.e10
Do ieros = 1, npts_Eros
   read  (TopEros_pts,*) (x_Eros(id,ieros),id=1,2), depth_Eros (ieros) 
   write (TopEros_chk,*) (x_Eros(id,ieros),id=1,2), depth_Eros (ieros)
   xe = x_Eros(1,ieros); ye = x_Eros(2,ieros) 
   if (xe.GE.xeros_max) xeros_max = xe
   if (ye.GE.yeros_max) yeros_max = ye
   if (xe.LE.xeros_min) xeros_min = xe
   if (ye.LE.yeros_min) yeros_min = ye
enddo

basal_props (2, itopo_zones) = xeros_min
basal_props (3, itopo_zones) = xeros_max
basal_props (4, itopo_zones) = yeros_min
basal_props (5, itopo_zones) = yeros_max

ixmin_eros = (xeros_min - xming)/deltxg + 1
ixmax_eros = (xeros_max - xming)/deltxg + 1
iymin_eros = (yeros_min - yming)/deltyg + 1
iymax_eros = (yeros_max - yming)/deltyg + 1

ndivSx = max (2.,(dx_S/deltxg))
ndivSy = max (2.,(dx_S/deltyg))

DO ip_Src = 1, npts_Eros
   xx = x_Eros(1,ip_Src)
   yy = x_Eros(2,ip_Src)
   ipoi2Dx    = (xx-xming)/deltxg + 1
   ipoi2Dy    = (yy-yming)/deltyg + 1

   inx0       = max(1,      ipoi2Dx - ndivSx )      ! max(1,      ipoi2Dx - ndivSx + 1 )
   iny0       = max(1,      ipoi2Dy - ndivSy     )
   inx1       = min(npoigx,ipoi2Dx + ndivSx + 1 )
   iny1       = min(npoigy,ipoi2Dy + ndivSy + 1 )   !   added +1

   inx0       = max(inx0 , ixmin_eros)
   inx1       = min(inx1 , ixmax_eros)
   iny0       = max(iny0 , iymin_eros)
   iny1       = min(iny1 , iymax_eros)

   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoigx + inx
      distx2  = ( xx - coorg(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - coorg(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = aux2D (1,ipoi2D) + 1 
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +        1.0/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) + depth_Eros(ip_Src)/dist      
   enddo
   enddo
ENDDO

!      ------  Obtain aux2D(4,ipoi2D) code at nodes    JUST in CASE
!                     4 ... interior node, connected to 4 rects. full dA
!                     3 ... boundary corner, convex, 3 rects..........dA/3
!                     2 ... boundary side, plane,    2 rects          dA/2
!                     1 ... corner node, 1 rect only..                dA/4


DO ipoi2Dy = 1, npoigy-1
DO ipoi2Dx = 1, npoigx-1
   ipoi2D  = (ipoi2Dy-1)*npoigx + ipoi2Dx
   n1      = ipoi2D
   n2      = n1 + 1 
   n3      = n1 + npoigx  
   n4      = n3 + 1
   if ( (aux2D(2,n1).ge.1.e-6).AND.(aux2D(2,n2).ge.1.e-6).AND.  &
        (aux2D(2,n3).ge.1.e-6).AND.(aux2D(2,n4).ge.1.e-6) ) THEN
         aux2D(4,n1) = aux2D(4,n1) + 1
         aux2D(4,n2) = aux2D(4,n2) + 1
         aux2D(4,n3) = aux2D(4,n3) + 1
         aux2D(4,n4) = aux2D(4,n4) + 1
   endif           
enddo
enddo


!      ------  Get values at topol nodes

ip   = 0 
DO ipoi2Dy = 1, npoigy
DO ipoi2Dx = 1, npoigx
   ip      = ip + 1
   if ( aux2D(2,ip).ge.1.e-6) then 
       topol (15,ip) = aux2D (3,ip)/aux2D (2,ip)  
       nplist = nplist + 1
       tmp_list (nplist) = ip 
       if (ip.LT.npmin) npmin = ip 
       if (ip.GT.npmax) npmax = ip       
   endif   
enddo
enddo

iplist = topozone_info (itopo_zones)%if_list   !  MP 3 March 2022
if (iplist.eq.1) then
    allocate (topozone_info (itopo_zones)%NodeList(nplist)) 
    topozone_info (itopo_zones)%if_list = 2
elseif(iplist.eq.2) then
    deallocate (topozone_info (itopo_zones)%NodeList)
    allocate (topozone_info (itopo_zones)%NodeList(nplist))
endif

do iplist = 1, nplist
   topozone_info (itopo_zones)%NodeList(iplist) = tmp_list(iplist)  !  MP 02 June
enddo
topozone_info (itopo_zones)%npmin  = npmin
topozone_info (itopo_zones)%npmax  = npmax
topozone_info (itopo_zones)%nplist = nplist

aux2D (1,:) = 0.0
if (type_Eros.eq.-1) then
    topol(1,:) = topol(1,:) + topol(15,:)
elseif (type_Eros.eq.-2) then
    aux2D(1,:) = topol(1,:)
    topol(1,:) = max ( topol(1,:) , topol(15,:) )
    topol(15,:) = topol(1,:) - aux2D(1,:)
elseif (type_Eros.eq.-3) then 
    topol(1,:) = topol(1,:) - topol(15,:)   !  MP CH 19th May 2019
endif

close ( TopEros_pts)
close ( TopEros_chk)
deallocate (x_Eros)
deallocate (depth_Eros, aux2D, tmp_list)

END SUBROUTINE Get_Eros_Pts



!--------------------------------------------------------

  subroutine get_cone
  
!--------------------------------------------------------  


!      Processes data for cone á la Stromboli

implicit none

integer(ink) ipoig

real   (irk) R_Top, R_Bottom, Z_Top, Z_bottom
real   (irk) xx, yy, rr 

character(60) text                    ! general purpose char string
1 format(a60) 

read (topDat_file,1) text
write(Topchk_file,1) text       
 
read  (topdat_file,*) R_Top, R_Bottom, Z_Top, Z_bottom
write (Topchk_file,*) R_Top, R_Bottom, Z_Top, Z_bottom

do ipoig = 1,npoig
   xx = coorg(1,ipoig)
   yy = coorg(2,ipoig)
   rr = (xx*xx + yy*yy)**0.5
   if ( rr.le.r_top ) then
      topol (1, ipoig) = Z_Top
   elseif ( rr.gt.R_Bottom) then
      topol (1,ipoig)  = Z_Bottom
   else
      topol (1,ipoig)  = Z_Bottom + (Z_Top - Z_Bottom)*(R_Bottom - rr)/(R_Bottom- R_Top)
   endif
enddo 


END SUBROUTINE Get_cone


!---------------------------------------------------------

  subroutine Get_DTM (ictop)
       
!---------------------------------------------------------


!      Processes data from a DTM file...( text, np..for all ip Xip Yip Zip ) 

!     for ictops 12=10 and 13=11 plus reading a rect. and its H
!                15 reads a window


implicit none

integer(ink) np_dtm, ictop   
integer(ink) ip    , jdimn, ipoig, ipoigx, ipoigy, inx0,  inx1,  iny0,  iny1,  inx,  iny
LOGICAL      discr

real   (irk), allocatable:: coor_dtm  (:,:), auxv(:,:)
real   (irk)  deltx
real   (irk)  xx, yy, zz, ztemp, xg, yg, distx2, disty2, dist, distg, Zp
real   (irk)  zmaxg, zming
real   (irk)  faknes
real   (irk)  xw0, xw1, yw0, yw1      !  coords of window corner

integer(ink) if_filter, filter_type
real   (irk) xxf0, yyf0, xxf1, yyf1, hhf

character(60) text                    ! general purpose char string

if_filter = 0
if (ictop.eq.12.or.ictop.eq.13) then  ! we initialize if_filter to 1  
    if_filter = 1		      ! and will read filter info
    ictop    = ictop - 2
endif    

read (topDat_file,1) text
write(Topchk_file,1) text
read (topDat_file,*) np_dtm, deltx
write(Topchk_file,*) np_dtm, deltx

if (.NOT.allocated (coor_dtm) ) allocate ( coor_dtm(3,np_dtm) )
    
xming    =  1.e9
yming    =  1.e9
xmaxg    = -1.e9
ymaxg    = -1.e9
zming    =  1.e9
zmaxg    = -1.e9
xw0      =  xming
xw1      =  xmaxg
yw0      =  yming
yw1      =  ymaxg

IF (ictop.eq.15.or.ictop.eq.16) then  !   ----- WINDOW
   read (topDat_file,1) text
   write(Topchk_file,1) text
   read (topDat_file,*) xw0, xw1, yw0, yw1
   write(Topchk_file,*) xw0, xw1, yw0, yw1
Endif
    
read (topDat_file,1) text
write(Topchk_file,1) text

DO ip = 1, np_dtm
   read (topDat_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 )
   write(Topchk_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 ) 
   if (coor_dtm(1,ip).le.xming) xming = coor_dtm(1,ip)
   if (coor_dtm(1,ip).ge.xmaxg) xmaxg = coor_dtm(1,ip)
   if (coor_dtm(2,ip).le.yming) yming = coor_dtm(2,ip)
   if (coor_dtm(2,ip).ge.ymaxg) ymaxg = coor_dtm(2,ip)
   if (coor_dtm(3,ip).le.zming) zming = coor_dtm(3,ip)
   if (coor_dtm(3,ip).ge.zmaxg) zmaxg = coor_dtm(3,ip)
ENDDO

npoigx  = (xmaxg-xming)/deltx + 1
npoigy  = (ymaxg-yming)/deltx + 1
npoig   = npoigx*npoigy 
deltxg  = (xmaxg-xming)/(npoigx-1)
deltyg  = (ymaxg-yming)/(npoigy-1)
xshiftg = 0.0
yshiftg = 0.0

IF (ictop.eq.11) then			!      ------  In this case we use relative coords
   xshiftg = xming
   yshiftg = yming
   DO ip = 1, np_dtm
      coor_dtm(1,ip) = coor_dtm(1,ip) - xming
      coor_dtm(2,ip) = coor_dtm(2,ip) - yming
      write(Topchk_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 )       
   ENDDO
   Xming_org      =  xming   ! Saeidmt 12Nov2023 Xming_org, Xmaxg_org, Yming_org, Ymaxg_org added to correctly read the deposit shape
   Xmaxg_org      =  xmaxg
   Yming_org      =  yming
   Ymaxg_org      =  ymaxg
   xmaxg = xmaxg - xming 
   ymaxg = ymaxg - yming 
   xming = 0.0
   yming = 0.0
ELSEIF (ictop.eq.15) then               
   xshiftg = xming
   yshiftg = yming
   xming   = xw0 - xshiftg
   xmaxg   = xw1 - xshiftg
   yming   = yw0 - yshiftg
   ymaxg   = yw1 - yshiftg
   npoigx  = (xmaxg-xming)/deltx + 1
   npoigy  = (ymaxg-yming)/deltx + 1
   npoig   = npoigx*npoigy 
   deltxg  = (xmaxg-xming)/(npoigx-1)
   deltyg  = (ymaxg-yming)/(npoigy-1)
   DO ip = 1, np_dtm
      coor_dtm(1,ip) = coor_dtm(1,ip) - xshiftg
      coor_dtm(2,ip) = coor_dtm(2,ip) - yshiftg
      write(Topchk_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 )       
   ENDDO
ELSEIF (ictop.eq.16) then   ! window in local coordinates             
   xshiftg = xming
   yshiftg = yming
   xming   = xw0  
   xmaxg   = xw1  
   yming   = yw0  
   ymaxg   = yw1  
   npoigx  = (xmaxg-xming)/deltx + 1
   npoigy  = (ymaxg-yming)/deltx + 1
   npoig   = npoigx*npoigy 
   deltxg  = (xmaxg-xming)/(npoigx-1)
   deltyg  = (ymaxg-yming)/(npoigy-1)
   DO ip = 1, np_dtm
      coor_dtm(1,ip) = coor_dtm(1,ip) - xshiftg
      coor_dtm(2,ip) = coor_dtm(2,ip) - yshiftg
      write(Topchk_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 )       
   ENDDO
   
ENDIF

if (.NOT.allocated (coorg) ) allocate (coorg(2,npoig),topol(ntopo,npoig),  poig_in_model(npoig) )

coorg         = 0.0
topol         = 0.0   
poig_in_model = 0

!  Get coordinates of structured mesh npoigx x npoigy (divs)

ipoig  = 0
DO ipoigy = 1,npoigy 
    yg     = yming + (ipoigy-1)*deltyg
    Do ipoigx = 1,npoigx 
       ipoig  = ipoig + 1
       xg     = xming + (ipoigx-1)*deltxg
       coorg(1,ipoig) = xg
       coorg(2,ipoig) = yg
    enddo
ENDDO

!  Get totpographic elevations at structured mesh npoigx x npoigy (divs)

if (.NOT.allocated (auxv) ) allocate (auxv(2,npoig))
auxv      = 0.0 
distg = (deltxg*deltxg + deltyg*deltyg)**0.5

DO ip = 1, np_dtm
   xx = coor_dtm(1,ip)
   yy = coor_dtm(2,ip)
   discr     = (xx.ge.xming).and.(xx.le.xmaxg).and.(yy.ge.yming).and.(yy.le.ymaxg)
   if (discr) then
      ipoigx    = (xx-xming)/deltxg + 1
      ipoigy    = (yy-yming)/deltyg + 1
      inx0      = max(1,ipoigx-1)
      iny0      = max(1,ipoigy-1)
      inx1      = min(npoigx,ipoigx+2)
      iny1      = min(npoigy,ipoigy+2)
      DO inx    = inx0,inx1
      DO iny    = iny0,iny1
         ipoig  = (iny-1)*npoigx + inx
         distx2 = ( xx - coorg(1,ipoig) )
         distx2 = distx2*distx2
         disty2 = ( yy - coorg(2,ipoig) )
         disty2 = disty2*disty2
         dist   = max( 1.e-6, (distx2+disty2)**0.5 )
         dist   = dist/distg
         auxv(1,ipoig) = auxv(1,ipoig) + coor_dtm(3,ip)/dist 
         auxv(2,ipoig) = auxv(2,ipoig) + 1.00/dist
      enddo
      enddo
   endif    
ENDDO

!      ------  Finally, we store in topol(1...) topo elevation

ipoig  = 0
DO ipoigy = 1,npoigy 
   Do ipoigx = 1,npoigx 
      ipoig  = ipoig + 1
      if (auxv(2,ipoig).ge.1.e-6) then
          topol(1,ipoig) = auxv(1,ipoig)/auxv(2,ipoig)
          poig_in_model(ipoig) = 1
      else
          !topol(1,ipoig) = zming
          !  2018.4.6 CHUAN to avoid mistake when create water particles           
          topol(1,ipoig) = zmaxg
          !  2018.4.6 CHUAN to avoid mistake when create water particles end
          poig_in_model(ipoig) = 0
      endif
   enddo
ENDDO 


!write(*,*) '  This is ICUNK 11 You can use filters '
!write(*,*) '  Type 0 to skip or 1 to apply filters '
! read (*,*)  icfilter


IF (if_filter.eq.1) THEN

   read (topDat_file,1) text
   write(Topchk_file,1) text
   read (topDat_file,*) filter_type
   write(Topchk_file,*) filter_type
   
   if (filter_type.eq.1) then
   
          read (topDat_file,1) text
          write(Topchk_file,1) text
          read (topDat_file,*) xxf0, yyf0, xxf1, yyf1, hhf
          write(Topchk_file,*) xxf0, yyf0, xxf1, yyf1, hhf
          
          do ipoig  = 1, npoig
             xx     = coorg(1,ipoig)
             yy     = coorg(2,ipoig)
             if ( (xx.ge.xxf0).and.(xx.le.xxf1).and.(yy.ge.yyf0).and.(yy.le.yyf1) ) then
                 topol (1,ipoig) = hhf 
             endif
          enddo
          
   elseif (filter_type.eq.2) then  !      ------  Filter for AKNES 
   				   !              Just add at the end factor faknes*Z
          ipoig  = 0
          DO ipoigy = 1,npoigy 
          Do ipoigx = 1,npoigx 
             ipoig  = ipoig + 1
             xx     = coorg(1,ipoig)
             yy     = coorg(2,ipoig)
             zz     = topol(1,ipoig)
             if (yy.ge.4075) then
                 ztemp = -400. + 500*(yy-4075)/(4475-4075)
             else if (yy.le.500) then
                 ztemp = 100.0 - 500*(yy-0.0)/(400.-0.)
             endif
             topol (1,ipoig) = max (zz,ztemp) 
          enddo
          ENDDO 
                    
	  read (topDat_file,1) text	 !   '  We can lower slopes by multiplying Z by a factor '
	  write(Topchk_file,1) text	 !   '  Type factor '
	  read (topDat_file,*) faknes
          write(Topchk_file,*) faknes
          
          if ( abs(faknes-1.0).gt.1.e-3 ) then
             do ipoig = 1, npoig
                ztemp = topol(1,ipoig)*faknes
                topol (1,ipoig) = ztemp 
             enddo
          endif 
   
   endif   
   
ENDIF

!   
!   write(*,*) '  This is for TEST 03 in Hong Kong set '
!   write(*,*) '  Type 0 or 1 to apply it '
!   read (*,*) iaknes
!   if (iaknes.eq.1) then
!       do ipoig = 1, npoig
!          xx     = coorg(1,ipoig)
!          yy     = coorg(2,ipoig)
!          if ( (9.4.le.xx).and.(xx.le.64).and.(8.4.le.yy).and.(yy.le.90.8) ) then
!             call Get_Z_Topo ( ipoig,xx,100.,Zp)
!             topol (1,ipoig) = Zp
!          endif
!       enddo
!   endif
!   
   

1 format(a60) 

deallocate (auxv, coor_dtm)


END SUBROUTINE  Get_DTM



!---------------------------------------------------------

  subroutine Get_DTM_FEM (ictop)
       
!---------------------------------------------------------


!      Processes data from a DTM file.  FEM format..
!      file is file.top, one file only
!     ( text, np..for all ip Xip Yip Zip ) 
!     ( text, nelem,  text, ielem, n1, n2, n3 )

implicit none

integer(ink) np_dtm, ne_dtm,  ictop		! nodes, elements and topo control(20 or 21)
real   (irk) deltx 				! topo structured grid size 
real   (irk), allocatable:: coor_dtm  (:,:)	! fem nodes coords (x,y,z)
integer(ink), allocatable:: intmat_dtm(:,:)	! fem element nodes

real   (irk)  xmin_e, xmax_e, ymin_e, ymax_e	! each element's limits
real   (irk)  xx, yy				! x and y of a grid point
real   (irk)  x1, x2, x3, y1, y2, y3, z1, z2, z3! x,y and zs of triangle nodes
real   (irk)  x31, x32, y31, y32, x30, y30, detx! x1-x3, x2-x3,..  used to get shape functions
real   (irk)  sh1, sh2, sh3			! values of shape functions at xx and yy
real   (irk)  xg, yg				! temp values for x and y at topo mesh

integer(ink) ip, ie, le, inode, lnode		! dummy nod, element, elem node
integer(ink) inx0, inx1, iny0, iny1, inx, iny	! limits in column and rows of grid for element
integer(ink) ipoig, jdimn, ipoigx, ipoigy	! point of topo grid and dummy index for dimension
integer(ink) node1, node2, node3		! nodes 1, 2 and 3 of element studied

character(60) text                    		! general purpose char string

integer(ink) ic_mask				! we will set a mask to provide a given Z_mask
real   (irk) Z_mask
real   (irk) Z_mask0				! This is the min Z at the fem mesh. 
						!      we set topol(1,:) = Zmask_0
Z_mask0  = +1.e9				! (b1) see also 
						
read (topDat_file,1) text
write(Topchk_file,1) text
read (topDat_file,*) np_dtm, deltx
write(Topchk_file,*) np_dtm, deltx

if (.NOT.allocated (coor_dtm) )  allocate ( coor_dtm(3,np_dtm) )
    
xming    =  1.e9
yming    =  1.e9
xmaxg    = -1.e9
ymaxg    = -1.e9
read (topDat_file,1) text
write(Topchk_file,1) text

DO ip = 1, np_dtm
   read (topDat_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 )
   if (coor_dtm(1,ip).le.xming)     xming = coor_dtm(1,ip)
   if (coor_dtm(1,ip).ge.xmaxg)     xmaxg = coor_dtm(1,ip)
   if (coor_dtm(2,ip).le.yming)     yming = coor_dtm(2,ip)
   if (coor_dtm(2,ip).ge.ymaxg)     ymaxg = coor_dtm(2,ip)
   if (coor_dtm(3,ip).lt.Z_mask0) Z_mask0 = coor_dtm(3,ip)   ! (b1) all exterior domain at Z_mask0
ENDDO							     ! use later mask to avoid escaping nodes

npoigx  = (xmaxg-xming)/deltx + 1
npoigy  = (ymaxg-yming)/deltx + 1
npoig   = npoigx*npoigy 
deltxg  = (xmaxg-xming)/(npoigx-1)
deltyg  = (ymaxg-yming)/(npoigy-1)
xshiftg = 0.0
yshiftg = 0.0

if (ictop.eq.21) then			!      ------  In this case we use relative coords
   xshiftg = xming
   yshiftg = yming
   DO ip = 1, np_dtm
      coor_dtm(1,ip) = coor_dtm(1,ip) - xming
      coor_dtm(2,ip) = coor_dtm(2,ip) - yming
      write(Topchk_file,*) ( coor_dtm(jdimn,ip), jdimn=1,3 )       
   ENDDO
   xmaxg = xmaxg - xming 
   ymaxg = ymaxg - yming 
   xming = 0.0
   yming = 0.0
endif

if (.NOT.allocated (coorg) ) allocate (coorg(2,npoig),topol(ntopo,npoig), poig_in_model(npoig) )
coorg    = 0.0
topol    = 0.0 
topol (1,:) = Z_mask0	!  (b1) this is done for all nodes not belonging to mesh 

!  Now read element info

read (topDat_file,1) text
write(Topchk_file,1) text
read (topDat_file,*) ne_dtm 
write(Topchk_file,*) ne_dtm 

if (.NOT.allocated (intmat_dtm) ) allocate ( intmat_dtm(3,ne_dtm) )
intmat_dtm = 0
    
read (topDat_file,1) text
write(Topchk_file,1) text
DO ie = 1, ne_dtm
   read (topDat_file,*) le, ( intmat_dtm(inode,le), inode=1,3 )
   write(Topchk_file,*) le, ( intmat_dtm(inode,le), inode=1,3 )
ENDDO

!  Get coordinates of structured mesh npoigx x npoigy (divs)

ipoig  = 0
DO ipoigy = 1,npoigy 
    yg     = yming + (ipoigy-1)*deltyg
    Do ipoigx = 1,npoigx 
       ipoig  = ipoig + 1
       xg     = xming + (ipoigx-1)*deltxg
       coorg(1,ipoig) = xg
       coorg(2,ipoig) = yg
    enddo
ENDDO

!  Loop in elements to find ipoigs belonging to them and interpolate

poig_in_model = 0   !  we initialize to zero, and set to 1 points within elements

DO ie = 1,ne_dtm    !  Loop on elements

   xmin_e = xmaxg; xmax_e = xming  !  we get min and max x and y of element
   ymin_e = ymaxg; ymax_e = yming  !  initialize to max and mins in grid
   DO inode = 1, 3
      lnode = intmat_dtm(inode,ie)
      if (coor_dtm(1,lnode).lt.xmin_e) xmin_e = coor_dtm(1,lnode) 
      if (coor_dtm(1,lnode).gt.xmax_e) xmax_e = coor_dtm(1,lnode)
      if (coor_dtm(2,lnode).lt.ymin_e) ymin_e = coor_dtm(2,lnode) 
      if (coor_dtm(2,lnode).gt.ymax_e) ymax_e = coor_dtm(2,lnode) 
   ENDDO
   
   node1 = intmat_dtm(1,ie)   
   node2 = intmat_dtm(2,ie)   
   node3 = intmat_dtm(3,ie)
   x1    = coor_dtm(1,node1) ; y1 = coor_dtm(2,node1) ; z1 = coor_dtm(3,node1) 
   x2    = coor_dtm(1,node2) ; y2 = coor_dtm(2,node2) ; z2 = coor_dtm(3,node2)
   x3    = coor_dtm(1,node3) ; y3 = coor_dtm(2,node3) ; z3 = coor_dtm(3,node3)
   x31   = x1 - x3 ; y31 = y1 - y3
   x32   = x2 - x3 ; y32 = y2 - y3 
   detx  = (x31*y32-y31*x32)
   
   inx0 = (xmin_e - xming)/deltxg +1
   iny0 = (ymin_e - yming)/deltyg +1
   if (inx0.lt.1) inx0 = 1
   if (iny0.lt.1) iny0 = 1

   inx1 = (xmax_e - xming)/deltxg + 2
   iny1 = (ymax_e - yming)/deltyg + 2
   if (inx1.gt.npoigx) inx1 = npoigx
   if (iny1.gt.npoigy) iny1 = npoigy
   
   DO inx   = inx0, inx1	!  We loop around points in the inscribing rectangle
   DO iny   = iny0, iny1
      ipoig = (iny-1)*npoigx + inx
      xx    = coorg (1,ipoig); yy = coorg(2,ipoig)
      x30   = xx - x3 ; y30 = yy - y3
      sh1   = (x30*y32-y30*x32)/detx
      sh2   = (y30*x31-x30*y31)/detx
      sh3   = 1.0 - sh1 - sh2
	  if ( (sh1.ge.0.and.sh1.le.1).AND.(sh2.ge.0.and.sh2.le.1).AND.(sh3.ge.0.and.sh3.le.1) ) then
         poig_in_model(ipoig) = 1
         topol (1, ipoig) = sh1*z1 + sh2*z2 + sh3*z3 ! interpolate Z using shape functions
      endif
   ENDDO
   ENDDO 

ENDDO 

write(*,*) '  ....  '
write(*,*) '  ....  '
write(*,*) ' do you want mask for heights upwind? (0,1) '
read (*,*)  ic_mask

if (ic_mask.eq.1) then
   read (topDat_file,1) text
   write(Topchk_file,1) text
   read (topDat_file,*) xmin_e, xmax_e, ymin_e, ymax_e, Z_mask
   write(Topchk_file,*) xmin_e, xmax_e, ymin_e, ymax_e, Z_mask   
   inx0 = (xmin_e - xming)/deltxg +1
   iny0 = (ymin_e - yming)/deltyg +1
   if (inx0.lt.1) inx0 = 1
   if (iny0.lt.1) iny0 = 1
   inx1 = (xmax_e - xming)/deltxg + 2
   iny1 = (ymax_e - yming)/deltyg + 2
   if (inx1.gt.npoigx) inx1 = npoigx
   if (iny1.gt.npoigy) iny1 = npoigy
   DO inx   = inx0, inx1	!  We loop around points in the inscribing rectangle
   DO iny   = iny0, iny1
      ipoig = (iny-1)*npoigx + inx
      if (poig_in_model(ipoig).eq.0) then
          topol (1, ipoig) = Z_mask
          poig_in_model(ipoig)=1
      endif
   ENDDO
   ENDDO 
endif   
         
 
1 format(a60) 

deallocate (intmat_dtm, coor_dtm)

END SUBROUTINE Get_DTM_FEM


!---------------------------------------------------------

  subroutine getgeom2D (neleg,intmag,geome) 
       
!---------------------------------------------------------


!      Obtains shape function derivatives. Note in our case only 2 diff elements, so it can be simplified!!

integer(ink) neleg, intmag(3,neleg)
integer(ink) ieleg, inodg, in 

real   (irk) geome (7,neleg)
real   (irk) x(3),y(3),nxi(3),net(3),elcod(2,3)     
real   (irk) x21, x31, y21, y31, rj, rj1, xix, xiy, etx, ety
real   (irk) rnxi, rnet

data  nxi/-1.0 , 1.0 , 0.0 / !derivadas de las funciones de forma respecto x
data  net/-1.0 , 0.0 , 1.0 / !derivadas de las funciones de forma respecto y

!     Shape functions are:
!		N1=1-xi-et
!		N2=xi
!		N3=et

DO  ieleg = 1,neleg

    do inodg = 1,3
       in       = intmag(inodg,ieleg)
       x(inodg) = coorg(1,in)
       y(inodg) = coorg(2,in)
       elcod(1,inodg) = x(inodg)
       elcod(2,inodg) = y(inodg)
    enddo

    x21 = x(2)-x(1)
    x31 = x(3)-x(1)
    y21 = y(2)-y(1)
    y31 = y(3)-y(1)
    rj  = x21*y31-x31*y21
    rj1 = 1./rj
 
    xix = y31*rj1
    xiy =-x31*rj1
    etx =-y21*rj1
    ety = x21*rj1

!     ----form n,x & n,y

    do in = 1,3
       rnxi = nxi(in)
       rnet = net(in)
       geome(in  ,ieleg) = xix*rnxi + etx*rnet
       geome(in+3,ieleg) = xiy*rnxi + ety*rnet
    enddo
 
    geome(7,ieleg) = rj    ! Area

ENDDO

END SUBROUTINE  getgeom2D 


!--------------------------------------------------------------------

       subroutine Get_Topol (neleg,intmag,geome)
       

!--------------------------------------------------------------------

!      Purpose : given  topol obtain its derivatives x and y
!                and second order derivatives   MP 01 04 2003  and for sph 28-05-2005

integer(ink) neleg, intmag(3,neleg) 
integer(ink) niter, iiter, id, in
integer(ink) ieleg, ipoig, jpoig, inodg, jnodg, idimn

real   (irk) geome (7,neleg), area3

real   (irk), allocatable:: topel   (:,:)   ! -- topo props at elements 
real   (irk), allocatable:: rhs     (:,:)   ! -- topo props at elements
real   (irk), allocatable:: rhelp   (:,:)   ! -- topo props at elements
real   (irk), allocatable:: mmatl   (:)     ! -- Lumped mass matrix

real   (irk)  mmat(3,3), derixe, deriye, coe, derixxe, derixye, deriyxe, deriyye 
real   (irk)  zx, zy, zxx, zxy, zyy, rt,  cm

niter = 5       !  ------  Nr. of Jacobi iterations 

if (.NOT.allocated (topel) ) allocate ( topel(ntopo,neleg) )
if (.NOT.allocated (rhs)   ) allocate ( rhs  (ntopo,npoig) )
if (.NOT.allocated (rhelp) ) allocate ( rhelp(ntopo,npoig) )
if (.NOT.allocated (mmatl) ) allocate ( mmatl(npoig) )
topel = 0.0
rhs   = 0.0
rhelp = 0.0
mmatl = 0.0

!      ------ A1  Get Mass matrix (consistent and lumped) 

mmat = 1.0
do inodg = 1,3
   mmat(inodg,inodg) = 2.0
enddo

mmatl = 0.0
do ieleg = 1, neleg
   area3 = geome(7,ieleg)/6.0
   do inodg = 1,3
      ipoig = intmag(inodg, ieleg)
      mmatl(ipoig) = mmatl(ipoig) + area3
   enddo
enddo

!      ------ A2   Obtain gradients at elements  

do ieleg     = 1,neleg
   derixe    = 0.0
   deriye    = 0.0
   do inodg  = 1,3
      ipoig  = intmag(inodg,ieleg)
      derixe = derixe + topol(1,ipoig)*geome(inodg,ieleg)
      deriye = deriye + topol(1,ipoig)*geome(inodg+3,ieleg)
   enddo
   topel(2,ieleg) = derixe
   topel(3,ieleg) = deriye
enddo

!      ------ A3  Build RHS Int(Ni.Grad) Store it in rhs

rhs = 0.0

do ieleg = 1,neleg
   coe = geome(7,ieleg)/6.
   do inodg = 1,3
      ipoig = intmag(inodg,ieleg)
      rhs (2,ipoig) = rhs(2,ipoig) + topel(2,ieleg)*coe
      rhs (3,ipoig) = rhs(3,ipoig) + topel(3,ieleg)*coe
   enddo
enddo

!      ------ A4   Get topol v(1)=v(0)+dt.inv(Ml)*(rhs-M.v(0)) v(0)=0

do ipoig = 1,npoig    
   cm    = mmatl(ipoig)
   topol(2,ipoig) = rhs(2,ipoig)/cm !+topol(2,ipoig) (EG)
   topol(3,ipoig) = rhs(3,ipoig)/cm !topol(3,ipoig)  (EG)
enddo

!      ------ A5  Iterate:       v(n+1)=v(n)+inv(Ml)*(rhs-M*V(n))
!                                M*v(n) is stored as rhelp
!       if (niter.ne.1) then

    do iiter = 2,niter
       rhelp = 0.0
       do ieleg = 1,neleg
          coe = geome(7,ieleg)/24.
          do inodg = 1,3
             ipoig = intmag(inodg,ieleg)
             do jnodg = 1,3
                jpoig = intmag(jnodg,ieleg)
                cm    = mmat(inodg,jnodg)
                do idimn = 1,2
                   rhelp(idimn+1,ipoig) = rhelp(idimn+1,ipoig) + cm*coe*topol(idimn+1,jpoig)
                enddo
             enddo
          enddo
       enddo

       do ipoig = 1,npoig
          cm = mmatl(ipoig)
          do idimn = 1,2
             topol(idimn+1,ipoig) = topol(idimn+1,ipoig) + (rhs(idimn+1,ipoig)-rhelp(idimn+1,ipoig))/cm
          enddo
       enddo
    enddo
!       ENDIF
       
!      ------ Second order derivatives

!      ------ Obtain gradients at elements  

do ieleg = 1,neleg
   derixxe = 0.0
   derixye = 0.0
   deriyxe = 0.0
   deriyye = 0.0
   do inodg = 1,3
      ipoig   = intmag(inodg,ieleg)
      derixxe = derixxe + topol(2,ipoig)*geome(inodg  ,ieleg)
      derixye = derixye + topol(2,ipoig)*geome(inodg+3,ieleg)
      deriyxe = deriyxe + topol(3,ipoig)*geome(inodg  ,ieleg)
      deriyye = deriyye + topol(3,ipoig)*geome(inodg+3,ieleg)
   enddo
   topel(4,ieleg) =  derixxe
   topel(5,ieleg) = (derixye + deriyxe)/2.
   topel(6,ieleg) =  deriyye
enddo

!      ------ B3  Build RHS Int(Ni.Grad) Store it in rhs

rhs = 0.0

do ieleg =1,neleg
   coe = geome(7,ieleg)/6.
   do inodg = 1,3
      ipoig = intmag(inodg,ieleg)
      rhs (1,ipoig) = rhs(1,ipoig) + topel(4,ieleg)*coe 
      rhs (2,ipoig) = rhs(2,ipoig) + topel(5,ieleg)*coe
      rhs (3,ipoig) = rhs(3,ipoig) + topel(6,ieleg)*coe
   enddo
enddo

!      ------  B4   Get topol v(1)=v(0)+dt.inv(Ml)*(rhs-M.v(0)) v(0)=0

do ipoig=1,npoig
   cm = mmatl(ipoig)
   topol(4,ipoig) =  rhs(1,ipoig)/cm !+ topol(4,ipoig) 
   topol(5,ipoig) =  rhs(2,ipoig)/cm !+ topol(5,ipoig) 
   topol(6,ipoig) =  rhs(3,ipoig)/cm !+ topol(6,ipoig) 
enddo

!      ------   B5  Iterate:       v(n+1)=v(n)+inv(Ml)*(rhs-M*V(n))
!                                  M*v(n) is stored as rhelp

!       if (niter.ne.1) then

do iiter = 2,niter
   rhelp = 0.0
   do ieleg = 1,neleg
      coe         = geome(7,ieleg)/24.
      do inodg    = 1,3
         ipoig    = intmag(inodg,ieleg)
         do jnodg = 1,3
            jpoig = intmag(jnodg,ieleg)
            cm    = mmat(inodg,jnodg)
            do id = 1,3
               rhelp(id,ipoig) = rhelp(id,ipoig) + cm*coe*topol(id+3,jpoig)
            enddo
         enddo
      enddo
   enddo

   do ipoig = 1,npoig
      cm    = mmatl(ipoig)
      do id = 1,3
         topol(id+3,ipoig) = topol(id+3,ipoig) +(rhs(id,ipoig)-rhelp(id,ipoig))/cm
      enddo
   enddo
enddo

!       ENDIF
       
!      ------ Store E=1+zx^2	F=zx.zy	G=1+zy^2 positions 4,5 y 6
!             L=zxx/root	M=zxy/root	N=zyy/root	positions 7, 8 y 9
!             root=(1+zx^2+zy^2)**0.5

Do ipoig = 1,npoig
   zx    = topol(2,ipoig)
   zy    = topol(3,ipoig)
   zxx   = topol(4,ipoig)
   zxy   = topol(5,ipoig)
   zyy   = topol(6,ipoig)
   rt    = (1.+zx*zx+zy*zy)**0.5
   topol( 4,ipoig) =(1.+zx*zx)
   topol( 5,ipoig) =(zx*zy)
   topol( 6,ipoig) =(1.+zy*zy)   
   topol( 7,ipoig) = zxx/rt
   topol( 8,ipoig) = zxy/rt
   topol( 9,ipoig) = zyy/rt
   topol(10,ipoig) = rt
   topol(11,ipoig) = zxx
   topol(12,ipoig) = zxy
   topol(13,ipoig) = zyy
Enddo

Deallocate (topel)
deallocate (  rhs)
deallocate (rhelp)
deallocate (mmatl)

END SUBROUTINE Get_Topol


!---------------------------------------------------------- 

subroutine ouplot (TOPO_problem_name, neleg,intmag)

!----------------------------------------------------------

implicit none

integer (ink) neleg, ieleg, ipoig, inodg  
integer (ink) intmag(3,neleg)
integer (ink) final_deposit_top, final_deposit_chk , my_iostat       ! we will read deposit  if a .fdeposit.top file exists 
integer (ink) npoid, ipoid, ipoigx, ipoigy, inx0, iny0, inx1, iny1, inx, iny
real    (irk) xx, yy, zp , time, distx2, disty2, dist , deltg 
real   (irk), allocatable::  coor_poid(:,:), ic10(:)                ! Saeidmt 12Nov2023              real   (irk), allocatable::  coor_poid(:,:), auxd1(:,:), ic10(:) 
real    (irk) Zx, Zy, Zgrad
character(60) Text
character(16) TOPO_problem_name
logical it_exists  
1 format(a60)

!      ---------------  new read and plot deposit extension (final)  ---------MP 30 april 18

!   ---------------  alternative 

final_deposit_top = 83
final_deposit_chk = 84
INQUIRE (final_deposit_top, EXIST = it_exists )
if (it_exists) then
   open( final_deposit_top,  STATUS = 'OLD', file=TOPO_problem_name(1:len1)//'.fdeposit.top', ERR=999, iostat = my_iostat)
   open( final_deposit_chk,                  file=TOPO_problem_name(1:len1)//'.fdeposit.chk')
endif

IF (it_exists)  then ! ------

   read (final_deposit_top,1, ERR=999, iostat = my_iostat) text
   write(final_deposit_chk,1) text
   read (final_deposit_top,*) npoid 
   write(final_deposit_chk,*) npoid 
   if (.NOT.allocated (coor_poid) ) then
      allocate (coor_poid(2,npoid), auxd1(2,npoig), ic10(npoid) )    ! bug in dim was ipoid!!
      coor_poid = 0.0; auxd1 = 0.0; ic10 = 0
   endif
   deltg = (deltxg*deltxg + deltyg*deltyg)**0.5
   read (final_deposit_top,1) text
   write(final_deposit_chk,1) text
   DO ipoid = 1, npoid
      read (final_deposit_top,*) coor_poid(1,ipoid), coor_poid(2,ipoid), ic10(ipoid) 
      write(final_deposit_chk,*) coor_poid(1,ipoid), coor_poid(2,ipoid), ic10(ipoid)
      xx = coor_poid(1,ipoid); yy =coor_poid(2,ipoid)
      ipoigx    = (xx-xming)/deltxg + 1
      ipoigy    = (yy-yming)/deltyg + 1
      inx0      = max(1,ipoigx)
      iny0      = max(1,ipoigy)
      inx1      = min(npoigx,ipoigx)  
      iny1      = min(npoigy,ipoigy+1)
      DO inx    = inx0,inx1
      DO iny    = iny0,iny1
         ipoig  = (iny-1)*npoigx + inx
         distx2 = ( xx - coorg(1,ipoig) )
         disty2 = ( yy - coorg(2,ipoig) )
         distx2 = distx2*distx2
         disty2 = disty2*disty2
         dist   = max( 1.e-6, (distx2+disty2)**0.5 )
         dist   = dist/deltg
         auxd1(1,ipoig) = auxd1(1,ipoig) + ic10(ipoid)/dist                
         auxd1(2,ipoig) = auxd1(2,ipoig) + 1/deltg               
      enddo
      enddo
   ENDDO

   DO ipoig = 1, npoig
      if (auxd1(1,ipoig).ge.1.e-6) then
          auxd1(1,ipoig) = auxd1(1,ipoig)/auxd1(2,ipoig)
      endif
   enddo

   close (final_deposit_top ) ; close (final_deposit_chk)
     
ENDIF   ! -------------------------------------------------------------

999 ipoig = 1

write (Topgid_msh,*) 'MESH SPH_tris    dimension   3    ElemType Triangle  Nnode 3 ' 
write (Topgid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)
   write(Topgid_msh,*) ipoig, xx, yy, zp, time
enddo  

write (Topgid_msh,*) ' End coordinates'

write (Topgid_msh,*) ' Elements'
do ieleg = 1,neleg
   write(Topgid_msh,*) ieleg, (intmag(inodg,ieleg),inodg=1,3),' 1'
enddo
write (Topgid_msh,*) ' End elements' 

close(Topgid_msh) 
       
time = 0.0
write (Topgid_res,*) 'GiD Post Results File 1.0'
      
write(Topgid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
enddo

write(Topgid_res,*) ' End Values '
      
write(Topgid_res,*) 'Result " Z GRAD " "   Zp   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '

do ipoig = 1, npoig 
   Zx = topol (2,ipoig)
   Zy = topol (3,ipoig)
   Zgrad = (Zx * Zx + Zy*Zy)
   if (Zgrad.GT.0.001) then
      Zgrad = Zgrad ** 0.5
   else
      Zgrad = 0.0
   endif
    write(Topgid_res,*)  ipoig,  Zx ,  Zy, Zgrad   ! ipoig, ' 0.   0.   ', Zgrad
enddo

write(Topgid_res,*) ' End Values '

if (my_iostat.eq.0) then         
   write(Topgid_res,*) 'Result " ic deposit " "   icdepos   " ',time,' Vector OnNodes "where" '
   write(Topgid_res,*) ' Values '
   do ipoig = 1, npoig 
       if (auxd1(1,ipoig).GT.1.0) then     !  changed MP May 2020
           auxd1(1,ipoig)=1.0
       endif
       write(Topgid_res,*) ipoig, ' 0.   0.   ', auxd1(1,ipoig)
   enddo
   write(Topgid_res,*) ' End Values '
endif

close(Topgid_res) 

END SUBROUTINE ouplot 


!-------------------------------------------------------------

subroutine ouplot_FDiffs (TOPO_problem_name)            ! Saeidmt 12Nov2023

!----------------------------------------------------------------

!      ------  OUT for topo fimnite diffs
!               

implicit none

integer  (ink) ipoig, ipoigx, ipoigy, ieleg, keleg, plteleg
integer (ink) final_deposit_top, final_deposit_chk , my_iostat       ! we will read deposit  if a .fdeposit.top file exists 
integer (ink) npoid, ipoid, inx0, iny0, inx1, iny1, inx, iny
real    (irk) xx, yy, zp , time, distx2, disty2, dist , deltg 
real   (irk), allocatable::  coor_poid(:,:), ic10(:)                ! Saeidmt 12Nov2023              real   (irk), allocatable::  coor_poid(:,:), auxd1(:,:), ic10(:)
real     (irk) Zpx, Zpy, ZpGrad
character(60) Text
character(16) TOPO_problem_name
1 format(a60)

! Saeidmt 12Nov2023   Whole !      ---------------  new read and plot deposit extension (final)  ---------MP 30 april 18
final_deposit_top = 83
final_deposit_chk = 84
open( final_deposit_top,  STATUS = 'OLD', file=TOPO_problem_name(1:len1)//'.fdeposit.top', ERR=999, iostat = my_iostat)
open( final_deposit_chk,                  file=TOPO_problem_name(1:len1)//'.fdeposit.chk')

IF (my_iostat.eq.0) then ! ------

   read (final_deposit_top,1, ERR=999, iostat = my_iostat) text
   write(final_deposit_chk,1) text
   read (final_deposit_top,*) npoid 
   write(final_deposit_chk,*) npoid 
   if (.NOT.allocated (coor_poid) ) then
      allocate (coor_poid(2,npoid), auxd1(2,npoig), ic10(npoid) )
   endif
   coor_poid = 0
   auxd1 = 0
   ic10 = 0
   deltg = (deltxg*deltxg + deltyg*deltyg)**0.5
   read (final_deposit_top,1) text
   write(final_deposit_chk,1) text
   DO ipoid = 1, npoid
      read (final_deposit_top,*) coor_poid(1,ipoid), coor_poid(2,ipoid), ic10(ipoid) 
      write(final_deposit_chk,*) coor_poid(1,ipoid), coor_poid(2,ipoid), ic10(ipoid)
      xx = coor_poid(1,ipoid); yy =coor_poid(2,ipoid)
      ipoigx    = (xx-Xming_org)/deltxg + 1
      ipoigy    = (yy-Yming_org)/deltyg + 1
      inx0      = max(1,ipoigx)
      iny0      = max(1,ipoigy)
      inx1      = min(npoigx,ipoigx+1)    ! Saeidmt 12Nov2023     Xming_org, Xmaxg_org, Yming_org, Ymaxg_org added to correct reading deposit shape
      iny1      = min(npoigy,ipoigy+1)
      DO inx    = inx0,inx1
      DO iny    = iny0,iny1
         ipoig  = (iny-1)*npoigx + inx
         distx2 = ( xx - coorg(1,ipoig) )
         disty2 = ( yy - coorg(2,ipoig) )
         distx2 = distx2*distx2
         disty2 = disty2*disty2
         dist   = max( 1.e-6, (distx2+disty2)**0.5 )
         dist   = dist/deltg
         auxd1(1,ipoig) = auxd1(1,ipoig) + ic10(ipoid)/dist                
         auxd1(2,ipoig) = auxd1(2,ipoig) + 1/dist               
      enddo
      enddo
   ENDDO

   DO ipoig = 1, npoig
      if (auxd1(1,ipoig).ge.1.e-6) then
          auxd1(1,ipoig) = auxd1(1,ipoig)/auxd1(2,ipoig)
      endif
   enddo

   close (final_deposit_top ) ; close (final_deposit_chk)
     
ENDIF   ! -------------------------------------------------------------

999 ipoig = 1

write (Topgid_msh,*) 'MESH SPH_quads    dimension   3    ElemType Quadrilateral  Nnode 4 ' 
write (Topgid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

time = 0.0

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)
   write(Topgid_msh,*) ipoig, xx, yy, zp, time
enddo  
write (Topgid_msh,*) ' End coordinates'

write (Topgid_msh,*) ' Elements'
do ipoigy = 2,npoigy-2
   do ipoigx  = 2, npoigx-2
      ieleg   = (ipoigy-1)*(npoigx-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoigx + ipoigx
      write (Topgid_msh,*) ieleg,ipoig,ipoig+1,ipoig+npoigx+1,ipoig+npoigx ,  '    1 '
   enddo
enddo
write (Topgid_msh,*) ' End elements' 

close(Topgid_msh) 

time = 0.0

write (Topgid_res,*) 'GiD Post Results File 1.0'
      
write(Topgid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
enddo
write(Topgid_res,*) ' End Values '
      
write(Topgid_res,*) 'Result " Zx topo " "   Zx   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(2,ipoig)
enddo
write(Topgid_res,*) ' End Values '
      
write(Topgid_res,*) 'Result " Zy topo " "   Zy   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(3,ipoig)
enddo
write(Topgid_res,*) ' End Values '
    
write(Topgid_res,*) 'Result " Grad topo " "   Grad  " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig  = 1, npoig 
   Zpx    = topol(2,ipoig)
   Zpy    = topol(2,ipoig)
   ZpGrad = (Zpx*Zpx + Zpy*Zpy)**0.5
   write(Topgid_res,*) ipoig, ' 0.   0.   ', ZpGrad
enddo
write(Topgid_res,*) ' End Values '
      
write(Topgid_res,*) 'Result " Zxx topo " "   Zxx   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(11,ipoig)
enddo
write(Topgid_res,*) ' End Values '
      
write(Topgid_res,*) 'Result " Zxy topo " "   Zxy   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(12,ipoig)
enddo
write(Topgid_res,*) ' End Values '

write(Topgid_res,*) 'Result " Zyy topo " "   Zyy   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(13,ipoig)
enddo
write(Topgid_res,*) ' End Values '

if (my_iostat.eq.0) then                       ! Saeidmt 12Nov2023         
   write(Topgid_res,*) 'Result " ic deposit " "   icdepos   " ',time,' Vector OnNodes "where" '
   write(Topgid_res,*) ' Values '
   do ipoig = 1, npoig 
       write(Topgid_res,*) ipoig, ' 0.   0.   ', auxd1(1,ipoig)
   enddo
   write(Topgid_res,*) ' End Values '
endif

close(Topgid_res)

End Subroutine ouplot_FDiffs

!-------------------------------------------------------------

subroutine ouplot_deposit ( TOPO_problem_name )

!----------------------------------------------------------------

!      ------  OUT for deposit
!               

implicit none

integer (ink) final_deposit_top, final_deposit_chk , my_iostat       ! read deposit  fdeposit.top file exists  
integer (ink) npoid, ipoid, ipoigx, ipoigy 
integer (ink)  inx0, iny0, inx1, iny1, inx, iny
integer (ink)  ieleg, ipoig, inodg 

real   (irk)  xx, yy, zp , deltg  
real    (irk) time, distx2, disty2, dist  
real   (irk), allocatable::  coor_poid(:,:), ic10(:)                ! Saeidmt 12Nov2023              real   (irk), allocatable::  coor_poid(:,:), auxd1(:,:), ic10(:)

character(60) Text
character(16) TOPO_problem_name
1 format(a60)


!      ---------------  new read and plot deposit extension (final)  ---------MP 30 april 18

final_deposit_top = 83
final_deposit_chk = 84
open( final_deposit_top,  STATUS = 'OLD', file=TOPO_problem_name(1:len1)//'.fdeposit.top', ERR=999, iostat = my_iostat)
open( final_deposit_chk,                  file=TOPO_problem_name(1:len1)//'.fdeposit.chk')

IF (my_iostat.eq.0) then ! ------

   read (final_deposit_top,1, ERR=999, iostat = my_iostat) text
   write(final_deposit_chk,1) text
   read (final_deposit_top,*) npoid 
   write(final_deposit_chk,*) npoid 
   if (.NOT.allocated (coor_poid) ) then 
       allocate (coor_poid(2,npoid), auxd1(2,npoig), ic10(npoid) )  ! Saeidmt 12Nov2023
   endif
   deltg = (deltxg*deltxg + deltyg*deltyg)**0.5
   read (final_deposit_top,1) text
   write(final_deposit_chk,1) text
   DO ipoid = 1, npoid
      read (final_deposit_top,*) coor_poid(1,ipoid), coor_poid(2,ipoid), ic10(ipoid) 
      write(final_deposit_chk,*) coor_poid(1,ipoid), coor_poid(2,ipoid), ic10(ipoid)
      xx = coor_poid(1,ipoid); yy =coor_poid(2,ipoid)
      ipoigx    = (xx-xming)/deltxg + 1
      ipoigy    = (yy-yming)/deltyg + 1
      inx0      = max(1,ipoigx)
      iny0      = max(1,ipoigy)
      inx1      = min(npoigx,ipoigx)  
      iny1      = min(npoigy,ipoigy+1)
      DO inx    = inx0,inx1
      DO iny    = iny0,iny1
         ipoig  = (iny-1)*npoigx + inx
         distx2 = ( xx - coorg(1,ipoig) )
         disty2 = ( yy - coorg(2,ipoig) )
         distx2 = distx2*distx2
         disty2 = disty2*disty2
         dist   = max( 1.e-6, (distx2+disty2)**0.5 )
         dist   = dist/deltg
         auxd1(1,ipoig) = auxd1(1,ipoig) + ic10(ipoid)/dist                
         auxd1(2,ipoig) = auxd1(2,ipoig) + 1/deltg               
      enddo
      enddo
   ENDDO

   DO ipoig = 1, npoig
      if (auxd1(1,ipoig).ge.1.e-6) then
          auxd1(1,ipoig) = auxd1(1,ipoig)/auxd1(2,ipoig)
      endif
   enddo

   close (final_deposit_top ) ; close (final_deposit_chk)
     
ENDIF   ! --

999  ipoig = 1

End Subroutine ouplot_deposit

!---------------------------------------------------------- 

subroutine ouplot_BK (neleg,intmag)

!----------------------------------------------------------


integer (ink) neleg, ieleg, ipoig, inodg  
integer (ink) intmag(3,neleg) 
real    (irk) xx, yy, zp , time 
write (Topgid_msh,*) 'MESH SPH_tris    dimension   3    ElemType Triangle  Nnode 3 ' 

write (Topgid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)
   write(Topgid_msh,*) ipoig, xx, yy, zp, time
enddo  

write (Topgid_msh,*) ' End coordinates'

write (Topgid_msh,*) ' Elements'
do ieleg = 1,neleg
   write(Topgid_msh,*) ieleg, (intmag(inodg,ieleg),inodg=1,3),' 1'
enddo
write (Topgid_msh,*) ' End elements' 

close(Topgid_msh) 
       
time = 0.0
write (Topgid_res,*) 'GiD Post Results File 1.0'
      
write(Topgid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
write(Topgid_res,*) ' Values '
do ipoig = 1, npoig 
    write(Topgid_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
enddo

write(Topgid_res,*) ' End Values '

close(Topgid_res)

END SUBROUTINE ouplot_BK 



! --------------------------------------------------------------------------------

SUBROUTINE Get_Z_Topo (ipoin,xx, yy, Zp, Zpx, Zpy, nfric, top4, top5, top6, top7, top8, top9)


! --------------------------------------------------------------------------------

!      Given xx and yy obtains the Z by interpolating in the topo grid

implicit none 

integer(ink) ipoig, ipoigx, ipoigy, ieleg, n1, n2, n3, n4, kpoin
real   (irk) xref, yref, xi, eta, sh1, sh2, sh3, sh4, den

integer(ink), intent(IN)           :: ipoin
real   (irk), intent(IN)           :: xx, yy
integer(ink), intent(IN), OPTIONAL :: nfric
real   (irk),             OPTIONAL :: Zp, Zpx, Zpy, top4, top5, top6, top7, top8, top9

logical   if_Gradient, if_nfric, if_V2 

if_Gradient = PRESENT (Zpx)
if_nfric    = PRESENT (nfric)
if_v2       = PRESENT (top4)

ipoigx = (xx-xming)/deltxg  + 1
ipoigy = (yy-yming)/deltyg  + 1

if (ipoigx.eq.npoigx) then
   if (abs(xmaxg-xx).le.1.e-3) then
      ipoigx = ipoigx -1
   else
      write(*,*) ' point x outside topo mesh', ipoin, xx, xmaxg
   endif
endif   

if (ipoigy.eq.npoigy) then
   if (abs(ymaxg-yy).le.1.e-3) then
      ipoigy = ipoigy -1
   else
      write(*,*) ' point y outside topo mesh', xx, xmaxg
   endif
endif
   
ieleg  = (ipoigy-1)*(npoigx-1) + ipoigx
ipoig  = (ipoigy-1)*npoigx + ipoigx
n1     = ipoig
n2     = ipoig + 1
n3     = ipoig + 1 + npoigx
n4     = n3 -1 
xref   = xming + (ipoigx-1)*deltxg
yref   = yming + (ipoigy-1)*deltyg
xi     = (xx-xref)/deltxg
eta    = (yy-yref)/deltyg
sh1    = (1.-xi) * (1.-eta)
sh2    =     xi  * (1.-eta)
sh3    =     xi  * eta
sh4    = (1.-xi) * eta

!if (xx.lt.0) then
!   write(*,*) 'Nano, ipoin=' , ipoin , 'xx=', xx
!endif
!if (n1.lt.1.OR.n2.lt.1.OR.n3.LT.1.OR.n4.LT.1) then
!   write (*,*) ' ^*******  BEWARE n1 2 3 or 4 less than 1'
!   write (*,*) ' ipoin=',ipoin, ' xx=', xx, 'yy =', yy, n1, n2, n3, n4
!   write (*,*) ' write a number to exit '
!   read  (*,*) kpoin
!endif 
   
Zp     = sh1*topol(1,n1) + sh2*topol(1,n2) + sh3*topol(1,n3) + sh4*topol(1,n4)
kpoin  = poig_in_model(n1) + poig_in_model(n2) + poig_in_model(n3) + poig_in_model(n4) 

if (kpoin.lt.4) THEN   ! MP 27 July 2019
     Zp     = poig_in_model(n1)*sh1*topol(1,n1) + poig_in_model(n2)*sh2*topol(1,n2) &
            + poig_in_model(n3)*sh3*topol(1,n3) + poig_in_model(n4)*sh4*topol(1,n4)
     den = (  poig_in_model(n1)*sh1 + poig_in_model(n2)*sh2 &
            + poig_in_model(n3)*sh3 + poig_in_model(n4)*sh4 )
     if (den.gt.1.e-3) Zp = Zp/den
ENDIF    
!        Zp = -1.e6  ipoin keeps the number of points of this grid element which belongs to computational domain

if (if_Gradient) then
    Zpx    = sh1*topol(2,n1) + sh2*topol(2,n2) + sh3*topol(2,n3) + sh4*topol(2,n4)
    Zpy    = sh1*topol(3,n1) + sh2*topol(3,n2) + sh3*topol(3,n3) + sh4*topol(3,n4) 
endif    

if (if_nfric) then							!      ------  Used for v2/R    
   if (if_V2) then		                    ! 	   nfric.eq.7.or.nfric.eq.8
      top4 = sh1*topol(4,n1) + sh2*topol(4,n2) + sh3*topol(4,n3) + sh4*topol(4,n4)
      top5 = sh1*topol(5,n1) + sh2*topol(5,n2) + sh3*topol(5,n3) + sh4*topol(5,n4)
      top6 = sh1*topol(6,n1) + sh2*topol(6,n2) + sh3*topol(6,n3) + sh4*topol(6,n4)
      top7 = sh1*topol(7,n1) + sh2*topol(7,n2) + sh3*topol(7,n3) + sh4*topol(7,n4)
      top8 = sh1*topol(8,n1) + sh2*topol(8,n2) + sh3*topol(8,n3) + sh4*topol(8,n4)
      top9 = sh1*topol(9,n1) + sh2*topol(9,n2) + sh3*topol(9,n3) + sh4*topol(9,n4)
   endif
endif   

 
END SUBROUTINE Get_Z_Topo


! --------------------------------------------------------------------------------

SUBROUTINE Get_Z_plus_hs (ipoin,xx, yy, Zp,hs_plus_Z_topo)


! --------------------------------------------------------------------------------

!      Given xx and yy obtains the Z by interpolating in the topo grid

implicit none 

integer(ink) ipoig, ipoigx, ipoigy, ieleg, n1, n2, n3, n4, kpoin
real   (irk) xref, yref, xi, eta, sh1, sh2, sh3, sh4

real   (irk) hs_plus_Z_topo(1,npoig)

integer(ink), intent(IN)           :: ipoin
real   (irk), intent(IN)           :: xx, yy
!integer(ink), intent(IN), OPTIONAL :: nfric
real   (irk),             OPTIONAL :: Zp   !, Zpx, Zpy, top4, top5, top6, top7, top8, top9



ipoigx = (xx-xming)/deltxg  + 1
ipoigy = (yy-yming)/deltyg  + 1

if (ipoigx.eq.npoigx) then
   if (abs(xmaxg-xx).le.1.e-3) then
      ipoigx = ipoigx -1
   else
      write(*,*) ' point x outside topo mesh', ipoin, xx, xmaxg
   endif
endif   

if (ipoigy.eq.npoigy) then
   if (abs(ymaxg-yy).le.1.e-3) then
      ipoigy = ipoigy -1
   else
      write(*,*) ' point y outside topo mesh', xx, xmaxg
   endif
endif
   
ieleg  = (ipoigy-1)*(npoigx-1) + ipoigx
ipoig  = (ipoigy-1)*npoigx + ipoigx
n1     = ipoig
n2     = ipoig + 1
n3     = ipoig + 1 + npoigx
n4     = n3 -1 
xref   = xming + (ipoigx-1)*deltxg
yref   = yming + (ipoigy-1)*deltyg
xi     = (xx-xref)/deltxg
eta    = (yy-yref)/deltyg
sh1    = (1.-xi) * (1.-eta)
sh2    =     xi  * (1.-eta)
sh3    =     xi  * eta
sh4    = (1.-xi) * eta

Zp     = sh1*hs_plus_Z_topo(1,n1) + sh2*hs_plus_Z_topo(1,n2) + sh3*hs_plus_Z_topo(1,n3) + sh4*hs_plus_Z_topo(1,n4)



 
END SUBROUTINE Get_Z_plus_hs


! --------------------------------------------------------------------------------

       SUBROUTINE Update_topol  

! --------------------------------------------------------------------------------

!      We have updated Z, but need to update derivatives

implicit none 


integer(ink) neleg, ipoig, ipoigx, ipoigy, ieleg, ieleg2
integer(ink), allocatable:: intmag  (:,:)   ! -- triangle mesh /1  (3*neleg)
real   (irk), allocatable:: geome   (:,:)   ! -- derivatives of shape functions 7*neleg 

!      ------   Now we will obtain the conectivities of the triangles /1 background mesh 

IF (ic_FDiff.eq.0) THEN
   neleg = 2*(npoigx-1)*(npoigy-1)

   if (.NOT.allocated (intmag) ) allocate ( intmag(3,neleg) )
   if (.NOT.allocated (geome)  ) allocate ( geome (7,neleg) )
   geome = 0
   do ipoigy = 1,npoigy-1
      do ipoigx = 1, npoigx-1
         ieleg  =   (ipoigy-1)*(npoigx-1) + ipoigx 
         ieleg2 = 2*ieleg
	     ipoig  = (ipoigy-1)*npoigx + ipoigx
         intmag(1,ieleg2-1) = ipoig   
         intmag(2,ieleg2-1) = ipoig  + 1
         intmag(3,ieleg2-1) = ipoig  + 1 + npoigx
         intmag(1,ieleg2)   = ipoig   
         intmag(2,ieleg2)   = ipoig  + 1 + npoigx
         intmag(3,ieleg2)   = ipoig      + npoigx
      enddo
   enddo

!      ------  Next, we obtain GEOME and TOPOL

   call GetGeom2D (neleg,intmag,geome)

   call Get_Topol (neleg,intmag,geome)

   deallocate (intmag, geome)

ELSEIF (ic_FDiff.eq.10)  THEN

   call FDiff_Derivs

ELSE
   write(*,*) ' !!!!!  wrong value of ic_Fdiff '
   read (*,*) neleg
ENDIF

END SUBROUTINE Update_topol


! --------------------------------------------------------------------------------

       SUBROUTINE Clean_topol  


! --------------------------------------------------------------------------------

!      Deallocate public variables

implicit none 

deallocate (topo_props, poig_in_model, coorg, topol)

END SUBROUTINE Clean_topol



!-------------------------------------------------------------------

       Subroutine Which_Topo_Zone(xx, yy, btype, index_topo_zone) 
        
!-------------------------------------------------------------------  
!
!    Given (xx,yy) of a point, 
!          with basal_type btype (1 friction, 2 erosion, 3 Pwp...)
!    Obtains the topo zone at which it belongs index_topo_zone
!            which is -999 if the points doesn't belong to any
!
!    uses: 
!       basal_type(itopo_zones) type of basal law 1..friction
!                                                 2..erosion
!                                                 3..PwpBCs
!       basal_props(ntopo_props,ntopo_zones)
!                   1.. type of geommetry of basal modif region
!                       1..delimited in Z by Zmin and Zmax (posns 2,3)
!                       2..delimited in X by Xmin and Xmax (posns 2,3)
!                       3..delimited in Y by Ymin and Ymax (posns 2,3)
!                       4..rectangle lim. corners  
!                                 XXmin, XXmax, YYmin, YYmax (2,3,4,5)
!                       5..tilted rectangle defined by  
!                          lower left corner (XXmin, YYmin) (2,3)
!                          sides a & b (lower and left)     (4,5)
!                          angle of a & X                   (6)
!                 

implicit none

integer(ink)  btype, index_topo_zone 
integer(ink)  iexit, itopo_zones, typ_geom, ipoin
integer(ink)  ipoigx, ipoigy, ipoig, n1, n2, n3, n4, npel   ! MP 02 June 
integer(ink)  npmin, npmax, nplist, ip, iplist              ! MP 02 June 

real   (irk)  xx, yy, Zp, Zpx, Zpy
real   (irk)  ZZmin, ZZmax, XXmin, XXmax, YYmin, YYmax
real   (irk)  ZZdiscr, XXdiscr, YYdiscr, ZZdiscr0
real   (irk)  a, b, thetadeg, theta, pi 
real   (irk)  avx, avy, bvx, bvy, Rvx, Rvy, Ra, Rb, costh, sinth
real   (irk)  xeros_min, xeros_max, yeros_min, yeros_max 

integer (ink) npoid, ipoid, inx0, iny0, inx1, iny1, inx, iny   ! Saeidmt 12Nov2023
real    (irk) time, distx2, disty2, dist , deltg               ! Saeidmt 12Nov2023

!      ------  Initialize

iexit = 0; itopo_zones = 0; index_topo_zone = -999; ipoin = 999
pi    = 4.*atan(1.0)

DO WHILE ((iexit.eq.0).and.(itopo_zones.lt.ntopo_zones))
   
   itopo_zones = itopo_zones + 1
   if (basal_type(itopo_zones).eq.btype) then
   
       typ_geom = basal_props(1,itopo_zones) + 0.1
       
       if (typ_geom.le.3) then   
           ZZmin = basal_props(2,itopo_zones)
           ZZmax = basal_props(3,itopo_zones)
           if  (typ_geom.eq.1) then
              call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
           elseif (typ_geom.eq.2) then   
              Zp = xx
           elseif (typ_geom.eq.3) then   
              Zp = yy
           endif 
           ZZdiscr = (ZZmax - Zp)*(Zp-ZZmin) 
       elseif (typ_geom.eq.4) then           ! square Pmin, Pmax
           XXmin = basal_props(2,itopo_zones)
           XXmax = basal_props(3,itopo_zones)
           YYmin = basal_props(4,itopo_zones)
           YYmax = basal_props(5,itopo_zones)
           XXdiscr = (XXmax - xx)*(xx-XXmin)
           YYdiscr = (YYmax - yy)*(yy-YYmin)
           ZZdiscr = -999               ! +ve value means it belongs
           if ((XXdiscr.ge.0).and.(YYdiscr.ge.0)) ZZdiscr = 999 
       elseif (typ_geom.eq.5) then  
           XXmin    = basal_props(2,itopo_zones)
           YYmin    = basal_props(3,itopo_zones)
           a        = basal_props(4,itopo_zones)
           b        = basal_props(5,itopo_zones)
           thetadeg = basal_props(6,itopo_zones)
           theta    = thetadeg*pi/180.
           costh    = cos (theta)
           sinth    = sin (theta)
           Rvx      = (xx - XXmin)
           Rvy      = (yy - YYmin)
           avx      =  costh
           avy      =  sinth
           bvx      = -sinth
           bvy      =  costh
           Ra       = avx*Rvx + avy*Rvy
           Rb       = bvx*Rvx + bvy*Rvy
           ZZdiscr = -999
           if ((Ra.ge.0.0).and.(Ra.le.a).and.(Rb.ge.0.0).and.(Rb.le.b)) ZZdiscr = 999
       elseif (typ_geom.eq.6) then          ! Saeidmt 12Nov2023
           if (.NOT.allocated (auxd1) ) then               
               allocate (auxd1(2,npoig))               
               auxd1 = 0.0           
           endif               
           ipoigx    = (xx-xming)/deltxg + 1
           ipoigy    = (yy-yming)/deltyg + 1
           inx0      = max(1,ipoigx)
           iny0      = max(1,ipoigy)
           ipoig  = (iny0-1)*npoigx + inx0
           ZZdiscr = -999
           if (auxd1(1,ipoig).ge.1.e-6) ZZdiscr = 999
                           
       elseif ( typ_geom.eq.999.or.typ_geom.eq.9999) then ! changed 2016
           ZZdiscr0 = - 999
           xeros_min = basal_props (2,itopo_zones)
           xeros_max = basal_props (3,itopo_zones)
           yeros_min = basal_props (4,itopo_zones)
           yeros_max = basal_props (5,itopo_zones)
           if ( (xx.ge.xeros_min).AND.(xx.le.xeros_max).AND.(yy.ge.yeros_min).AND.(yy.le.yeros_max) ) THEN
               ZZdiscr0 = 999
           endif
           
        elseif ( typ_geom.eq.101.or.typ_geom.eq.102 ) then          ! MP June 2019
           ZZdiscr = - 999
           ipoigx = (xx-xming)/deltxg  + 1
           ipoigy = (yy-yming)/deltyg  + 1

           if (ipoigx.eq.npoigx) then
               if (abs(xmaxg-xx).le.1.e-3) then
                  ipoigx = ipoigx -1
               else
                  write(*,*) ' point x outside topo mesh', ipoin, xx, xmaxg
               endif
           endif   

           if (ipoigy.eq.npoigy) then
              if (abs(ymaxg-yy).le.1.e-3) then
                 ipoigy = ipoigy -1
              else
                 write(*,*) ' point y outside topo mesh', xx, xmaxg
              endif
           endif

           ipoig  = (ipoigy-1)*npoigx + ipoigx
           n1     = ipoig
           n2     = ipoig + 1
           n3     = ipoig + 1 + npoigx
           n4     = n3 -1
           
           npel   = 0
           npmin  = TopoZone_Info (itopo_zones)%npmin
           npmax  = TopoZone_Info (itopo_zones)%npmax
           nplist = TopoZone_Info (itopo_zones)%nplist

           do ip = 1, nplist
              iplist = topozone_info (itopo_zones)%NodeList(ip)  
              if (iplist.eq.n1) npel = npel + 1
              if (iplist.eq.n2) npel = npel + 1
              if (iplist.eq.n3) npel = npel + 1
              if (iplist.eq.n4) npel = npel + 1
           enddo
           
           if (npel.GT.0) then
               ZZdiscr = 999
           endif     
       else
           write(*,*) 'code stopped at Which_Topo_Zone'
           write(*,*) 'because typ_geom not valid'
           PAUSE
       endif
       
       if (ZZdiscr.ge.0) then
          index_topo_zone = itopo_zones
          iexit = 1
       endif
   endif    
ENDDO           
           
END SUBROUTINE   Which_Topo_Zone

!---------------------------------------------------------- 

SUBROUTINE Fdiff_derivs   

!----------------------------------------------------------


!      ------  Compute Zx Zy Zxx Zxy Zyy etc using FDIff formulas
!              domain boundary excluded

implicit none

integer(ink) ipoig, ipoigx, ipoigy, ipoigy0
integer(ink) iA, iAL, iAR, iAT, iAB, iATR, iATL, iABL, iABR
real   (irK) del2x, del2y, delx2, dely2, delxy
real   (irk) Z, Zx, Zy, Zxx, Zxy, Zyy, rt

del2x = 2.*deltxg 
del2y = 2.*deltyg
delx2 = deltxg*deltxg
dely2 = deltyg*deltyg
delxy = deltxg*deltyg   
   
DO ipoigy = 2, npoigy-1

   ipoigy0 = (ipoigy-1)*npoigx
   Do ipoigx = 2, npoigx-1
      
      iA    = ipoigy0 + ipoigx
      ipoig = iA
      
      if (poig_in_model(iA).eq.1) then
      
         iAL  = iA  - 1
         iAR  = iA  + 1
         iAT  = iA  + npoigx
         iATL = iAT - 1
         iATR = iAT + 1
         iAB  = iA  - npoigx
         iABL = iAB - 1
         iABR = iAB + 1
      
         Zx   = (topol(1,iAR)-topol(1,iAL))/del2x
         Zy   = (topol(1,iAT)-topol(1,iAB))/del2y
         rt   = (1.+zx*zx+zy*zy)**0.5
         Zxx  = (topol(1,iAL)-2.*topol(1,iA)+topol(1,iAR))/delx2
         Zyy  = (topol(1,iAB)-2.*topol(1,iA)+topol(1,iAT))/dely2
         Zxy  = (topol(1,iATR)-topol(1,iABR)-topol(1,iATL) +topol(1,iABL))/(4.*delxy)
      
         topol( 2,ipoig)  =  Zx
         topol( 3,ipoig)  =  Zy 
         topol( 4,ipoig)  = (1.+zx*zx)
         topol( 5,ipoig)  = (zx*zy)
         topol( 6,ipoig)  = (1.+zy*zy)   
         topol( 7,ipoig)  = zxx/rt
         topol( 8,ipoig)  = zxy/rt
         topol( 9,ipoig)  = zyy/rt
         topol(10,ipoig)  = rt
         topol(11,ipoig)  = zxx
         topol(12,ipoig)  = zxy
         topol(13,ipoig)  = zyy
      else
         topol (2:13,ipoig) = 0.0
      endif
   enddo   
enddo  

END SUBROUTINE Fdiff_derivs  



!---------------------------------------------------------- 

SUBROUTINE Get_Topo_restart  (ic_erosion) 

!----------------------------------------------------------

!      ------  read from file restart.top
!              topo 1, 14 and 15
!              characterizing current topo, eroded and available

implicit none

integer (ink) In_Source, Chk_Source, ic_erosion, i, ip
real    (irk) time_aux
character(60) text 

1005 format (a60)

! time_sph_restart  is time at which computation is restarted 

In_Source    = 86
Chk_Source   = 88
open (  In_Source ,  file = 'restart.top'   )
open (  Chk_Source,  file = 'restart.top.chk'   )

read ( In_Source,  1005) text 
write( Chk_Source, 1005) text
read ( In_Source,  *)  npoig, time_aux     
write( Chk_Source, *)  npoig, time_aux
read ( In_Source,  1005) text
write( Chk_Source, 1005) text 
DO i = 1, npoig 
   read ( In_Source,*)  ip, Topol (1,i) , Topol (14,i) , Topol (15,i) 
   write( Chk_Source,*) ip, Topol (1,i) , Topol (14,i) , Topol (15,i) 
ENDDO

close ( In_Source  )
close ( Chk_Source )

END SUBROUTINE Get_Topo_restart



!---------------------------------------------------------- 

SUBROUTINE Get_obstacles_old  (ic_obst, ictop) 

!----------------------------------------------------------

!      ------  adds obstacles to topol(1,:)

implicit none

integer (ink)  ic_obst, ictop
integer (ink)  ipx0, ipx1, ipy0, ipy1, ipx, ipy, ipoig  
real    (irk)  X0_Tx, Xf_Tx, Lambda_Tx, A0_Tx
real    (irk)  x0,x1,y0,y1, wlength, a0, th_deg, Th_rad
real    (irk)  xi, eta, xx, yy, txepa
character(60) text 

if (ictop.eq.11.OR.ictop.EQ.21) then
    x0 = xshiftg
endif

if (ic_obst.EQ.1) then   
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    read  (topdat_file,*) X0_Tx, Xf_Tx, Lambda_Tx, A0_Tx
    write (Topchk_file,*) X0_Tx, Xf_Tx, Lambda_Tx, A0_Tx
    do ipoig = 1,npoig
       xx = coorg(1,ipoig)
       txepa = 0.0
       if ( xx.ge.X0_Tx.AND.xx.le.Xf_Tx) then
          txepa = A0_Tx*sin( 2*3.141592*(xx-X0_Tx)/lambda_Tx)
          if (txepa.lt.0.0) txepa = 0.
       endif
       topol(1,ipoig) = topol(1,ipoig) + txepa
    enddo
elseif (ic_obst.EQ.2) then
   read  (topdat_file,1) text
   write (Topchk_file,1) text
   read  (topdat_file,*) X0, X1, Y0, Y1, wlength, a0, Th_deg
   write (Topchk_file,*) X0, X1, Y0, Y1, wlength, a0, Th_deg
   Th_rad = Th_deg * 3.141592/180.
   ipx0   = (x0-xming)/deltxg 
   ipx1   = (x1-xming)/deltxg + 1
   ipy0   = (y0-yming)/deltyg
   ipy1   = (y1-yming)/deltyg + 1
   if (ipx0.le.0) ipx0 = 1
   if (ipy0.le.0) ipy0 = 1
   if (ipx1.gt.npoigx) ipx1 = npoigx
   if (ipy1.gt.npoigy) ipy1 = npoigy
   DO ipy = ipy0, ipy1
      yy  = yming + (ipy-1)*deltyg
      DO ipx = ipx0, ipx1
         xx  = xming + (ipx-1)*deltxg
         xi  =  xx*cos(th_rad) + yy*sin(th_rad)
         eta = -xx*sin(th_rad) + yy*cos(th_rad)
         txepa = a0* sin((2*3.141592*xi)/wlength)*sin((2*3.141592*eta)/wlength)
         if (txepa.LT.0.0) txepa = 0.0
         ipoig = (ipy-1)*npoigx + ipx
         topol(1,ipoig) = topol(1,ipoig) + txepa
      enddo
   enddo

endif
1 format (a60) 

END SUBROUTINE Get_obstacles_old

!---------------------------------------------------------- 

SUBROUTINE Get_obstacles  (ic_obst, ictop) 

!----------------------------------------------------------

!      ------  adds obstacles to topol(1,:)

implicit none

integer (ink)  ic_obst, ictop
integer (ink)  ipx0, ipx1, ipy0, ipy1, ipx, ipy, ipoig 
integer (ink)  icx, icy, ihx, ihy, ihh
real    (irk)  cx, hx, ex, cy, hy, ey, remx, remy
real    (irk)  X0_Tx, Xf_Tx, Lambda_Tx, A0_Tx
real    (irk)  x0,x1,y0,y1, wlength, a0, th_deg, Th_rad
real    (irk)  xi, eta, xx, yy, txepa
real    (irk)  xQL, yQL, xQR, yQR, cosTh, sinTh, aL, bL
character(60) text 

if (ictop.eq.11.OR.ictop.EQ.21) then
    x0 = xshiftg
endif

if (ic_obst.EQ.1) then   
    read  (topdat_file,1) text
    write (Topchk_file,1) text
    read  (topdat_file,*) X0_Tx, Xf_Tx, Lambda_Tx, A0_Tx
    write (Topchk_file,*) X0_Tx, Xf_Tx, Lambda_Tx, A0_Tx
    do ipoig = 1,npoig
       xx = coorg(1,ipoig)
       txepa = 0.0
       if ( xx.ge.X0_Tx.AND.xx.le.Xf_Tx) then
          txepa = A0_Tx*sin( 2*3.141592*(xx-X0_Tx)/lambda_Tx)
          if (txepa.lt.0.0) txepa = 0.
       endif
       topol(1,ipoig) = topol(1,ipoig) + txepa
    enddo
elseif (ic_obst.EQ.2) then
   read  (topdat_file,1) text
   write (Topchk_file,1) text
   read  (topdat_file,*) X0, X1, Y0, Y1, wlength, a0, Th_deg
   write (Topchk_file,*) X0, X1, Y0, Y1, wlength, a0, Th_deg
   Th_rad = Th_deg * 3.141592/180.
   ipx0   = (x0-xming)/deltxg 
   ipx1   = (x1-xming)/deltxg + 1
   ipy0   = (y0-yming)/deltyg
   ipy1   = (y1-yming)/deltyg + 1
   if (ipx0.le.0) ipx0 = 1
   if (ipy0.le.0) ipy0 = 1
   if (ipx1.gt.npoigx) ipx1 = npoigx
   if (ipy1.gt.npoigy) ipy1 = npoigy
   DO ipy = ipy0, ipy1
      yy  = yming + (ipy-1)*deltyg
      DO ipx = ipx0, ipx1
         xx  = xming + (ipx-1)*deltxg
         xi  =  xx*cos(th_rad) + yy*sin(th_rad)
         eta = -xx*sin(th_rad) + yy*cos(th_rad)
         txepa = a0* sin((2*3.141592*xi)/wlength)*sin((2*3.141592*eta)/wlength)
         if (txepa.LT.0.0) txepa = 0.0
         ipoig = (ipy-1)*npoigx + ipx
         topol(1,ipoig) = topol(1,ipoig) + txepa
      enddo
   enddo
elseif (ic_obst.EQ.3) then
   read  (topdat_file,1) text
   write (Topchk_file,1) text
   read  (topdat_file,*) X0, X1, Y0, Y1, hx, ex, hy, ey, a0, Th_deg
   write (Topchk_file,*) X0, X1, Y0, Y1, hx, ex, hy, ey, a0, Th_deg
   cx = hx + ex; cy = hy +  ey 
   Th_rad = Th_deg * 3.141592/180.
   ipx0   = (x0-xming)/deltxg 
   ipx1   = (x1-xming)/deltxg + 1
   ipy0   = (y0-yming)/deltyg
   ipy1   = (y1-yming)/deltyg + 1
   if (ipx0.le.0) ipx0 = 1
   if (ipy0.le.0) ipy0 = 1
   if (ipx1.gt.npoigx) ipx1 = npoigx
   if (ipy1.gt.npoigy) ipy1 = npoigy
   DO ipy = ipy0, ipy1
      yy  = yming + (ipy-1)*deltyg
      eta = -xx*sin(th_rad) + yy*cos(th_rad)
      icy  = (eta-Y0)/cy
      remy = (eta-Y0) - icy*cy
      ihy  = remy/hy
      DO ipx = ipx0, ipx1
         xx  = xming + (ipx-1)*deltxg
         xi  =  xx*cos(th_rad) + yy*sin(th_rad)  ! only Theta =0 
         icx  = (xi-x0)/cx
         remx = (xi-x0) - icx*cx
         ihx  = remx/hx
         ihh  = ihx + ihy
         if (ihh.eq.0) then
            txepa = a0  
            ipoig = (ipy-1)*npoigx + ipx
            topol(1,ipoig) = topol(1,ipoig) + txepa
         endif
      enddo
   enddo

elseif (ic_obst.EQ.4) then
   read  (topdat_file,1) text
   write (Topchk_file,1) text
   read  (topdat_file,*) X0, Y0, aL, bL, hx, ex, hy, ey, a0, Th_deg
   write (Topchk_file,*) X0, Y0, aL, bL, hx, ex, hy, ey, a0, Th_deg
   cx = hx + ex; cy = hy +  ey 
   Th_rad = Th_deg * 3.141592/180.
   cosTh = cos (Th_rad); sinTh = sin (Th_rad)
   
   !      ------  Task (i) obtain left down and rif¡ght upper corners of the houses domain 
  
   xQL = X0 - bL*sinTh
   yQL = Y0  
   xQR = X0 + aL*cosTh
   yQR = Y0 + aL*sinTh + bL*cosTh
                      !      associated indexes in topo grid
   ipx0   = (xQL-xming)/deltxg 
   ipy0   = (yQL-yming)/deltyg
   ipx1   = (xQR-xming)/deltxg + 1
   ipy1   = (yQR-yming)/deltyg + 1
   if (ipx0.le.0) ipx0 = 1
   if (ipy0.le.0) ipy0 = 1
   if (ipx1.gt.npoigx) ipx1 = npoigx
   if (ipy1.gt.npoigy) ipy1 = npoigy
   
   !      ------  Task (ii) Loop around x,y topo grid nodes
   
   DO ipy = ipy0, ipy1
      yy  = yming + (ipy-1)*deltyg
      DO ipx = ipx0, ipx1
         xx  = xming + (ipx-1)*deltxg  

   !      ------  Task (iii) Obtain relative coords of all xx, yy -> xi, eta

         xi  =  (xx-X0) * cosTh + (yy-y0) * sinTh
         eta = -(xx-X0) * sinTh + (yy-y0) * cosTh
   !             skip points outside rectangle
         if ( .NOT. ( (xi.GE.0).AND.(xi.LE.aL).and.(eta.GE.0).and.(eta.LE.bL) )  ) CYCLE
         icx  = xi / cx
         remx = xi - icx*cx
         ihx  = remx/hx
         icy  = eta / cy
         remy = eta - icy*cy
         ihy  = remy/hy
         ihh  = ihx + ihy
         if (ihh.eq.0) then
            txepa = a0  
            ipoig = (ipy-1)*npoigx + ipx
            topol(1,ipoig) = topol(1,ipoig) + txepa
         endif
      ENDDO
   ENDDO

ENDIF   !  end of case 4
1 format (a60) 

END SUBROUTINE Get_obstacles 




!---------------------------------------------------------- 

SUBROUTINE Get_Geo_Zones

!----------------------------------------------------------

!      ------  Obtains geo zones (only convex)

implicit none

! integer (ink) n_GZones   !  temporal for dbuggin will come from driver vars
integer (ink) np_max, np_zone 
integer (ink) inx0, inx1,iny0, iny1 , ip1, ip2
integer (ink) i, izone, ix, iy, ipz, ipoig
real    (irk) xgzmin, xgzmax, ygzmin, ygzmax  
real    (irk) x1, y1, x2, y2, xx, yy, discr
real   (irk), allocatable:: xp_Zone  (:), yp_Zone(:)     ! -- x's and y's for boundary
character(60) text 

1 format (a60) 

topol (16,:) = 1    ! Sets default to zone 1
np_max       = 100  ! max nr of points defining each zone         
allocate ( xp_Zone(np_max), yp_Zone(np_max) )
! n_GZones = 1  ACTIVATE FOR TEST Gz

DO izone = 1, n_GZones-1 
   read  (topdat_file,1) text
   write (Topchk_file,1) text 
   read  (topdat_file,*) np_zone
   write (Topchk_file,*) np_zone
   read  (topdat_file,1) text
   write (Topchk_file,1) text 
   read  (topdat_file,*) (xp_Zone(i), i=1, np_zone)
   read  (topdat_file,*) (yp_Zone(i), i=1, np_zone)
   xgzmin = 1.e10; xgzmax = -1.e10;  ygzmin = 1.e10;  ygzmax = -1.e10
   do i  = 1, np_zone
      if (xp_Zone(i).LT.xgzmin) xgzmin = xp_Zone(i)
      if (xp_Zone(i).GT.xgzmax) xgzmax = xp_Zone(i)
      if (yp_Zone(i).LT.ygzmin) ygzmin = yp_Zone(i)
      if (yp_Zone(i).GT.ygzmax) ygzmax = yp_Zone(i)  !  bug mp 9may 2022
   enddo
   inx0 = int ( (xgzmin -  xming) / deltxg + 0.001) + 1 
   iny0 = int ( (ygzmin -  yming) / deltyg + 0.001) + 1 
   inx1 = int ( (xgzmax -  xming) / deltxg + 0.001) + 1 ; if (inx1.GT.npoigx) inx1 = npoigx
   iny1 = int ( (ygzmax -  yming) / deltyg + 0.001) + 1 ; if (iny1.GT.npoigx) iny1 = npoigy     
   DO ix = inx0, inx1
   Do iy = iny0, iny1
      ipoig = (iy-1)*npoigx + ix
      xx    = coorg (1, ipoig) ; yy    = coorg (2, ipoig)
      Do ipz = 1, np_Zone
         ip1 = ipz; ip2 = ip1 +1 ; if (ip2.GT.np_Zone) ip2 = 1
         x1  = xp_Zone (ip1); x2 = xp_Zone (ip2)
         y1  = yp_Zone (ip1); y2 = yp_Zone (ip2)
         discr = (x1 - xx)*(y2 - yy) -(x2 - xx)*(y1 - yy)  ! Area*2, positive if xx inside
         if (discr.LT.-0.001) EXIT
      enddo
      if (discr.GE.-0.001) then
         topol (16, ipoig) = izone +1 
      endif
   enddo
   enddo
ENDDO


END SUBROUTINE Get_Geo_Zones


!---------------------------------------------------------- 

SUBROUTINE Get_Geo_Zones_NX

!----------------------------------------------------------

!      ------  adds obstacles to topol(1,:)

implicit none

integer (ink) np_max, np_zone 
integer(ink)  ip, ip1    
integer (ink) i, izone, ix, iy, ipz, ipoig
integer (ink) inx0, iny0, inx1, iny1

real   (irk), allocatable:: xp_Zone  (:), yp_Zone(:)     ! -- x's and y's for boundary
real   (irk)  Pi, x0,y0, x1, y1, x2, y2
real   (irk)  x01, y01, x02, y02, vprod, v01, v02, Vesc 
real   (irk)  DTh, Th, sindTh, cosdTh  

real    (irk) xgzmin, xgzmax, ygzmin, ygzmax 

character(60) text 

1 format (a60) 

Pi = 4.0*atan(1.00) 

topol (16,:) = 1    ! Sets default to zone 1
np_max       = 100  ! max nr of points defining each zone         
if (.NOT.allocated(xp_Zone)) then
    allocate ( xp_Zone(np_max), yp_Zone(np_max) )
endif

DO izone = 1, n_GZones-1 

   xgzmin =  1.e10 ; ygzmin=  1.e10  !  Saeid April 2024
   xgzmax = -1.e10 ; ygzmax= -1.e10  

   read  (topdat_file,1) text
   write (Topchk_file,1) text 
   read  (topdat_file,*) np_zone
   write (Topchk_file,*) np_zone
   read  (topdat_file,1) text
   write (Topchk_file,1) text 
   if (ic_topo_props.eq.-101) then
      do i = 1, np_zone                                !   MP & Saeid April 2024
         read  (topdat_file,*)  xp_Zone(i), yp_Zone(i) 
         write (Topchk_file,*)  xp_Zone(i), yp_Zone(i) 
      enddo
   elseif (ic_topo_props.eq.101) then
     read  (topdat_file,*) (xp_Zone(i), i=1, np_zone)
     read  (topdat_file,*) (yp_Zone(i), i=1, np_zone)
   endif
   do i  = 1, np_zone
      if (xp_Zone(i).LT.xgzmin) xgzmin = xp_Zone(i)
      if (xp_Zone(i).GT.xgzmax) xgzmax = xp_Zone(i)
      if (yp_Zone(i).LT.ygzmin) ygzmin = yp_Zone(i)
      if (yp_Zone(i).GT.ygzmax) ygzmax = yp_Zone(i)  !  bug mp 9may 2022
   enddo

   inx0 = int ( (xgzmin -  xming) / deltxg + 0.001) + 1   ! rectangle indexes containing region
   iny0 = int ( (ygzmin -  yming) / deltyg + 0.001) + 1 
   inx1 = int ( (xgzmax -  xming) / deltxg + 0.001) + 1 ; if (inx1.GT.npoigx) inx1 = npoigx
   iny1 = int ( (ygzmax -  yming) / deltyg + 0.001) + 1 ; if (iny1.GT.npoigx) iny1 = npoigy     
   DO ix = inx0, inx1
   Do iy = iny0, iny1
      ipoig  = (iy-1)*npoigx + ix              ! candidate point to be INside
      x0     = coorg (1, ipoig) ; y0    = coorg (2, ipoig)
      Th     = 0.0;  dTh  =  0.0  
      DO ip  = 1, np_Zone  !  loop about points of the boundary, Vectors 01 and 02
         ip1 = ip + 1
         if (ip1.gt.np_Zone) ip1 = 1
         x01 = xp_Zone(ip)  - x0
         y01 = Yp_Zone(ip)  - y0
         x02 = xp_Zone(ip1) - x0
         y02 = Yp_Zone(ip1) - y0
         v01   = sqrt( x01*x01 + y01*y01)
         v02   = sqrt( x02*x02 + y02*y02)
         vprod = x01*y02 - x02*y01
         Vesc = x01*x02 + y01*y02
         sindTh = Vprod/(v01*v02)
         cosdTh = Vesc/(v01*v02)
         dTh    = ATAN2 (sindTh,cosdTh)
         Th = Th + Dth
      enddo 
      if (Abs(Th).LT.0.01) then
!         write(*,*) ' x is OUTside'
      elseif ( (Abs(Th)-2*pi).LT.0.1) then
         topol (16, ipoig) = izone +1
      else
         write(*,*) ' x is not OUTside nor INside  x,y, ipoig, Th', x0, y0,ipoig, Th
         write(*,*) ' check data ...code will pause '
         write(*,*) ' type 1 to exit '
         STOP
      endif
   enddo
   enddo  ! ...... enddo candidate ipoigs in zone

ENDDO  !....  enddo zones

deallocate ( xp_Zone, yp_Zone )

END SUBROUTINE Get_Geo_Zones_NX

END MODULE SPH_FEM_topo_2019