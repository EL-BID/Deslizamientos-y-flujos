       
      
!      -------------------------------------------------------

                MODULE GFL_MAIN_2019

!
!               Geoflow_Main  October 2008
!               Geoflow 2005 aiguablava!      
!               2005  MP for aiguablava. Add ictop=-5, and icunk=-6
!       
!
!      
!
!      -------------------------------------------------------

!      changed npoin to npoinF    and ndimn to ndimnF (sph compatibility)           

USE SPH_FEM_variable_types_2019   !  Variable precision defined here.

USE SPH_FEM_Driver_time_vars_2019 !  Time variables common to sph and geoflow

USE GFL_Global_Vars_2019          !  Main SPH variables are defined here

USE SPH_FEM_topo_2019             !  getinput uses this module to import a DTM

USE SPH_GFL_SW_Interactions_2019, & 

    ONLY: hs_w,              &  !  we use hs_w to modify HH and pres term 
          gradhs_w,          &  !  we use only gradhs in nodes and elems
          gradhs_w_elem,     &  !     in sbroutines 1stStep, 2ndStep and sources
          WIR_Interact,      &  !  We skip interactions when .FALSE:
          n_ggates,          &  !  nr of gates to provide input to sph
          mps_ggates,        &
          nps_ggates,        &
          list_pts_ggates,   &
          
          type_of_fem_sph_interaction
                 
                                !  This module provides interaction whith sph
                                !  Uses: gradhs_w (npoinF) grad of soil at fe nodes
                                !        WIR_Interact .TRUE. when there is interaction
                                !  Routines are: 1st step, 2nd step and sources
     
PRIVATE                         !      ------  No variable will be accesible from outside 
                                !              the module unless declared public (none indeed)
                                

public:: Init_GFL               !  called from Driver 
public:: Time_Integration_GFL
public:: Adaptive_dt_GFL

public:: Out_Plot_GFL

public:: Out_print_GFL

! public:: clean_up_sph  ** closes files, deallocates all 


CONTAINS
  

!-------------------------------------------------------------------

       SUBROUTINE Init_GFL              !  REVISED

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none

integer (ink) cont       !   characters counter  
        
nsour = GFL_t_Integ_Alg
if(nsour.ne.0.and.nsour.ne.2) then
   write(*,*)'Error selecting  integration of the source term'
   PAUSE
endif

!write(*,*) 'Enter 1 for advected pwp, 0 otherwise'
!read (*,*) icpwpF
       
!      ------  Open files
         
cont = len_trim(GFL_SW_problem_name)

datg_file  =  5
chkg_file  =  6 
resg_file1 =  8  
gidg_msh   = 15  
gidg_res   = 16  
gid3g_msh  = 17 
gid3g_res  = 18  

open (unit= datg_file,  file = GFL_SW_problem_name(1:cont)//'.dat'         )
open (unit= chkg_file,  file = GFL_SW_problem_name(1:cont)//'.echo'        )
open (unit= resg_file1, file = GFL_SW_problem_name(1:cont)//'.res'         )
open (unit= gidg_msh,   file = GFL_SW_problem_name(1:cont)//'.post.msh'  )
open (unit= gidg_res,   file = GFL_SW_problem_name(1:cont)//'.post.res'  )
!open (unit= gid3g_msh, file = GFL_SW_problem_name(1:cont)//'3d.flavia.bon')
!open (unit= gid3g_res, file = GFL_SW_problem_name(1:cont)//'3d.flavia.res')

!      ------  Read basic parameters

call in0tg1d  


!      ------  Allocate space for arrays and matrices

allocate ( coord (ndimnF,npoinF))
allocate ( intmat (nnode,nelem))
allocate ( geome (ngeom,nelem) )
allocate ( topol (ntopo,npoinF))
allocate ( topel (ntopo,nelem) )

allocate ( intmab (nnode,nelemb))
allocate ( geomb  (ngeob,nelemb))

allocate ( unkno (namat,npoinF))
allocate ( unkne (namat,ntotg) )
allocate ( delun (namat,npoinF))

allocate ( fluxp (namat,npoinF))
allocate ( fluyp (namat,npoinF))
allocate ( gpoin (namat,npoinF))
allocate ( gelem (namat,nelem) )
allocate ( rhs0  (namat,npoinF))
allocate ( rhs1  (namat,npoinF))
allocate ( rhs2  (namat,npoinF))
allocate ( mmatl ( npoinF)     )
allocate ( deltel (nelem)      )
allocate ( iactnd (npoinF)     )
allocate ( iactel (nelem)      )

allocate ( unpre  (nunpre,nprer) )
allocate ( rtanr  (ndimnF,nbour) )
allocate ( rtang  (ndimnF,npoinF))
allocate ( btanr  (nbour)        )
allocate ( bconl  (npoinF)       )
allocate ( ltang  (npoinF)       )
allocate ( bprer  (2,nprer)      )
allocate ( babso  (nabso)        )
allocate ( babhz  (nabshz)       )
allocate ( bload  (nload)        )
allocate ( lprpoi (nprpoi)       )

allocate ( bindex (nelemb)      )
allocate ( denom  (nelem)       )

allocate ( constF (nconstF)      )
allocate ( patma  (npatma,npoinF))
allocate ( patmel (npatma,nelem) )
allocate ( twind  (ntwind)       )
allocate ( datm   (ndatm )       )
allocate ( wave   (nwave )       )
allocate ( lweir  (nweir)        )
allocate ( zweir  (nweir)        )

coord  = 0.0
topol  = 0.0
topel  = 0.0
geome  = 0.0
mmatl  = 0.0
geomb  = 0.0
rtanr  = 0.0
rtang  = 0.0
unkno  = 0.0
unkne  = 0.0
unpre  = 0.0
delun  = 0.0
rhs0   = 0.0
rhs1   = 0.0
rhs2   = 0.0
gelem  = 0.0
gpoin  = 0.0
fluxp  = 0.0
fluyp  = 0.0
patma  = 0.0
patmel = 0.0
constF = 0.0
twind  = 0.0
datm   = 0.0
wave   = 0.0
denom  = 0.0

iactnd = 1              !     --  We start with all nodes set to active mode
iactel = 1              !     --  and all elements activated too
bindex = 0

call getinput           !     --  read  all input data

call getgeom2D          !     --  obtain the geometrical parameters needed

call getmmatl           !     --  fill the lumped mass mx : mmatl

call Get_Topol          !     --  Obtain gradients of topography at nodes

call get_boun_nodes_2D  !     --  find out bound. nodes and normalized tangents

call bounel             !     --  find boundary elements and form conectivities

call output             !     --  Output initial state (file prn)

call ouplot_mesh        !     --  Plot mesh and heading of res  
call ouplot_res         !         Plot Initial state
call Out_print_GFL              !         Print initial state
         
END  SUBROUTINE Init_GFL



!-------------------------------------------------------------------

       SUBROUTINE Time_Integration_GFL

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none

integer(ink) n1_RK              !      ------  used as index of calls at RK  
real   (irk) gratm

delun = 0.0
rhs1  = 0.0
n1_RK    = 0

call  delun0       !      ------  get delun for the prescribed nodes

!      ------  evaluate atm.pressur gradient

gratm = datm(4)
if(abs(gratm).ge.0.0001) then
    call getpatm
endif

!   --------- sources term (Runge-Kutta)    

if(ntype.eq.3.and.nsour.eq.2)then
   call sources_RK (n1_RK)
endif

!   ---------

call first_step  

call second_step 

call solve

if (csmooth.gt.0)  call smooth  !      ------   smooth shocks

call  addthem                   !     ------  add the incr to  old solution

!   --------- sources term (Runge-Kutta)

if(ntype.eq.3.and.nsour.eq.2)then
   call sources_RK (n1_RK) 
endif

!      ------  check there is not any h+eta' < 0.0, if so : h+eta=0.0 !!

call check  !MP new check for -ve pwp

!      ------   apply absorbent and incident wave conditions

IF (nabso.ne.0)   call abshz    !      ------  Abs BCs
          
IF (nload.ne.0)   call loadbc   !      ------  Incident wave 

IF (nweir.ne.0)   call weir     !      ------  apply weir conditions

!      ------  check again for any h+eta' < 0.0, if so : h+eta=0.0 !!

call check 


END  SUBROUTINE Time_Integration_GFL


!-------------------------------------------------------------------

       SUBROUTINE Adaptive_dt_GFL

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none

call       check_delt 

END SUBROUTINE Adaptive_dt_GFL

!-------------------------------------------------------------------

       SUBROUTINE Out_PlotMesh_GFL

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none

call  ouplot_mesh
             
END SUBROUTINE Out_PlotMesh_GFL


!-------------------------------------------------------------------

       SUBROUTINE Out_Plot_GFL

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none

call  ouplot_res
             
END SUBROUTINE Out_Plot_GFL


!-------------------------------------------------------------------

       SUBROUTINE Out_Print_GFL

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none
 
call output                             
             
END SUBROUTINE Out_Print_GFL


!======================      NEXT ROUTINES ARE CALLED FROM ABOVE ONES   ======================



!----------------------------------------------------------

        subroutine in0tg1d               
                    
!----------------------------------------------------------


implicit none

integer(ink) nline, iline
character*72 text

!      Read number of text lines and lines

       read (datg_file,*) nline
       write(chkg_file,*) nline
       do iline = 1,nline
          read (datg_file,1) text
          write(chkg_file,1) text
       enddo
 1     format (a72)

!      ------  nprpoi is the nr of control nodes for output
!              icoutp  is 1 for classical SW (elevation over meanWL)

       read (datg_file,1) text
       write(chkg_file,1) text
       read (datg_file,*) nelem, npoinF, nbour, nprer, nload, nabso, ntype, &
                  nprpoi,icoutp,nweir,nabshz,ngaus, nacti
       write(chkg_file,*) nelem, npoinF, nbour, nprer, nload, nabso, ntype, &
                  nprpoi,icoutp,nweir,nabshz,ngaus, nacti

       if(ngaus.ne.1.and.ngaus.ne.3)then
          write(*,*) 'Error writing ngaus: just 1 or 3 allowed'
          stop
       end if

!      1. Zcomp When h<Zcomp we make h=0 

       read (datg_file,1) text
       write(chkg_file,1) text
       read (datg_file,*) icPwpf, comp    ! Icpwp for Pwp (advected)    
       write(chkg_file,*) icPwpf, comp    !
       if(comp.lt.0.)then
          write(*,*)' Error entering Zcomp: must be greater than 0.'
          stop
       endif
        
       nnode     = 3
       ndimnF    = 2
       ntopo     = 9       ! MP2 z zx zy L(zxx) M(zxy) N(zyy) E F G
       nconstF   = 13      ! MP
       namat     = 3
       if (icpwpF.eq.1) then
           namat  = 4             ! MP
           nconstF= 14
       endif
       ngeom     = 7
       ngeob     = 3
       nelemb    = nelem + 2
       npatma    = 3
       ntwind    = 3
       ndatm     = 10
       nwave     = 10
       nunpre    = 5
       ntotg     = nelem*ngaus

        
END  SUBROUTINE in0tg1d

! ---------------------------------------------------------

        subroutine getinput

! --------------------------------------------------------


implicit none

integer(ink) i, j, k, ielem, ipoin, ip, iptsZ, iconst
integer(ink) nptsZ, ictop, icunkno, nfric
integer(ink) nptsh, iptsh
integer(ink) idum3,idum4,idum5,idum6
integer(ink) icdry, icpwpo, iwv

real   (irk) Zconst, twopi, dum1,dum2
real   (irk) xx, x1, f0, f1, yy, zz, Zp          
real   (irk) x0, y0, xs, ys
real   (irk) Rout,Depthout,Rint,Depthint, Xl0, depth2, TolerR
real   (irk) SWL0, Ref_depth, x, y, rx, ry, r, vnx, vny, vmod
real   (irk) dpth ,dot, dot2, dZ, eta0
                                         
real   (irk) xcenter,xlong,deltah, h,xhmin, xhmax       
real   (irk) a0
real   (irk) h0, dx, dy, rr
real   (irk) y1,b,d0, xc, yc, a, vnxi, vneta
real   (irk) vvxi, vveta, xi, eta, dh, hh 
real   (irk) r0, cota, r2, xd, yd, dist2
real   (irk) ph0, ph1, phh, xn0, xn1
real   (irk) yn0, yn1, h1

real   (irk), allocatable:: xptsZ  (:)      ! -- x's for Z segments
real   (irk), allocatable:: zptsZ  (:)      ! -- z's
real   (irk), allocatable:: xptsh  (:)      ! -- x's for def. of h
real   (irk), allocatable:: hptsh  (:)      ! -- h's

! 
       character*80 text
! 
! 
!      ----this subroutine reads all the input data
!           data is read in the following manner :
! 
!      1. Zcomp from keyboard. When h<Zcomp we make h=0
!      2. Interconnectivity matrix intmat
!      3. Topography or Bathymetry  plus nodal coordinates
!      4. Initial values of unknowns
!      5. Initial pore pressures pwp0, if we have activated this option
!      6. Boundary Conditions:
!               6.a Tangent vectors
!               6.b Prescribed unknowns
!               6.c Absorbing boundaries
!               6.d Loaded BCs (Incident wave)
!      7. Constants or material properties
!      8. Niter and csmooth
!      9. Output nodes
!     10. Wind stress
!     11. Atmospheric pressure info (if any)

! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

!      1. Zcomp When h<Zcomp we make h=0 moved up

!read (datg_file,1) text
!write(chkg_file,1) text
!read (datg_file,*) icPwpf, comp    ! Icpwp for Pwp (advected)    
!write(chkg_file,*) icPwpf, comp    !
!if(comp.lt.0.)then
   !write(*,*)' Error entering Zcomp: must be greater than 0.'
   !stop
!endif

!      2. Interconnectivity matrix intmat           

read (datg_file,1) text
write(chkg_file,1) text
do  i = 1,nelem
    read (datg_file,*) ielem,(intmat(j,i),j=1,nnode)
    write(chkg_file,5) ielem,(intmat(j,i),j=1,nnode)
enddo
    
!      3. Topography or Bathymetry plus nodal coordinates    


SWL0 = 0.0

read (datg_file,1) text         !  --  ictop is a switch for topo routines
write(chkg_file,1) text         !      we will use +ve for h avalanches
read (datg_file,*) ictop        !                  -ve for SW problems
write(chkg_file,*) ictop        !

if (ictop.eq.0) then            ! --- Bottom elev. is constant

    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) Zconst
    write(chkg_file,*) Zconst
    read (datg_file,1) text
    write(chkg_file,1) text
    do  i = 1,npoinF
        read (datg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
        write(chkg_file,3) ipoin,(coord(j,ipoin),j=1,ndimnF)
        topol(1,ipoin) = Zconst
    enddo

elseif (ictop.eq.1) then        ! --- Read x,y,z

    read (datg_file,1) text
    write(chkg_file,1) text
    do  i = 1,npoinF
        read (datg_file,*) ip,(coord(j,ip),j=1,ndimnF),topol(1,ip)
        write(chkg_file,3) ip,(coord(j,ip),j=1,ndimnF),topol(1,ip)
    enddo

elseif (ictop.eq.-3) then       ! --- Multilinear law in X

    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) nptsZ
    write(chkg_file,*) nptsz
    
    allocate ( xptsZ(nptsZ) )
    allocate ( ZptsZ(nptsZ) )
    
    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) (xptsZ(i),i=1,nptsZ)
    write(chkg_file,*) (xptsZ(i),i=1,nptsZ)
    read (datg_file,*) (ZptsZ(i),i=1,nptsZ)
    write(chkg_file,*) (ZptsZ(i),i=1,nptsZ)
    read (datg_file,1) text
    write(chkg_file,1) text

    do ip = 1,npoinF
       read (datg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
       write(chkg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
       xx = coord(1,ipoin)
       do iptsZ = 1,nptsZ-1
          x0 = xptsZ(iptsZ)
          x1 = xptsZ(iptsZ+1)
          if (xx.ge.x0.and.xx.le.x1) then
              f0 = (xx-x1)/(x0-x1)
              f1 = 1.-f0
              zz = f0*ZptsZ(iptsZ) + f1*ZptsZ(iptsZ+1)
              topol(1,ipoin) = zz
              exit
          endif
       enddo
    enddo
    deallocate(xptsZ)
    deallocate(ZptsZ) 
    

elseif (ictop.eq.-5) then       ! --- Aiguablava. 
                                !     Ref_depth is depth from swl of ref level. 
                                !     It is also SWL
     read (datg_file,1) text
     write(chkg_file,1) text
     read (datg_file,*) x0,y0,xs,ys
     write(chkg_file,*) x0,y0,xs,ys
     read (datg_file,1) text
     write(chkg_file,1) text
     read (datg_file,*) Rout,Depthout,Rint,Depthint, Xl0, depth2, TolerR
     write(chkg_file,*) Rout,Depthout,Rint,Depthint, Xl0, depth2, TolerR
     read (datg_file,1) text
     write(chkg_file,1) text
     read (datg_file,*) Ref_depth
     write(chkg_file,*) Ref_depth
     SWL0 = Ref_depth
     read (datg_file,1) text
     write(chkg_file,1) text
     do i = 1,npoinF
        read (datg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
        write(chkg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
        x    = coord (1,ipoin)
        y    = coord (2,ipoin)
        rx   = x-x0
        ry   = y-y0
        r    = (rx*rx + ry*ry)**0.5
        vnx  = -(ys-y0)
        vny  =  (xs-x0)
        vmod = (vnx*vnx + vny*vny)**0.5
        dot  = (rx*vnx + ry*vny)
        if (dot.ge.0.or.y.lt.0.0) then
           if (r.gt.rout) then
               dpth = depthout
           elseif (r.ge.rint) then
               dpth = ((r-rint)/(rout-rint))*(depthout-depthint) + depthint
           else
               dpth = depthint
           endif
        else
           dot2 = -dot/vmod
           if (dot2.le.xl0) then
               dpth = (dot2/xl0)*(depth2-depthint) + depthint
           else
               dpth = depth2
           endif
        endif
        topol(1,ipoin) = Ref_depth - dpth
     enddo
     
elseif (ictop.eq.-10) then       ! --- We obtain h's from TOPO module  

    read (datg_file,1) text
    write(chkg_file,1) text
    do ip = 1,npoinF
       read (datg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
       write(chkg_file,*) ipoin,(coord(j,ipoin),j=1,ndimnF)
    enddo 
    
    if (npoig.eq.0) then        !  -- This means we haven't called before TOPO
        call TOPO_MAIN
    endif
    
    Do ip = 1, npoinF
       xx = coord (1,ip)
       yy = coord (2,ip)
       call Get_Z_Topo (ip,xx,yy,Zp)    !  -- Zp interpolated by TOPO in geom grid 
       topol(1,ip) = Zp
    Enddo
     
ELSE
     write(*,*) ' Error: topography not implemented'
     stop
     
ENDIF    


!      4. Initial values of unknowns    -----------------------


read (datg_file,1) text
write(chkg_file,1) text
read (datg_file,*) icunkno
write(chkg_file,*) icunkno


if     (icunkno.eq.1) then      !  --  h, hu and hv are zero everywhere

        unkno = 0.0

elseif (icunkno.eq.-1) then     !  --  water level is ct, v are zero

        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) swl  !      still water level
        write(chkg_file,*) swl
        do ipoin=1,npoinF
           h = swl - topol(1,ipoin)
           if(h.lt.comp)  h=0.0
           unkno(1,ipoin) = h
        enddo

elseif (icunkno.eq.-2) then ! --- swl with hump in X

        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) swl  ! still water level
        write(chkg_file,*) swl
        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) xcenter,xlong,deltah
        write(chkg_file,*) xcenter,xlong,deltah
        xhmin = xcenter - xlong/2.
        xhmax = xcenter + xlong/2.
        do ipoin=1,npoinF
           xx = coord(1,ipoin)
           h  = swl - topol(1,ipoin)
           if (h.lt.comp)  h = 0.
           if (xx.le.xhmax.and.xx.ge.xhmin) h = deltah + h
           unkno(1,ipoin) = h
        enddo

elseif (icunkno.eq.-3) then ! --- swl with hump in X

        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) swl  ! still water level
        write(chkg_file,*) swl
        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) xcenter,xlong,deltah
        write(chkg_file,*) xcenter,xlong,deltah
        xhmin = xcenter - xlong/2.
        xhmax = xcenter + xlong/2.
        do ipoin=1,npoinF
           xx = coord(2,ipoin)
           h  = swl - topol(1,ipoin)
           if (h.lt.comp)  h = 0.
           if (xx.le.xhmax.and.xx.ge.xhmin) h = deltah + h
           unkno(1,ipoin) = h
        enddo
                

else if (icunkno.eq.-4) then  ! --- swl with cosinus ---

        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) swl , a0, xlong  ! still water level
        write(chkg_file,*) swl , a0, xlong
        do ipoin = 1,npoinF
           x     = coord(1,ipoin)
           h     = swl - topol(1,ipoin) + a0*cos(3.1416*x/xlong)
           if(h.lt.comp)  h = 0.0
           unkno(1,ipoin) = h
        enddo
        
else if (icunkno.eq.-5) then    ! --- cone ---

        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) swl , x0, y0, R, h0
        write(chkg_file,*) swl , x0, y0, R, h0
        do ipoin = 1,npoinF
           dx    = coord(1,ipoin)-x0
           dy    = coord(2,ipoin)-y0
           rr    = sqrt(dx*dx+dy*dy)
           h     = swl - topol(1,ipoin)
           if(h.lt.comp)  h=0.0
           if (rr.lt.R) h = h + h0*(1.-rr/R)
           unkno(1,ipoin) = h
        enddo
        
elseif (icunkno.eq.-6) then     !  cosinus +1 wave x0->x1

      read (datg_file,1) text
      write(chkg_file,1) text
      read (datg_file,*) x0,y0,x1,y1,b,d0,swl
      write(chkg_file,*) x0,y0,x1,y1,b,d0,swl
      xc    = (x0+x1)/2.
      yc    = (y0+y1)/2.
      a     = ( (x1-xc)**2.+(y1-yc)**2 )**0.5
      vnxi  = (x1-xc)/a
      vneta = (y1-yc)/a
      vvxi  = - vneta
      vveta =   vnxi
      do ipoin = 1,npoinF
         xi    = (coord(1,ipoin)-xc)*vnxi + (coord(2,ipoin)-yc)*vneta
         eta   = (coord(1,ipoin)-xc)*vvxi + (coord(2,ipoin)-yc)*vveta
         dh    = 0.0
         if (abs(eta).le.b) then
             dh = (d0/2.)*(1.+ cos(3.1416*eta/b))
         endif
         h  = swl - topol(1,ipoin) + dh
         unkno(1,ipoin) = h  
     enddo

elseif (icunkno.eq. 3) then ! --- h with hump in Y

    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) y0,h0,y1,h1
    write(chkg_file,*) y0,h0,y1,h1
    do ipoin = 1,npoinF
       yy    = coord(2,ipoin)
       hh    = 0.0
       if (yy.le.y1.and.yy.ge.y0) then
           yn0 = (yy-y1)/(y0-y1)
           yn1 = 1.-yn0
           hh  = yn0*h0+yn1*h1
       endif
       if(hh.lt.comp)  hh=0.0
       unkno(1,ipoin) = hh
    enddo
        
elseif (icunkno.eq.6) then      ! --- Multilinear law in X

    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) nptsh
    write(chkg_file,*) nptsh
    allocate ( xptsh(nptsh) )
    allocate ( hptsh(nptsh) )
    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) (xptsh(i),i=1,nptsh)
    write(chkg_file,*) (xptsh(i),i=1,nptsh)
    read (datg_file,*) (hptsh(i),i=1,nptsh)
    write(chkg_file,*) (hptsh(i),i=1,nptsh)
    do ipoin = 1,npoinF
       xx    = coord(1,ipoin)
       do iptsh = 1,nptsh-1
          x0 = xptsh(iptsh)
          x1 = xptsh(iptsh+1)
          if (xx.ge.x0.and.xx.le.x1) then
              f0 = (xx-x1)/(x0-x1)
              f1 = 1.-f0
              hh = f0*hptsh(iptsh)+f1*hptsh(iptsh+1)
              if(hh.lt.comp)  hh = 0.0
              unkno(1,ipoin) = hh
          endif
       enddo
    enddo
    deallocate(xptsh)
    deallocate(hptsh)
    
elseif (icunkno.eq.7) then      ! --- cylindrical deposit: x0,y0,r,cota_surf

      read (datg_file,1) text
      write(chkg_file,1) text
      read (datg_file,*) x0,y0,r0,cota
      write(chkg_file,*) x0,y0,r0,cota
      r2 = r0*r0
      do ipoin = 1,npoinF
         xd = coord(1,ipoin)-x0
         yd = coord(2,ipoin)-y0
         dist2 = xd*xd + yd*yd
         if (dist2.le.r2) then
             hh = cota - topol(1,ipoin)
             if (hh.lt.comp) then
                 hh = 0.0
                 write(*,*) ' DEPOSIT hh negative at pt ',ipoin
             endif
             unkno(1,ipoin) = hh
         endif
      enddo

    
elseif (icunkno.eq.20) then     ! --- Reservoir

      call Read_reservoir_SW_GFL
      
else
       write(*,*) 'Error: icunkno = ',icunkno, ' not implemented'
       stop
endif

!      ------  check in SW computations for dry points

icdry = 0

if (ntype.eq.2) then

   DRY: do ipoin = 1, npoinF
           if (topol(1,ipoin).ge. swl) then
                       icdry = 1
                        exit DRY
           endif
   enddo DRY

   if (icdry.eq.1) then
      write(*,*) '            *****  '
      write(*,*) ' there are dry points. Please low their Z by dZ'
      write(*,*) ' Input dZ '
      read (*,*) dZ
      do ipoin = 1,npoinF
         if (  (swl - topol(1,ipoin).lt.dZ)  ) then                
                topol(1,ipoin) = swl - dZ
                unkno(1,ipoin) = dZ
                eta0 = eta0 !  just to pause it while debugging
         endif
      enddo
   endif
   
endif   


!      5. Initial pore pressures pwp0, if we have activated this option

IF (icpwpF.eq.1) THEN

    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) icpwpo 
    write(chkg_file,*) icpwpo

    If (icpwpo.eq. 2) then      ! --- trapezoidal law in x
    
        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) x0,ph0,x1,ph1
        write(chkg_file,*) x0,ph0,x1,ph1
        do ipoin = 1,npoinF
           xx    = coord(1,ipoin)
           phh   = 0.0
           if (xx.ge.x0.and.xx.le.x1) then
               xn0 = (xx-x1)/(x0-x1)
               xn1 = 1.-xn0
               phh  = xn0*ph0 + xn1*ph1
           endif
           unkno(4,ipoin) = phh
        enddo
    else if (icpwpo.eq.1) then  ! --- constant pressure
    
        read (datg_file,1) text
        write(chkg_file,1) text
        read (datg_file,*) ph0
        write(chkg_file,*) ph0
        do ipoin = 1,npoinF
             phh=0.
             xx=unkno(1,ipoin)
             if (xx.gt.comp) then
               phh  = ph0
             endif
             unkno(4,ipoin) = phh
        enddo
    else
       write(*,*) ' unknown type of pwp introduced. STOP'
       STOP
    endif
    
ENDIF

!      6. Boundary Conditions:


!         6.a Tangent vectors  
!                     btanr : node nr.
!                     rtanr : tangential vector

if(nbour.ne.0) then

   read (datg_file,1) text
   write(chkg_file,1) text
   do  i = 1,nbour
       read (datg_file,*)  btanr(i),k,(rtanr(j,i),j=1,ndimnF)
       if(btanr(i).eq.0.or.btanr(i).gt.npoinF.or.nbour.gt.npoinF) then
           write(*,*)'Error writing tangential vectors'
           stop
       endif
       write(chkg_file,19) btanr(i),k,(rtanr(j,i),j=1,ndimnF)
   enddo
   
endif


  
!       6.b prescribed values : bprer : node nr. & variable
!                               unpre : (a0+a1*sin(wt-fi))*(1-exp(-t/tf))
!                                        a0,a1,w,fi,tf
!           if bprer(2,ipoin)= iunkn = -2 or -3 we prescribe fluxes
!           
!           if time curves are used: 
!              v3 is given a negative value  (frequency never is)
!              v1 is the time curve being used
!              v2 is the factor by which we multiply all values f(t)
!              time will be read later, in this section
!              f    is a f(t). Value of unknown is v2*f(t)



if (nprer.ne.0) then


    k = 0
    twopi = acos(0.0)*4.0/360.
    read (datg_file,1) text
    write(chkg_file,1) text
    do i = 1,nprer
       read (datg_file,*) (bprer(j,i),j=1,2),(unpre(j,i),j=1,5)
       if(bprer(1,i).eq.0.or.bprer(1,i).gt.npoinF.or.nprer.gt.npoinF  &   
                    .or.bprer(2,i).gt.3) then
             write(*,*)'Error writing prescribed values'
             stop
       endif
       write(chkg_file,4) (bprer(j,i),j=1,2),(unpre(j,i),j=1,5)
       unpre(4,i)=unpre(4,i)*twopi
       
       if ( unpre(3,i).lt.0 ) then
           k = - unpre(3,i)
           k = int(k+0.01)  !   time curve nr
       endif
       
    enddo
    
     
    if (k.gt.0) then                    !    read time curves
                                        !    variables belong to sph_gfl_time_vars
       read (datg_file,1) text
       write(chkg_file,1) text
       read (datg_file,*) ntcurves, mptstcurves ! number of time curves
       write(chkg_file,*) ntcurves, mptstcurves ! and max nr of pts in all
       
       allocate ( nptstcurves(ntcurves) )               !
       allocate ( ttcurves(ntcurves,mptstcurves) )      !
       allocate ( ftcurves(ntcurves,mptstcurves) )      !
       
       do i = 1, ntcurves
          read (datg_file,1) text
          write(chkg_file,1) text
          read (datg_file,*) nptstcurves(i)     !    nr of pts in each time curve
          write(chkg_file,*) nptstcurves(i)
          read (datg_file,1) text
          write(chkg_file,1) text
          read (datg_file,*) (ttcurves(i,j), j=1,nptstcurves(i))        !  t values
          read (datg_file,*) (ftcurves(i,j), j=1,nptstcurves(i))        !  f values
          write(chkg_file,*) (ttcurves(i,j), j=1,nptstcurves(i))  
          write(chkg_file,*) (ftcurves(i,j), j=1,nptstcurves(i)) 
       enddo
          
     endif
          

endif

 
  
!      6.c Absorbing nodes (take care, it is the super critical version)
  

if(nabshz.ne.0) then

   read (datg_file,1) text
   write(chkg_file,1) text
   read (datg_file,*) (babhz(j),j=1,nabshz)
   do j = 1,nabshz
      if(babhz(j).eq.0.or.babhz(j).gt.npoinF.or.nabshz.gt.npoinF)then
         write(*,*)'Error writing aborbing nodes with vn>0'
         stop
      endif
   enddo
   write(chkg_file,5) (babhz(j),j = 1,nabshz)
endif

  
  
!       6.d loaded nodes    : bload : node nr.
!                                thetax: angle of x with wave ray
!                                omega : freq.(rad.s-1)
!                                amplit: wave amplitude


if (nload.ne.0) then

    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) (bload(j),j=1,nload),ntypbc
    write(chkg_file,*) (bload(j),j=1,nload),ntypbc
    do  j = 1,nload
        if(bload(j).eq.0.or.bload(j).gt.npoinF.or.nload.gt.npoinF)then
           write(*,*)'Error writing loaded nodes'
           stop
        endif
    enddo
    read (datg_file,1)  text
    write(chkg_file,1)  text
    read (datg_file,*) (wave(iwv),iwv=1,nwave)
    write(chkg_file,*) (wave(iwv),iwv=1,nwave)
    
endif

  
!      6.e absorb.nodes    : babso : node nr.
  
if(nabso.ne.0) then
      read (datg_file,1) text
      write(chkg_file,1) text
      read (datg_file,*) (babso(j),j=1,nabso)
      do j = 1,nabso
         if(babso(j).eq.0.or.babso(j).gt.npoinF.or.nabso.gt.npoinF)then
         write(*,*)'Error writing aborbing nodes'
         stop
         endif
      enddo
      write(chkg_file,5) (babso(j),j=1,nabso)
endif

!     6.f  Gates for interaction with sph

if (type_of_fem_sph_interaction.eq.2) then

      read (datg_file,1) text
      write(chkg_file,1) text
      read (datg_file,*) n_ggates, mps_ggates   ! nr of gates and 
      write(chkg_file,*) n_ggates, mps_ggates   ! max nr of nodes in all gates
      do j = 1, n_ggates
         allocate (     nps_ggates(n_ggates))
         allocate (list_pts_ggates(n_ggates,2,mps_ggates))
         read (datg_file,1) text
         write(chkg_file,1) text
         read (datg_file,*) nps_ggates(j)
         write(chkg_file,5) nps_ggates(j)
         i = nps_ggates(j)
         read (datg_file,*) (list_pts_ggates(j,1,k), i=1,i)
         write(chkg_file,5) (list_pts_ggates(j,1,k), i=1,i)
         list_pts_ggates(j,2,:) = 0             !  2nd row sph nodes not known=0
      enddo 
           
endif

!      7. Constants or material properties

 
read (datg_file,1) text
write(chkg_file,1) text
read (datg_file,*) (constF(i),i=1,nconstF)    ! MP   note that constant 14 is hrelpwp = hsat/h
do iconst = 1,nconstF
   constF(iconst) = abs(constF(iconst))
enddo
write(chkg_file,7) (constF(i),i = 1,nconstF)
nfric = constF(5)+0.01
if (nfric.eq.2.and.constF(2).le.0.001) then
  write(*,*)' Error: density should not be zero for Julien frict'
  stop
endif
if(constF(10).le.1.e-32)then
  write(*,*)' Error: hfrict must be greater than 1.e-32'
  stop
endif


!      8. Niter and csmooth


read (datg_file,1) text   !  NOTE we only use niter and csmooth.
write(chkg_file,1) text   !  the rest have been read in DRIVER
read (datg_file,*) niter, csmooth 
write(chkg_file,*) niter, csmooth

csmooth = abs(csmooth)
if (csmooth.gt.2) then
   write(*,*)'Error writing csmooth: must be smaller than 2'
   stop
endif


!      9. Output nodes


if   (nprpoi.ne.0) then
     read (datg_file,1) text
     write(chkg_file,1) text
     read (datg_file,*) (lprpoi(i),i = 1,nprpoi)
     write(chkg_file,*) (lprpoi(i),i = 1,nprpoi)
endif

  
!       ------  wind stress
  
read (datg_file,1) text
write(chkg_file,1) text
read (datg_file,*) (twind(i),i = 1,ntwind)
write(chkg_file,*) (twind(i),i = 1,ntwind)
write(chkg_file,*) twind(1),twind(2)
  
!       ------  atmospheric pressure
  
 read (datg_file,1) text
 write(chkg_file,1) text
 read (datg_file,*) (datm(i),i = 1,ndatm)
 write(chkg_file,*) (datm(i),i = 1,ndatm)

!       ------  weir info (cts, list, heights)

if (nweir.ne.0) then
    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) cweir(1), cweir(2)
    write(chkg_file,*) cweir(1), cweir(2)
    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) (lweir(i),i = 1,nweir)
    do j = 1,nweir
       if(lweir(j).eq.0.or.lweir(j).gt.npoinF.or.nweir.gt.npoinF)then
          write(*,*)'Error writing weir nodes'
          stop
       endif
    enddo
    write(chkg_file,*) (lweir(i),i=1,nweir)
    read (datg_file,1) text
    write(chkg_file,1) text
    read (datg_file,*) (zweir(i),i=1,nweir)
    write(chkg_file,*) (zweir(i),i=1,nweir)
endif
 
!      ------  converts pwp to hw values multiplying it by  h  MP
 
if (icpwpF.eq.1) then 
    do ipoin = 1, npoinF
       unkno(4,ipoin) = unkno(4,ipoin)*unkno(1,ipoin) 
    enddo
endif

 
!     -----formats.....
 
 1 format(a80)
 3 format(i10,10f15.4)
 4 format(2i8,5f15.7)
 5 format(10i8)
 7 format(e15.6,2x,e15.6,2x,e15.6,2x,e15.6,2x,e15.6)
 9 format(2f10.4,6i10)
19 format(2(i10,2x),4f15.7)
!
 
END  SUBROUTINE getinput

! ----------------------------------------------------------------------   

      subroutine Read_reservoir_SW_GFL 

! ----------------------------------------------------------------------   

!     This routine gets reservoir water nodes info from topo routine... 
!     Reservoir is defined by a dam (we give x0,y0 and xf,yf and Z of water surface 
!     nx and ny_dam are components of vector normal to dam, pointing inwards.
!     to check that a point is in the reservoir, its Z topo should be below water surface Z, Z_swl, 
!     and the product (x-x0)*v >0

implicit none  

integer(ink) ipoin

real   (irk) Z_swl, dgx, dgy, facthsml, lgx, lgy, xx, yy, RR 
real   (irk) Zp, dArea
real   (irk) x0_dam, y0_dam, xf_dam, yf_dam, nx_dam, ny_dam             !   initial and final point of dam, and normal vector
real   (irk) xmin_rsv, xmax_rsv, ymin_rsv, ymax_rsv                     !   the reservoir is within rectangle defined by xmin xmax, ymin and ymax_dam

read (datg_file,1005) text 
write(chkg_file,1005) text
read (datg_file,*)    xmin_rsv, xmax_rsv, ymin_rsv, ymax_rsv, dgx, dgy  ! dgx is deltax for spacing nodes in X 
write(chkg_file,*)    xmin_rsv, xmax_rsv, ymin_rsv, ymax_rsv, dgx, dgy  
read (datg_file,1005) text 
write(chkg_file,1005) text
read (datg_file,*)    Z_swl, x0_dam, y0_dam, xf_dam, yf_dam  
write(chkg_file,*)    Z_swl, x0_dam, y0_dam, xf_dam, yf_dam 

1005     format(a60) 

if (ndimnF.eq.1) then
    write(*,*) ' routine read_reservoir only for 2D problems '
    write(*,*) ' hit any key '
    read (*,*) text
    STOP
endif

nx_dam = -(yf_dam-y0_dam)       ! this is the normal pointing inwards reservoir
ny_dam =   xf_dam-x0_dam

lgx   = (xmax_rsv-xmin_rsv-dgx)
lgy   = (ymax_rsv-ymin_rsv-dgy)

DO ipoin = 1, npoinF

   xx     = coord (1,ipoin)     
   yy     = coord (2,ipoin)
   Zp     = topol (1,ipoin)
   
   RR = (xx-x0_dam)*nx_dam + (yy-y0_dam)*ny_dam
   
   if (Zp.lt.Z_swl.and.RR.ge.0.0.and.Zp.gt.-1.e6) then     
      unkno(1,ipoin) = Z_swl - Zp
   endif
   
ENDDO 
 
 
END SUBROUTINE  Read_reservoir_SW_GFL


!---------------------------------------------------------

        subroutine getgeom2D 
     
!--------------------------------------------------------

!      La subrutina getgeom2D calcula las derivadas de las funciones de forma respecto x e y y el área de los elementos
!      almacenando todo en la variable geome
!      Se han quitado las variable ntotg y ngaus

implicit none

integer(ink) ielem, inode, in 

real   (irk) x21, x31, y21, y31, rj, rj1, xix, xiy, etx, ety
real   (irk) rnxi, rnet 
real   (irk) x(3),y(3),nxi(3),net(3)!,posgp(2,3),shapef(4)
real   (irk) elcod(2,3)
 

data  nxi/-1.0 , 1.0 , 0.0 / !derivadas de las funciones de forma respecto x
data  net/-1.0 , 0.0 , 1.0 / !derivadas de las funciones de forma respecto y

!     Las funciones de forma son:
!               N1=1-xi-et
!               N2=xi
!               N3=et
!
!            This sub evaluates n,x & n,y for each element
!            shape function (linear triangles) and the
!            jacobian (=2*area)

do  ielem = 1,nelem   !     -----  loop over the elements

    do inode = 1,nnode
       in = intmat(inode,ielem)
       x(inode) = coord(1,in)
       y(inode) = coord(2,in)
       elcod(1,inode) = x(inode)
       elcod(2,inode) = y(inode)
    enddo

! ------      evaluate the geometrical quantities needed

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
       geome(in  ,ielem)=xix*rnxi+etx*rnet
       geome(in+3,ielem)=xiy*rnxi+ety*rnet
    enddo
!
    geome(7,ielem) = rj    ! Area
    
ENDDO           !     -----  end of loop over the elements

END subroutine getgeom2D

! --------------------------------------------------------------------

        subroutine Get_Topol 

! --------------------------------------------------------------------

!      Purpose : given  topol obtain its derivatives x and y
!                and second order derivatives   MP 01 04 2003

!                       Topol(1,ipoin)...Z
!                       Topol(2,ipoin)...Zx
!                       Topol(3,ipoin)...Zy
!                       Topol(4,ipoin)...E  (1.+zx*zx) 
!                       Topol(5,ipoin)...F  (zx*zy)
!                       Topol(6,ipoin)...G  (1.+zy*zy)
!                       Topol(7,ipoin)...L  zxx/rt 
!                       Topol(8,ipoin)...M  zxy/rt
!                       Topol(9,ipoin)...N  zyy/rt
!                                        where rt is (1.+zx*zx+zy*zy)**0.5

implicit none

real   (irk)  mmat(3,3)
real   (irk), allocatable:: rhs  (:,:)      ! --  Zxx Zxy Zyy
real   (irk), allocatable:: rhelp(:,:)      ! -- auxiliar


integer(ink) id, idimn, ielem, iiter, inode, ipoin, jnode, jpoin
integer(ink) niter_t

real   (irk) coe, cm, derixe, deriye, derixxe, derixye, deriyxe, deriyye
real   (irk) zx, zy, zxx, zxy, zyy, rt



niter_t = 5                     !      perform 5 iters to solve the system
allocate ( rhs(3,npoinF), rhelp (ndimnF+1,npoinF ) )                                


!      ------ A1  Get Mass matrix  

mmat = 1.0
do inode = 1,nnode
   mmat(inode,inode)=2.0
enddo

!      ------ A2   Obtain gradients at elements  

do ielem = 1,nelem
   derixe = 0.0
   deriye = 0.0
   do inode = 1,nnode
      ipoin = intmat(inode,ielem)
      derixe = derixe + topol(1,ipoin)*geome(inode,  ielem)
      deriye = deriye + topol(1,ipoin)*geome(inode+3,ielem)
   enddo
   topel(2,ielem) = derixe
   topel(3,ielem) = deriye
enddo

!      ------  A3  Build RHS Int(Ni.Grad) Store it in rhs

rhs = 0.0

do ielem = 1,nelem
   coe = geome(7,ielem)/6.
   do inode = 1,nnode
      ipoin = intmat(inode,ielem)
      rhs (2,ipoin) = rhs(2,ipoin) + topel(2,ielem)*coe
      rhs (3,ipoin) = rhs(3,ipoin) + topel(3,ielem)*coe
   enddo
enddo

!      ------ A4   Get topol v(1)=v(0)+dt.inv(Ml)*(rhs-M.v(0)) v(0)=0

do ipoin = 1,npoinF
   cm = mmatl(ipoin)
   topol(2,ipoin) = rhs(2,ipoin)/cm !+topol(2,ipoin) (EG)
   topol(3,ipoin) = rhs(3,ipoin)/cm !topol(3,ipoin) (EG)
enddo

!      ------ A5  Iterate:       v(n+1) = v(n) + inv(Ml)*(rhs-M*V(n))
!                                M*v(n) is stored as rhelp

do iiter = 2,niter_t
   rhelp = 0.0
   do ielem = 1,nelem
      coe = geome(7,ielem)/24.
      do inode = 1,nnode
         ipoin = intmat(inode,ielem)
         do jnode = 1,nnode
            jpoin = intmat(jnode,ielem)
            cm    = mmat(inode,jnode)
            do idimn = 1,ndimnF
               rhelp(idimn+1,ipoin) =  rhelp(idimn+1,ipoin) & 
                                       + cm*coe*topol(idimn+1,jpoin)
            enddo
         enddo
      enddo
   enddo

   do ipoin = 1,npoinF
      cm = mmatl(ipoin)
      do idimn = 1,ndimnF
         topol(idimn+1,ipoin) = topol(idimn+1,ipoin) &
              + (rhs(idimn+1,ipoin)-rhelp(idimn+1,ipoin))/cm
      enddo
   enddo
enddo
       
!      ------  Second order derivatives

!      ------ B2   Obtain gradients at elements  

do ielem = 1,nelem
   derixxe = 0.0
   derixye = 0.0
   deriyxe = 0.0
   deriyye = 0.0
   do inode = 1,nnode
      ipoin   = intmat(inode,ielem)
      derixxe = derixxe + topol(2,ipoin)*geome(inode  ,ielem)
      derixye = derixye + topol(2,ipoin)*geome(inode+3,ielem)
      deriyxe = deriyxe + topol(3,ipoin)*geome(inode  ,ielem)
      deriyye = deriyye + topol(3,ipoin)*geome(inode+3,ielem)
   enddo
   topel(4,ielem) =  derixxe
   topel(5,ielem) = (derixye + deriyxe)/2.
   topel(6,ielem) =  deriyye
enddo

!      ------  B3  Build RHS Int(Ni.Grad) Store it in rhs

rhs = 0.0

do ielem =1,nelem
   coe = geome(7,ielem)/6.
   do inode = 1,nnode
      ipoin = intmat(inode,ielem)
      rhs (1,ipoin) = rhs(1,ipoin) + topel(4,ielem)*coe 
      rhs (2,ipoin) = rhs(2,ipoin) + topel(5,ielem)*coe
      rhs (3,ipoin) = rhs(3,ipoin) + topel(6,ielem)*coe
   enddo
enddo

!      ------ B4   Get topol v(1)=v(0)+dt.inv(Ml)*(rhs-M.v(0)) v(0)=0

do ipoin = 1,npoinF
   cm = mmatl(ipoin)
   topol(4,ipoin) =  rhs(1,ipoin)/cm    !  ** + topol(4,ipoin)  CHECK This prunning
   topol(5,ipoin) =  rhs(2,ipoin)/cm    !  ** + topol(5,ipoin) 
   topol(6,ipoin) =  rhs(3,ipoin)/cm    !  ** + topol(6,ipoin) 
enddo

!      ------ B5  Iterate:       v(n+1)=v(n)+inv(Ml)*(rhs-M*V(n))
!                                M*v(n) is stored as rhelp

do iiter = 2,niter_t

   rhelp = 0.0
   do ielem = 1,nelem
      coe   = geome(7,ielem)/24.
      do inode = 1,nnode
         ipoin = intmat(inode,ielem)
         do jnode = 1,nnode
            jpoin = intmat(jnode,ielem)
            cm    = mmat(inode,jnode)
            do id = 1,3
               rhelp(id,ipoin) = rhelp(id,ipoin) + cm*coe*topol(id+3,jpoin)
            enddo
         enddo
      enddo
   enddo

   do ipoin=1,npoinF
      cm = mmatl(ipoin)
      do id =1,3
         topol(id+3,ipoin) = topol(id+3,ipoin)  + (rhs(id,ipoin)-rhelp(id,ipoin))/cm
      enddo
   enddo
   
enddo

       
!      ------  Store E=1+zx^2    F=zx.zy G=1+zy^2 positions 4,5 y 6
!              L=zxx/root       M=zxy/root      N=zyy/root      positions 7, 8 y 9
!              root=(1+zx^2+zy^2)**0.5

Do ipoin = 1,npoinF
   zx   = topol(2,ipoin)
   zy   = topol(3,ipoin)
   zxx  = topol(4,ipoin)
   zxy  = topol(5,ipoin)
   zyy  = topol(6,ipoin)
   rt   = (1.+zx*zx+zy*zy)**0.5
   topol(4,ipoin) = (1.+zx*zx)
   topol(5,ipoin) = (zx*zy)
   topol(6,ipoin) = (1.+zy*zy)
   topol(7,ipoin) = zxx/rt
   topol(8,ipoin) = zxy/rt
   topol(9,ipoin) = zyy/rt
Enddo

deallocate (rhs, rhelp)

 END subroutine Get_Topol

! ---------------------------------------------------------------

        subroutine get_boun_nodes_2D   

! ---------------------------------------------------------------

!      ------  obtain tangents in the global mesh

!               nboun           number of boundary nodes found
!               rtang(2,nboun)  tangent vectors at boundary nodes
!               bconl(nboun)    list of boundary nodes
!
!               we use topel as an auxiliar matrix. 
!               Remember, after it will be re defined for elems topo

implicit none

integer(ink)  ib, iboun, id, ie, ielem, in, in1, in2, ip, ip1, ip2, it
integer(ink)  ipoin, inode
real   (irk)  amodl, c13

!      ------  initialization

ltang = 0
bconl = 0
topel = 0.0

c13=1./3.

do ie  = 1,nelem                !  --  Algorithm to identify B nodes
do in  = 1,nnode                !      they have ltang .ne. 0
   in1 = in+1
   in2 = in+2
   if(in1.gt.3) in1 = in1-3
   if(in2.gt.3) in2 = in2-3
   ip  = intmat(in,ie)
   ip1 = intmat(in1,ie)
   ip2 = intmat(in2,ie)
   ltang(ip) = ltang(ip)+ip2-ip1
   do it = 1,ntopo              !  --  evaluate mean topological data
      topel(it,ie) = topel(it,ie)+c13*topol(it,ip)
   enddo
enddo
enddo

!      ------  find the tangent vector

nboun = 0

DO ipoin = 1,npoinF

   if(ltang(ipoin).eq.0) CYCLE  !  --  because it is not a boundary node

   nboun = nboun+1
   bconl(nboun)= ipoin

   do id = 1,ndimnF
      rtang(id,nboun) = 0.0
   enddo

   do ielem = 1,nelem   !  -- loop over the elements surrounding that point
      ip  = intmat(1,ielem)
      ip1 = intmat(2,ielem)
      ip2 = intmat(3,ielem)
      if(ip.ne.ipoin.and.ip1.ne.ipoin.and.ip2.ne.ipoin) CYCLE
      do inode = 1,nnode
         in = inode
         if(intmat(in,ielem).eq.ipoin) EXIT
      enddo
      rtang(1,nboun) = rtang(1,nboun) - geome(in+3,ielem)*geome(7,ielem)
      rtang(2,nboun) = rtang(2,nboun) + geome(in  ,ielem)*geome(7,ielem)
    enddo 

    do ib = 1,nbour             !  --  add the prescribed tangents
       ip = btanr(ib)
       if(ip.ne.ipoin) CYCLE
       do id = 1,ndimnF
          rtang(id,nboun)=rtanr(id,ib)
       enddo
    enddo

ENDDO 

!      ------  normalize and output the tangent vectors

write(chkg_file,10)

do iboun = 1,nboun
   amodl = sqrt(rtang(1,iboun)**2+rtang(2,iboun)**2)
   if(amodl.lt.1.0e-06) CYCLE
   do id = 1,ndimnF
      rtang(id,iboun) = rtang(id,iboun)/amodl
   enddo
   write(chkg_file,12) bconl(iboun),(rtang(id,iboun),id = 1,2)
enddo

10     format(' surrounding domain and prescribed tangents ')
12     format(i10,2f15.6)


END subroutine get_boun_nodes_2D

! -------------------------------------------------------------

        subroutine bounel

! -------------------------------------------------------------

!      ------  this subroutine finds out the boundary elements and
!              their conectivities. Geommetrical parameters are also
!              obtained

!               neleb             number of boundary sides of wall type
!               intmab(1:2,ieleb) nodes
!               intmab(3,  ieleb) element to which side belongs
!               bindex(neleb)     1 if  boundary side is of wall type

implicit none


integer(ink) i,ie, in, in1, ip, ip1, it, it1, je, jn, jn1
integer(ink) ieleb, ipoin1, ipoin2, j 
real   (irk) aleng, cosx, cosy, x1, x2, y1, y2

neleb =0

IE100: do ie = 1,nelem
IN101: do in = 1,nnode

    in1 = in+1
    if(in1.eq.4) in1 = 1
    ip = intmat(in ,ie)
    ip1= intmat(in1,ie)

    JE200: do je = 1,nelem
    JN201: do jn = 1,nnode
        it = intmat(jn ,je)
        if(ip1.ne.it) CYCLE
        jn1 = jn+1
        if(jn1.eq.4) jn1 = 1
        it1 = intmat(jn1,je)
        if(ip.eq.it1) CYCLE IN101  
    enddo JN201 
    enddo JE200 

    neleb = neleb + 1           !  --  conectivity
    intmab(1,neleb)= ip
    intmab(2,neleb)= ip1
    intmab(3,neleb)= ie

    x1 = coord(1,ip )           !  --  geometrical parameters
    y1 = coord(2,ip )
    x2 = coord(1,ip1)
    y2 = coord(2,ip1)

    aleng = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
    cosx  = (y2-y1)/aleng
    cosy  =-(x2-x1)/aleng
    geomb(1,neleb) = cosx
    geomb(2,neleb) = cosy
    geomb(3,neleb) = aleng

enddo IN101 
enddo IE100

!   Wall segments: if both nodes are abs, presc, loaded they are not wall-like

IELEB3000: do ieleb = 1,neleb

     ipoin1 = intmab(1,ieleb)
     ipoin2 = intmab(2,ieleb)

     do i = 1,nabso
       if(babso(i).eq.ipoin1)then
          do j = 1, nabso
             if(babso(j).eq.ipoin2)    CYCLE IELEB3000
          end do
       end if
     enddo

     do i = 1,nprer
        if(bprer(1,i).eq.ipoin1)then
           do j = 1,nprer
              if(bprer(1,j).eq.ipoin2) CYCLE IELEB3000
           end do
        end if
     end do

     do i = 1,nabshz
        if(babhz(i).eq.ipoin1)then
           do j = 1,nabshz
              if(babhz(j).eq.ipoin2)   CYCLE IELEB3000
           end do
        end if
     end do

     do i = 1,nload
        if(bload(i).eq.ipoin1)then
           do j = 1,nload
              if(bload(j).eq.ipoin2)   CYCLE IELEB3000
           end do
        end if
     enddo
   
     bindex(ieleb) = 1
     
enddo IELEB3000   

END subroutine bounel


! --------------------------------------------------------------

        subroutine delun0

! --------------------------------------------------------------


!      ------  Apply Dirichlet BC's
!
!
!               We use -2 and -3 for fluxes, i.e. q=hv
!               See also changes in solve routine

implicit none

integer(ink) iprer, ipoin, iunkn, i0
real   (irk) a0,a1,w,fi,tt,fact,argum
real   (irk) hinv

integer(ink) it_curves, ipts_tcurves
real   (irk) unk_pres, tt0, tt1, xi

IF(nprer.ne.0) then

   DO iprer = 1,nprer  
   
      ipoin = bprer(1,iprer)
      iunkn = bprer(2,iprer)
      a0    = unpre(1,iprer)    !        factor to multiply f(t)  
      a1    = unpre(2,iprer)    !         
      w     = unpre(3,iprer)    !        if -ve, is a the time curve index
      fi    = unpre(4,iprer)
      tt    = unpre(5,iprer)
      if(tt.gt.0.001) then
         fact = 1.0-exp(-time_gfl/tt)
      else
         fact = 1.00
      endif
      
      if (w.ge.-0.001) then
         
         argum = w*time_gfl - fi
         unk_pres = (a0+a1*sin(argum))*fact
         
      else
      
         it_curves = int (a0)

         do ipts_tcurves = 1, nptstcurves(it_curves)-1
            tt0 = ttcurves(it_curves,ipts_tcurves)
            tt1 = ttcurves(it_curves,ipts_tcurves+1)
            if ( time_gfl.ge.tt0.AND.time_gfl.le.tt1 ) EXIT
            tt0 = -1000.
         enddo
         
         if (tt0.ge.0.0) then
            xi = (time_gfl - tt0) /(tt1-tt0)
            unk_pres = (1.-xi) * ftcurves(it_curves , ipts_tcurves)  &
                         + xi  * ftcurves(it_curves , ipts_tcurves+1) 
         else
            unk_pres = 0.0
         endif  
         unk_pres = unk_pres * a1         


      endif    
      
      if (iunkn.eq.1) then
         if (ntype.eq.2) then
             delun(1,ipoin) = swl - topol(1,ipoin) & 
                             + unk_pres - unkno(1,ipoin)
         elseif (ntype.eq.3) then
             delun(1,ipoin) =  unk_pres - unkno(1,ipoin)
         endif        
         
      elseif(iunkn.eq.2.or.iunkn.eq.3) then     ! note that in ntype=3 
         i0   = iabs(iunkn)                     ! we should prescribe q and not v
         hinv = 0.0
         if (unkno(1,ipoin).ge.comp)   hinv = 1/unkno(1,ipoin)
         delun(i0,ipoin) =  unk_pres - unkno(i0,ipoin)*hinv
         
      elseif( (iunkn.eq.-2).or.(iunkn.eq.-3) ) then
         i0 = iabs(iunkn)
         delun(i0,ipoin) =  unk_pres - unkno(i0,ipoin)
         
      elseif(iunkn.eq.4) then !se ha prescrito una presión (EG)
         hinv = 0.0
         if (unkno(1,ipoin).lt.comp)   hinv = 0.0
         delun(4,ipoin) =   unk_pres - unkno(4,ipoin)*hinv
      endif
      
   ENDDO  
   
ENDIF

END subroutine delun0


! -----------------------------------------------------------------

        subroutine first_step 

! -----------------------------------------------------------------

!      ------  this subroutine calculates the first step in TG algorithm


implicit none


integer(ink) ia, igaus, ipoin, ielem, inode, ktotg
 
real   (irk) cgra,  cori, Bfact, hrelpw2, alpha 

real   (irk) h, hu, hv, pwp, HH, factw, tww
real   (irk) frictx, fricty, slopx, slopy, correcx, correcy

real   (irk) hinv, pres, correx, correy, c13

real   (irk) delt2,  anx, any
real   (irk) posgp(2,3), shapef(3)
real   (irk) source(4),  variab(4), divef2(4) 

real   (irk) mini_topo(ntopo), mini_patm(npatma)   !  -- compatibility call to frictlaws  


!      ------  Get fluxes and sources at nodes

cgra     = constF(1)
cori     = constF(4)
Bfact    = constF(13)
if (icpwpF.eq.1) hrelpw2  = constF(14)*constF(14)     !MP

tww      = twind(3)             !  -- wind contribution to sources
if(tww.lt.0.0001) tww = 0.0001
factw    = time_gfl/tww
if(factw.gt.1.00) factw = 1.00

DO ipoin = 1,npoinF              ! --------------------    DO IPOIN

   if (iactnd(ipoin).eq.0) CYCLE
   
   h  = unkno(1,ipoin)
   hu = unkno(2,ipoin)
   hv = unkno(3,ipoin)
   if (icpwpF.eq.1) then
       pwp = unkno(4,ipoin)   
   else
       pwp = 0.0
   endif

   hinv  = 0.0
   if (h.ge.comp) hinv = 1./h
   
   if (ntype.eq.2) then
      HH   = swl - topol(1,ipoin)
      if (WIR_Interact) then
         HH = HH - hs_w(ipoin) 
      endif
      pres = 0.5*cgra*(h*h-HH*HH)
   elseif (ntype.eq.3) then
      HH   = 0.0
      pres = 0.5 * cgra*h*h
   else
      write(*,*) ' Error: ntype not implemented'
      stop
   endif
   
   fluxp(1,ipoin) = hu
   fluxp(2,ipoin) = hu*hu*hinv + pres
   fluxp(3,ipoin) = hu*hv*hinv
   fluyp(1,ipoin) = hv
   fluyp(2,ipoin) = hu*hv*hinv
   fluyp(3,ipoin) = hv*hv*hinv + pres
   
   if (icpwpF.eq.1) then              
       fluxp(4,ipoin) = pwp*hinv*hu
       fluyp(4,ipoin) = pwp*hinv*hv
   endif
   
   mini_topo(:) = topol(:,ipoin)
   mini_patm(:) = patma(:,ipoin)
   
!      ------  interaction waves in reservoirs with sph
!              we add to minitopol 2,3  gradhs in fem nodes
!              if there is interaction, WIR_Interact = .TRUE.
!              in frictlaws, slopx = -g.(h-HH)*(Zx+gradhs)

   if (WIR_Interact) then
      mini_topo(2) = mini_topo(2) + gradhs_w(1,ipoin)
      mini_topo(3) = mini_topo(3) + gradhs_w(2,ipoin)
   endif
    
   call frictlaws (mini_topo, ntopo, mini_patm, npatma, &
                   h, hu, hv, pwp, HH, factw,           &
                   frictx, fricty,                      &
                   slopx, slopy, correx, correy            ) 

!   -----  Sources

   if(ntype.eq.3.and.nsour.eq.1) then
      gpoin(1,ipoin) =  0.0
      gpoin(2,ipoin) =  cori*hv - frictx + correx
      gpoin(3,ipoin) = -cori*hu - fricty + correy
   else
      gpoin(1,ipoin) =  0.0
      gpoin(2,ipoin) =  cori*hv - frictx + slopx + correx
      gpoin(3,ipoin) = -cori*hu - fricty + slopy + correy
   endif
   if (icpwpF.eq.1) gpoin(4,ipoin) =  -Bfact*hinv*hinv*unkno(4,ipoin)/hrelpw2   

ENDDO            

!      ------  If sources are integrated separately then:

if(ntype.eq.3.and.nsour.eq.2)then
   gpoin = 0.
endif

!      ------  Perform first step  time level n+alpha*n

unkne = 0.0
alpha  = 0.5
c13    = 1./3.

if (ngaus.eq.3) then
    ktotg      = 0
    posgp(1,1) = 1.0/6.0
    posgp(2,1) = 1.0/6.0
    posgp(1,2) = 4.0/6.0
    posgp(2,2) = 1.0/6.0
    posgp(1,3) = 1.0/6.0
    posgp(2,3) = 4.0/6.0
elseif (ngaus.eq.1) then
    ktotg      = 0
    posgp(1,1) = 1.0/3.0
    posgp(2,1) = 1.0/3.0
endif

DO ielem =1,nelem               !  --  loop over the elements  IELEM

   if (iactel(ielem).eq.0) CYCLE  
   
    if (ic_adapt_dt_elem_gfl.eq.0) then 
        delt2  = alpha*dt_gfl
    elseif (ic_adapt_dt_elem_gfl.eq.1) then
        delt2  = min(alpha*deltel(ielem),dt_gfl)
    endif
   
    divef2 = 0.0

    do inode = 1,nnode
       anx   = geome(inode  ,ielem) !  --  form n,x & n,y
       any   = geome(inode+3,ielem)
       ipoin = intmat(inode,ielem)

       do ia = 1,namat
          divef2(ia) = divef2(ia) + anx*fluxp(ia,ipoin)+any*fluyp(ia,ipoin)
       enddo

    enddo 
                
    ktotg = (ielem-1)*ngaus     !  --  As we consider only active elements...
    
    do igaus = 1,ngaus
         
       ktotg     = ktotg+1
       shapef(1) = posgp(1,igaus)
       shapef(2) = posgp(2,igaus)
       shapef(3) = 1.0 - shapef(1) - shapef(2)
       source    = 0.
       variab    = 0.
       do inode  = 1,nnode
          ipoin  = intmat(inode,ielem)
          do ia  = 1,namat
             source(ia) =  source(ia) + shapef(inode)*gpoin(ia,ipoin)
             variab(ia) =  variab(ia) + shapef(inode)*unkno(ia,ipoin)
          enddo
       enddo

       do ia = 1,namat
          unkne(ia,ktotg) = variab(ia) + delt2*(source(ia) - divef2(ia))
       end do
       
    end do
      
            
ENDDO  !                      -------------  end loop over elements
 
END subroutine first_step

! ------------------------------------------------------------------

        subroutine second_step 

! -----------------------------------------------------------------

implicit none

integer(ink) ktotg, ielem, igaus, inode, ipoin, in, ia, nfrict
integer(ink) ieleb, ipoin1, ipoin2, iaux, inodb, jnodb, jn, ip
integer(ink) iload, iabso, ipres 

real   (irk) cgra, dens, cmanning, cori, tauy0, constK, visco
integer(ink) nfric, itens
real   (irk) tanfi8, hfrict0, tanfi0, Bfact, c73, hrelpwp2, tanfi
real   (irk) rhog, gtanfi, gtanfi8, cgra05 
real   (irk) tww, factw, frictx, fricty, dvolu, h, hu, hv, pwp
real   (irk) hinv, HH, pres  
real   (irk) dudx, dudy, dvdx, dvdy, dxz, dyz
real   (irk) tauy, hinv2, ahat, gi,visco2, hfrict, hvmod, vmod 
real   (irk) hvmodinv 
real   (irk) slopx, slopy, correx, correy
real   (irk) anx, any, xni
real   (irk) raiz3
real   (irk) cosx, cosy, aleng, coe, cm, aux
real   (irk) hrelpw2, fi

real   (irk) posgp(2,3), shapef(3)
real   (irk) mmatb(2,2)
real   (irk) fluxe2(4),    fluye2(4),  gelem2(4), divef2(4)
real   (irk) bflux(4,2),   bfluy(4,2),  fin05(4)

real   (irk) mini_topo(ntopo), mini_patm(npatma)   !  -- compatibility call to frictlaws

rhs1     = 0.0

fluxe2   = 0.0
fluye2   = 0.0

denom    = 0.0

cgra     = constF(1)  
cgra05   = cgra*0.5
dens     = constF(2)
cmanning = constF(3)
cori     = constF(4)
nfric    = constF(5) + 0.01
tauy0    = constF(6)
constK   = constF(7)
visco    = constF(8)
tanfi8   = constF(9)
hfrict0  = constF(10)
itens    = constF(11)
tanfi0   = constF(12)
Bfact    = constF(13)
c73      = 7./3.

if (icpwpF.eq.1) hrelpw2  = constF(14)*constF(14)     !MP

if (icpwpF.eq.1) hrelpw2  = constF(14)*constF(14)     !MP

if (Bfact.ge.1.e-6) then
    tanfi  = tanfi8 + (tanfi0-tanfi8)*exp(-time_gfl/Bfact)
else
    tanfi  = tanfi8
endif

fi       = atan(tanfi)
gtanfi   = cgra*tanfi
gtanfi8  = cgra*tanfi8

if (icpwpF.eq.1) hrelpw2  = constF(14)*constF(14)    
rhog     = dens*cgra

c73      = 7./3.

tww      = twind(3)
if(tww.lt.0.0001) tww=0.0001
factw      = time_gfl/tww
if(factw.gt.1.00) factw=1.00

if (ngaus.eq.3) then
    ktotg      = 0
    posgp(1,1) = 1.0/6.0
    posgp(2,1) = 1.0/6.0
    posgp(1,2) = 4.0/6.0
    posgp(2,2) = 1.0/6.0
    posgp(1,3) = 1.0/6.0
    posgp(2,3) = 4.0/6.0
elseif (ngaus.eq.1) then
    ktotg      = 0
    posgp(1,1) = 1.0/3.0
    posgp(2,1) = 1.0/3.0
endif

         
DO IELEM = 1,NELEM
       
    if (iactel(ielem).eq.0)    CYCLE  
    ktotg = (ielem-1)*ngaus  
 
    DO igaus    = 1,ngaus
    
       frictx   = 0.
       fricty   = 0.
       ktotg    = ktotg+1
       shapef(1)= posgp(1,igaus)
       shapef(2)= posgp(2,igaus)
       shapef(3)= 1.0-shapef(1)-shapef(2)
       dvolu    = geome(7,ielem)/(2*ngaus)      ! weigth (area assoc to GP)
       h        = unkne(1,ktotg)
       hu       = unkne(2,ktotg)
       hv       = unkne(3,ktotg)
       pwp     = 0.0
       if (icpwpF.eq.1) pwp = unkne(4,ktotg)

       hinv    = 0.0
       if (h.ge.comp) hinv=1./h

       if(ntype.eq.2)  then
          HH = 0.0
          do inode = 1,nnode
             ipoin = intmat(inode,ielem)
             HH = HH + shapef(inode)*topol(1,ipoin)
             if (WIR_Interact) HH = HH + shapef(inode)*hs_w(ipoin)
          end do
          HH   = swl-HH
          pres = cgra05*(h*h-HH*HH)
       elseif (ntype.eq.3)then
          HH=0.0
          pres = cgra05*h*h
       end if

       fluxe2(1)  = hu          !    ------  Obtain fluxes
       fluxe2(2)  = hu*hu*hinv + pres
       fluxe2(3)  = hu*hv*hinv
       fluye2(1)  = hv
       fluye2(2)  = fluxe2(3)
       fluye2(3)  = hv*hv*hinv + pres
       if (icpwpF.eq.1) then               
           fluxe2(4) = pwp*hinv*hu   
           fluye2(4) = pwp*hinv*hv   
       endif

       if(abs(itens).gt.0)then  ! -- Viscous contribution
       
           dudx = 0.0
           dudy = 0.0
           dvdx = 0.0
           dvdy = 0.0
           dxz  = 0.0
           dyz  = 0.0
          
           do inode = 1,nnode
              in    = intmat(inode,ielem)
              hinv  = 0.0
              h     = unkno(1,in)
              if(h.ge.comp) hinv = 1./h
              dudx  = dudx + geome(inode  ,ielem)*unkno(2,in)*hinv
              dudy  = dudy + geome(inode+3,ielem)*unkno(2,in)*hinv
              dvdx  = dvdx + geome(inode  ,ielem)*unkno(3,in)*hinv
              dvdy  = dvdy + geome(inode+3,ielem)*unkno(3,in)*hinv
           enddo
        
           if (nfric.le.8) then

               h     = unkne(1,ktotg)
               hinv  = 0.0
               if (h.ge.comp) hinv=1./h
               tauy  = tauy0*cos(fi)+(pres*hinv)*sin(fi)
               hinv2 = hinv*hinv

               if (tauy.ge.1.e-3) then
                   hvmod = sqrt(hu*hu+hv*hv)
                   ahat  = 6.0*constK*hvmod*hinv2/tauy
                   gi    = 114.+32.*ahat
                   gi    =(gi-sqrt(gi*gi-(4.*48.*65.)))/(2.*48.)
               else
                 gi=0.0
               endif

               dxz = 3./(2.+gi)*hu*hinv2
               dyz = 3./(2.+gi)*hv*hinv2

               denom(ielem) = sqrt(dudx**2.+dvdy**2.+dudx*dvdy + &
                            ((dudy+dvdx)*0.5)**2.+dxz**2.+dyz**2.)
                            
!                       Neglect Di3 terms:
!                       denom(ielem) = sqrt(0.5*dudx**2.+0.5*dvdy**2.+
!                                     ((dudy+dvdx)*0.5)**2.)

               if (denom(ielem).ge.1.e-6)then
                   visco2=(tauy/denom(ielem)+2.*visco)*h/dens
                   fluxe2(2) = fluxe2(2) - visco2*dudx
                   fluxe2(3) = fluxe2(3) - visco2*(dudy+dvdx)*0.5
                   fluye2(2) = fluxe2(3)
                   fluye2(3) = fluye2(3) - visco2*dvdy
               endif
               
           endif
          
       endif                    !  ----end of viscous contribution

!   -- Source term contribution

       h              = unkne(1,ktotg)  
       hfrict         = max(h,hfrict0)
       hvmod          = sqrt(hu*hu+hv*hv)
       vmod           = hvmod*hinv
       hvmodinv       = 0.0
       if (hvmod.ge.1.e-5) hvmodinv = 1./hvmod

        
       mini_topo(:) = topel (:,ielem)
       mini_patm(:) = patmel(:,ielem)
   
!      ------  interaction waves in reservoirs with sph
!              we add to minitopol 2,3  gradhs in fem nodes
!              if there is interaction, WIR_Interact = .TRUE.
!              in frictlaws, slopx = -g.(h-HH)*(Zx+gradhs)

       if (WIR_Interact) then
           mini_topo(2) = mini_topo(2) + gradhs_w_elem(1,ielem)
           mini_topo(3) = mini_topo(3) + gradhs_w_elem(2,ielem)
       endif
           
       
       call frictlaws (mini_topo, ntopo, mini_patm, npatma, &
                       h, hu, hv, pwp, HH, factw,           &
                       frictx, fricty,                      &
                       slopx, slopy, correx, correy            )  

!   -----------Sources

       if(ntype.eq.3.and.nsour.eq.1) then
          gelem2(1) =  0.0
          gelem2(2) =  cori*hv - frictx + correx
          gelem2(3) = -cori*hu - fricty + correy
       else
          gelem2(1) =  0.0
          gelem2(2) =  cori*hv - frictx + slopx + correx
          gelem2(3) = -cori*hu - fricty + slopy + correy
       endif
       if (icpwpF.eq.1) then
          gelem2(4) = - Bfact*hinv*hinv*unkne(4,ielem)/hrelpw2   
       endif

!   If sources are integrated separately:

       if (ntype.eq.3.and.nsour.eq.2) gelem2=0.0

!   --- Obtain integrals rhs= integral of sources + integral of fluxes
!       rhs1=rhs+ dt_gfl*gelem* Weigth gausspt pg + dt_gfl...

          do inode = 1,nnode
             ipoin = intmat(inode,ielem)
             anx   = geome(inode,ielem)
             any   = geome(inode+3,ielem)
             xni   = shapef(inode)
             do ia = 1,namat
                rhs1(ia,ipoin) = rhs1(ia,ipoin)  + (anx*fluxe2(ia) + any*fluye2(ia) &
                                 +  xni*gelem2(ia))*dvolu*dt_gfl
             end do
          end do
          
       ENDDO                    !  --  Gauss points
       
    ENDDO                       !  --  Elements

!      ------  Boundary integral

   mmatb      = 1.0
   mmatb(1,1) = 2.0
   mmatb(2,2) = 2.0

   raiz3    = sqrt(3.0)
   DO ieleb = 1,neleb
      cosx     = geomb(1,ieleb)       !  ----  length and cosines
      cosy     = geomb(2,ieleb)
!     aleng    = geomb(3,ieleb)*0.5   !  Para calcular la integral con varios
!                                        ptos de gauss
      aleng    = geomb(3,ieleb)
      coe      = aleng/6.
      ipoin1   = intmab(1,ieleb)
      ipoin2   = intmab(2,ieleb)
      ielem    = intmab(3,ieleb)
      iaux     = bindex(ieleb)

!    Boundary integral only 1 gauss point

     do inodb = 1,2
        in    = intmab(inodb,ieleb)
        do jnodb = 1,2
           jn    = intmab(jnodb,ieleb)
           cm    = mmatb (inodb,jnodb)*coe*dt_gfl
           do ia = 1,namat
              rhs1(ia,in) = rhs1(ia,in) - cm*(fluxp(ia,jn)*cosx & 
                                            + fluyp(ia,jn)*cosy)
              aux=rhs1(ia,in)
           enddo
        enddo
     enddo

END DO

END subroutine second_step

!-------------------------------------------------------------------------

        subroutine sources_RK (n1_RK)

!   -------------------------------------------------------------------------

implicit none

integer(ink) n1_RK, ipoin,  nfric, ia, ipres, i, iload, in, iabso 

real   (irk) cgra, dens, cmanning, cori, tauy0  
real   (irk) constK, tanfi8, hfrict0, tanfi0, Bfact
real   (irk) hrelpw2, rhog, c73, tanfi,  gtanfi, gtanfi8
real   (irk) tww, factw 
real   (irk) h, hu, hv, pwp_new,  HH 
real   (irk) ux, uy, pwp_old, hinv, hfrict, frictx, fricty  
real   (irk) hvmod, vmod, hvmodinv, amodl
real   (irk) anx, any, vx, vy, vmd, sqhg, rn
real   (irk) slopx, slopy, correx, correy, pwp 
 
real   (irk) coef(3)
real   (irk) mini_topo(ntopo), mini_patm(npatma)   !  -- compatibility call to frictlaws

data coef/4.0,3.0,2.0/

!      ------  Get fluxes and sources at nodes

cgra     = constF(1)
dens     = constF(2)
cmanning = constF(3)
cori     = constF(4)
nfric    = constF(5) + 0.01
tauy0    = constF(6)
constK   = constF(7)
tanfi8   = constF(9)
hfrict0  = constF(10)
tanfi0   = constF(12)
Bfact    = constF(13)

if (icpwpF.eq.1) hrelpw2  = constF(14)*constF(14) 
rhog     = dens*cgra
c73      = 7./3.

gpoin = 0.

n1_RK    = n1_RK + 1

if (Bfact.ge.1.e-6) then
  tanfi  = tanfi8 + (tanfi0-tanfi8)*exp(-time/Bfact)
else
  tanfi  = tanfi8
endif
gtanfi   = cgra*tanfi
gtanfi8  = cgra*tanfi8

tww        = twind(3)
if(tww.lt.0.0001) tww = 0.0001
factw      = time_gfl/tww
if(factw.gt.1.00) factw=1.00

DO IPOIN = 1,npoinF  

   h           = unkno(1,ipoin)
   if(h.lt.comp) CYCLE
   hu      = unkno(2,ipoin)
   hv      = unkno(3,ipoin)
   pwp_new = 0.0

   Ux = hu
   Uy = hv

   if (icpwpF.eq.1) then          ! MP
       pwp_new = unkno(4,ipoin)
       pwp_old = unkno(4,ipoin)
   endif
    
   hinv=1./h

   
   HH      = 0.0
   if (ntype.eq.2)   then               ! Interaction SW SPH
       HH = swl-topol(1,ipoin)
       if (WIR_Interact) HH = HH - hs_w(ipoin)
   endif

   hfrict      = max(h,hfrict0)

   do ia = 1,3
      frictx = 0.
      fricty = 0.
      hvmod       = sqrt(hu*hu+hv*hv)
      vmod        = hvmod*hinv
      hvmodinv    = 0.0
      if (hvmod.ge.1.e-5) hvmodinv = 1./hvmod
      
      mini_topo(:) = topol(:,ipoin)
      mini_patm(:) = patma(:,ipoin)   
   
!      ------  interaction waves in reservoirs with sph
!              we add to minitopol 2,3  gradhs in fem nodes
!              if there is interaction, WIR_Interact = .TRUE.
!              in frictlaws, slopx = -g.(h-HH)*(Zx+gradhs)

      if (WIR_Interact) then
         mini_topo(2) = mini_topo(2) + gradhs_w(1,ipoin)
         mini_topo(3) = mini_topo(3) + gradhs_w(2,ipoin)
      endif
          
      call frictlaws ( mini_topo, ntopo, mini_patm, npatma, &
                       h, hu, hv, pwp, HH, factw,           &
                       frictx, fricty,                      &
                       slopx, slopy, correx, correy            )  
                       
      gpoin(1,ipoin) =  0.0
      gpoin(2,ipoin) =  cori*hv - frictx + slopx + correx
      gpoin(3,ipoin) = -cori*hu - fricty + slopy + correy
      hu      = Ux      + (dt_gfl/2.)*gpoin(2,ipoin)/coef(ia)
      hv      = Uy      + (dt_gfl/2.)*gpoin(3,ipoin)/coef(ia)
      if (icpwpF.eq.1) then
          gpoin(4,ipoin) = -Bfact*hinv*hinv*pwp_new/hrelpw2         
          pwp_new = pwp_old + (dt_gfl/2.)*gpoin(4,ipoin)/coef(ia)    
      endif
   enddo

   frictx = 0.
   fricty = 0.
   hvmod       = sqrt(hu*hu+hv*hv)
   vmod        = hvmod*hinv
   hvmodinv    = 0.0
   if (hvmod.ge.1.e-5) hvmodinv = 1./hvmod
 
         
   mini_topo(:) = topol(:,ipoin)
   mini_patm(:) = patma(:,ipoin)
      
   !      ------  interaction waves in reservoirs with sph
   !              we add to minitopol 2,3  gradhs in fem nodes
   !              if there is interaction, WIR_Interact = .TRUE.
   !              in frictlaws, slopx = -g.(h-HH)*(Zx+gradhs)
   
   if (WIR_Interact) then
      mini_topo(2) = mini_topo(2) + gradhs_w(1,ipoin)
      mini_topo(3) = mini_topo(3) + gradhs_w(2,ipoin)
   endif
    
    
   call frictlaws (mini_topo, ntopo, mini_patm, npatma, &
                   h, hu, hv, pwp, HH, factw,           &
                   frictx, fricty,                      &
                   slopx, slopy, correx, correy            ) 

   gpoin(1,ipoin) =  0.0
   gpoin(2,ipoin) =  cori*hv - frictx + slopx + correx
   gpoin(3,ipoin) = -cori*hu - fricty + slopy + correy
   if (icpwpF.eq.1)  gpoin(4,ipoin) = -Bfact*hinv*hinv*pwp_new/hrelpw2        

   IF (nprer.ne.0) THEN         !  --  Do not add nodes with presc velocities
       do ipres = 1,nprer
          in    = bprer(1,ipres)
          if(in.eq.ipoin.and.iabs(bprer(2,ipres)).eq.2)then
             gpoin(2,ipoin) = 0.
          endif
          if(in.eq.ipoin.and.iabs(bprer(2,ipres)).eq.3) then
             gpoin(3,ipoin) = 0.
          endif
          if(in.eq.ipoin.and.bprer(2,ipres).eq.4) then
             gpoin(4,ipoin)=0.
          endif
       enddo
   ENDIF

   unkno(2,ipoin) = Ux      + (dt_gfl/2.)* gpoin(2,ipoin)
   unkno(3,ipoin) = Uy      + (dt_gfl/2.)* gpoin(3,ipoin)
   if (icpwpF.eq.1) unkno(4,ipoin) = pwp_old + (dt_gfl/2.)* gpoin(4,ipoin)  

ENDDO

if(n1_RK.eq.2)then

   I3000: DO i = 1,nboun

       amodl = abs(rtang(1,i))+abs(rtang(2,i))
       if(amodl.lt.1.e-06) CYCLE I3000

       ipoin = bconl(i)
       
       ILOAD3500: do iload = 1,nload    !  --  skip abs. & loaded points
                in = bload(iload)
                if(in.eq.ipoin)          CYCLE I3000
       ENDDO ILOAD3500

       do iabso = 1,nabso
          in    = babso(iabso)
          if(in.eq.ipoin)                CYCLE I3000
       enddo

       do iabso = 1,nabshz
       
          in    = babhz(iabso)
          
          if (in.eq.ipoin) then
          
             anx  =  rtang(2,i)
             any  = -rtang(1,i)
             Vx   =  unkno(2,ipoin)
             Vy   =  unkno(3,ipoin)
             vmd  =  sqrt(vx**2.+vy**2.)
             sqhg =  sqrt(cgra*unkno(1,ipoin))
             rn   =  vx*anx+vy*any
             if(rn.lt.0.0.and.vmd.gt.sqhg) then !  --  Only supercritical
                  unkno(2,ipoin) = 0.
                  unkno(3,ipoin) = 0.
             end if
             CYCLE I3000
             
          end if
       enddo


       IF (nprer.ne.0) THEN
           do ipres = 1,nprer
              in    = bprer(1,ipres)
              if(in.eq.ipoin)            CYCLE I3000
           enddo
       endif

!   -- Here we have a wall node!

       anx =  rtang(2,i)
       any = -rtang(1,i)
       rn  =  unkno(2,ipoin)*anx + unkno(3,ipoin)*any
       unkno(2,ipoin) = unkno(2,ipoin)-rn*anx
       unkno(3,ipoin) = unkno(3,ipoin)-rn*any


   enddo I3000 
endif

END subroutine sources_RK



! -------------------------------------------------------------------

        subroutine solve 

! -------------------------------------------------------------------

implicit none


integer(ink) i, ia,iabso, iamat, iiter, iload, in, inode, iunkn 
integer(ink)  ipoin, ip, ipres,  jnode, jpoin, ielem, i0 

real   (irk) h, coe, cgra

real   (irk) amodl,anx, any, gra, cm, rn, Vx, Vy, vmd, sqhg
real   (irk) mmat(3,3)

cgra = constF(1)

!      ------  Get Mass matrix

mmat = 1.0
do inode = 1,nnode
   mmat(inode,inode) = 2.0
enddo

do ipres = 1,nprer
   ipoin = bprer(1,ipres)
   iunkn = bprer(2,ipres)
   i0    = iabs(iunkn)
   rhs1(i0,ipoin) = 0.0
enddo

I13: DO i = 1,nboun

   ipoin = bconl(i)
   do iabso = 1,nabso
      in = babso(iabso)
      if(in.eq.ipoin)           CYCLE I13
   end do
   do iabso = 1,nabshz
      in = babhz(iabso)
      if(in.eq.ipoin)           CYCLE I13
   end do
   do iload = 1,nload
      in = bload(iload)
      if(in.eq.ipoin)           CYCLE I13
   end do
   do ipres=1,nprer
      in=bprer(1,ipres)
      if(in.eq.ipoin)           CYCLE I13
   end do  

   anx  =  rtang(2,i)          !   -- It is a wall node
   any  = -rtang(1,i)
   rn   = rhs1(2,ipoin)*anx+rhs1(3,ipoin)*any
   rhs1(2,ipoin) = rhs1(2,ipoin) - rn*anx
   rhs1(3,ipoin) = rhs1(3,ipoin) - rn*any
   
 enddo I13

!      ------  u(n+1)=u(n)+inv(ML)*(Rhs-M*u(n))
!              where delun2 -> u & u(0)=0
!              and   M*u(n) -> rhs2

rhs0 = 0.0

DO iiter = 1,niter

   IF (iiter.eq.1) THEN
       do ipoin = 1,npoinF
          cm    = mmatl(ipoin)
          do iamat = 1,namat
             rhs0(iamat,ipoin) = rhs1(iamat,ipoin)/cm
          enddo
       enddo
   ELSE
       rhs2=0.0
       
       do ielem = 1,nelem
          coe   = geome(7,ielem)/24.
          do inode = 1,nnode
             ipoin = intmat(inode,ielem)
             do jnode = 1,nnode
                jpoin = intmat(jnode,ielem)
                cm    = mmat(inode,jnode)
                do iamat = 1,namat
                   rhs2(iamat,ipoin) =  rhs2(iamat,ipoin) & 
                                       + cm*coe*rhs0(iamat,jpoin)
                enddo
             enddo
          enddo
       enddo

       do ipoin = 1,npoinF
          cm = mmatl(ipoin)
          do iamat = 1,namat
             rhs0(iamat,ipoin) = rhs0(iamat,ipoin) & 
                                + (rhs1(iamat,ipoin)-rhs2(iamat,ipoin))/cm
          enddo
       enddo
   ENDIF

!      ------  Apply BC's : (I) prescribed values

   IF (nprer.ne.0) THEN
      do ipres = 1,nprer
         ipoin = bprer(1,ipres)
         iunkn = bprer(2,ipres)
         h     = unkno(1,ipoin)+rhs0(1,ipoin)
         if(iunkn.eq.2)then
            rhs0(2,ipoin) = delun(2,ipoin)     !    *h
         elseif(iunkn.eq.3)then
            rhs0(3,ipoin) = delun(3,ipoin)     !    *h
         elseif ((iunkn.eq.-2).or.(iunkn.eq.-3)) then
            i0 = iabs(iunkn)
            rhs0(i0, ipoin) = delun(i0, ipoin)                                       
         else
            rhs0(iunkn,ipoin) = delun(iunkn,ipoin)
         end if
      enddo
   ENDIF

!      ------  filter out the prescribed vnormal = 0 - conditions

!                       **  check iadda = 1
   I3000: DO i = 1,nboun

      amodl = abs(rtang(1,i))+abs(rtang(2,i))
      if(amodl.lt.1.e-06)       CYCLE I3000  

      ipoin = bconl(i)

      do  iload = 1,nload             !      ------  skip abs. & loaded points
         in = bload(iload)
         if(in.eq.ipoin)        CYCLE I3000 
      enddo

      do iabso = 1,nabso
         in = babso(iabso)
         if(in.eq.ipoin)        CYCLE I3000 
      enddo

     if(ntype.eq.3)then
        do iabso = 1,nabshz
           in = babhz(iabso)
           if (in.eq.ipoin) then
              anx  = rtang(2,i)
              any  =-rtang(1,i)
              Vx   = unkno(2,ipoin) + rhs0(2,ipoin)
              Vy   = unkno(3,ipoin) + rhs0(3,ipoin)
              vmd  = sqrt(vx**2.+vy**2.)
              sqhg = sqrt(cgra*unkno(1,ipoin))
              rn   = vx*anx+vy*any
              if(rn.lt.0.0.and.vmd.gt.sqhg) then ! Only supercritical conds
                 rhs0(2,ipoin) = -unkno(2,ipoin)
                 rhs0(3,ipoin) = -unkno(3,ipoin)
              end if
              CYCLE I3000 
           end if
        enddo
     endif

     do ipres = 1,nprer
        in    = bprer(1,ipres)
        if(in.eq.ipoin)         CYCLE I3000 
     enddo

!   -- If code exec arrives here, it is a wall node

     anx =  rtang(2,i)
     any = -rtang(1,i)
     rn  = rhs0(2,ipoin)*anx + rhs0(3,ipoin)*any
     rhs0(2,ipoin) = rhs0(2,ipoin) - rn*anx
     rhs0(3,ipoin) = rhs0(3,ipoin) - rn*any

   ENDDO I3000

ENDDO  ! ------------------------------  ends loop in iterations

!      ------  transcribe rhelp onto delun

do  ip=1,npoinF
    do  ia=1,namat
        delun(ia,ip)= rhs0(ia,ip)
    enddo
enddo

 END subroutine solve

! ------------------------------------------------------------

        subroutine getmmatl 
                    
! ------------------------------------------------------------

implicit none 

integer(ink) ielem, inode, in 
real   (irk) area3

 mmatl=0.0   ! --- Initializes ML

 do ielem = 1,nelem
    area3 = geome(7,ielem)/6.           ! --- area of element
    do inode = 1,nnode
       in = intmat(inode,ielem)
       mmatl(in) = mmatl(in) + area3
    enddo
 enddo

END subroutine getmmatl


!-------------------------------------------------------

        subroutine ouplot_mesh

!------------------------------------------------------

!      ------  Mesh output ------------ 

implicit none

integer(ink)  idimn, ielem, ipoin, inode 

write (gidg_msh,*) 'MESH geoflow   dimension   3    ElemType Triangle  Nnode 3 '

write (gidg_msh,*) ' coordinates'
do ipoin = 1,npoinF
   write(gidg_msh,*) ipoin,(coord(idimn,ipoin),idimn=1,ndimnF), topol(1,ipoin)
enddo
write(gidg_msh,*) ' End coordinates' 

write (gidg_msh,*) ' Elements'
do ielem=1,nelem
   write(gidg_msh,*) ielem, (intmat(inode,ielem),inode=1,nnode),' 1'
enddo 
write (gidg_msh,*) ' End elements'   
      
close(gidg_msh)

write (gidg_res,*) 'GiD Post Results File 1.0'
   
END subroutine ouplot_mesh

  

!-------------------------------------------------------

        subroutine ouplot_res

!------------------------------------------------------

!      ------  res output ------------ 

implicit none

integer(ink)  idim, ielem, ipoin  
integer(ink)  inode, idimn 
real   (irk)   grav, froude, h, hinv, u, v, vmod, zxx, zxy, zyy
real   (irk)  valor, depth, pwp


!      ------  Write h or eta (icoutp=2 or 1)

write(gidg_res,*) 'Result "height" "Height" ',time_gfl,' Vector OnNodes "where" '
write(gidg_res,*) ' Values '
if (icoutp.eq.1) then
    do ipoin = 1,npoinF
       h = unkno(1,ipoin)
       h =(h+topol(1,ipoin))-swl
       write(gidg_res,*) ipoin, ' 0.   0.   ', h
    enddo
else
    do ipoin = 1,npoinF
       h = unkno(1,ipoin)
       write(gidg_res,*) ipoin, ' 0.   0.   ', h
    enddo
endif
write(gidg_res,*) ' End Values ' 


!      ------  Write velocities

write(gidg_res,*) 'Result "vel" "veloc" ',time_gfl,' Vector OnNodes "where" '
write(gidg_res,*) ' Values '

do ipoin = 1,npoinF
   h     = unkno(1,ipoin)
   if (h.ge.comp) then
       u = unkno(2,ipoin)/h
       v = unkno(3,ipoin)/h
   else
       u = 0.0
       v = 0.0
   endif
   write(gidg_res,7)  ipoin,u,v, '  0.0  '
enddo
write(gidg_res,*) ' End Values ' 

!      ------  And PWPs  -----------------------------

if (icpwpF.eq.1) then
   write(gidg_res,*) 'Result "Pwp" "pwpress" ',time_gfl,' Scalar OnNodes "where" '
   write(gidg_res,*) ' Values '
   do ipoin = 1, npoinF
      h     = unkno(1,ipoin)
      if (h.ge.comp) then
         pwp = unkno(4,ipoin)/h
      else
         pwp = 0.0
      endif 
      write(gidg_res,*) ipoin, pwp
   enddo    
   write(gidg_res,*) ' End Values '
endif 

!      ------  depths  -----------------------------


write(gidg_res,*) 'Result "Dpth Z" "depth Z" ',time_gfl,' Scalar OnNodes "where" '
write(gidg_res,*) ' Values '
if (icoutp.eq.1) then
    do ipoin = 1,npoinF
       depth = swl - topol(1,ipoin)
       write(gidg_res,*) ipoin,depth
    enddo
else
   do ipoin = 1,npoinF
      write(gidg_res,*) ipoin, topol(1,ipoin)
   enddo
endif
write(gidg_res,*) ' End Values ' 

 7 format(i5,2x,6(g11.4,1x))
91 format(a12,i5,(1x,e15.8),3i5)

 
END subroutine ouplot_res

!-----------------------------------------

        subroutine output

!------------------------------------------

implicit none

integer(ink) ic, ipoin, in 
real   (irk) h, u, v, vmod 
real   (irk) aux(nprpoi,3)

do  ic =1,nprpoi
    if(lprpoi(ic).eq.0)     lprpoi(ic) = 1
    if(lprpoi(ic).gt.npoinF) lprpoi(ic) = npoinF
    ipoin   = lprpoi(ic)
    h       = unkno(1,ipoin)
    if (h.ge.comp) then
       u = unkno(2,ipoin)/h
       v = unkno(3,ipoin)/h
    else
       u = 0.0
       v = 0.0
    endif
    vmod      = (u*u+v*v)**0.5
    
    if (icoutp.eq.1) h = (h+topol(1,ipoin)) - swl

    aux (ic,1) = h
    aux (ic,2) = u
    aux (ic,3) = unkno(2,ipoin)

enddo

if (time_gfl.eq.0) then
    write(resg_file1,*) 'nodo(s):',(lprpoi(in),in=1,nprpoi) 
endif

write(resg_file1,10) time_gfl,(aux(in,1),aux(in,2),aux(in,3),in = 1,nprpoi)

 10    format(21(2x,f10.4)) !solo podemos escribir los diez
                            !primeros nodos de control
 11    format(2x,a8,20(2x,i10,a1,2x,i10,a4))

END subroutine output

!-------------------------------------------------------

        subroutine ouplot

!------------------------------------------------------

!                  OLD VERSION NOT USED HERE
implicit none

integer(ink)  idim, ielem, ipoin 
integer(ink)  inode, idimn 
real   (irk)   grav, froude, h, hinv, u, v, vmod, zxx, zxy, zyy
real   (irk)  valor, depth, pwp



       !      ------  salida para gid ---------------------------------

if(itimestep_gfl.eq.0) then

    write(gidg_msh,*) ' line 1: output of GEOFLOW2000'
    write(gidg_msh,*) ' line 2: output of GEOFLOW2000'
    write(gidg_msh,*) ' line 3: output of GEOFLOW2000'
    write(gidg_msh,*) ' line 4: output of GEOFLOW2000'
    write(gidg_msh,*) ' line 5: output of GEOFLOW2000'
    write(gidg_msh,*) ' nelem   npoin   nelem_type '
    write(gidg_msh,*)   nelem, npoinF, nnode
    write(gidg_msh,*) ' ----  coordenadas ----'
    do ipoin = 1,npoinF
       write(gidg_msh,*) ipoin,(coord(idim,ipoin),idim=1,ndimnF)
    enddo
    write(gidg_msh,*) ' --- ielem   n1   n2   n3   imat -----'
    
    do ielem=1,nelem
       write(gidg_msh,*) ielem, (intmat(inode,ielem),inode=1,nnode),' 1'
    enddo
    
end if

write(gidg_res,*) '   deri2       1 ', time_gfl, ' 3   1   0'   !  ---- Curvatures

do ipoin = 1,npoinF
   zxx = topol(4,ipoin)
   zxy = topol(5,ipoin)
   zyy = topol(6,ipoin)
   write(gidg_res,7)  ipoin,zxx,zxy,zyy
enddo
    
write(gidg_res,*) '   veloc       1 ', time_gfl, ' 2   1   0'   !  ---- Velocity

    do ipoin = 1,npoinF
       h     = unkno(1,ipoin)
       if (h.ge.comp) then
           u = unkno(2,ipoin)/h
           v = unkno(3,ipoin)/h
       else
           u = 0.0
           v = 0.0
       endif
       write(gidg_res,7)  ipoin,u,v
    enddo

write(gidg_res,*) '   froud       1 ', time_gfl, '  1  1    0'   !---- Velocity
    grav      = constF(1)
    do ipoin  = 1,npoinF
       h      = unkno(1,ipoin)
       hinv   = 0.0
       if (h.ge.comp) hinv = 1./h
       u      = unkno(2,ipoin)*hinv
       v      = unkno(3,ipoin)*hinv
       vmod   = (u*u+v*v)**0.5
       froude = vmod*(hinv/grav)**0.5
       write(gidg_res,*)  ipoin,froude
    enddo

write (gidg_res,*)'   eta         1 ', time_gfl, '  1  1    0'  ! ---- ETA ----
if (icoutp.eq.1) then
    do ipoin = 1,npoinF
       h = unkno(1,ipoin)
       h =(h+topol(1,ipoin))-swl
       write(gidg_res,*) ipoin, h
    enddo
else
    do ipoin = 1,npoinF
       h = unkno(1,ipoin)
       write(gidg_res,*) ipoin, h
    enddo
endif

if (icpwpF.eq.1) then
   write (gidg_res,*)'   pwp         1 ', time_gfl, '  1  1    0'  ! ---- pwp ----
      do ipoin = 1,npoinF
         h   = unkno(1,ipoin)
         if (h.ge.comp) then
            pwp = unkno(4,ipoin)/h
         else
            pwp = 0.0
         endif 
         write(gidg_res,*) ipoin, pwp
      enddo
endif

write (gidg_res,*)'   dpth        1 ', time_gfl, '  1  1    0'  ! ---- TOPOLOGY
if (icoutp.eq.1) then
    do ipoin = 1,npoinF
       depth = swl - topol(1,ipoin)
       write(gidg_res,*) ipoin,depth
    enddo
else
   do ipoin = 1,npoinF
      write(gidg_res,*) ipoin, topol(1,ipoin)
   enddo
endif

!      ------  salida para gid 3d---------------------------------

rewind(gid3g_msh)

write (gid3g_msh,*) ' line 1: output of GEOFLOW2000'
write (gid3g_msh,*) ' line 2: output of GEOFLOW2000'
write (gid3g_msh,*) ' line 3: output of GEOFLOW2000'
write (gid3g_msh,*) ' line 4: output of GEOFLOW2000'
write (gid3g_msh,*) ' line 5: output of GEOFLOW2000'
write (gid3g_msh,*) 'nelem   npoin   nelem_type'
write (gid3g_msh,*) nelem, npoinF,7
write (gid3g_msh,*) ' ----  coordenadas ----'
do ipoin = 1,npoinF
   if (icoutp.lt.3) then
       write(gid3g_msh,*) ipoin,(coord(idim,ipoin),idim=1,ndimnF),0.
   elseif (icoutp.eq.3) then
       write(gid3g_msh,*) ipoin,(coord(idim,ipoin),idim=1,ndimnF), &
                         topol(1,ipoin)
   endif
enddo
write(gid3g_msh,*) ' --- ielem   n1   n2   n3   imat -----'
do ielem = 1,nelem
   write(gid3g_msh,*) ielem, (intmat(inode,ielem),inode=1,nnode),' 1'
enddo

write (gid3g_res,91)'   eta 3d       ',1, time_gfl,2,1,0  ! ---- 3d ----
if (icoutp.eq.1) then
    do ipoin=1,npoinF
       h = unkno(1,ipoin)
       h =(h+topol(1,ipoin))-swl
       write(gid3g_res,*) ipoin,0.,0., h
    enddo
elseif (icoutp.eq.2) then
    do ipoin=1,npoinF
       h = unkno(1,ipoin)
       write(gid3g_res,*) ipoin,0.,0., h
    enddo
elseif (icoutp.eq.3) then
    do ipoin=1,npoinF
       h = unkno(1,ipoin)
       write(gid3g_res,*) ipoin,0.,0., h
    enddo
endif


if (icpwpF.eq.1) then      
        write (gid3g_res,91)'   pwp 3d       ',1, time_gfl,2,1,0  ! ---- 3d ----
      if (icoutp.eq.1) then
        do ipoin=1,npoinF
                 h=unkno(1,ipoin)
                         pwp = unkno(4,ipoin)/h
           write(gid3g_res,*) ipoin,0.,0., pwp
        enddo
    elseif (icoutp.eq.2) then
        do ipoin=1,npoinF
                     h=unkno(1,ipoin)
                         pwp = unkno(4,ipoin)/h
           write(gid3g_res,*) ipoin,0.,0., pwp
        enddo
    elseif (icoutp.eq.3) then
        do ipoin=1,npoinF
                     h=unkno(1,ipoin)
                         pwp = unkno(4,ipoin)/h
           write(gid3g_res,*) ipoin,0.,0., pwp
        enddo
    endif
endif               

if (nacti.eq.1) then
     write (gid3g_res,91)'   nodo activo       ',1, time_gfl,2,1,0  ! ---- 3d ----
     if (icoutp.eq.1) then
         do ipoin = 1,npoinF
            valor = iactnd(ipoin)
            write(gid3g_res,*) ipoin,0.,0., valor
         enddo
     elseif (icoutp.eq.2) then
         do ipoin = 1,npoinF
            valor = iactnd(ipoin)
            write(gid3g_res,*) ipoin,0.,0., valor
         enddo
     elseif (icoutp.eq.3) then
         do ipoin = 1,npoinF
            valor = iactnd(ipoin)
            write(gid3g_res,*) ipoin,0.,0., valor
         enddo
     endif
endif

 7 format(i5,2x,6(g11.4,1x))
91 format(a12,i5,(1x,e15.8),3i5)

 
END subroutine ouplot

! ----------------------------------------------------------------

       subroutine addthem  
       
! ----------------------------------------------------------------

implicit none

integer(ink) i, ia, iprer

!      --------adding the increments

DO i = 1, npoinF

   DO ia = 1,namat
      unkno(ia,i) = unkno(ia,i) + delun(ia,i)       
   ENDDO   
   
ENDDO   

 
END subroutine addthem

! --------------------------------------------------------

       subroutine getpatm 
       
! --------------------------------------------------------

implicit none 

integer(ink) icatm, ipoin, ie, ip, it, it1, in
real   (irk) x0, y0, grad, rend, vxatm, vyatm, tlimit, factt
real   (irk) vatm, anvx, anvy, xc, yc, fact
real   (irk) x,y,rx, ry, r, d, c13, anx, any 

icatm = datm(1)+0.1
x0    = datm(2)
y0    = datm(3)
grad  = datm(4)
rend  = datm(5)
vxatm = datm(6)
vyatm = datm(7)
tlimit= datm(8)
if(tlimit.le.0.0001) tlimit = 0.0001
factt = time_gfl/tlimit
if(factt.gt.1.) factt = 1.0
if(icatm.eq.1.or.icatm.eq.2) factt = 1.0
vatm  = (vxatm*vxatm+vyatm*vyatm)**0.5
if(vatm.le.0.0001) vatm=0.0001
anvx = vxatm/vatm
anvy = vyatm/vatm

if(icatm.lt.3) then
   xc = x0 + vxatm*time_gfl
   yc = y0 + vyatm*time_gfl
else
   xc = x0
   yc = y0
endif

DO ipoin = 1,npoinF
   x  = coord(1,ipoin)
   y  = coord(2,ipoin)
   rx = (x-xc)
   ry = (y-yc)
   if (icatm.eq.1.or.icatm.eq.3) then
      r = (rx*rx+ry*ry)**0.5
      if(r.gt.0.0001) then
          anx = rx/r
          any = ry/r
      else
          anx = anvx
          any = anvy
      endif
      fact = 1.0
      if(r.gt.rend) fact = 0.0
         patma(2,ipoin) = grad*anx*fact*factt
         patma(3,ipoin) = grad*any*fact*factt
      endif
   if (icatm.eq.2.or.icatm.eq.4) then
      d    = rx*anvx+ry*anvy
      fact = 0.0
      if(d.lt.0..and.abs(d).le.rend) fact=1.0
      patma(2,ipoin)=grad*anvx*fact*factt
      patma(3,ipoin)=grad*anvy*fact*factt
   endif
ENDDO

!      ------  get p,x     & p,y   averaged on elements

c13    = 1./3.
patmel = 0.

do ie  = 1,nelem
   do in = 1,nnode
      ip = intmat(in,ie)
      do it = 1,2
         it1 = it+1
         patmel(it1,ie) = patmel(it1,ie)+c13*patma(it1,ip)
      enddo
   enddo
enddo


END subroutine getpatm

! ----------------------------------------------------------------

       subroutine check    
       
! ----------------------------------------------------------------

implicit none

integer(ink) ipoin, iprer, iunk0
real   (irk) h, pwpmax, pwpmin
real   (irk) ff(namat)

pwpmin = comp        !MP

do ipoin = 1,npoinF

   h    = unkno(1,ipoin)
   
   if (icpwpF.eq.1) then                                   !MP
       pwpmax = h !constF(14)*h
       if (unkno(4,ipoin).lt.pwpmin) unkno(4,ipoin)  = 0.0
       if (unkno(4,ipoin).gt.pwpmax)  unkno(4,ipoin) = pwpmax
   endif
   
   if(h.lt.comp) then                   !    if node is prescribed  
      ff = 0.0                          !    we make sure we don't zero
      do iprer = 1, nprer               !    prescribed variables
         if (bprer(1,iprer).eq.ipoin) then
             iunk0 = iabs (bprer(2,iprer))
             ff (iunk0) = 1.0
         endif
      enddo
          unkno(1,ipoin) = ff(1)*unkno(1,ipoin)
      unkno(2,ipoin) = ff(2)*unkno(2,ipoin)
      unkno(3,ipoin) = ff(3)*unkno(3,ipoin)
   endif
   
enddo


END subroutine check


! ------------------------------------------------------------

       subroutine abshz 

! ------------------------------------------------------------

!      ------  this subroutine applies absorbing bc's of the type:
!
!               (i) subcritical
!                     r1 = 2.c+u   r2=2.c-u=2.(g.H)**0.5
!                     h=(r1+r0)^2/(16g)  u=(r1-r0)/2
!               (ii) supercritical
!                     everything is free.
!

implicit none

integer(ink) iabso, ipoin, iboun
real   (irk) Z, depth0, h, hinv, vx, vy, v, g
real   (irk) anx, any, vn, vtx, vty, c, c0, r10, r20, hnew, vnew


g = constF(1)

DO iabso  = 1,nabso
   ipoin  = babso(iabso)
   Z      = topol(1,ipoin)
   if (ntype.eq.2) then
       depth0 = swl-Z
   elseif (ntype.eq.3) then
       depth0 = 0.0
   endif
   h      = unkno(1,ipoin)
   hinv   = 0.0
   if (h.ge.comp) hinv=1./h
   vx     = unkno(2,ipoin)*hinv
   vy     = unkno(3,ipoin)*hinv
   v      = (vx*vx+vy*vy)**0.5
   do  iboun = 1,nboun
       if(bconl(iboun).eq.ipoin) goto 1120
   enddo
1120      anx  = rtang(2,iboun)
   any  =-rtang(1,iboun)
   vn   = vx*anx + vy*any
   vtx  = vx - anx*vn
   vty  = vy - any*vn
   c    = 0.0
   if (h.gt.1.e-5)  c  = (g*h)**0.5
   c0   = 0.0
   if (depth0.gt.1.e-5) c0 = (g*depth0)**0.5
   r10  = 2.*c + vn
   r20  = 2.*c0
   if (vn.gt.0.and.v.gt.c) goto 1130
   hnew = ((r10+r20)**2.)/(16.*g)
   vnew = (r10-r20)/2.
   unkno(1,ipoin) = hnew
   unkno(2,ipoin) = (vtx+vnew*anx)*hnew
   unkno(3,ipoin) = (vty+vnew*any)*hnew
1130      continue

ENDDO


END subroutine abshz

! ------------------------------------------------------------

        subroutine loadbc 

! ------------------------------------------------------------

!      ------  this subroutine applies incoming wave bc's 
!              incident wave (linear & nonlinear)
!                 r2   = (g/h)**0.5 n0 cos(kr-wt)1-cosb
!                 r1   = (g/h)**0.5 n + vn
!
!                 n    = ( r1 + r2 ) / 2.*(g/h)**0.5
!                 vn   = ( r2 - r1 ) / 2.
!

implicit none

integer(ink)  iload, ipoin, iboun  
 
real   (irk)  g, thetax, omega, amplit, bspeed, x0, y0 
real   (irk)  h, gh, eta, vx, vy, pi2, c, tt
real   (irk)  alength, akx, aky, phase, anx, any, cosb, r2, r1
real   (irk)  un  

g       = constF(1)
thetax  = wave (1)
omega   = wave (2)
amplit  = wave (3)
bspeed  = wave (4)
x0      = wave (5)
y0      = wave (6)


DO iload = 1,nload

   ipoin = bload(iload)
   h     = swl - topol(1,ipoin)
   gh    = (g/h)**0.5

   if (ntype.eq.1.or.ntype.eq.3) then
       write(*,*) ' no ntype=1 nor 3'
       write(*,*) ' program stopped at load bcs'
       pause
       eta = unkno(1,ipoin)
       vx  = unkno(2,ipoin)
       vy  = unkno(3,ipoin)
   endif
   if(ntype.eq.2) then
       eta = unkno(1,ipoin)-h
       vx  = unkno(2,ipoin)/unkno(1,ipoin)    
       vy  = unkno(3,ipoin)/unkno(1,ipoin)
   endif

   pi2     = 3.1415*2.
   c       = (g*h)**0.5
   tt      = pi2/omega
   alength = c*tt
   akx     = pi2*cos(thetax)/alength
   aky     = pi2*sin(thetax)/alength
   phase   = coord(1,ipoin)*akx+coord(2,ipoin)*aky-omega*time_gfl

   iboun2010: do  iboun=1,nboun
      if(bconl(iboun).eq.ipoin) EXIT iboun2010
   enddo iboun2010

   anx     =  rtang(2,iboun)
   any     = -rtang(1,iboun)
   cosb    =  anx*cos(thetax) + any*sin(thetax)
   
   r1 = gh * eta + vx*anx + vy*any
   if(cosb.gt.0.) then
      r2 = 0.0      
   else
      r2 = amplit * sin(phase) * (1.-cosb) * gh
   endif   
   eta = (r1 + r2) / (2.*gh)
   un  = (r1 - r2) / 2.

   if (ntype.eq.1) then
      unkno(1,ipoin) = eta
      unkno(2,ipoin) = un*anx
      unkno(3,ipoin) = un*any
   endif
   if(ntype.eq.2) then
      unkno(1,ipoin) = eta + h
      unkno(2,ipoin) = un * anx *unkno(1,ipoin)
      unkno(3,ipoin) = un * any *unkno(1,ipoin)
   endif

ENDDO
       
END subroutine loadbc

! --------------------------------------------------------------------

        subroutine check_delt 

! -------------------------------------------------------------------

!      ------  This routine obtains dt min in nodes.
!              Uses characteristic length 1/L^2= (N,x^2+N,y^2)
!              if dt_gfl<dtcrit stops, and write unknowns for restart

implicit none

integer(ink) i, ielem, inode, ipoin, nfric
real   (irk) auxdt,dt_crit_gfl, h
real   (irk) cgra, dens, visco, hfrict0
real   (irk) rkx, c43, c6, alfa, dtelm
real   (irk) hfrict, hinv, vx, vy, vmod, vmodinv
real   (irk) anx, any, xle, celert, vcomp, dtnod, slope, Pe
real   (irk) auxc, rkxy


if (nacti.eq.1) then            !  --  only when option to use active
                                !      nodes and elements is active
    iactnd = 0                  !  --  switch off nodes
    iactel = 0                  !  --  switch off elements

    Do i = 1,3                  !  --  Do it 3 times, to create a safety zone     

       Do ielem = 1, nelem      !  --  Activate elements with h>comp
          do inode = 1,3
             ipoin = intmat(inode,ielem)
             h = unkno(1,ipoin)
             if (h.ge.comp) iactel(ielem)=1         
          enddo
       Enddo

       Do ielem = 1,nelem
          if (iactel(ielem).eq.1) then
              do inode = 1,3
                 ipoin = intmat(inode,ielem)
                 iactnd(ipoin)=1
              enddo
          endif
       Enddo
   
    Enddo
    
endif    


auxdt          = 1.0e10
dt_crit_gfl    = 1.0e10
deltel         = 0.0  

cgra    = constF(1)
dens    = constF(2)
visco   = constF(8)
hfrict0 = constF(10)

if(constF(3).gt.0.0) auxdt = 2.0/(constF(3)**2.*constF(1))
rkxy  = constF(7)/16.0/constF(2)
nfric = constF(5) + 0.01
c43   = 4.0/3.0
c6    = 1.0/6.0
alfa  = 1.0/sqrt(3.0)
if(niter.eq.1)then
   alfa = 1.0
   c6   = 0.5
end if

DO ielem = 1,nelem
               
   if (iactel(ielem).eq.0) CYCLE
   dtelm = 0.0
   DO inode = 1,nnode
      ipoin = intmat(inode,ielem)
      h     = unkno(1,ipoin)
      if (h.ge.comp) then       !  --  Activate elems with one active node
          iactel(ielem) = 1
          hfrict        = amax1(hfrict0,h)
          hinv          = 1./hfrict
          vx            = unkno(2,ipoin)/h
          vy            = unkno(3,ipoin)/h
          vmod          = sqrt(vx*vx+vy*vy)
          vmodinv       = 1.0e10
          if(vmod.ge.1.0e-5) vmodinv = 1.0/vmod
          anx           = geome(inode  ,ielem)
          any           = geome(inode+3,ielem)
          xle           = 1./sqrt(anx*anx+any*any)
          celert        = sqrt(cgra*h)
          vcomp         = vmod + celert
          dtnod         = 0.9*xle/vcomp*alfa
          slope         = sqrt(topol(2,ipoin)**2+topol(3,ipoin)**2)

          if (abs(visco).gt.1.0e-6) then         !  --  limit by visco
              Pe        = (vcomp*xle*dens)/(2.*visco)
              if(Pe.gt.0.1)then
                 dtnod  = min(dtnod,0.9*xle/vcomp*(sqrt(1./Pe**2.+alfa**2.)-1./Pe))
              else
                 dtnod  = min(dtnod,0.9*xle*xle*c6/visco)
              end if
          end if
           
          if (h.ge.hfrict0) then                !  --  Limit by source terms
              auxc = abs(2.*celert-vmod)
              if (auxc.le.1.e-4) auxc = 2.*celert+vmod
              if ((nfric.eq.1.or.nfric.eq.5.or.nfric.eq.6.or.nfric.eq.7        &
                   .or.nfric.eq.8).and.constF(3).ne.0.0)                        &
                  dtnod = min(dtnod,0.9*2.0*hfrict**c43*auxc/cgra/             &
                    ((vmod*constF(3))**2+hfrict**c43*slope))      !  -- with nfric.eq.5:
                                                                 !     sometimes fails
              if(nfric.eq.4)  dtnod = min(dtnod,0.9*1.0                        & 
                                 /(vmod/(auxdt*hfrict**c43)+rkxy*hinv**2))
          end if

          if (dtnod.lt.dt_crit_gfl) then
              dt_crit_gfl   = dtnod
          endif

          dtelm = dtelm + dtnod
      
      endif
    
  ENDDO
  
  deltel(ielem) = dtelm/3.
  
ENDDO

write(*,*) ' time_gfl = ',time_gfl,' dt = ',dt_gfl,' dtcrit= ',dt_crit_gfl
  

if (ic_adapt_dt_gfl.eq.1) then
    if (dt_crit_gfl.lt.dt_gfl) then
        dt_gfl = dt_crit_gfl
    else
        dt_gfl = amin1(1.1*dt_gfl,dt_crit_gfl)
    endif
endif

if (dt_gfl.gt.dt_crit_gfl) then
    write(*,*)          ' Time increment fixed larger than critical ' 
    write(chkg_file,*) 'Time increment fixed larger than critical '
    STOP
endif


END subroutine check_delt

! ----------------------------------------------------------------

        subroutine weir  

! ----------------------------------------------------------------

implicit none    

integer(ink) iweir, ipoin, iboun
real   (irk) h, Z, ZW,d, hvx, hvy, c1, c2  
real   (irk) anx, any, hvn, hvtx, hvty, hvnweir

DO iweir = 1,nweir
   ipoin = lweir(iweir)
   h     = unkno(1,ipoin)
   Z     = topol(1,ipoin)
   ZW    = zweir(iweir)
   d     = (Z + h - ZW)
   if (d.lt.comp) CYCLE
   hvx   = unkno(2,ipoin)
   hvy   = unkno(3,ipoin)
   c1    = cweir(1)
   c2    = cweir(2)
   do iboun = 1,nboun
      if (bconl(iboun).eq.ipoin) goto 1000
   enddo
1000      continue
   anx     =   rtang(2,iboun)
   any     = - rtang(1,iboun)
   hvn     =   hvx*anx + hvy*any
   hvtx    =   hvx - anx*hvn
   hvty    =   hvy - any*hvn
   hvnweir =   c1*(d**c2)
   unkno(2,ipoin) = hvtx + hvnweir*anx
   unkno(3,ipoin) = hvty + hvnweir*any
ENDDO

END subroutine weir


! -------------------------------------------------------------------

        subroutine smooth 

! -------------------------------------------------------------------

!      Smoothing after OCZs Bible Vol. Fluids

implicit none

integer(ink) ia, ieleb, ielem, ip,ip1, ip2, ip3

real   (irk) anx, any
real   (irk) ce, delt2, fsw, ps1, ps2, ps3, p11, p22, p33 
real   (irk) Swel, xle1, xle2, xle3, xps, xpd

ce   = csmooth
rhs0 = 0.0
rhs1 = 0.0

Do ielem = 1,nelem
   ip1 = intmat(1,ielem)
   ip2 = intmat(2,ielem)
   ip3 = intmat(3,ielem)
   ps1 = unkno (1,ip1)
   ps2 = unkno (1,ip2)
   ps3 = unkno (1,ip3)
   p11 = (ps1-ps2)+(ps1-ps3)
   p22 = (ps2-ps3)+(ps2-ps1)
   p33 = (ps3-ps1)+(ps3-ps2)
   rhs0(1,ip1) = rhs0(1,ip1) + p11
   rhs0(1,ip2) = rhs0(1,ip2) + p22
   rhs0(1,ip3) = rhs0(1,ip3) + p33
   rhs1(1,ip1) = rhs1(1,ip1) + abs(ps1-ps2) + abs(ps1-ps3)
   rhs1(1,ip2) = rhs1(1,ip2) + abs(ps2-ps3) + abs(ps2-ps1)
   rhs1(1,ip3) = rhs1(1,ip3) + abs(ps3-ps1) + abs(ps3-ps2)
Enddo

Do ieleb = 1,neleb
   ip1 = intmab(1,ieleb)
   ip2 = intmab(2,ieleb)
   ps1 = unkno (1,ip1)
   ps2 = unkno (1,ip2)
   !          xps = ps1 + ps2
   xpd = ps1 - ps2
   rhs0(1,ip1) = rhs0(1,ip1) + xpd
   rhs0(1,ip2) = rhs0(1,ip2) - xpd
   rhs1(1,ip1) = rhs1(1,ip1) + abs(xpd)
   rhs1(1,ip2) = rhs1(1,ip2) + abs(xpd)
Enddo

Do ip = 1,npoinF
   if (rhs1(1,ip).lt.0.1*unkno(1,ip)) rhs1(1,ip)=unkno(1,ip)
   if (rhs1(1,ip).lt.1.e-4) rhs1(1,ip)=1.e-4
   rhs0(1,ip) = abs(rhs0(1,ip))/rhs1(1,ip)
Enddo

rhs1 = 0.0

Do ielem = 1, nelem

   delt2 = min(0.5*deltel(ielem),dt_gfl)
   if (delt2.eq.0.) delt2 = dt_gfl
      
   anx   = geome(1  ,ielem)
   any   = geome(4,ielem)
   xle1   = 1./sqrt(anx*anx+any*any)
   anx   = geome(2  ,ielem)
   any   = geome(5,ielem)
   xle2   = 1./sqrt(anx*anx+any*any)
   anx   = geome(3  ,ielem)
   any   = geome(6,ielem)
   xle3   = 1./sqrt(anx*anx+any*any)
         
   ip1   = intmat(1,ielem)
   ip2   = intmat(2,ielem)
   ip3   = intmat(3,ielem)
   Swel  = (rhs0(1,ip1)+rhs0(1,ip2)+rhs0(1,ip3))/ 3.0
   fsw   = Swel*ce/4.0
   
   Do ia = 1, namat
      rhs1(ia,ip1) = rhs1(ia,ip1) & 
                     + fsw*(-2.*unkno(ia,ip1)+unkno(ia,ip2)+unkno(ia,ip3))      &
                          * xle1**2.*dt_gfl/delt2
      rhs1(ia,ip2) = rhs1(ia,ip2) &
                     + fsw*(-2.*unkno(ia,ip2)+unkno(ia,ip3)+unkno(ia,ip1))      &
                          * xle2**2.*dt_gfl/delt2
      rhs1(ia,ip3) = rhs1(ia,ip3) &
                    + fsw*(-2.*unkno(ia,ip3)+unkno(ia,ip1)+unkno(ia,ip2))       &
                          * xle3**2.*dt_gfl/delt2
   Enddo
Enddo

Do ip = 1, npoinF

   Do ia = 1, namat
      unkno(ia,ip) = unkno(ia,ip) + rhs1(ia,ip)
   Enddo
   
Enddo

END subroutine smooth



! -------------------------------------------------------------------

        subroutine frictlaws (mini_topo, ntopo, mini_patm, npatma, &
                              h, hu, hv, pwp, HH, factw,           &
                              frictx, fricty,                      &
                              slopx, slopy, correx, correy            ) 

! -------------------------------------------------------------------
!       

!       Obtain friction, slope and wind stress contributions to sources
!              both for 1st and 2n step

!       Use local vars 
!           mini_topo(ntopo)  <- topol (ntopo,ipoin) or topel  (ntopo,ielem)
!           mini_patm         <- patma(npatma,ipoin) or patmel(npatma,ielem)
!       
     
implicit none

integer(ink) ntopo, npatma

real   (irk) h, hu, hv, pwp, HH, factw                          !  input
real   (irk) frictx, fricty, slopx, slopy, correx, correy       !  output

real   (irk) cgra, dens, cmanning, cori, tauy0, constK, visco
integer(ink) nfric
real   (irk) tanfi8, hfrict0, c73, tanfi0, Bfact, hrelpwp2, tanfi
real   (irk) rhog, gtanfi, gtanfi8, hinv0, hfrict, hvmod, hvmodinv
real   (irk) cf1,  cf2a, cf3a, cf3b, a_hat, a2, b2, c2, xh, Tau0
real   (irk) cf4b, factpw, cf5, cf8, hinv, hrelpw2 
real   (irk) vx, vy, vmod, cv2byR, cosa, sina, cos2, sin2, cbs
real   (irk) xnum, xden, rinv, factaux 

real   (irk) mini_topo(ntopo), mini_patm(npatma)   !  -- compatibility call to frictlaws

cgra     = constF(1) !gravedad
dens     = constF(2)
cmanning = constF(3)
cori     = constF(4)
nfric    = constF(5) + 0.01
tauy0    = constF(6)
constK   = constF(7)
visco    = constF(8)
tanfi8   = constF(9)
hfrict0  = constF(10)
c73      = 7./3.
tanfi0   = constF(12)
Bfact    = constF(13)

if (icpwpF.eq.1) hrelpw2  = constF(14)*constF(14)     !MP

if (Bfact.ge.1.e-6) then
    tanfi  = tanfi8 + (tanfi0-tanfi8)*exp(-time_gfl/Bfact)
else
    tanfi  = tanfi8
endif

rhog     = dens*cgra
gtanfi   = cgra*tanfi
gtanfi8  = cgra*tanfi8

hinv     = 0.0
if (h.gt.0.) hinv = 1./h
hfrict   = max(h,hfrict0)
hvmod    = sqrt(hu*hu+hv*hv)
vmod     = hvmod*hinv
hvmodinv = 0.0
if(hvmod.ge.1.e-5) hvmodinv = 1./hvmod

frictx = 0.0
fricty = 0.0 

! ..... Obtain friction terms:

if(nfric.eq.1) then                     ! (1) Newtonian fluid turbulent
   cf1     = cgra*cmanning*cmanning
   frictx  = cf1*hu*hvmod/hfrict**c73 
   fricty  = cf1*hv*hvmod/hfrict**c73 
                     
elseif (nfric.eq.2) then                ! (2) Newtonian fluid viscous 
   cf2a    = constK/(8.*dens)
   frictx  = cf2a*hu/hfrict**2.
   fricty  = cf2a*hv/hfrict**2.

elseif (nfric.eq.3) then                ! (3) Coussot rheological law 
   cf3a    = 1.93*(constK**0.9)*(tauy0**0.1)
   cf3b    = (tauy0 + cf3a*(vmod/hfrict)**0.3)/dens
   frictx  = cf3b*hu*hvmodinv
   fricty  = cf3b*hv*hvmodinv        
   
elseif (nfric.eq.4) then                ! (4)  Bingham via Tchebichef approximation
   if (tauy0.lt.1.e-32)then
       write(*,*)' Error entering Yield: must be greater than 1.e-32'
       stop
   endif
   a_hat  = (6.*constK*hvmod)/(tauy0*hfrict**2.)
   a2     = 48.
   b2     = -32.*a_hat-114.
   c2     = 65.
   xh     = (-b2-(b2*b2-4.*a2*c2)**0.5)/(2.*a2)
   Tau0   = tauy0/xh
   cf4b   = Tau0/dens
   frictx = cf4b*hu*hvmodinv
   fricty = cf4b*hv*hvmodinv

elseif (nfric.eq.5) then                ! (5) Frictional fluid 
   cf5    = gtanfi*h
   frictx = cf5*hu*hvmodinv
   fricty = cf5*hv*hvmodinv

elseif(nfric.eq.6) then                 ! (6) Frictional with pwp
   factpw = 1. - pwp*hinv 
   if (factpw.gt.1.0) factpw = 1.0
   if (factpw.lt.0.0) factpw = 0.0
   cf5    = gtanfi8*factpw*h
   frictx = cf5*hu*hvmodinv
   fricty = cf5*hv*hvmodinv

elseif (nfric.eq.7) then                ! (7) Frictional fluid, v2/R term       
   vx     = hu*hinv
   vy     = hv*hinv
   vmod   = (vx*vx+vy*vy)**0.5
   cv2byR = 0.0
   if(vmod.ge.1e-5) then
      cosa   = vx/vmod
      sina   = vy/vmod
      cos2   = cosa*cosa
      sin2   = sina*sina
      cbs    = sina*cosa
      xnum   = mini_topo(7)*cos2+mini_topo(8)*cbs + mini_topo(9)*sin2
      xden   = mini_topo(4)*cos2+mini_topo(5)*cbs + mini_topo(6)*sin2
      rinv   = xnum/xden
      cv2byR =vmod*vmod*rinv
   endif
   factaux = cgra+cv2byR
   if (factaux.lt.0.0) factaux=0.0
   cf5    = factaux*tanfi8*h
   frictx = frictx+cf5*hu*hvmodinv
   fricty = fricty+cf5*hv*hvmodinv

elseif (nfric.eq.8) then        
   vx     = hu*hinv
   vy     = hv*hinv
   vmod   = (vx*vx+vy*vy)**0.5
   cv2byR = 0.0
   if(vmod.ge.1e-5) then
       cosa   = vx/vmod
       sina   = vy/vmod
       cos2   = cosa*cosa
       sin2   = sina*sina
       cbs    = sina*cosa
       xnum   = mini_topo(7)*cos2+mini_topo(8)*cbs +mini_topo(9)*sin2
       xden   = mini_topo(4)*cos2+mini_topo(5)*cbs +mini_topo(6)*sin2
       rinv   = xnum/xden
       cv2byR = vmod*vmod*rinv
   endif
   factaux = cgra+cv2byR
   if (factaux.lt.0.0) factaux=0.0 
   factpw = 1. - pwp*hinv 
   if (factpw.gt.1.0) factpw = 1.0
   if (factpw.lt.0.0) factpw = 0.0
   cf8     = factaux*tanfi8*h*factpw
   frictx  = cf8*hu*hvmodinv
   fricty  = cf8*hv*hvmodinv

elseif(nfric.gt.8)then
   write(*,*) ' Error: nfric not implemented'
   stop
endif

!      ------  We add  manning if necessary
        
if (cmanning.ne.0.and.nfric.ne.1) then
    cf1      = cgra*cmanning*cmanning
    frictx   = frictx + cf1*hu*hvmod/hfrict**c73 
    fricty   = fricty + cf1*hv*hvmod/hfrict**c73 
endif
        
!      ------  We obtain slope terms        
!              In case       ntype = 2 we use HH as obtained in sw,
!              otherwise, if ntype = 3 (DF)  HH=0

slopx   = - cgra*(h-HH)*mini_topo(2)
slopy   = - cgra*(h-HH)*mini_topo(3)
      
!       ------  Finally, we obtain atmospheric gradient and wind terms

correx      = twind(1)*factw - h*mini_patm(2)/1000.
correy      = twind(2)*factw - h*mini_patm(3)/1000.


END subroutine frictlaws
        
! *****************************************************************+++
        
END MODULE GFL_MAIN_2019        