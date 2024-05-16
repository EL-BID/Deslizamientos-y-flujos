!---------------------------------------------------------------------


  MODULE SPH_GFL_SW_Interactions_2019


!---------------------------------------------------------------------

!   Provides interactions between SPH and Geoflow SW codes
!   Interactions can be:
!         (1) waves in reservoirs: geoflow fem water and sph avalanche
!         (2) gates for histogram input: geoflow channel with histogram
!   Routines are called from SPH_SW and Geoflow modules
!   Output routine for WIR called from Main SPH_GFL driver

USE SPH_FEM_variable_types_2019

USE SPH_FEM_Driver_time_vars_2019
USE SPH_SW_Vars_2019
USE GFL_Global_Vars_2019 
USE SPH_Main_Vars_2019

USE SPH_FEM_topo_2019

implicit none
private


integer(ink), public             :: type_of_fem_sph_interaction   
                                                       !  0...no interaction
                                                       !  1   WIR type
                                                       !  2   gate histogram 
                                        
logical,      public             :: WIR_Interact       ! FALSE  no interaction w-s 

!      ------  WIR information computed here used by geoflow & sph-sw                                                       ! TRUE      interaction w-s 

real   (irk), public, allocatable:: hs_w           (:) ! grad h  soil at fem nodes
real   (irk), public, allocatable:: gradhs_w     (:,:) ! grad h  soil at fem nodes
real   (irk), public, allocatable:: gradhs_w_elem(:,:) ! grad h  soil at fem elms
real   (irk), public, allocatable:: gradhw_s     (:,:) ! grad h fluid at fem nodes

!      ------  Gate information read from geoflow

integer(ink), public             :: n_ggates           ! number of interaction gates
integer(ink), public             :: mps_ggates         ! max nr of nodes in all gates
integer(ink), public, allocatable:: nps_ggates(:)      ! n of nodes in each gate
integer(ink), public, allocatable:: list_pts_ggates(:,:,:) 
                                                       ! 1st index for gates
                                                       ! 2nd for SPH and Geoflow node nums
                                                       ! 3rd for max nr of nodes in gates
                                                       ! Read in Geoflow_Main at the end
                                                       ! when type_of_fem_sph_interaction=2
                                                      
integer(ink), public, allocatable::   which_element(:) ! for all sph nods the elem they belong 

public::  get_pics_in_element   !
public::  get_gradhs_w          !
public::  get_gradhw_s          !
public::  plot_wir_msh          !
public::  plot_wir_res          !
public::  Init_gfl_sph_interaction


integer (ink) gidsg_msh, gidsg_res
character(20) problem_gfl_sph_name                     ! name of problem being solved.

CONTAINS

!-------------------------------------------------------------------

       SUBROUTINE Init_gfl_sph_interaction

!-------------------------------------------------------------------

 
!      ------  Initializes output files

implicit none

integer (ink) cont       !   characters counter  
 
write(*,*)' ...........................' 
write(*,*)'    Interaction module ' 
write(*,*)' ...........................' 
write(*,*)'  Please input name for gid joint output (w/o extension) '
read (*,*)   problem_gfl_sph_name
            
!      ------  Open files
         
cont = len_trim(problem_gfl_sph_name)


gidsg_msh   = 35  
gidsg_res   = 36 

open (unit= gidsg_msh,   file = problem_gfl_sph_name(1:cont)//'.post.msh'  )
open (unit= gidsg_res,   file = problem_gfl_sph_name(1:cont)//'.post.res'  )

call plot_wir_msh       !     --  Plot mesh and heading of res  
call plot_wir_res       !         Plot Initial state 
         
END  SUBROUTINE Init_gfl_sph_interaction


!-------------------------------------------------------------------

       SUBROUTINE  get_pics_in_element

!-------------------------------------------------------------------

 
!  we assume same ndimn in sph and ndimnF in geoflow     

implicit none

integer(ink)  nel_full                   !   total pics in elements
                                         !   if zero, we set WIR_Interact .FALSE.
!      ------  vars Task A   ---- 

integer(ink)  idimn, itotal, ielem, inode 
real   (irk)  hminF, anx, any, xle, length, length_new
real   (irk)  xminS (ndimn), xmaxS (ndimn), deltS(ndimn)  
real   (irk)  xminF(ndimnF), xmaxF(ndimnF)

integer(ink), SAVE :: ndivt, mdivt      !  -- save to avoid deallocates and allocates
integer(ink), SAVE, allocatable:: which_cell(:), info_picell(:,:), list_picell(:)
integer(ink), SAVE             :: ndivx(3)
real   (irk), SAVE             :: xmin (3), xmax (3), deltx (3)

!      ------  vars Task B   ------ 

integer(ink) idivt, ipost, idest

!      ------  Vars Task C   ------ 

type picbox                                              
     integer(ink) ipic
     type(picbox), pointer::next
end type picbox

type list_of_picboxes
     integer(ink) npics
     type (picbox), pointer:: current, last
end type list_of_picboxes

type (list_of_picboxes),SAVE, allocatable::sph_in_elems(:)   !  we save it, do not deallocate

integer(ink) lnode, node1, node2, node3, inx,iny, inx0, inx1, iny0, iny1
integer(ink) npics, idest0, ipics, ipoin
integer(ink) idivx(ndimn)

real   (irk) x1, x2, x3, y1, y2, y3, x31, x32, y31, y32, det, detx
real   (irk) xx, yy, x30, y30, sh1, sh2, sh3
real   (irk) xmin_e, xmax_e, ymin_e, ymax_e

!      ---------  General vars    ----- 

character(60) text                      ! general purpose char string
integer(ink), allocatable:: Fst_Time(:) ! to know wheter is the first time we
                                        ! check whether it is allocated
  
 
!      ------  Task A obtain SPH min and max and min spacing     

xminS = 1.e+10; xmaxS = -xminS; deltS = 0.0

Do idimn = 1,ndimn
   Do itotal = 1, npoin
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      if (x(idimn,itotal).lt.xminS(idimn)) xminS(idimn) = x(idimn,itotal)
      if (x(idimn,itotal).gt.xmaxS(idimn)) xmaxS(idimn) = x(idimn,itotal)
      if (deltS(idimn).lt.hsml(itotal))    deltS(idimn) = hsml(itotal)
   enddo  
Enddo

deltS = deltS*2./3.
 
!      ------  obtain for FE mesh min and max and min spacing 

xminF = 1.e+10; xmaxF = -xminF; hminF = 1.e+10

DO ielem = 1,nelem
   DO inode = 1, nnode
      ipoin = intmat(inode,ielem)
      anx   = geome (inode,ielem)
      any   = geome (inode+3,ielem)
      xle   = 1./sqrt(anx*anx + any*any)
      if (xle.lt.hminF) hminF = xle
      DO idimn = 1,ndimnF
         if ( coord(idimn,ipoin).lt.xminF(idimn)) xminF(idimn)=coord(idimn,ipoin)
         if ( coord(idimn,ipoin).gt.xmaxF(idimn)) xmaxF(idimn)=coord(idimn,ipoin)
      enddo
   enddo
enddo

!      ------  Obtain intersection of sph and fe domain

WIR_Interact = .TRUE.
DO idimn = 1, ndimn
   if ( (xminF(idimn).gt.xmaxS(idimn)).OR.(xminS(idimn).gt.xmaxF(idimn)) ) then  !
      WIR_Interact = .FALSE.
   endif
enddo

IF (WIR_Interact) then   ! ***** follows up to end of routine. Thing of exit
   
   ndivx = 1
   
   do idimn = 1,ndimn   !  ----  only valid for ndimn = ndimnF
      xmin(idimn) = max (xminS(idimn),xminF(idimn))
      xmax(idimn) = min (xmaxS(idimn),xmaxF(idimn))
      deltx(idimn) = (deltS(idimn) + hminF )/2.
   enddo
   
   do idimn = 1, ndimn          !  ---  cosmetics, we enlarge region a bit
                                !       only valid for ndimn = ndimnF
      length       = xmax(idimn) - xmin(idimn)
      ndivx(idimn) = (length/deltx(idimn)) + 1
      length_new   = ndivx(idimn)*deltx(idimn)
      xmin (idimn) = xmin(idimn) - (length_new-length)/2 - 0.001*length
      xmax (idimn) = xmax(idimn) + (length_new-length)/2 + 0.001*length
      deltx(idimn) = (xmax (idimn)-xmin(idimn))/ndivx(idimn)
   Enddo
   
   ndivt = ndivx(1)*ndivx(2)*ndivx(3)  !  this is the number of cells we need
   
   !      ------  Now we can allocate or reallocate some control arrays
   
   if (.not.allocated(which_cell) ) then
       mdivt = ndivt
       allocate ( which_cell(npoin), info_picell(3,mdivt),  list_picell(npoin) )
   else if (ndivt.gt.mdivt) then
       deallocate(info_picell)
       mdivt = max (2*mdivt, ndivt)
       allocate  (info_picell(3,mdivt) )
   endif
   
   which_cell = 0; info_picell=0; list_picell=0
      
!      ------  Task B: Fill grid structure with sph nodes (Comes from Grid_Find_New)

   DO itotal = 1, npoin
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      DO idimn = 1, ndimn       !       only valid for ndimn = ndimnF
         idivx(idimn) = (x(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
         if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
      enddo
      idivt =  ndivx(1)*(idivx(2)-1) + idivx(1)   !  in 3d add ndivx(1)*ndivx(2)*(idivx(3)-1) + 
      which_cell(itotal)   = idivt
      info_picell(1,idivt) = info_picell(1,idivt) +1
   Enddo

!      ------  in 2nd row we store, for each cell, the position in array list_picell(ntotal)
!              where we start storing particles in this cell

   ipost = 1
   do idivt = 1,ndivt
      info_picell(2,idivt) = ipost
      ipost = ipost + info_picell(1,idivt)
   enddo

   info_picell(3,:) = -1

   do itotal = 1, npoin
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      idivt  = which_cell  (itotal)
      info_picell(3,idivt) = info_picell(3,idivt) + 1 
      idest  = info_picell(2,idivt) + info_picell(3,idivt)
      list_picell(idest) = itotal
   enddo      
       
!      ------  Task C: obtain to which element belongs each sph node

  if (.NOT.allocated (which_element) ) then
     allocate ( which_element(npoin) )
     allocate ( sph_in_elems(nelem)  )
  else                                  !  we clean all links created previously
     do ielem = 1, nelem
        sph_in_elems(ielem)%current => sph_in_elems(ielem)%last
        do while (associated (sph_in_elems(ielem)%last) )
           sph_in_elems(ielem)%last => sph_in_elems(ielem)%current%next
           deallocate ( sph_in_elems(ielem)%current )
           sph_in_elems(ielem)%current => sph_in_elems(ielem)%last
        enddo
     enddo
  endif

  do ielem = 1,nelem                    !   Inititalizes
     nullify ( sph_in_elems(ielem)%last )
     sph_in_elems(ielem)%npics = 0
  enddo
  which_element = 0
  nel_full = 0                          !  this is total nr of sph nodes in Felems
  
  DO ielem = 1, nelem
     
     xmin_e = xmaxF(1); ymin_e = xmaxF(2)
     xmax_e = xminF(1); ymax_e = xminF(2)
     do inode = 1, nnode
        lnode = intmat(inode,ielem)
        if (coord(1,lnode).lt.xmin_e) xmin_e = coord(1,lnode)
        if (coord(2,lnode).lt.ymin_e) ymin_e = coord(2,lnode)
        if (coord(1,lnode).gt.xmax_e) xmax_e = coord(1,lnode)
        if (coord(2,lnode).gt.ymax_e) ymax_e = coord(2,lnode)
     enddo
                                ! intersection element sph region void if
     if ( xmin_e.gt.xmax(1).OR.xmin(1).gt.xmax_e) CYCLE
     if ( ymin_e.gt.xmax(2).OR.xmin(2).gt.ymax_e) CYCLE

     node1 = intmat (1,ielem)
     node2 = intmat (2,ielem)
     node3 = intmat (3,ielem)
     x1    = coord (1,node1); y1 = coord (2,node1)
     x2    = coord (1,node2); y2 = coord (2,node2)
     x3    = coord (1,node3); y3 = coord (2,node3)
     x31   = x1 - x3; y31 = y1 - y3
     x32   = x2 - x3; y32 = y2 - y3
     detx  = ( x31*y32 - y31*x32 )
     inx0  = ( xmin_e - xmin(1) )/deltx(1) 
     iny0  = ( ymin_e - xmin(2) )/deltx(2)  
     if (inx0.lt.1) inx0 = 1
     if (iny0.lt.1) iny0 = 1
     inx1  = (xmax_e - xmin(1))/deltx(1) 
     iny1  = (ymax_e - xmin(2))/deltx(2)  
     if (inx1.gt.ndivx(1)-1) inx1 = ndivx(1)-1  
     if (iny1.gt.ndivx(2)-1) iny1 = ndivx(2)-1  
     
     DO inx = inx0, inx1
     DO iny = iny0, iny1
        idivt  = inx + 1 + ndivx(1)*iny   
        npics  = info_picell (1,idivt)
        idest0 = info_picell (2,idivt)
        DO ipics = 1, npics
           idest = idest0 + ipics - 1
           ipoin = list_picell (idest)
           xx    = x(1,ipoin)
           yy    = x(2,ipoin)
           x30   = xx - x3 ; y30  = yy - y3
           sh1   = (x30*y32-y30*x32)/detx
           sh2   = (y30*x31-x30*y31)/detx
           sh3   = 1.0 - sh1 - sh2
           if ( (sh1.ge.0.and.sh1.le.1).AND.(sh2.ge.0.and.sh2.le.1).AND.(sh3.ge.0.and.sh3.le.1) ) then
              which_element (ipoin) = ielem
              nel_full = nel_full +1
              allocate (sph_in_elems(ielem)%current)
              sph_in_elems(ielem)%npics = sph_in_elems(ielem)%npics + 1
              sph_in_elems(ielem)%current%ipic = ipoin
              sph_in_elems(ielem)%current%next => sph_in_elems(ielem)%last 
              sph_in_elems(ielem)%last => sph_in_elems(ielem)%last
           endif
        ENDDO
     ENDDO
     ENDDO
     
  ENDDO
   
   
ENDIF   !  endif for existing interaction WIR_Interact = .TRUE.

if ( nel_full.eq.0) WIR_Interact = .FALSE.

END SUBROUTINE  get_pics_in_element


!-------------------------------------------------------------------

       SUBROUTINE  get_gradhs_w

!-------------------------------------------------------------------


implicit none

integer(ink) ielem, lpoiF, ipoiF
integer(ink) id, idimn, iiter, inode, jnode,ipoin, jpoin
integer(ink) niter_t

real   (irk) xn, yn, x0, y0, rdist, anx, any
real   (irk)  mmat(3,3)
real   (irk) coe, cm, derixe, deriye 

real   (irk), allocatable:: aux  (:,:)
real   (irk), allocatable:: rhs  (:,:)       
real   (irk), allocatable:: rhelp(:,:)

! real   (irk), allocatable:: areatmp(:), areagrad(:,:)  !  used to smooth the gradient

integer(ink) active_ipoin(npoinF)   ! 1 if ipoin contributes to grad hs
integer(ink) active_ielem(nelem)    ! 1 if ielem contributes to grad hs

character(60) text                    ! general purpose char string

real   (irk) temp_hs_w(npoinF)


!      ------ Initialize and allocate

niter_t = 5                           !      perform 5 iters to solve the system

if (.NOT.allocated (hs_w) ) then             ! allocation instead of direct dimens 
   allocate ( hs_w(npoinF)                )  ! saves memory when not used
   allocate ( gradhs_w(ndimnF,npoinF)     )
   allocate ( gradhs_w_elem(ndimnF,nelem) )
endif 

allocate ( aux (ndimnF, npoinF ) )
allocate ( rhs(ndimnF,npoinF)    )
allocate ( rhelp (ndimnF,npoinF ))

!allocate ( areatmp(npoinF), areagrad(ndimnF, npoinF) )

hs_w = 0.0 ; gradhs_w = 0.0 ; gradhs_w_elem = 0.0
aux  = 0.0 ; rhs     = 0.0 ; rhelp         = 0.0
active_ipoin = 0; active_ielem = 0

!     ------- obtain active nbodes and elements

do ipoin = 1, npoin             !  activates elems with sph nodes inside
   ielem = which_element(ipoin)
   if (ielem.ne.0) then
       active_ielem(ielem) = 1
   endif
enddo

do ielem = 1, nelem             !  activate nodes of active elements
   if ( active_ielem(ielem).eq.1 ) then
      do inode = 1,nnode
         ipoin = intmat (inode, ielem)
         active_ipoin(ipoin) = 1
      enddo
   endif
enddo

do ielem = 1, nelem             !  activate elements with active nodes
   if ( active_ielem(ielem).eq.0 ) then 
      do inode = 1,nnode
         ipoin = intmat (inode, ielem)
         if ( active_ipoin(ipoin).eq.1 ) then
            active_ielem(ielem) = 1
            EXIT
         endif
      enddo
   endif
enddo

!      ------ A1   Obtain hs_w at fe mesh nodes  

do ipoin = 1, npoin     !   we obtain h = sum (h/r) / sum(1/r)  
   ielem = which_element(ipoin)
   if (ielem.eq.0) CYCLE
   do inode = 1, nnode
      lpoiF = intmat(inode,ielem)
      xn    = coord (1,lpoiF) ; yn = coord (2,lpoiF)
      x0    = x     (1,ipoin) ; y0 = x     (2,ipoin)
      rdist = sqrt( (xn-x0)*(xn-x0) + (yn-y0)*(yn-y0) )
      aux(1,lpoiF) = aux(1,lpoiF) + rho (ipoin)/rdist
      aux(2,lpoiF) = aux(2,lpoiF) +          1./rdist
   enddo
enddo

do ipoiF = 1, npoinF
   if (aux(2,ipoiF).lt.1.e-6) CYCLE
   hs_w(ipoiF) = aux(1,ipoiF)/aux(2,ipoiF)
enddo

temp_hs_w = hs_w

!      ------ B1   Obtain gradients at elements  

do ielem = 1,nelem
   derixe = 0.0
   deriye = 0.0
   do inode = 1,nnode
      ipoin = intmat(inode,ielem)
      derixe = derixe + hs_w(ipoin)*geome(inode,  ielem)
      deriye = deriye + hs_w(ipoin)*geome(inode+3,ielem)
   enddo
   gradhs_w_elem(1,ielem) = derixe
   gradhs_w_elem(2,ielem) = deriye
enddo

!      ------ B2  Get Mass matrix  

mmat = 1.0
do inode = 1,nnode
   mmat(inode,inode)=2.0
enddo


!      ------  B3  Build RHS Int(Ni.Grad) Store it in rhs

rhs = 0.0

do ielem = 1,nelem
   coe = geome(7,ielem)/6.    !  this is A/3 = geome/2.3
   do inode = 1,nnode
      ipoin = intmat(inode,ielem)
      rhs (1,ipoin) = rhs(1,ipoin) + gradhs_w_elem(1,ielem)*coe
      rhs (2,ipoin) = rhs(2,ipoin) + gradhs_w_elem(2,ielem)*coe
   enddo
enddo

!      ------ B4   Get  v(1)=v(0)+dt.inv(Ml)*(rhs-M.v(0)) v(0)=0

do ipoin = 1,npoinF
   cm = mmatl(ipoin)
   gradhs_w (1,ipoin) = gradhs_w (1,ipoin)/cm  
   gradhs_w (2,ipoin) = gradhs_w (2,ipoin)/cm  
enddo

!      ------ B5  Iterate: v(n+1) = v(n) + inv(Ml)*(rhs-M*V(n))
!                          M*v(n) is stored as rhelp

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
               rhelp(idimn,ipoin) =  rhelp(idimn,ipoin) & 
                                    + cm*coe*gradhs_w(idimn,jpoin)
            enddo
         enddo
      enddo
   enddo

   do ipoin = 1,npoinF
      cm = mmatl(ipoin)
      do idimn = 1,ndimnF
         gradhs_w(idimn,ipoin) = gradhs_w(idimn,ipoin) &
              + (rhs(idimn,ipoin)-rhelp(idimn,ipoin))/cm
      enddo
   enddo
   
enddo

!areatmp = 0.0 ; areagrad = 0.0        !       smooth gradient    CHANGE 28 dic
!do ielem = 1, nelem
!   if (active_ielem(ielem).eq.1) then
!      coe = geome(7,ielem)
!      do inode = 1,nnode
!         ipoin = intmat (inode,ielem)
!         areatmp   (ipoin) = areatmp   (ipoin) + coe
!         areagrad(:,ipoin) = areagrad(:,ipoin) + coe*gradhs_w(:,ipoin)
!      enddo
!   endif
!enddo
!
!do ipoin = 1,npoinF
!   if (active_ipoin(ipoin).eq.1) then
!      gradhs_w(:,ipoin) = areagrad(:,ipoin)/areatmp(ipoin)
!   endif
!enddo

 
deallocate (aux, rhs, rhelp )
! deallocate (areatmp,areagrad)
      
END SUBROUTINE get_gradhs_w 

!-------------------------------------------------------------------

       SUBROUTINE get_gradhw_s 

!-------------------------------------------------------------------

 
!      Obtains grad hw at sph nodes
!      knowing the elem where ipoin is

implicit none  


integer(ink) ipoin, ipoinF, ielem, inode
real   (irk) anx, any

if (.NOT.allocated (gradhw_s) ) then
   allocate ( gradhw_s(ndimnF,npoinF) )
endif
gradhw_s = 0.0

do ipoin = 1, npoin                     !    loop in all sph nodes
   ielem = which_element (ipoin)        !    this is the element 
   if (ielem.eq.0) CYCLE
   do inode = 1,nnode                   !    hw is height in SW: unkno(1,.)  
      ipoinF = intmat (inode  , ielem)  !    
      anx    = geome  (inode  , ielem)
      any    = geome  (inode+3, ielem)
      gradhw_s (1,ipoinF) = gradhw_s (1,ipoinF) + anx*unkno(1,ipoinF)
      gradhw_s (2,ipoinF) = gradhw_s (1,ipoinF) + anx*unkno(1,ipoinF)
   enddo
enddo

END SUBROUTINE  get_gradhw_s

 


!-------------------------------------------------------------------

       SUBROUTINE  plot_wir_msh

!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  for geoflow and sph codes --------

implicit none

integer  (ink) ipoin, i,is_end, idimn, inodg, inode
integer  (ink) ipoig, ipoigx, ipoigy, ieleg, neleg, ielem
integer  (ink) keleg, plteleg 

integer  (ink), allocatable:: intmag4(:,:)

real     (irk) xx, yy, zp

!      ------  generates background grid, plotting purposes only

neleg =  (npoigx-1)*(npoigy-1)
allocate ( intmag4(4,neleg) )
intmag4 = 0
keleg   = 0   ! chg
do ipoigy = 1,npoigy-1
   do ipoigx  = 1, npoigx-1
      ieleg   = (ipoigy-1)*(npoigx-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoigx + ipoigx
      plteleg = poig_in_model(ipoig)*poig_in_model(ipoig+1)*poig_in_model(ipoig+npoigx+1)*poig_in_model(ipoig+npoigx)  ! **** chg
      if (plteleg.eq.1) then! **** chg
          keleg = keleg+1
          intmag4(1,keleg) = ipoig   
          intmag4(2,keleg) = ipoig  + 1
          intmag4(3,keleg) = ipoig  + npoigx + 1
          intmag4(4,keleg) = ipoig  + npoigx 
      endif
   enddo
enddo

      
write (gidsg_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gidsg_msh,*) ' coordinates'      !    - COORDS  --------------------

do ipoig = 1, npoig                     !    --  Q1 This is first set of g nodes 1:npoig
   xx    = coorg (1,ipoig)
   yy    = coorg (2,ipoig)
   zp    = topol(1,ipoig)
   write(gidsg_msh,*) ipoig, xx, yy, zp
enddo

do ipoin = 1, npoinF                    !   --  T1 This is 2nd set of F nodes npoig+1:2npoig+npoiF
   i     = npoig + ipoin
   xx    = coord (1,ipoin)
   yy    = coord (2,ipoin)
   zp    = topol (1,ipoin)
   write(gidsg_msh,*) i, xx, yy, zp
enddo

Call Check_Out_Domain_SWInt

do ipoin  = 1, npoin                    !    --  3rd set: coordinates of SPH nodes
   i      = ipoin + npoig + npoinF
   xx     = x(1,ipoin)
   if (ndimn.eq.1) then
      yy  = (yming + ymaxg)/2.
   else
      yy  = x(2,ipoin)
   endif 
   call Get_Z_Topo ( ipoin, xx, yy, Zp)

   Z00(ipoin) = Zp                      !    ---  Here we store the initial Z of particles. Used for plotting
   write(gidsg_msh,*) i, xx, yy, zp 
enddo
write (gidsg_msh,*) ' End coordinates'      
                                        !    ---  ELEMENTS  ------ 
                                                
write (gidsg_msh,*) 'MESH SPH_Q1   dimension   3    ElemType Quadrilateral  Nnode 4 '  !-------    
write (gidsg_msh,*) ' Elements'
do ieleg = 1, keleg               
   write (gidsg_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gidsg_msh,*) ' End elements' 

      
write (gidsg_msh,*) 'MESH GFL_T1   dimension   3    ElemType Triangle  Nnode 3 '  !------- 
write (gidsg_msh,*) ' Elements'
do ielem = 1, nelem
   write (gidsg_msh,*) ielem, (intmat(inode,ielem)+npoig, inode = 1,3), '    3 '
enddo
write (gidsg_msh,*) ' End elements'


write (gidsg_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  SOIL ---- 
write (gidsg_msh,*) ' Elements'
do i = 1, npoin                                          
   is_end = i
   if (iabs(itype(i)).eq.2) EXIT
   write(gidsg_msh,*) i,  npoinF + npoig+i, iabs(itype(i))
enddo
write (gidsg_msh,*) ' End elements' 

     
write (gidsg_res,*) 'GiD Post Results File 1.0'
      
close(gidsg_msh)
deallocate (intmag4)


END SUBROUTINE  plot_wir_msh

 


!-------------------------------------------------------------------

       SUBROUTINE  plot_wir_res

!-------------------------------------------------------------------

 
!      ------  OUT for waves in reservoirs --------

implicit none

integer  (ink)  i, ipoin, ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer  (ink)  inx0, inx1, inx, iny0, iny1, iny
integer  (ink)  idest
real     (irk)  xx, yy, zz, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp 
real     (irk)  xdisp, ydisp, zdisp, vxx, vyy
         
integer  (ink), allocatable, SAVE:: ictrack(:) 

real     (irk), allocatable::disp(:)
real     (irk), allocatable::auxv(:,:)

real     (irk)  distx2, disty2, dist, distg 
real     (irk)  dArea 

allocate ( disp(ndimn) ) 
allocate ( auxv(5,npoig))

Call Check_Out_Domain_SWInt


!      ------   First of all we get values at topo grid. Add for each cell value*(1/dist) and divide by sum (1/dist)

auxv    = 0.0

if(.not.allocated(ictrack)) then                        !      ------  Added
   allocate (ictrack(npoig))
   ictrack = 0
endif
                                                        !      ------  Added
distg   = (deltxg*deltxg + deltyg*deltyg)**0.5

DO ipoin = 1,npoin                                      !      ------  ch
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   idest = 0                                            !  5*(iabs(itype(ipoin))-1)    0 for soil and 5 for water
   xx = x(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
   inx0      = max(1,ipoigx-1)
   iny0      = max(1,ipoigy-1)
   inx1      = min(npoigx,ipoigx+2)
   iny1      = min(npoigy,ipoigy+2)
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      distx2 = ( x(1,ipoin) - coorg(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x(2,ipoin) - coorg(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      dist   = dist/distg
      auxv(1,ipoig) = auxv(1,ipoig) + rho(ipoin)/dist                   ! h        soil  and water
      auxv(2,ipoig) = auxv(2,ipoig) + hsml(ipoin)/dist                  ! temporal hsml s and w
      auxv(4,ipoig) = auxv(4,ipoig) + (u(ipoin)/rho(ipoin))/dist        ! pwp
      auxv(5,ipoig) = auxv(5,ipoig) + 1.00/dist                         ! dist for soil and water
   enddo
   enddo
ENDDO



DO ipoig = 1, npoig
   if (auxv(5,ipoig).ge.1.e-6) then
       auxv(1:4,ipoig) = auxv(1:4,ipoig)/auxv(5,ipoig) 
   endif
   if (auxv(1,ipoig).ge.1.e-6) ictrack(ipoig)=1                 !    ADDED and changed to 0.001 from 0.1
enddo

!      ------   And now plot:
 
!      ------  Erosion  --------- 


If (constF(4).gt.1.e-6) then
   write(gidsg_res,*) 'Result "erosion" " erosion " ',time_sph,' Vector OnNodes "where" '
   write(gidsg_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots erosioned soil
      if ( topol(14,ipoig).le.1.e-6) CYCLE
      write(gidsg_res,*) ipoig, ' 0.   0.   ', -topol(14,ipoig) 
   enddo
   write(gidsg_res,*) ' End Values ' 
endif   

!      ------  h soil  ------------- 

write(gidsg_res,*) 'Result "height soil" "Height soil" ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values '
do ipoig = 1, npoig                                    !      ------  plots h SOIL at quad nodes  ADDED
   if ( ictrack(ipoig).ne.1) CYCLE
   write(gidsg_res,*) ipoig, ' 0.   0.   ', auxv(1,ipoig)             ! max (0.1,auxv(1,ipoig))
enddo
write(gidsg_res,*) ' End Values ' 

!      ------  Topo Z  --------------- 
write(gidsg_res,*) 'Result " Z topo " "   Z   " ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values '
do ipoig = 1, npoig 
   write(gidsg_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
enddo
write(gidsg_res,*) ' End Values ' 

!      ------  Displacements at particles   --- 

write(gidsg_res,*) 'Result "dis" "Disp" ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values  '  
DO i = 1,npoin                                        !      ------  plots disp at particles 
   IF ( if_Out_Domain(i).eq.1) CYCLE
      disp(:) = x(:,i) - x00(:,i)
      xdisp   = disp(1)
      if (ndimn.eq.1) then
          ydisp = 0.0
      else
          ydisp = disp(2)
      endif
      xx     =  x(1,i)
      yy     = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,i) 
      
      call Get_Z_Topo (ipoin, xx,yy,Zp)

      zdisp  = Zp - Z00(i)
      write(gidsg_res,*) i + npoinF + npoig, xdisp, ydisp, zdisp
ENDDO
write(gidsg_res,*) ' End Values ' 


!      ------  And PWPs  -----------------------------

if (icpwpF.eq.1) then
   write(gidsg_res,*) 'Result "Pwp" "pwpress" ',time_sph,' Scalar OnNodes "where" '
   write(gidsg_res,*) ' Values '
   do ipoig = 1, npoig
      if ( auxv(1,ipoig).le.1.e-6) CYCLE
      write(gidsg_res,*) ipoig, auxv(4,ipoig)
   enddo
   do i = 1,npoin 
      if ( if_Out_Domain(i).eq.1) CYCLE
      write(gidsg_res,*) i , u(i)
   enddo
   write(gidsg_res,*) ' End Values '
endif  

!      ------  hw in the finite element T1 mesh   --- 

write(gidsg_res,*) 'Result "hw" "h_w" ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values  '  
DO i = 1,npoinF                                        !      ------  plots hw at fem nodes 
      xx     =  coord(1,i)
      yy     =  coord(2,i)
      zz     =  unkno(1,i)
      write(gidsg_res,*) i + npoig , '    0.   .0  ' , zz
ENDDO
write(gidsg_res,*) ' End Values ' 


!      ------  eta in the finite element T1 mesh: only for ntype=2   --- 

if (ntype.eq.2) then 
   write(gidsg_res,*) 'Result "eta" "eta w" ',time_sph,' Vector OnNodes "where" '
   write(gidsg_res,*) ' Values  '  
   if (WIR_Interact) then
      do i = 1,npoinF          !      ------  plots hw at fem nodes  when interaction exists
         xx     =  coord(1,i)
         yy     =  coord(2,i)
         Zp     =  topol(1,i) + unkno(1,i) + hs_w(i) - SWL
         write(gidsg_res,*) i + npoig , '     0.   .0  ' , Zp
      enddo
   else
      do i = 1,npoinF          !      ------  same  when interaction does not exists
         xx     =  coord(1,i)
         yy     =  coord(2,i)
         Zp     =  topol(1,i) + unkno(1,i) - SWL
         write(gidsg_res,*) i + npoig , '     0.   .0  ' , Zp
      enddo
   endif
   write(gidsg_res,*) ' End Values ' 
endif

!      ------  hs_w + hw in the finite element T1 mesh: only for ntype=2   --- 

write(gidsg_res,*) 'Result "hsplushw" "hshw" ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values  '  
 
if (WIR_Interact) then
    do i = 1,npoinF          !      ------  plots hw at fem nodes  when interaction exists
       xx     =  coord(1,i)
       yy     =  coord(2,i)
       Zp     =  unkno(1,i) + hs_w(i)  
       write(gidsg_res,*) i + npoig , '     0.     .0  ' , Zp
    enddo
else
    do i = 1,npoinF          !      ------  same  when interaction does not exists
       xx     =  coord(1,i)
       yy     =  coord(2,i)
       Zp     =  unkno(1,i)    
       write(gidsg_res,*) i + npoig , '     0.      .0     ' , Zp
    enddo
endif
write(gidsg_res,*) ' End Values ' 


!      ------  hs_w   in the finite element T1 mesh: only for ntype=2   --- 

write(gidsg_res,*) 'Result "hs w" "hs_w" ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values  '  
 
if (WIR_Interact) then
    do i = 1,npoinF          !      ------  plots hw at fem nodes  when interaction exists
       xx     =  coord(1,i)
       yy     =  coord(2,i)
       Zp     =    hs_w(i)  
       write(gidsg_res,*) i + npoig , '     0.     .0  ' , Zp
    enddo
else
    do i = 1,npoinF          !      ------  same  when interaction does not exists
       xx     =  coord(1,i)
       yy     =  coord(2,i)
       Zp     =  unkno(1,i)    
       write(gidsg_res,*) i + npoig , '     0.      .0     ' , Zp
    enddo
endif
write(gidsg_res,*) ' End Values ' 

!      ------  Gradhs_w in the finite element T1 mesh: only for ntype=2   --- 

write(gidsg_res,*) 'Result "gradhsw" "gradhs_hw" ',time_sph,' Vector OnNodes "where" '
write(gidsg_res,*) ' Values  '  
 
if (WIR_Interact) then
    do i = 1,npoinF          !      ------  plots hw at fem nodes  when interaction exists
       xx     =  gradhs_w(1,i)
       yy     =  gradhs_w(2,i)
       Zp     =  0  
       write(gidsg_res,*) i + npoig , xx, yy , Zp
    enddo
endif
write(gidsg_res,*) ' End Values ' 
deallocate ( disp )
deallocate ( auxv )



END SUBROUTINE  plot_wir_res


! ----------------------------------------------------------------------   

      subroutine  Check_Out_Domain_SWInt  

! ----------------------------------------------------------------------   

!   Subroutine to check wheter a given point is outside the computational domain  
 
implicit none

integer(ink)  ipoin, idimn
real   (irk)  dxx 

!      ------  Initialization 

!   If_Out_Domain = 0    CHANGED 19 JULY 2007

!      ------  

DO ipoin = 1,npoin
   Do idimn = 1, ndimn
      dxx   = (x(idimn,ipoin)-Xmin_Domain(idimn))*(x(idimn,ipoin)-Xmax_Domain(idimn))
      if ( dxx.gt.0.0 ) If_Out_Domain(ipoin) = 1
   Enddo
Enddo

END  SUBROUTINE Check_Out_Domain_SWInt

END MODULE SPH_GFL_SW_Interactions_2019