MODULE My_SPH_Tools_2019
!
USE SPH_FEM_variable_types_2019
USE SPH_FEM_Driver_time_vars_2019
USE SPH_Main_Vars_2019 
USE SPH_SW_Vars_2019

implicit none

PRIVATE    ! MP reactivated it
! PRIVATE

Public :: grid_find_NEW_tools
Public :: grid_find_NEW_tools2
Public :: grid_find_NEW_tools3
Public :: grid_find_NEW_Tools_TEST  !   ***** MP AUG 2022 to test increm Mb
Public :: Pint_Update_SW_Tools
Public :: sum_density_Tools
Public :: con_density_Tools
Public :: Solid_Slab
Public :: kernel_Tools
Public :: get_wir_density_TOOLS
Public :: get_gradient_TOOLS 
Public :: get_div_TOOLS
Public :: get_grad_TOOLS
Public :: deallocate_list_TOOLS

contains


!---------------------------------------------------------------

     subroutine get_div_TOOLS (Field, ic_norm_div, div_V )

!---------------------------------------------------------------

implicit none

real    (irk) Field(ndimn, npoin), div_V (npoin)
integer (ink) ic_norm_div
real    (irk), allocatable:: A_norm         (:,:)
real    (irk), allocatable:: div_V_tmp(:,:,:)

integer(ink) i, j, idimn, ipoin
real   (irk) hf1i, hf1j, hf2i, hf2j, gradW

if (.not.allocated(div_V_tmp)) then                                                 
    allocate ( div_V_tmp(ndimn,ndimn,npoin) );     div_V_tmp = 0.0
    allocate ( A_norm                 (5,npoin) ); A_norm    = 0.0
    div_V     = 0.0
endif

!!!!!!!!!!!!! For testing gradV* ==> V*(x) = x ==> gradV*=1

Do ipoin= 1, npoin 
   Field = x
Enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

current => last
DO WHILE (associated(current))
   i = current%pair_i 
   j = current%pair_j 
   if (current%Pint_type.eq.0) then     !  i..node j..node
       gradW = current%dwdx(1)
       hf1i  =   mass(j)*gradW/rho(j)
       hf1j  = - mass(i)*gradW/rho(i)
       DO idimn = 1, ndimn
	      div_V_tmp(idimn,1,i) =  div_V_tmp(idimn,1,i) + (Field(idimn,j)- Field(idimn,i)) * hf1i
          div_V_tmp(idimn,1,j) =  div_V_tmp(idimn,1,j) + (Field(idimn,i)- Field(idimn,j)) * hf1j
	      If (ndimn.eq.2) then
	         hf2i =   mass(j)*current%dwdx(2)/rho(j)
             hf2j = - mass(i)*current%dwdx(2)/rho(i)
             div_V_tmp(idimn,2,i) =  div_V_tmp(idimn,2,i) & 
                                   + (Field(idimn,j)- Field(idimn,i)) * hf2i
              div_V_tmp(idimn,2,j) =  div_V_tmp(idimn,2,j) & 
                                   + (Field(idimn,i)- Field(idimn,j)) * hf2j                
	      Endif		

		  If (ic_norm_div.eq.1) then  ! FS MP 25th April 17
    		  If (idimn.eq.1) then
	    	      A_norm(1,i)=A_norm(1,i)+(x(1,j)-x(1,i))*hf1i		    ! A11 for node i
	    	      A_norm(1,j)=A_norm(1,j)+(x(1,i)-x(1,j))*hf1j            ! A11 for node j
		          If (ndimn.eq.2) then
				      A_norm(2,i)=A_norm(2,i)+(x(2,j)-x(2,i))*hf1i		    ! A12 for node i
			          A_norm(3,i)=A_norm(3,i)+(x(1,j)-x(1,i))*hf2i		    ! A21 for node i
				      A_norm(4,i)=A_norm(4,i)+(x(2,j)-x(2,i))*hf2i		    ! A22 for node i
				      A_norm(2,j)=A_norm(2,j)+(x(2,i)-x(2,j))*hf1j		    ! A12 for node j
			          A_norm(3,j)=A_norm(3,j)+(x(1,i)-x(1,j))*hf2j		    ! A21 for node j
				      A_norm(4,j)=A_norm(4,j)+(x(2,i)-x(2,j))*hf2j		    ! A22 for node j
		          Endif
		      endif
          Endif
       Enddo
   Endif
   current=>current%next
ENDDO

If (ic_norm_div.eq.1) then
	If (ndimn.eq.1) then
		A_norm (2,:) = 0
		A_norm (3,:) = 0
		A_norm (4,:) = 1
	Endif
	do ipoin = 1,npoin 
	    A_norm(5,ipoin)=(A_norm(1,ipoin)*A_norm(4,ipoin))-(A_norm(2,ipoin)*A_norm(3,ipoin))
		If (abs(A_norm(5,ipoin)).lt.1.e-9) then
			A_norm(5,ipoin)=0
	    else
			A_norm(5,ipoin)=1/A_norm(5,ipoin)
		endif
	enddo
Endif

If (ndimn.eq.1) then
	if (ic_norm_div.eq.1) then
	   Do i=1,npoin 
			div_V(i) = A_norm(5,i)*div_V_tmp(1,1,i)
	   enddo
	elseif (ic_norm_div.eq.0) then
	   Do i=1,npoin 
			div_V(i) = div_V_tmp(1,1,i)
	   Enddo
	endif

Elseif (ndimn.eq.2) then
	if (ic_norm_div.eq.1) then
	   Do i=1,npoin 
			div_V(i) = A_norm(5,i)*((A_norm(4,i)*div_V_tmp(1,1,i)) 	& 
                                  - (A_norm(2,i)*div_V_tmp(1,2,i))  &  
						          - (A_norm(3,i)*div_V_tmp(2,1,i))  &
                                  + (A_norm(1,i)*div_V_tmp(2,2,i)))
	   enddo
	elseif (ic_norm_div.eq.0) then
	   Do i=1,npoin 
			div_V(i) = (div_V_tmp(1,1,i) + div_V_tmp(2,2,i))
	   Enddo
	endif
Endif
 

 end subroutine get_div_TOOLS


!---------------------------------------------------------------

     subroutine get_grad_TOOLS (Fi, ic_norm_grad,  grad_Fi )

!---------------------------------------------------------------
 

!     Subroutine to calculate the gradient grad_Field(nd,np) of a scalar field  Field(np)
!
!     Fi    : scalar field of dimension nd defined at np points (sph)
!     grad_Fi (ndimn,npoin) gradient of Field (npoin)              [out]   
!     ic_norm_grad: if set to 1 normalize            
implicit none

real   (irk) Fi (npoin), grad_Fi (ndimn,npoin)
integer(ink) ic_norm_grad

real    (irk), allocatable:: A_norm      (:,:)
real    (irk), allocatable:: grad_Fi_tmp (:,:)

integer(ink) i, j, idimn, ipoin 
real   (irk) hf1i, hf1j, hf2i, hf2j, gradW 

if (.not.allocated ( grad_Fi_tmp)) then
   allocate (A_norm        (5,npoin) )
   allocate (grad_Fi_tmp(ndimn,npoin)) 
endif

grad_Fi          = 0.0
A_norm           = 0.0
grad_Fi_tmp      = 0.0 

 !Fi = 0.0
! Do ipoin=1,npoin                          !!  ** MP test 
 !   do idimn=1,ndimn
 !     Fi  (ipoin)  = Fi  (ipoin) + x(idimn,ipoin)  
 !   enddo                                                           
 !Enddo                                                                

!      ------   !st, we obtain gradient of Fi at time n+1

current => last
DO WHILE (associated(current))
     i = current%pair_i 
     j = current%pair_j 
     if (current%Pint_type.eq.0) then     !  i..node j..node
         DO idimn = 1, ndimn
            gradW = current%dwdx(idimn)
            hf1i  =   mass(j)*gradW/rho(j)
            hf1j  = - mass(i)*gradW/rho(i)
            grad_Fi_tmp(idimn,i) =  grad_Fi_tmp(idimn,i) + ( Fi (j)- Fi(i)) * hf1i
            grad_Fi_tmp(idimn,j) =  grad_Fi_tmp(idimn,j) + ( Fi (i)- Fi(j)) * hf1j                 
            If (ic_norm_grad.eq.1) then
                If (idimn.eq.1) then
                    A_norm(1,i)=A_norm(1,i)+(x(1,j)-x(1,i))*hf1i        ! A11 for node i
                    A_norm(1,j)=A_norm(1,j)+(x(1,i)-x(1,j))*hf1j        ! A11 for node j
                    If (ndimn.eq.2) then
                        hf2i =   mass(j)*current%dwdx(2)/rho(j)
                        hf2j = - mass(i)*current%dwdx(2)/rho(i)
                        A_norm(2,i) = A_norm(2,i)+(x(2,j)-x(2,i))*hf1i   ! A12 for node i
                        A_norm(3,i) = A_norm(3,i)+(x(1,j)-x(1,i))*hf2i   ! A21 for node i
                        A_norm(4,i) = A_norm(4,i)+(x(2,j)-x(2,i))*hf2i   ! A22 for node i
                        A_norm(2,j) = A_norm(2,j)+(x(2,i)-x(2,j))*hf1j   ! A12 for node j
                        A_norm(3,j) = A_norm(3,j)+(x(1,i)-x(1,j))*hf2j   ! A21 for node j
                        A_norm(4,j) = A_norm(4,j)+(x(2,i)-x(2,j))*hf2j   ! A22 for node j
                    Endif
                endif
             Endif
         Enddo
     Endif
     current=>current%next
ENDDO

If (ic_norm_grad.eq.1) then
    If (ndimn.eq.1) then
        A_norm (2,:) = 0
        A_norm (3,:) = 0
        A_norm (4,:) = 1
    Endif
    do ipoin=1,npoin
        A_norm(5,ipoin)=(A_norm(1,ipoin)*A_norm(4,ipoin))-(A_norm(2,ipoin)*A_norm(3,ipoin))
        If (abs(A_norm(5,ipoin)).lt.1.e-9) then
            A_norm(5,ipoin)=0
        else
            A_norm(5,ipoin)=1/A_norm(5,ipoin)
        endif
    enddo
Endif

If (ndimn.eq.1) then
    if (ic_norm_grad.eq.1) then
        Do i=1,npoin
            grad_Fi (1,i) = A_norm(5,i)*grad_Fi_tmp(1,i)
        enddo
    elseif (ic_norm_grad.eq.0) then
         Do i=1,npoin
            grad_Fi (1,i) = grad_Fi_tmp(1,i)
         Enddo
    endif
Elseif (ndimn.eq.2) then
     if (ic_norm_grad.eq.1) then
        Do i=1,npoin
            grad_Fi (1,i) = A_norm(5,i)*( A_norm(4,i)*grad_Fi_tmp(1,i) - A_norm(2,i)*grad_Fi_tmp(2,i))
            grad_Fi (2,i) = A_norm(5,i)*(-A_norm(3,i)*grad_Fi_tmp(1,i) + A_norm(1,i)*grad_Fi_tmp(2,i))
         enddo
      elseif (ic_norm_grad.eq.0) then
         Do i=1,npoin
             grad_Fi (1,i) = grad_Fi_tmp(1,i)
             grad_Fi (2,i) = grad_Fi_tmp(2,i)
         Enddo
      endif
Endif 

END subroutine get_grad_TOOLS



!---------------------------------------------------------------

     subroutine get_gradient_TOOLS (Field, nd, np, grad_Field )

!---------------------------------------------------------------
 

!     Subroutine to calculate the gradient grad_Field(nd,np) of a scalar field  Field(np)
!
!     Field  : vector field of dimension nd defined at np points (sph)
!     np     : Number of particles (s+w)                            [in]
!     nd     : dimensions                                           [in]
!     grad_Field (nd,np) gradient of Field (np)                     [out]   

implicit none

integer(ink) nd, np
integer(ink) I, J, id, jd, ip
  
real   (irk) Field( np), grad_Field (nd,np) 
real   (irk) FJI, rJI (nd), OmegaI, OmegaJ
real   (irk) rhoIJ, h
real   (irk), allocatable:: xMat (:,:,:)
real   (irk) a11, a12, a21, a22, InvDet, b1, b2
if (.NOT.allocated (xMat) ) allocate ( xMat(nd,nd,np) )

grad_Field = 0.0 ; xMat = 0.0; FJI=0.0; rJI = 0.0

Field = 0.0
do ip=1,np
   field(ip) = x(1,ip)
enddo    

current=> last
Do while (associated(current))
   I = current%pair_i 
   J = current%pair_j 
   OmegaI = (mass(I)/rho(I)) 
   OmegaJ = (mass(J)/rho(J))
   if (current%Pint_type.eq.0) then
      FJI = Field(J) - Field(I) 
      do id = 1,nd
         rJI(id) = x(id,J) - x(id,I)
         grad_Field(id,I) = grad_Field(id,I) + OmegaJ * FJI * current%dwdx(id)   ! lcassical grad approx
         grad_Field(id,J) = grad_Field(id,J) + OmegaI * FJI * current%dwdx(id)
         do jd = 1,nd
            xMat(id,jd,I) = xMat(id,jd,I) + OmegaJ * current%dwdx(id) * rJI (jd)   
            xMat(id,jd,J) = xMat(id,jd,J) + OmegaI * current%dwdx(id) * rJI (jd) 
         enddo
      enddo 
   endif
   current=>current%next
enddo
!      ------  normalize. Solve M * grad_norm = grad Mab = RJI(b)*grad(a)
DO ip = 1,np
   if (nd.eq.1) then
      grad_Field(1,ip) = grad_Field(1,ip)/Xmat(1,1,ip)
   elseif(nd.eq.2) then
      a11 = xMat(1,1,ip)
      a12 = xMat(1,2,ip)    ! a12
      a21 = xMat(2,1,ip)    ! a21
      a22 = xMat(2,2,ip)
      Invdet = 1. / (a11*a22 - a12*a21)
      b1  = grad_Field(1,ip)
      b2  = grad_Field(2,ip)
      grad_Field (1,ip) = InvDet * (  a22*b1 - a12*b2)
      grad_Field (2,ip) = InvDet * ( -a21*b1 + a11*b2)
   else
     STOP 
   endif 
enddo

grad_Field = 0.0
current=> last          !    Initializes linked list to last element 
DO WHILE (associated(current))

   i = current%pair_i 
   j = current%pair_j 

   IF ( current%Pint_type.eq.0 ) THEN
 
      if(pa_sph.eq.1) then                  !   ----  For SPH algorithm 1  
         rhoij = 1.e0/(rho(i)*rho(j))                      
         do id = 1,nd
            h = (Field(i) + Field(j))*current%dwdx(id)        
            h = h*rhoij
            grad_Field(id,i) = grad_Field(id,i) + mass(j)*h
            grad_Field(id,j) = grad_Field(id,j) - mass(i)*h
         enddo  
             
      else if (pa_sph.eq.2) then            !      ----- For SPH algorithm 2     
         do id = 1,nd 
            h = (Field(i)/rho(i)**2 + Field(j)/rho(j)**2)*current%dwdx(id)  
            grad_Field(id,i) = grad_Field(id,i) + mass(j) * h
            grad_Field(id,j) = grad_Field(id,j) - mass(i) * h
         enddo
      endif

   ENDIF
   
   current=>current%next 
   
ENDDO

deallocate ( Xmat )

END SUBROUTINE get_gradient_TOOLS


!---------------------------------------------------------------

     subroutine get_divergence_TOOLS (Field, nd, np, div_Field)

!---------------------------------------------------------------
 

!     Subroutine to calculate the divergence of a vector field F
!
!     Field  : vector field of dimension nd defined at np points (sph)
!     np     : Number of particles (s+w)                            [in]
!     nd     : dimensions                                           [in]
!     dwdx   : derivation of Kernel for all interaction pairs       [in]
!     div_Field : divergence                                        [out]   

implicit none

integer(ink) nd, np, i ,j ,id, ip
  
real   (irk) Field(nd, np), div_Field (np), div_Corr(np)
real   (irk) dF (nd), dR (nd), Fcc, Rcc
real   (irk), allocatable:: dField(:)

if (.NOT.allocated (dField) ) allocate ( dField(nd) )

div_Field = 0.0 ; div_Corr = 0.0; dF=0.0; dR = 0.0

do i = 1,np
   Field(1,i) = x(1,i) ! 0.0         
   if (nd.eq.2) Field(2,i) =  0.0   ! x(2,i)      
enddo

current=> last
Do while (associated(current))
   i = current%pair_i 
   j = current%pair_j 
   if (current%Pint_type.eq.0) then
      do id = 1,nd
         dF(id) = Field(id,J) - Field(id,I) 
         dr(id) =     x(id,J) -     x(id,I)
      enddo        
      Fcc = dF(1)*current%dwdx(1)   
      Rcc = dr(1)*current%dwdx(1)     
      do id  = 2,nd 
         Fcc = Fcc + dF(id )*current%dwdx(id)
         Rcc = Rcc + dr(id) *current%dwdx(id)
      enddo     
      div_Field(i) = div_Field(i) + mass(j)*Fcc/rho(J)
      div_Field(j) = div_Field(j) + mass(i)*Fcc/rho(I)    
      div_Corr (i) = div_Corr (i) + mass(j)*Rcc/rho(J)
      div_Corr (j) = div_Corr (j) + mass(i)*Rcc/rho(I)
   endif
   current=>current%next
enddo
!      ------  normalize
DO ip = 1, np
   div_Field(ip) = div_Field(ip)/div_Corr(ip)
ENDDO 

END SUBROUTINE get_divergence_TOOLS



!-------------------------------------------------------------------

       Subroutine Solid_Slab
       
!-------------------------------------------------------------------

!      ----   18th may 2019  solid slab with accels

implicit none

integer (ink) ipoin, idimn

real    (irk) Mtot, rgi2, Ig, m_rgia, alphaG, fact_s2w
real    (irk), allocatable, SAVE::  sum_ma (:), sum_mx(:), aG (:), XG(:), rgi (:,:), acci(:)    !

If (.NOT.allocated ( sum_ma) ) then
   allocate ( sum_ma (ndimn), sum_mx(ndimn), aG (ndimn), XG(ndimn), rgi (ndimn,npoin_s), acci(ndimn) )
endif

Mtot=0.0;    rgi2  =  0.0; Ig  = 0.0; m_rgia = 0.0; alphaG = 0.0
sum_ma = 0.0; sum_mx = 0.0; aG  = 0.0; XG     = 0.0; rgi    = 0.0; acci = 0.0

if (ic_end_solid.eq.1) then        !  solid behaviour via decreasing p
    fact_s2w = 1.0                 ! default  MP 29 Jan 2034
    if (law_end_solid.eq.0) then
        if (time_sph.lt.t1_end_solid) fact_s2w = 0.0
    elseif ( law_end_solid.eq.1 ) then
       if  (time_sph.lt.t1_end_solid) then
           fact_s2w = (time_sph/t1_end_solid)**fct_exp_end_solid 
       endif
    elseif (law_end_solid.eq.2) then
       fact_s2w = 0.0
       if  (time_sph.gt.t1_end_solid) then
          fact_s2w = 1.0 - exp( -fct_exp_end_solid *  (time_sph-t1_end_solid) )
       endif 
    endif
endif
 
DO ipoin = 1, npoin_s
   Mtot       = Mtot + mass(ipoin)
   sum_ma (:) = sum_ma(:) + mass(ipoin)*exdvxdt(:,ipoin)
   sum_mx (:) = sum_mx(:) + mass(ipoin)*x      (:,ipoin)
enddo

aG = sum_ma/Mtot
xG = sum_mx/Mtot
 
if (ndimn.eq.2) then

   Do ipoin = 1, npoin_s
      rgi(:,ipoin) = x(:,ipoin) - XG(:)
      rgi2         = ( rgi(1,ipoin)*  rgi(1,ipoin) +  rgi(2,ipoin)* rgi(2,ipoin) )
      IG           = IG + mass(ipoin) * rgi2
      m_rgia       = m_rgia + mass(ipoin)*(rgi(1,ipoin)*exdvxdt(2,ipoin)-rgi(2,ipoin)*exdvxdt(1,ipoin))
   enddo
   
   alphaG = m_rgia/IG

   do ipoin = 1,npoin_s
      acci(1) = aG(1) - alphaG * rgi(2,ipoin)
      acci(2) = aG(2) + alphaG * rgi(1,ipoin)
      exdvxdt(1:ndimn,ipoin) = (1-fact_s2w)*(acci(1:ndimn)) + fact_s2w * exdvxdt(1:ndimn,ipoin)
   enddo

elseif (ndimn.eq.1) then
   acci(1) = aG(1)
   do ipoin = 1,npoin_s
      exdvxdt(1,ipoin) = (1-fact_s2w)*(acci(1)) + fact_s2w * exdvxdt(1,ipoin)
   enddo
endif

End Subroutine Solid_Slab

!----------------------------------------------------------------------
 
       subroutine grid_find_NEW_tools2  (itimestep_sph)  
 
!----------------------------------------------------------------------

!      ------  This is the new routine
 
implicit none
 
integer(ink), SAVE:: mdivt, ndivt, m_pairs, n_pairs 

integer(ink), SAVE, allocatable:: which_cell(:), info_picell(:,:), list_picell(:)
real   (irk), SAVE, allocatable:: xmin(:), xmax(:), deltx(:)
real   (irk), SAVE, allocatable:: dxiac(:), tdwdx(:)

real   (irk), SAVE, allocatable:: xR(:,:)         !   rotated coordinates to minimize nr of cells (ndimn, ntotal)
real   (irk)  IXX, IXY, IYY, Thetam, cTh, sTh     ! inertia moments. rotate coords along main inertia axes CH MP april 2020

!! type (pairs), SAVE, pointer::keepit0, keepit1 26 Aug   !!!! ! *****   MP sept 2022  , current     
                                                         ! to global main vars module 
 
integer(ink)  ndivx(3),idivx(3),jdivx(3)
integer(ink)  itimestep_sph
integer(ink)  idimn, itotal, jtotal
integer(ink)  idivt,   jdivt
integer(ink)  idest, ipost,     i0,   idelt
integer(ink)  idx, idy, idz, jdx, jdy, jdz, idx0, idy0, idz0, idx1, idy1, idz1
integer(ink)  sumiac, maxiac, miniac, noiac, maxp, minp 
integer(ink)  npicsi, npicsj, jj, jdest, j0
integer(ink)    i,   j,   k
integer(ink)  scale_k 
integer(ink)  countiac_real(ntotal)

 
real   (irk)  length, length_new 
real   (irk)  driac, r, mhsml, w 

countiac_real = 0
!      ------  The first time we have not allocated anything. Do it

if (.NOT.allocated(dxiac) ) then                             ! ** MP FS
   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
   allocate ( dxiac(ndimn), tdwdx(ndimn) )
   allocate ( xR(ndimn, ntotal))                            ! CHGD MP April2020
endif

!if (itimestep_sph.eq.1) then
!   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
!   allocate ( dxiac(ndimn), tdwdx(ndimn) )
!endif
xmin = 1.e+10; xmax = -xmin; deltx = 0.0; ndivx = 1; idivx = 1; jdivx = 1
 
if (skf.eq.1.or.skf.eq.4) then    ! Saeidmt 24Oct2020 Wendland kernel (skf=4)
   scale_k = 2  
else if (skf.eq.2) then 
   scale_k = 3 
else if (skf.eq.3) then 
   scale_k = 3 
endif 

niac =0; countiac = 0

!      ------  Task 1a: find ppal inertia axes and rotate coords x ->x R

if (ndimn.eq.2 ) then 
   
   IXX = 0.0; IXY = 0.0; IYY = 0.0
   Do itotal=1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      IXX = IXX + x(2,itotal)*x(2,itotal)
      IXY = IXY + x(1,itotal)*x(2,itotal)
      IYY = IYY + x(1,itotal)*x(1,itotal)
   enddo

   Thetam = -0.5*atan2 (2*IXY, (IXX-IYY))
   cTh    = cos (Thetam)
   sTh    = sin (Thetam)

   Do itotal=1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      xR(1,itotal) =  cTH*x(1,itotal) +  sTh *x(2,itotal) 
      xR(2,itotal) = -sTH*x(2,itotal) +  cTh *x(2,itotal) 
   enddo
ELSE
   xR = x
ENDIF

!      ------  Task 1b: Setup Grid parameters
   
Do idimn = 1,ndimn
   Do itotal = 1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      if (xR(idimn,itotal).lt.xmin(idimn)) xmin(idimn)  = xR(idimn,itotal)
      if (xR(idimn,itotal).gt.xmax(idimn)) xmax(idimn)  = xR(idimn,itotal)
      if (deltx(idimn).lt.hsml(itotal))   deltx(idimn) = hsml(itotal)
   enddo
   deltx(idimn) = deltx(idimn)*2       
   length       = xmax(idimn) - xmin(idimn)
   ndivx(idimn) = (length/deltx(idimn)) + 1
   length_new   = ndivx(idimn)*deltx(idimn)
   xmin(idimn)  = xmin(idimn) - (length_new-length)/2 - 0.001*length
   xmax(idimn)  = xmax(idimn) + (length_new-length)/2 + 0.001*length
Enddo

ndivt = ndivx(1)*ndivx(2)*ndivx(3) 

!      ------  Now we can allocate or reallocate some control arrays

if (.NOT.allocated ( which_cell)  ) then   !  ** MP FS
    mdivt = ndivt
    allocate ( which_cell(ntotal), info_picell(3,mdivt),  list_picell(ntotal) )
else if (ndivt.gt.mdivt) then
    deallocate(info_picell)
    mdivt = max (2*mdivt, ndivt)
    allocate  (info_picell(3,mdivt) )
endif

which_cell = 0; info_picell=0; list_picell=0

!      ------  Task 2: Fill Grid structure with particles

DO itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   DO idimn = 1, ndimn
      idivx(idimn) = (xR(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
      if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
   enddo
   idivt =  ndivx(1)*ndivx(2)*(idivx(3)-1) + ndivx(1)*(idivx(2)-1) + idivx(1)
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

do itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   idivt  = which_cell  (itotal)
   info_picell(3,idivt) = info_picell(3,idivt) + 1 
   idest  = info_picell(2,idivt) + info_picell(3,idivt)
   list_picell(idest) = itotal
enddo

!      ------  Task 3: Get list of interactions. As nodes are in cells, we use x again


!      ------  Task 3: Get list of interactions

if (itimestep_sph.LE.1.AND.nproblems_index.EQ.1) then   !  MP 26 Aug 2022          
   m_pairs = 0
   n_pairs = 0
   nullify (last)
   allocate (keepit00, keepit11)
else                                   ! *****   MP aug 2022
   if (.NOT.associated(current)) then  ! *****   MP aug 2022
      current => last                  ! *****   MP aug 2022
   endif                               ! *****   MP aug 2022
   if (n_pairs.lt.m_pairs) then
       keepit00%next => keepit11     
       current => last
!       call cross_link 
   endif
endif      

n_pairs = 0
DO idz = 1, ndivx(3)
DO idy = 1, ndivx(2)
DO idx = 1, ndivx(1)
   idivt  =  ndivx(1)*ndivx(2)*(idz-1) + ndivx(1)*(idy-1) + idx 
   npicsi = info_picell(1,idivt)
   DO i = 1, npicsi
      i0     = info_picell(2,idivt)
      idest  = i0 + i -1 
      itotal = list_picell(idest)
      idx0   = max (       1,idx-1); idy0   = max (       1,idy-1); idz0   = max (       1,idz-1)
      idx1   = min (ndivx(1),idx+1); idy1   = min (ndivx(2),idy+1); idz1   = min (ndivx(3),idz+1)
      DO jdz = idz, idz1        !  before it was idz0, idz1
      DO jdy = idy, idy1        !  before it was idy0, idy1
      DO jdx = idx0,idx1 
         jdivt = ndivx(1)*ndivx(2)*(jdz-1) + ndivx(1)*(jdy-1) + jdx
         if (jdivt.lt.idivt) CYCLE
         jj = 1
         if (idivt.eq.jdivt) jj = i + 1
         npicsj = info_picell(1,jdivt)
         DO j = jj,npicsj
            j0     = info_picell(2,jdivt)
            jdest  = j0 + j -1
            jtotal = list_picell(jdest)
            dxiac(1) = x(1,itotal)-x(1,jtotal)          !          note jdivt is always greater than idivt
            driac    = dxiac(1)*dxiac(1)
        	do idimn = 2, ndimn
               dxiac(idimn) = x(idimn,itotal)-x(idimn,jtotal)
               driac = driac + dxiac(idimn)*dxiac(idimn)
            enddo
            mhsml =  (hsml(itotal)+hsml(jtotal))/2.
            !mhsml = min(hsml(itotal),hsml(jtotal))
            r = sqrt(driac)
            if (r.lt.scale_k*mhsml) then
               n_pairs = n_pairs + 1
               countiac(itotal) = countiac(itotal) + 1
               countiac(jtotal) = countiac(jtotal) + 1
               If (itype(itotal).ne.itype(jtotal)) then
			   countiac_real(itotal) = countiac_real(itotal) +1
			   countiac_real(jtotal) = countiac_real(jtotal) +1 
			   Endif
			   if (n_pairs.gt.m_pairs) then
                   allocate(current)                    !      create next pair structure
                   allocate (current%dwdx(ndimn))       !      allocate dwdx
                   current%next => last                 !      
                   last         => current              !      
                   m_pairs      =  m_pairs + 1          ! 
               endif                                    !    
               
               call kernel_TOOLS (r, dxiac, mhsml, w, tdwdx)   !  chgd MP April 2020
               current%niac      = n_pairs             !      this is the pair number 
               current%pair_i    = itotal              !      first particle
               current%pair_j    = jtotal              !      and second particle
               current%Pint_type = 0                   !      Sets interaction type to default (same material)
               current%w         = w                   !      Weigth
               current%dwdx(:)   = tdwdx(:)            !      Gradient of weight
               
               if (n_pairs.lt.m_pairs) then            !
                   keepit00  => current                !   keepit00 is global  Sept 2022 MP                 
                   current   => current%next           !   
               endif
               
            endif
         ENDDO
         
      ENDDO
      ENDDO
      ENDDO
      
   ENDDO
   
ENDDO
ENDDO
ENDDO


If (n_pairs.lt.m_pairs) then   ! ****      
   keepit11 => current  
   nullify (keepit00%next)
endif   

current=> last
 
!      ------  Statistics for the interaction

sumiac = 0
maxiac = 0
miniac = 1000
noiac  = 0
DO i = 1, ntotal
   sumiac = sumiac + countiac(i)
   if (countiac(i).gt.maxiac) then
      maxiac = countiac(i)
      maxp = i
   endif
   if (countiac(i).lt.miniac) then 
      miniac = countiac(i)
      minp = i
   endif
   if (countiac(i).eq.0)      noiac  = noiac + 1
enddo
 
!if (mod(itimestep_sph,print_step).eq.0) then

if (t_print_reset.ge.time_print) then
       write(*,*) '   '
       write(*,*) '  statistics from sph main grid_find  '
       write(*,*) '  time= ',time, ' dt= ', dt_sph
       write(*,*) '  max interactions= ', maxiac, ' at pic ', maxp
       write(*,*) '  min interactions= ', miniac, ' at pic ', minp
       write(*,*) '  average ', real(sumiac)/real(ntotal)
       write(*,*) '  Total pairs : ',niac 
       write(*,*) '  Pics w/o interaction: ', noiac
       write(*,*) '  maxboxes, nboxes   ', mdivt, ndivt      ! MP Jan 2022
       write(*,*) '  Mpairs, npairs     ', m_pairs, n_pairs
endif    

END SUBROUTINE  grid_find_NEW_Tools2 

!----------------------------------------------------------------------
 
       subroutine grid_find_NEW_tools3  (itimestep_sph)  
 
!----------------------------------------------------------------------

!      ------  This is the new routine
 
implicit none
 
integer(ink), SAVE:: mdivt, ndivt, m_pairs, n_pairs 

integer(ink), SAVE, allocatable:: which_cell(:), info_picell(:,:), list_picell(:)
real   (irk), SAVE, allocatable:: xmin(:), xmax(:), deltx(:)
real   (irk), SAVE, allocatable:: dxiac(:), tdwdx(:)

real   (irk), SAVE, allocatable:: xR(:,:)         !   rotated coordinates to minimize nr of cells (ndimn, ntotal)
real   (irk)  IXX, IXY, IYY, Thetam, cTh, sTh     ! inertia moments. rotate coords along main inertia axes CH MP april 2020
real   (irk), SAVE, allocatable:: xG(:)           ! grav centre coords
real   (irk)  xrelG, yrelG
integer(ink)  itotal_sum                          ! nr active nodes (not out of domain)

!! type (pairs), SAVE, pointer::keepit0, keepit1 26 Aug   !!!! ! *****   MP aug 2022  , current     
                                                         ! to global main vars module 
integer(ink)  ndivx(3),idivx(3),jdivx(3)
integer(ink)  itimestep_sph
integer(ink)  idimn, itotal, jtotal
integer(ink)  idivt,   jdivt
integer(ink)  idest, ipost,     i0,   idelt
integer(ink)  idx, idy, idz, jdx, jdy, jdz, idx0, idy0, idz0, idx1, idy1, idz1
integer(ink)  sumiac, maxiac, miniac, noiac, maxp, minp 
integer(ink)  npicsi, npicsj, jj, jdest, j0
integer(ink)    i,   j,   k
integer(ink)  scale_k 
integer(ink)  countiac_real(ntotal)

 
real   (irk)  length, length_new 
real   (irk)  driac, r, mhsml, w 

countiac_real = 0
!      ------  The first time we have not allocated anything. Do it

if (.NOT.allocated(dxiac) ) then                             ! ** MP FS
   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
   allocate ( dxiac(ndimn), tdwdx(ndimn) )
   allocate ( xR(ndimn, ntotal))                            ! CHGD MP April2020
   allocate ( xG(ndimn) )
endif

!if (itimestep_sph.eq.1) then
!   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
!   allocate ( dxiac(ndimn), tdwdx(ndimn) )
!endif
xmin = 1.e+10; xmax = -xmin; deltx = 0.0; ndivx = 1; idivx = 1; jdivx = 1
 
if (skf.eq.1.or.skf.eq.4) then    ! Saeidmt 24Oct2020 Wendland kernel (skf=4)
   scale_k = 2 
else if (skf.eq.2) then 
   scale_k = 3 
else if (skf.eq.3) then 
   scale_k = 3 
endif 

niac =0; countiac = 0

!      ------  Task 1a: find ppal inertia axes and rotate coords x ->x R

if (ndimn.eq.2 ) then 
   
   xG = 0.0 ; itotal_sum = 0
   Do itotal=1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      itotal_sum = itotal_sum +1
      xG (:)     = xG(:) + x(:,itotal)  
   enddo
   xG (:) = xG (:)/itotal_sum
      
   IXX = 0.0; IXY = 0.0; IYY = 0.0
   Do itotal = 1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      xrelG = x(1,itotal)-xG(1)  
      yrelG = x(2,itotal)-xG(2)
      IXX = IXX + yrelG * yrelG
      IXY = IXY + xrelG * yrelG
      IYY = IYY + xrelG * xrelG
   enddo

 !  Thetam = -0.5*atan2 (2*IXY, (IXX-IYY))
   Thetam =  - 0.5*atan2 (2*IXY, (IYY-IXX))
   cTh    = cos (Thetam)
   sTh    = sin (Thetam)

   Do itotal=1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      xrelG = x(1,itotal)-xG(1)  
      yrelG = x(2,itotal)-xG(2)
      xR(1,itotal) =  cTH*xrelG +  sTh *yrelG
      xR(2,itotal) = -sTH*xrelG +  cTh *yrelG
   enddo
ELSE
   xR = x
ENDIF

!      ------  Task 1b: Setup Grid parameters
   
Do idimn = 1,ndimn
   Do itotal = 1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      if (xR(idimn,itotal).lt.xmin(idimn)) xmin(idimn)  = xR(idimn,itotal)
      if (xR(idimn,itotal).gt.xmax(idimn)) xmax(idimn)  = xR(idimn,itotal)
      if (deltx(idimn).lt.hsml(itotal))   deltx(idimn) = hsml(itotal)
   enddo
   deltx(idimn) = deltx(idimn)*2       
   length       = xmax(idimn) - xmin(idimn)
   ndivx(idimn) = (length/deltx(idimn)) + 1
   length_new   = ndivx(idimn)*deltx(idimn)
   xmin(idimn)  = xmin(idimn) - (length_new-length)/2 - 0.001*length
   xmax(idimn)  = xmax(idimn) + (length_new-length)/2 + 0.001*length
Enddo

ndivt = ndivx(1)*ndivx(2)*ndivx(3) 

!      ------  Now we can allocate or reallocate some control arrays

if (.NOT.allocated ( which_cell)  ) then   !  ** MP FS
    mdivt = ndivt
    allocate ( which_cell(ntotal), info_picell(3,mdivt),  list_picell(ntotal) )
else if (ndivt.gt.mdivt) then
    deallocate(info_picell)
    mdivt = max (2*mdivt, ndivt)
    allocate  (info_picell(3,mdivt) )
endif

which_cell = 0; info_picell=0; list_picell=0

!      ------  Task 2: Fill Grid structure with particles

DO itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   DO idimn = 1, ndimn
      idivx(idimn) = (xR(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
      if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
   enddo
   idivt =  ndivx(1)*ndivx(2)*(idivx(3)-1) + ndivx(1)*(idivx(2)-1) + idivx(1)
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

do itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   idivt  = which_cell  (itotal)
   info_picell(3,idivt) = info_picell(3,idivt) + 1 
   idest  = info_picell(2,idivt) + info_picell(3,idivt)
   list_picell(idest) = itotal
enddo

!      ------  Task 3: Get list of interactions

if (itimestep_sph.LE.1.AND.nproblems_index.EQ.1) then   !  MP 7 sept 2022          
   m_pairs = 0
   n_pairs = 0
   nullify (last)
   allocate (keepit00, keepit11)
else                                   ! *****   MP sept 2022
   if (.NOT.associated(current)) then  ! *****   MP aug 2022
      current => last                  ! *****   MP aug 2022
   endif                               ! *****   MP aug 2022
   if (n_pairs.lt.m_pairs) then
       keepit00%next => keepit11     
       current => last
!       call cross_link 
   endif
endif  

n_pairs = 0
DO idz = 1, ndivx(3)
DO idy = 1, ndivx(2)
DO idx = 1, ndivx(1)
   idivt  =  ndivx(1)*ndivx(2)*(idz-1) + ndivx(1)*(idy-1) + idx 
   npicsi = info_picell(1,idivt)
   DO i = 1, npicsi
      i0     = info_picell(2,idivt)
      idest  = i0 + i -1 
      itotal = list_picell(idest)
      idx0   = max (       1,idx-1); idy0   = max (       1,idy-1); idz0   = max (       1,idz-1)
      idx1   = min (ndivx(1),idx+1); idy1   = min (ndivx(2),idy+1); idz1   = min (ndivx(3),idz+1)
      DO jdz = idz, idz1        !  before it was idz0, idz1
      DO jdy = idy, idy1        !  before it was idy0, idy1
      DO jdx = idx0,idx1 
         jdivt = ndivx(1)*ndivx(2)*(jdz-1) + ndivx(1)*(jdy-1) + jdx
         if (jdivt.lt.idivt) CYCLE
         jj = 1
         if (idivt.eq.jdivt) jj = i + 1
         npicsj = info_picell(1,jdivt)
         DO j = jj,npicsj
            j0     = info_picell(2,jdivt)
            jdest  = j0 + j -1
            jtotal = list_picell(jdest)
            dxiac(1) = x(1,itotal)-x(1,jtotal)          !          note jdivt is always greater than idivt
            driac    = dxiac(1)*dxiac(1)
        	do idimn = 2, ndimn
               dxiac(idimn) = x(idimn,itotal)-x(idimn,jtotal)
               driac = driac + dxiac(idimn)*dxiac(idimn)
            enddo
            mhsml =  (hsml(itotal)+hsml(jtotal))/2.
            !mhsml = min(hsml(itotal),hsml(jtotal))
            r = sqrt(driac)
            if (r.lt.scale_k*mhsml) then
               n_pairs = n_pairs + 1
               countiac(itotal) = countiac(itotal) + 1
               countiac(jtotal) = countiac(jtotal) + 1
               If (itype(itotal).ne.itype(jtotal)) then
			   countiac_real(itotal) = countiac_real(itotal) +1
			   countiac_real(jtotal) = countiac_real(jtotal) +1 
			   Endif
			   if (n_pairs.gt.m_pairs) then
                   allocate(current)                    !      create next pair structure
                   allocate (current%dwdx(ndimn))       !      allocate dwdx
                   current%next => last                 !      
                   last         => current              !      
                   m_pairs      =  m_pairs + 1          ! 
               endif                                    !    
               
               call kernel_TOOLS (r, dxiac, mhsml, w, tdwdx)   !  chgd MP April 2020
               current%niac      = n_pairs             !      this is the pair number 
               current%pair_i    = itotal              !      first particle
               current%pair_j    = jtotal              !      and second particle
               current%Pint_type = 0                   !      Sets interaction type to default (same material)
               current%w         = w                   !      Weigth
               current%dwdx(:)   = tdwdx(:)            !      Gradient of weight
               
               if (n_pairs.lt.m_pairs) then            !
                   keepit00  => current                !   keepit00 is global   sept 2022 MP              
                   current   => current%next           !   
               endif
               
               
            endif
         ENDDO
         
      ENDDO
      ENDDO
      ENDDO
      
   ENDDO
   
ENDDO
ENDDO
ENDDO

If (n_pairs.lt.m_pairs) then   ! ****      
   keepit11 => current  
   nullify (keepit00%next)
endif

current=> last
 
!      ------  Statistics for the interaction

sumiac = 0
maxiac = 0
miniac = 1000
noiac  = 0
DO i = 1, ntotal
   sumiac = sumiac + countiac(i)
   if (countiac(i).gt.maxiac) then
      maxiac = countiac(i)
      maxp = i
   endif
   if (countiac(i).lt.miniac) then 
      miniac = countiac(i)
      minp = i
   endif
   if (countiac(i).eq.0)      noiac  = noiac + 1
enddo
 
!if (mod(itimestep_sph,print_step).eq.0) then

if (t_print_reset.ge.time_print) then
       write(*,*) '  -------  ALG IXY   3  --------------- '
       write(*,*) '  statistics from sph main grid_find  '
       write(*,*) '  time= ',time, ' dt= ', dt_sph
       write(*,*) '  max interactions= ', maxiac, ' at pic ', maxp
       write(*,*) '  min interactions= ', miniac, ' at pic ', minp
       write(*,*) '  average ', real(sumiac)/real(ntotal)
       write(*,*) '  Total pairs : ',niac 
       write(*,*) '  Pics w/o interaction: ', noiac
       write(*,*) '  maxboxes, nboxes   ', mdivt, ndivt      ! MP Jan 2022
       write(*,*) '  Mpairs, npairs     ', m_pairs, n_pairs
endif    


END SUBROUTINE  grid_find_NEW_Tools3 


!----------------------------------------------------------------------
 
       subroutine grid_find_NEW_tools  (itimestep_sph)  

!----------------------------------------------------------------------

!      ------  This is the memory corrected routine Aug 2022 MP
 
implicit none
 
integer(ink), SAVE:: mdivt, ndivt, m_pairs, n_pairs 

integer(ink), SAVE, allocatable:: which_cell(:), info_picell(:,:), list_picell(:)
real   (irk), SAVE, allocatable:: xmin(:), xmax(:), deltx(:)
real   (irk), SAVE, allocatable:: dxiac(:), tdwdx(:)

!! type (pairs), SAVE, pointer::keepit0, keepit1 26 Aug   !!!! ! *****   MP aug 2022  , current     
                                                         ! to global main vars module
integer(ink)  ndivx(3),idivx(3),jdivx(3)
integer(ink)  itimestep_sph
integer(ink)  idimn, itotal, jtotal
integer(ink)  idivt,   jdivt
integer(ink)  idest, ipost,     i0,   idelt
integer(ink)  idx, idy, idz, jdx, jdy, jdz, idx0, idy0, idz0, idx1, idy1, idz1
integer(ink)  sumiac, maxiac, miniac, noiac, maxp, minp 
integer(ink)  npicsi, npicsj, jj, jdest, j0
integer(ink)    i,   j,   k
integer(ink)  scale_k 
integer(ink)  countiac_real(ntotal)

 
real   (irk)  length, length_new 
real   (irk)  driac, r, mhsml, w 

countiac_real = 0
!      ------  The first time we have not allocated anything. Do it

if (.NOT.allocated(dxiac) ) then                             ! ** MP FS
   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
   allocate ( dxiac(ndimn), tdwdx(ndimn) )
endif

!if (itimestep_sph.eq.1) then
!   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
!   allocate ( dxiac(ndimn), tdwdx(ndimn) )
!endif
xmin = 1.e+10; xmax = -xmin; deltx = 0.0; ndivx = 1; idivx = 1; jdivx = 1
 
if (skf.eq.1.or.skf.eq.4) then    ! Saeidmt 24Oct2020 Wendland kernel (skf=4)
   scale_k = 2 
else if (skf.eq.2) then 
   scale_k = 3 
else if (skf.eq.3) then 
   scale_k = 3 
endif 

niac =0; countiac = 0

!      ------  Task 1: Setup Grid parameters

  
Do idimn = 1,ndimn
   Do itotal = 1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      if (x(idimn,itotal).lt.xmin(idimn)) xmin(idimn)  = x(idimn,itotal)
      if (x(idimn,itotal).gt.xmax(idimn)) xmax(idimn)  = x(idimn,itotal)
      if (deltx(idimn).lt.hsml(itotal))   deltx(idimn) = hsml(itotal)
   enddo
   deltx(idimn) = deltx(idimn)*2       
   length       = xmax(idimn) - xmin(idimn)
   ndivx(idimn) = (length/deltx(idimn)) + 1
   length_new   = ndivx(idimn)*deltx(idimn)
   xmin(idimn)  = xmin(idimn) - (length_new-length)/2 - 0.001*length
   xmax(idimn)  = xmax(idimn) + (length_new-length)/2 + 0.001*length
Enddo

ndivt = ndivx(1)*ndivx(2)*ndivx(3) 

!      ------  Now we can allocate or reallocate some control arrays

if (.NOT.allocated ( which_cell)  ) then   !  ** MP FS
    mdivt = ndivt
    allocate ( which_cell(ntotal), info_picell(3,mdivt),  list_picell(ntotal) )
else if (ndivt.gt.mdivt) then
    deallocate(info_picell)
    mdivt = max (2*mdivt, ndivt)
    allocate  (info_picell(3,mdivt) )
endif

which_cell = 0; info_picell=0; list_picell=0

!      ------  Task 2: Fill Grid structure with particles

DO itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   DO idimn = 1, ndimn
      idivx(idimn) = (x(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
      if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
   enddo
   idivt =  ndivx(1)*ndivx(2)*(idivx(3)-1) + ndivx(1)*(idivx(2)-1) + idivx(1)
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

do itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   idivt  = which_cell  (itotal)
   info_picell(3,idivt) = info_picell(3,idivt) + 1 
   idest  = info_picell(2,idivt) + info_picell(3,idivt)
   list_picell(idest) = itotal
enddo

!      ------  Task 3: Get list of interactions

if (itimestep_sph.LE.1.AND.nproblems_index.EQ.1) then   !  MP 26 Aug 2022          
   m_pairs = 0
   n_pairs = 0
   nullify (last)
   allocate (keepit00, keepit11)
else                                   ! *****   MP aug 2022
   if (.NOT.associated(current)) then  ! *****   MP aug 2022
      current => last                  ! *****   MP aug 2022
   endif                               ! *****   MP aug 2022
   if (n_pairs.lt.m_pairs) then
       keepit00%next => keepit11     
       current => last
!       call cross_link 
   endif
endif   

n_pairs = 0
DO idz = 1, ndivx(3)
DO idy = 1, ndivx(2)
DO idx = 1, ndivx(1)
   idivt  =  ndivx(1)*ndivx(2)*(idz-1) + ndivx(1)*(idy-1) + idx 
   npicsi = info_picell(1,idivt)
   DO i = 1, npicsi
      i0     = info_picell(2,idivt)
      idest  = i0 + i -1 
      itotal = list_picell(idest)
      idx0   = max (       1,idx-1); idy0   = max (       1,idy-1); idz0   = max (       1,idz-1)
      idx1   = min (ndivx(1),idx+1); idy1   = min (ndivx(2),idy+1); idz1   = min (ndivx(3),idz+1)
      DO jdz = idz, idz1        !  before it was idz0, idz1
      DO jdy = idy, idy1        !  before it was idy0, idy1
      DO jdx = idx0,idx1 
         jdivt = ndivx(1)*ndivx(2)*(jdz-1) + ndivx(1)*(jdy-1) + jdx
         if (jdivt.lt.idivt) CYCLE
         jj = 1
         if (idivt.eq.jdivt) jj = i + 1
         npicsj = info_picell(1,jdivt)
         DO j = jj,npicsj
            j0     = info_picell(2,jdivt)
            jdest  = j0 + j -1
            jtotal = list_picell(jdest)
            dxiac(1) = x(1,itotal)-x(1,jtotal)          !          note jdivt is always greater than idivt
            driac    = dxiac(1)*dxiac(1)
        	do idimn = 2, ndimn
               dxiac(idimn) = x(idimn,itotal)-x(idimn,jtotal)
               driac = driac + dxiac(idimn)*dxiac(idimn)
            enddo
            mhsml =  (hsml(itotal)+hsml(jtotal))/2.
            !mhsml = min(hsml(itotal),hsml(jtotal))
            r = sqrt(driac)
            if (r.lt.scale_k*mhsml) then
               n_pairs = n_pairs + 1
               countiac(itotal) = countiac(itotal) + 1
               countiac(jtotal) = countiac(jtotal) + 1
               If (itype(itotal).ne.itype(jtotal)) then
			   countiac_real(itotal) = countiac_real(itotal) +1
			   countiac_real(jtotal) = countiac_real(jtotal) +1 
			   Endif
			   if (n_pairs.gt.m_pairs) then
                   allocate(current)                    !      create next pair structure
                   allocate (current%dwdx(ndimn))       !      allocate dwdx
                   current%next => last                 !      
                   last         => current              !      
                   m_pairs      =  m_pairs + 1          ! 
               endif                                    !    
               
!               call kernel(r, dxiac, mhsml, w, tdwdx)  !      using r, dxiac and mshml obtains w and its gradient
               call kernel_TOOLS (r, dxiac, mhsml, w, tdwdx)   !  chgd MP April 2020
               current%niac      = n_pairs             !      this is the pair number 
               current%pair_i    = itotal              !      first particle
               current%pair_j    = jtotal              !      and second particle
               current%Pint_type = 0                   !      Sets interaction type to default (same material)
               current%w         = w                   !      Weigth
               current%dwdx(:)   = tdwdx(:)            !      Gradient of weight
               
               if (n_pairs.lt.m_pairs) then            !
                   keepit00  => current                !   keepit00 is global                 
                   current   => current%next           !   
               endif
               
            endif
         ENDDO
         
      ENDDO
      ENDDO
      ENDDO
      
   ENDDO
   
ENDDO
ENDDO
ENDDO

If (n_pairs.lt.m_pairs) then   ! ****      
   keepit11 => current  
   nullify (keepit00%next)
endif

current=> last
 
!      ------  Statistics for the interaction

sumiac = 0
maxiac = 0
miniac = 1000
noiac  = 0
DO i = 1, ntotal
   sumiac = sumiac + countiac(i)
   if (countiac(i).gt.maxiac) then
      maxiac = countiac(i)
      maxp = i
   endif
   if (countiac(i).lt.miniac) then 
      miniac = countiac(i)
      minp = i
   endif
   if (countiac(i).eq.0)      noiac  = noiac + 1
enddo
 
!if (mod(itimestep_sph,print_step).eq.0) then

if (t_print_reset.ge.time_print) then 
       write(*,*) '   '
       write(*,*) '  statistics from sph main grid_find  '
       write(*,*) '  time= ',time, ' dt= ', dt_sph
       write(*,*) '  max interactions= ', maxiac, ' at pic ', maxp
       write(*,*) '  min interactions= ', miniac, ' at pic ', minp
       write(*,*) '  average ', real(sumiac)/real(ntotal)
       write(*,*) '  Total pairs : ',niac 
       write(*,*) '  Pics w/o interaction: ', noiac
       write(*,*) '  maxboxes, nboxes   ', mdivt, ndivt      ! MP Jan 2022
       write(*,*) '  Mpairs, npairs     ', m_pairs, n_pairs
endif    


END SUBROUTINE  grid_find_NEW_Tools


!----------------------------------------------------------------------
 
       subroutine grid_find_NEW_tools_TEST  (itimestep_sph)  
 
!----------------------------------------------------------------------

!      ------  This is the OLD routine which spilled memory
 
implicit none
 
integer(ink), SAVE:: mdivt, ndivt, m_pairs, n_pairs 

integer(ink), SAVE, allocatable:: which_cell(:), info_picell(:,:), list_picell(:)
real   (irk), SAVE, allocatable:: xmin(:), xmax(:), deltx(:)
real   (irk), SAVE, allocatable:: dxiac(:), tdwdx(:)

type (pairs), SAVE, pointer::keepit0, keepit1, current   
 
integer(ink)  ndivx(3),idivx(3),jdivx(3)
integer(ink)  itimestep_sph
integer(ink)  idimn, itotal, jtotal
integer(ink)  idivt,   jdivt
integer(ink)  idest, ipost,     i0,   idelt
integer(ink)  idx, idy, idz, jdx, jdy, jdz, idx0, idy0, idz0, idx1, idy1, idz1
integer(ink)  sumiac, maxiac, miniac, noiac, maxp, minp 
integer(ink)  npicsi, npicsj, jj, jdest, j0
integer(ink)    i,   j,   k
integer(ink)  scale_k 
integer(ink)  countiac_real(ntotal)

 
real   (irk)  length, length_new 
real   (irk)  driac, r, mhsml, w 

countiac_real = 0
!      ------  The first time we have not allocated anything. Do it

if (.NOT.allocated(dxiac) ) then                             ! ** MP FS
   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
   allocate ( dxiac(ndimn), tdwdx(ndimn) )
endif

!if (itimestep_sph.eq.1) then
!   allocate ( xmin(ndimn),  xmax(ndimn),  deltx(ndimn) )
!   allocate ( dxiac(ndimn), tdwdx(ndimn) )
!endif
xmin = 1.e+10; xmax = -xmin; deltx = 0.0; ndivx = 1; idivx = 1; jdivx = 1
 
if (skf.eq.1.or.skf.eq.4) then    ! Saeidmt 24Oct2020 Wendland kernel (skf=4)
   scale_k = 2 
else if (skf.eq.2) then 
   scale_k = 3 
else if (skf.eq.3) then 
   scale_k = 3 
endif 

niac =0; countiac = 0

!      ------  Task 1: Setup Grid parameters

  
Do idimn = 1,ndimn
   Do itotal = 1, ntotal
      if ( If_Out_Domain(itotal).eq.1 ) CYCLE
      if (x(idimn,itotal).lt.xmin(idimn)) xmin(idimn)  = x(idimn,itotal)
      if (x(idimn,itotal).gt.xmax(idimn)) xmax(idimn)  = x(idimn,itotal)
      if (deltx(idimn).lt.hsml(itotal))   deltx(idimn) = hsml(itotal)
   enddo
   deltx(idimn) = deltx(idimn)*2       
   length       = xmax(idimn) - xmin(idimn)
   ndivx(idimn) = (length/deltx(idimn)) + 1
   length_new   = ndivx(idimn)*deltx(idimn)
   xmin(idimn)  = xmin(idimn) - (length_new-length)/2 - 0.001*length
   xmax(idimn)  = xmax(idimn) + (length_new-length)/2 + 0.001*length
Enddo

ndivt = ndivx(1)*ndivx(2)*ndivx(3) 

!      ------  Now we can allocate or reallocate some control arrays

if (.NOT.allocated ( which_cell)  ) then   !  ** MP FS
    mdivt = ndivt
    allocate ( which_cell(ntotal), info_picell(3,mdivt),  list_picell(ntotal) )
else if (ndivt.gt.mdivt) then
    deallocate(info_picell)
    mdivt = max (2*mdivt, ndivt)
    allocate  (info_picell(3,mdivt) )
endif

which_cell = 0; info_picell=0; list_picell=0

!      ------  Task 2: Fill Grid structure with particles

DO itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   DO idimn = 1, ndimn
      idivx(idimn) = (x(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
      if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
   enddo
   idivt =  ndivx(1)*ndivx(2)*(idivx(3)-1) + ndivx(1)*(idivx(2)-1) + idivx(1)
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

do itotal = 1, ntotal
   if ( If_Out_Domain(itotal).eq.1 ) CYCLE
   idivt  = which_cell  (itotal)
   info_picell(3,idivt) = info_picell(3,idivt) + 1 
   idest  = info_picell(2,idivt) + info_picell(3,idivt)
   list_picell(idest) = itotal
enddo

!      ------  Task 3: Get list of interactions

! if (itimestep_sph.eq.1) then
if (itimestep_sph.LE.1) then             ! tried MP 2020 
   m_pairs = 0
   n_pairs = 0
   nullify (last)
   allocate (keepit0, keepit1)
else
   if (n_pairs.lt.m_pairs) then
       keepit0%next => keepit1
       current => last
!       call cross_link 
   endif
endif   

n_pairs = 0
DO idz = 1, ndivx(3)
DO idy = 1, ndivx(2)
DO idx = 1, ndivx(1)
   idivt  =  ndivx(1)*ndivx(2)*(idz-1) + ndivx(1)*(idy-1) + idx 
   npicsi = info_picell(1,idivt)
   DO i = 1, npicsi
      i0     = info_picell(2,idivt)
      idest  = i0 + i -1 
      itotal = list_picell(idest)
      idx0   = max (       1,idx-1); idy0   = max (       1,idy-1); idz0   = max (       1,idz-1)
      idx1   = min (ndivx(1),idx+1); idy1   = min (ndivx(2),idy+1); idz1   = min (ndivx(3),idz+1)
      DO jdz = idz, idz1        !  before it was idz0, idz1
      DO jdy = idy, idy1        !  before it was idy0, idy1
      DO jdx = idx0,idx1 
         jdivt = ndivx(1)*ndivx(2)*(jdz-1) + ndivx(1)*(jdy-1) + jdx
         if (jdivt.lt.idivt) CYCLE
         jj = 1
         if (idivt.eq.jdivt) jj = i + 1
         npicsj = info_picell(1,jdivt)
         DO j = jj,npicsj
            j0     = info_picell(2,jdivt)
            jdest  = j0 + j -1
            jtotal = list_picell(jdest)
            dxiac(1) = x(1,itotal)-x(1,jtotal)          !          note jdivt is always greater than idivt
            driac    = dxiac(1)*dxiac(1)
        	do idimn = 2, ndimn
               dxiac(idimn) = x(idimn,itotal)-x(idimn,jtotal)
               driac = driac + dxiac(idimn)*dxiac(idimn)
            enddo
            mhsml =  (hsml(itotal)+hsml(jtotal))/2.
            !mhsml = min(hsml(itotal),hsml(jtotal))
            r = sqrt(driac)
            if (r.lt.scale_k*mhsml) then
               n_pairs = n_pairs + 1
               countiac(itotal) = countiac(itotal) + 1
               countiac(jtotal) = countiac(jtotal) + 1
               If (itype(itotal).ne.itype(jtotal)) then
			   countiac_real(itotal) = countiac_real(itotal) +1
			   countiac_real(jtotal) = countiac_real(jtotal) +1 
			   Endif
			   if (n_pairs.gt.m_pairs) then
                   allocate(current)                    !      create next pair structure
                   allocate (current%dwdx(ndimn))       !      allocate dwdx
                   current%next => last                 !      
                   last         => current              !      
                   m_pairs      =  m_pairs + 1          ! 
               endif                                    !    
               
!               call kernel(r, dxiac, mhsml, w, tdwdx)  !      using r, dxiac and mshml obtains w and its gradient
               call kernel_TOOLS (r, dxiac, mhsml, w, tdwdx)   !  chgd MP April 2020
               current%niac      = n_pairs             !      this is the pair number 
               current%pair_i    = itotal              !      first particle
               current%pair_j    = jtotal              !      and second particle
               current%Pint_type = 0                   !      Sets interaction type to default (same material)
               current%w         = w                   !      Weigth
               current%dwdx(:)   = tdwdx(:)            !      Gradient of weight
               
               if (n_pairs.lt.m_pairs) then            !
                   keepit0      => current             !               
                   current      => current%next        ! 
               endif
               
            endif
         ENDDO
         
      ENDDO
      ENDDO
      ENDDO
      
   ENDDO
   
ENDDO
ENDDO
ENDDO

If (n_pairs.lt.m_pairs) then
   keepit1 => current
   nullify (keepit0%next)
endif

current=> last
 
!      ------  Statistics for the interaction

sumiac = 0
maxiac = 0
miniac = 1000
noiac  = 0
DO i = 1, ntotal
   sumiac = sumiac + countiac(i)
   if (countiac(i).gt.maxiac) then
      maxiac = countiac(i)
      maxp = i
   endif
   if (countiac(i).lt.miniac) then 
      miniac = countiac(i)
      minp = i
   endif
   if (countiac(i).eq.0)      noiac  = noiac + 1
enddo
 
!if (mod(itimestep_sph,print_step).eq.0) then

if (t_print_reset.ge.time_print) then
       write(*,*) '   '
       write(*,*) '  statistics from sph main grid_find  '
       write(*,*) '  time= ',time, ' dt= ', dt_sph
       write(*,*) '  max interactions= ', maxiac, ' at pic ', maxp
       write(*,*) '  min interactions= ', miniac, ' at pic ', minp
       write(*,*) '  average ', real(sumiac)/real(ntotal)
       write(*,*) '  Total pairs : ',niac 
       write(*,*) '  Pics w/o interaction: ', noiac
       write(*,*) '  maxboxes, nboxes   ', mdivt, ndivt      ! MP Jan 2022
       write(*,*) '  Mpairs, npairs     ', m_pairs, n_pairs
endif    


END SUBROUTINE  grid_find_NEW_Tools_TEST 



!---------------------------------------------------------------

     subroutine deallocate_list_TOOLS

!---------------------------------------------------------------

implicit none                         !  MP august 2022 to avoid increase of memory in multi-problems

current => last
do while (associated(current))
   last => current%next
   deallocate (current)
   current =>last
enddo



end subroutine deallocate_list_TOOLS


! -----------------------------------------------------------------
 
       SUBROUTINE Pint_Update_SW_tools
 
! -----------------------------------------------------------------

!   Subroutine to update soil-reservoir water interactions to 1
!   Soil=1  Water=2  VPs=-3 Abs=-12 Vn0=-14
 
implicit none

integer(ink) i, j, ii, jj, isumm, idiff 
integer(ink) i_ss_niac

ss_niac = 0

current=>last

DO WHILE (associated(current))
   
      i     = current%pair_i
      j     = current%pair_j
      ii    = iabs(itype(i))
      jj    = iabs(itype(j))
      isumm = itype(i) + itype(j)
      idiff = itype(i) - itype(j)
      
      if (jj.lt.ii) then                 !   First check. Valid combinations 
          current%pair_i    = j          !   (1,2);(1,-3);(1,-12);(2,-3),(2,-12)
          current%pair_j    = i
          current%dwdx(:)   = - current%dwdx(:)
      endif       
            
      i     = current%pair_i             ! make sure indexes are well ordered
      j     = current%pair_j    

      if (isumm.gt.0) then
      
          if (idiff.eq.0) then
              current%Pint_Type = 0                     !------  Soil-Soil (1,1) or water-water (2,2)
              ss_niac = ss_niac +1          ! MP ** FS  Chuan 24 april 2017
          else
              current%Pint_Type = 1                     !------  Soil-water interaction (1,2)
          endif
          
      elseif ((isumm.eq.-1).or.(isumm.eq.-2)) THEN      !------  This is soil-VP (1,-3) or water-VP (2,-3)
          current%Pint_Type = -1
          
      elseif ((isumm.eq.-10).or.(isumm.eq.-11)) THEN    !------  Soil-Abs(1,-12) or Water-Abs (2,-12) 
          current%Pint_Type = -11
          
      elseif ((isumm.eq.-12).or.(isumm.eq.-13)) THEN    !------  Soil-Vn0(1,-14) or Water-Vn0 (2,-14) 
          current%Pint_Type = -21
          
      else                                              !------  Abs-Abs(-12,-12) V-Abs(-3,-12) or V-V(-3,-3)...
          current%Pint_Type = -999
          
!      else
!          write(*,*) ' ====== pair wrong  i, type, j, type ======'
!          write(*,*) i, itype(i), j, itype(j)
!          write(*,*) ' stroke a number to finish '
!          read (*,*) isumm
      
      endif    
      
      current=>current%next
      
ENDDO
!      new ** MP FS  allocation of list of ss interactions for CG routine
!                    similar to Intmat: 

IF (ic_semi_implicit.eq.1) then   !  (SW Vars) alternative  SPH_t_Integ_alg=6 (MAIN vars)
   
   if (ss_niac.gt.ss_niac0) then                 !  if we have more interactions, dealloc & alloc
      if (allocated (list_ss_interacts)) then    !  only deallocate if present
          deallocate (list_ss_interacts)
      endif
      allocate ( list_ss_interacts (2,ss_niac) )  ! allocate the newly dimensioned  list
   else
      if(.NOT.allocated(list_ss_interacts)) then
         allocate ( list_ss_interacts (2,ss_niac) )
      endif    
   endif

   list_ss_interacts = 0
   
   i_ss_niac = 0    ! counter

   current=>last
   DO WHILE (associated(current))  
      i     = current%pair_i
      j     = current%pair_j
      ii    = iabs(itype(i))
      jj    = iabs(itype(j))
      isumm = itype(i) + itype(j)
      idiff = itype(i) - itype(j)
      if (isumm.gt.0) then 
         if (idiff.eq.0) then    !------  Soil-Soil (1,1) or water-water (2,2)
             i_ss_niac = i_ss_niac + 1
             list_ss_interacts (1, i_ss_niac) = i           
             list_ss_interacts (2, i_ss_niac) = j   
          endif              
      endif
      current=>current%next
   ENDDO

ENDIF

END subroutine Pint_Update_SW_Tools


!---------------------------------------------------------------

     subroutine sum_density_Tools

!---------------------------------------------------------------
 
!   Subroutine to calculate the density with SPH summation algorithm.

!     ntotal    : Number of particles                                  [in]
!     hsml      : Smoothing Length                                     [in]
!     mass      : Particle masses                                      [in]
!     niac      : Number of interaction pairs                          [in]
!     pair_i    : List of first partner of interaction pair            [in]
!     pair_j    : List of second partner of interaction pair           [in]
!     w         : Kernel for all interaction pairs                     [in]
!     itype     : type of particles                                   [in]
!     x         : Coordinates of all particles                        [in]
!     rho       : Density                                             [out]
!     wi(ntotal): integration of the kernel itself 
!     wii       : Self density of each particle: Wii (Kernel for distance 0)
!                 and take contribution of particle itself:
    
implicit none
integer(ink)  i, j     
real   (irk)  selfdens, r 
real   (irk), allocatable:: hv(:), wi(:)

allocate ( hv(ndimn), wi(ntotal) )
hv = 0.0; wi = 0.0
 
hv = 0.e0
r  = 0.
      
!      ------  First  calculate the integral of the kernel over the space

do i = 1,ntotal
   if (if_out_domain(i).eq.1) cycle
   ! call kernel(r,hv,hsml(i),selfdens,hv)
   call kernel_TOOLS (r,hv,hsml(i),selfdens,hv)  ! chgd MP April 2020
   wi(i)=selfdens*mass(i)/rho(i)
enddo

!      ------  use the linked list. Note we do not modify last, so we will can re-use it...

current=> last
Do while (associated(current))
   if (current%Pint_type.eq.0) then
      i = current%pair_i 
      j = current%pair_j 
      wi(i) = wi(i) + mass(j)/rho(j)*current%w 
      wi(j) = wi(j) + mass(i)/rho(i)*current%w
   endif   
   current=>current%next
enddo

!      ------  Compute the rho integral over the space

do i = 1,ntotal
   if (if_out_domain(i).eq.1) cycle
   !call kernel(r,hv,hsml(i),selfdens,hv)
   call kernel_TOOLS (r,hv,hsml(i),selfdens,hv)    ! CHGD MP April 2020
   rho(i) = selfdens*mass(i)
enddo

!      ------  Calculate SPH sum for rho:

current=> last
Do while (associated(current))
   if (current%Pint_type.eq.0) then
      i = current%pair_i 
      j = current%pair_j 
      rho(i) = rho(i) + mass(j)*current%w 
      rho(j) = rho(j) + mass(i)*current%w
   endif      
   current=>current%next
enddo

!      ------  Compute the normalized rho, rho=sum(rho)/sum(w)
!              at vn=0 nodes we normalize density 

if (nor_density) then 
   do i = 1, ntotal
   if (if_out_domain(i).eq.1) cycle                                             
      rho(i)=rho(i)/wi(i)
   enddo
endif  

if (ic2_vn0.eq.1.AND..not.nor_density) then   !   Bug Chuan found
   do i = 1, npoin
      if (if_out_domain(i).eq.1) cycle 
      if ( if_correct(i).eq.1) then
         rho(i)=rho(i)/wi(i)
      endif
   enddo
endif   
 
deallocate (hv,wi)

END SUBROUTINE sum_density_Tools


!---------------------------------------------------------------

     subroutine con_density_TOOLS

!---------------------------------------------------------------
 

!     Subroutine to calculate the density with SPH continuiity approach.

!     ntotal : Number of particles                                  [in]
!     mass   : Particle masses                                      [in]
!     niac   : Number of interaction pairs                          [in]
!     pair_i : List of first partner of interaction pair            [in]
!     pair_j : List of second partner of interaction pair           [in]
!     dwdx   : derivation of Kernel for all interaction pairs       [in]
!     vx     : Velocities of all particles                          [in]
!     itype   : type of particles                                   [in]
!     x      : Coordinates of all particles                         [in]
!     rho    : Density                                              [in]
!     drhodt : Density change rate of each particle                [out]   

implicit none
  

integer(ink) i ,j ,idimn    
real   (irk) vcc  
real   (irk), allocatable:: dvx(:) 
  
real   (irk)  selfdens, r 
real   (irk), allocatable:: drho_tmp(:), hv(:), wi(:)

if (.NOT.allocated (dvx) ) allocate ( dvx(ndimn)  )

drho = 0.     

if (pa_sph.eq.4) then                   !      Valid only for 1D approx of cons gradient
    current=> last
    Do while (associated(current))
       i = current%pair_i 
       j = current%pair_j 
       if (current%Pint_type.eq.0) then
          do idimn = 1,ndimn
             dvx(idimn) = vx(idimn,J) - vx(idimn,I) 
          enddo        
          vcc = dvx(1)*current%dwdx(1)        
          do idimn = 2,ndimn
             vcc = vcc + dvx(idimn)*current%dwdx(idimn)
          enddo     
          drho(i) = drho(i) + mass(j)*vcc/rho(J)
          drho(j) = drho(j) + mass(i)*vcc/rho(I)
       endif
       current=>current%next
    enddo
    do i = 1,ntotal
       drho(i) = drho(i)*rho(i)/alphaJB(i)
    enddo
else                                    !      Valid for standar cont density and JB alpha (pa_sph=3)
    current=> last
    Do while (associated(current))
       i = current%pair_i 
       j = current%pair_j 
       if (current%Pint_type.eq.0) then
          do idimn = 1,ndimn
             dvx(idimn) = vx(idimn,i) - vx(idimn,j) 
          enddo        
          vcc = dvx(1)*current%dwdx(1)        
          do idimn = 2,ndimn
             vcc = vcc + dvx(idimn)*current%dwdx(idimn)
          enddo     
          drho(i) = drho(i) + mass(j)*vcc
          drho(j) = drho(j) + mass(i)*vcc 
       endif
       current=>current%next
    enddo

    if (pa_sph.eq.3) then                       !      This is the alpha correction factor from JB
        do i = 1,ntotal
           drho(i) = drho(i)/alphaJB(i)
        enddo
    endif

endif

!      ------ normalize drhoI = ( sumJ drhoJ*Wij*massJ/rhoJ ) / ( sumJ Wij*massJ/rhoJ )

if (nor_density) then

    if (.NOT.allocated (hv) ) then 
       allocate ( hv(ndimn), wi(ntotal), drho_tmp(ntotal)  )
    endif
    drho_tmp = 0.0
    hv       = 0.e0
    r        = 0.
!                                       ------  First we get sumJ Wij*massJ/rhoJ
    do i = 1,npoin
   if (if_out_domain(i).eq.1) cycle
     ! call kernel       (r,hv,hsml(i),selfdens,hv)
       call kernel_TOOLS (r,hv,hsml(i),selfdens,hv)   ! CHGD MP april 2020
       wi(i) = selfdens*mass(i)/rho(i)
    enddo

    current=> last
    Do while (associated(current))
       if (current%Pint_type.eq.0) then
           i     = current%pair_i 
           j     = current%pair_j 
           wi(i) = wi(i) + mass(j)/rho(j)*current%w 
           wi(j) = wi(j) + mass(i)/rho(i)*current%w
       endif   
       current=>current%next
    enddo
    
!                                       ------  Then, sumJ drhoJ*Wij*massJ/rhoJ
    do i = 1,ntotal
   if (if_out_domain(i).eq.1) cycle
    !  call kernel(r,hv,hsml(i),selfdens,hv)
       call kernel_TOOLS (r,hv,hsml(i),selfdens,hv)   ! CHGD MP april 2020
       drho_tmp (i) = selfdens * mass(i)*drho(i)/rho(i)
    enddo

    current=> last
    Do while (associated(current))
       if (current%Pint_type.eq.0) then
           i = current%pair_i 
           j = current%pair_j 
           drho_tmp(i) = drho_tmp(i) + mass(j)*current%w * drho(j)/rho(j) 
           drho_tmp(j) = drho_tmp(j) + mass(i)*current%w * drho(i)/rho(i)
       endif      
       current=>current%next
     enddo
 
     do i = 1, npoin
        drho(i) = drho_tmp(i)/wi(i)
     enddo

     deallocate ( drho_tmp, wi, hv )
     
endif   

deallocate ( dvx  ) 

 END SUBROUTINE con_density_TOOLS


!---------------------------------------------------------------

     subroutine get_wir_density_TOOLS

!---------------------------------------------------------------
 
implicit none

integer(ink) i ,j ,idimn 

real   (irk) vcc, selfdens, r  
real   (irk), allocatable:: dvx(:), hv(:) 

allocate ( dvx(ndimn), hv(ndimn) )

drho = 0.
hv   = 0.
r    = 0.

do i = 1,ntotal                                   !  This is for SOIL summ density  (first step)
   if (if_out_domain(i).eq.1) cycle
   if (itype(i).eq.1) then
   !   call kernel(r,hv,hsml(i),selfdens,hv)       ! CHGD MP april 2020
       call kernel_TOOLS (r,hv,hsml(i),selfdens,hv)
       rho(i) = selfdens*mass(i) 
   endif
enddo

current=> last
Do while (associated(current))
   if (current%Pint_type.eq.0) then
      i = current%pair_i 
      j = current%pair_j 
      if ( itype(i).eq.1) then                  !  This is for SOIL summ density  (2nd step)
         rho(i) = rho(i) + mass(j)*current%w 
         rho(j) = rho(j) + mass(i)*current%w
      elseif ( itype(i).eq.2) then              !  This is for WATER con density   
         do idimn = 1,ndimn
            dvx(idimn) = vx(idimn,i) - vx(idimn,j) 
         enddo
         vcc = 0.0
         do idimn=1,ndimn
            vcc = vcc + dvx(idimn)*current%dwdx(idimn)
         enddo
         drho(i) = drho(i) + mass(j)*vcc
         drho(j) = drho(j) + mass(i)*vcc 
      endif 
   endif 
   current=>current%next
enddo

deallocate ( dvx, hv  )

END    subroutine get_wir_density_TOOLS
! -------------------------------------------------------------

     subroutine kernel_TOOLS (r,dx,hsml,w,dwdx)   

! -------------------------------------------------------------

!   Subroutine to calculate the smoothing kernel wij and its 
!              derivatives dwdxij.
!   if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!          = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!          = 3, Quintic kernel (Morris 1997)
!           =4, Wendland kernel  ! Saeidmt 24Oct2020

!   r    : Distance between particles i and j                     [in]
!   dx   : x-, y- and z-distance between i and j                  [in]  
!   hsml : Smoothing length                                       [in]
!   w    : Kernel for all interaction pairs                      [out]
!   dwdx : Derivative of kernel with respect to x, y and z       [out]

implicit none
      
real   (irk) r, dx(ndimn), hsml, w, dwdx(ndimn)    
real   (irk) q, factor
integer(ink) idimn  

q = r/hsml 
w       = 0. 
dwdx    = 0.

if (skf.eq.1) then 
    if (ndimn.eq.1) then
       factor = 1.e0/hsml
    elseif (ndimn.eq.2) then
       factor = 15.e0/(7.e0*pi*hsml*hsml)
    elseif (ndimn.eq.3) then
       factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
    else
       print *,' >>> Error <<< : Wrong dimension: Dim =',ndimn
       stop
    endif                                           
    if (q.ge.0.and.q.le.1.e0) then          
       w = factor * (2./3. - q*q + q**3 / 2.)
       do idimn = 1, ndimn
          dwdx(idimn) = factor * (-2.+3./2.*q)/hsml**2 * dx(idimn)       
       enddo   
    else if (q.gt.1.e0.and.q.le.2) then          
       w = factor * 1.e0/6.e0 * (2.-q)**3 
       do idimn = 1, ndimn
          dwdx(idimn) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx(idimn)/r)        
       enddo              
    else
       w    = 0.0
       dwdx = 0.0           
    endif  
    
else if (skf.eq.2) then      
    factor = 1.e0 / (hsml**ndimn * pi**(ndimn/2.))      
    if(q.ge.0.and.q.le.3) then
       w = factor * exp(-q*q)
       do idimn = 1, ndimn
          dwdx(idimn) = w * ( -2.* dx(idimn)/hsml/hsml)
       enddo 
    else
    w    = 0.
    dwdx = 0.0     
endif 

else if (skf.eq.3) then      
    if (ndimn.eq.1) then
        factor = 1.e0 / (120.e0*hsml)
    elseif (ndimn.eq.2) then
        factor = 7.e0 / (478.e0*pi*hsml*hsml)
    elseif (ndimn.eq.3) then
        factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
    else
        print *,' >>> Error <<< : Wrong dimension: Dim =',ndimn
        stop
    endif              
    if(q.ge.0.and.q.le.1) then
       w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
       do idimn= 1, ndimn
         dwdx(idimn) = factor * ( (-120 + 120*q - 50*q**2) / hsml**2 * dx(idimn) ) 
       enddo 
    else if(q.gt.1.and.q.le.2) then
       w = factor * ( (3-q)**5 - 6*(2-q)**5 )
       do idimn= 1, ndimn
          dwdx(idimn) = factor * (-5*(3-q)**4 + 30*(2-q)**4) / hsml * (dx(idimn)/r) 
       enddo 
    else if(q.gt.2.and.q.le.3) then
       w = factor * (3-q)**5 
       do idimn= 1, ndimn
          dwdx(idimn) = factor * (-5*(3-q)**4) / hsml * (dx(idimn)/r) 
       enddo 
    else   
       w    = 0.
       dwdx = 0.0
    endif 
else if (skf.eq.4) then           ! Saeidmt 24Oct2020 Wendland kernel
    if (ndimn.eq.1) then
        factor = 5.e0 / (8.e0*hsml)
    elseif (ndimn.eq.2) then
        factor = 7.e0 / (4.e0*pi*hsml*hsml)
    elseif (ndimn.eq.3) then
        factor = 21.e0 / (16.e0*pi*hsml*hsml*hsml)
    else
        print *,' >>> Error <<< : Wrong dimension: Dim =',ndimn
        stop
    endif
    if(q.ge.0.and.q.le.2) then
       w = factor * ( ((1-q/2)**4)*(2*q+1) )    ! MP25Oct2020 was 2q :-)
       do idimn= 1, ndimn
!         dwdx(idimn) = factor * ( (-5*q*(1-(q/2))**3 ) ) * ( (1/hsml)*(dx(idimn)/r) )
          dwdx(idimn) = factor * ( (-5*(1-(q/2))**3 ) ) * ( (1/hsml/hsml)*(dx(idimn)) )  ! MP 25Oct2020
       enddo 
    else   
       w    = 0.
       dwdx = 0.0
    endif          
endif 
        
END SUBROUTINE kernel_TOOLS 

END MODULE My_SPH_Tools_2019