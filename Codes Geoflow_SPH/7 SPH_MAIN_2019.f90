       
      
!      -------------------------------------------------------
!
!       
              MODULE SPH_MAIN_2019
!      
!
!      -------------------------------------------------------

!     Semiimplicit search for string MP s_implicit 17.03.18
 
!     ---Changes in MODULE SPH_MAIN_2014 MArch 21 2014---
! 
! (a) In both RK4 and RK4_WIR
!
! (a.1) Moved average_velocity before specific SW operations
! 
! (a.2) Within the RK loop call ExtForcesSW (0)
! 
! (a.3) After RK loop, call external forces (1)
!       and do not call erosion anymore
!
!      SPH90_Main  Octobre 2008
!      from Tom  coupled 19 01 2011

!     *****  saved from 9th April version. Change: deactiv call to gets at w (done in INTforces)

! --------------------------------------------              

USE SPH_FEM_variable_types_2019  	!  Variable precision defined here.
USE SPH_FEM_Driver_time_vars_2019   !  Time variables common to sph and geoflow
USE SPH_Main_Vars_2019        		!  Main SPH variables are defined here
USE SPH_SW_2019				        !  For SPH SW problems
USE GCsph_solvers_2019              !  CG for semi.implicit SW from TB ** MP 2017/03/25
USE My_SPH_Tools_2019               !  migrate to find routines in module tools 

! USE SPH_SW_Vars_2014          	!  for SW  (pwp, curvatures, etc)

! USE SPH_TGSolid_2008            !  For TGSolid problems
! USE SPH_NS_2008                 !  for NS  
! USE SPH_DF_2008                 !  for DF 

PRIVATE       !      ------  No variable will be accesible from outside the module
              !              unless declared public (none indeed)

public:: Init_SPH             !  called from Driver
public:: time_integration_sph
public:: adaptive_dt_sph
public:: Out_plot_sph
public:: Out_Plot_Trigger_SPH 
public:: Out_print_sph
public:: Out_save_sph
public:: clean_up_sph 
public:: Update_Const_SPH
public:: Stop_case_SPH
public:: Close_files_SPH
public:: Get_Prob_Propag_SPH     ! MP June 2022
public:: deallocate_list_SPH    !  MP Aug 2022 deallocates linked list

public:: Manage_pixel_Hazard_SPH     ! MP J24 uly 2023 

!  Added in RK4 and RKWir a call to a filtering routine which avoids average veolcity over nodes behind gates
!  Routine is get_npoin0_SW
!  Drum 2013
 


CONTAINS


!-------------------------------------------------------------------

       SUBROUTINE Init_sph 

!-------------------------------------------------------------------

 
!      ------  Gets  data for SW, allocates and  initializes

implicit none

integer(ink)  ipoin 

real   (irk)  total_mass  

! write(*,*) ' ------------------------------------------ '
! write(*,*) '                                            '
! write(*,*) '    NOTICES '
! write(*,*) ' ------------------------------------------ '
! write(*,*) '     '
! write(*,*) ' ** central integration routine removed  '
! write(*,*) '                                            '
! write(*,*) ' ** Added in RK4 and RKWir a call to a filtering routine...'
! write(*,*) '    which avoids average veolcity over nodes behind gates...'
! write(*,*) '    Routine is get_npoin0_SW'
! write(*,*) '   '
! write(*,*) ' ** Can use icadaptdt= 0, 1, 0r 2 v+c'
! write(*,*) '    labs dimensioned to 1000'
! write(*,*) '    routine check abs just to check, useless  '
! write(*,*) '   '
! write(*,*) ' ** Call to get Z which is in TOPO module '
! write(*,*) '    ramp ONLY FOR AKNES FILE to avoid BCs in TOPO'
! write(*,*) '    smoothing factor for slope in topo '
! write(*,*) '    modified int forces to put in p(i) h*Z,i'
! write(*,*) '   '
! write(*,*) '  **  NORMALIZED RHO also normalizes drho in con dens '
! write(*,*) '   '
! write(*,*) '  ** RK10  changed for prescribed velocity - study drho '
! write(*,*) '      we have introduced a new routine which substitutes con and sum density'
! write(*,*) '      water is con and soil sum'
! write(*,*) '   '
! write(*,*) '   '
! pause 

!    ------ Set Pw and Theta stress integration indicators

if_coupled_pw = 0
ic_stress_integ = 0

!    ------- SPH_Problem has been read in driver and stored in driver vars

If (SPH_ProblemType.eq.11) then
    if_coupled_pw = 1
    SPH_ProblemType = 10
Elseif (SPH_ProblemType.eq.20) then
    ic_stress_integ = 1
    SPH_ProblemType = 10    
Elseif (SPH_ProblemType.eq.21) then
    ic_stress_integ = 1
    SPH_ProblemType = 10   
    if_coupled_pw = 1
Endif    
    
If ( SPH_ProblemType.eq.1) then        
     
     Call Get_Init_SPH_SW               !    ------  Get Coorg, Topol and Xming,Xmaxg,..npoig, deltas... 

!* ElseIf (SPH_ProblemType.eq.2) then
!*      call Input_NS 
!* ElseIf (SPH_ProblemType.eq.3) then
!*      call Topo_main                      
!*      call Input_DF 
!* ElseIf (SPH_ProblemType.eq.10) then        
!*      Call Get_Init_SPH_TGSolid         !    ------  

else
    write(*,*) ' wrong type of problem= ', SPH_ProblemType
    pause
Endif

!* If ( if_coupled_pw.eq.1) then
!*      Call Get_Init_SPHTG_Solid_FrStep_Pw
!* endif

If (SPH_t_Integ_Alg.eq.4.or.SPH_t_Integ_Alg.eq.5   &
                        .or.SPH_t_Integ_Alg.eq.100 &
                        .or.SPH_t_Integ_Alg.eq.110 &
                        .or.SPH_t_Integ_Alg.eq.10  &
                        .or.SPH_t_Integ_Alg.eq.6  &
                        .or.SPH_t_Integ_Alg.eq.7  &
                        .or.SPH_t_Integ_Alg.eq.8) then       ! MP s_implicit 17.03.18
   if (.NOT.allocated (Dx_RK) ) then 
      allocate ( Dx_RK(ndimn,ntotal),Dvx_RK(ndimn,ntotal),Du_RK(ntotal),Drho_RK(ntotal) )
   endif
   Dx_RK = 0.0; Dvx_RK = 0.0 ;Du_RK = 0.0; Drho_RK = 0.0
Endif

If (pa_sph.eq.3.or.pa_sph.eq.4) then
   if (.NOT.allocated (AlphaJB) ) allocate ( AlphaJB(ntotal) ) !  ------ Alpha correction factor from JB
endif

if (SPH_ProblemType.eq.1.AND.ic_TRIGGER.Eq.0) then   !     MP March 2022 ------  Plot initial mesh ( post.msh file )
    total_mass = 0.0
    do ipoin = 1,npoin 
       total_mass = total_mass + mass(ipoin)
    enddo
    write(*,*)'.........................................'
    write(*,*)'   '
    write(*,*)' Initial volume (SW) or mass =', total_mass
    write(*,*)'   '  
    write(*,*)'.........................................'
endif 

x00 = x0; dvx = 0.0; drho = 0.0; du  = 0.0; av = 0.0    !      ------  Initialize 

if (SPH_ProblemType.eq.1) then                          !      ------  Plot initial mesh ( post.msh file )
       
    if (ic_trigger.eq.0) then                           ! MP19th March 2022
       call OutputMesh_SW
       if (ic_Cases_mat.EQ.0 ) then                        !      ------  FOR CASES, AVOID WRITTING THE FIRST
          call OutputRes_SW 
       endif
       call output_points_SW
    else
       if (nproblems_index.eq.1)then
          call OutputMesh_SW
       endif
    endif
    
!elseif (SPH_ProblemType.eq.10) then              ! TGSolid    
!       
!    call OutputMesh_TGSolid
!    call OutputRes_TGSolid 
!    call output_points_TGSolid
!    call OutputMesh_GMESH
!    call OutputRes_GMESH
       
else
    write(*,*) ' wrong option for type of problem = ', SPH_ProblemType
    pause
endif

END SUBROUTINE Init_sph    



!-------------------------------------------------------------------

       SUBROUTINE time_integration_sph

!-------------------------------------------------------------------

 
!      ------  Called from SPH_GFL_2008: perform time integration for sph
!              at time time_sph and itime_step_sph

implicit none

if (SPH_ProblemType.eq.1.AND.Ic_TRIGGER.EQ.0 ) then  !      MP March 2022
   call Inject_histogram_SW (time_sph)               !      ------  Inject particles if necessary   
   call Check_Out_Domain                             !      ------  Set up If_Out_Domain for 
endif                                                !              all real particles
      
!      ------       Time integration algorithms

If     (SPH_t_Integ_Alg.eq.0) then
    call RK45TRI  (itimestep_sph)
!   call Central(itimestep_sph )
ElseIf (SPH_t_Integ_Alg.eq.1) then
!   call Central1(itimestep_sph )
ElseIf (SPH_t_Integ_Alg.eq.4) then
   call RK45 (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK )   ! changed MP 20 Jan 2019 from RK4
ElseIf (SPH_t_Integ_Alg.eq.5) then
   call RK45 (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK )   ! changed MP 20 Jan 2019 from RK4_WIR
!ElseIf (SPH_t_Integ_Alg.eq.6) then                             ! **MP s_implicit 17.03.18                   
!     call FS_semi_implicit_heat (itimestep_sph)  
!ElseIf (SPH_t_Integ_Alg.eq.60) then                     
!     call FS_semi_implicit (itimestep_sph)                     ! **MP s_implicit 17.03.18                   
ElseIf (SPH_t_Integ_Alg.eq.7) then                             ! **MP s_implicit fem 17.09.14                   
     call FS_semi_implicit_fem (itimestep_sph)                   
ElseIf (SPH_t_Integ_Alg.eq.8) then                             ! **MP s_implicit fem 17.09.14                   
      call FS_semi_implicit_SFM (itimestep_sph)    
ElseIf (SPH_t_Integ_Alg.eq.10) then
      call RK10_slow (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) !   *** SLOW
ElseIf (SPH_t_Integ_Alg.eq.100.OR.SPH_t_Integ_Alg.eq.110) then      ! new call MP August 2022
     call RK00   (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) !   stochastic block
ElseIf (SPH_t_Integ_Alg.eq.200) then
    call Trigger00 
Else
   write(*,*) ' Check SPH_t_Integ_Alg variable. Should be 0, 1 or 4 '
   pause
   stop
Endif

!      ------  Check If_Out_Domain for all real particles

if (ic_Trigger.NE.1) then   ! added 25th April 2022 MP  
   call Check_Out_Domain     
   if (time_SPH.ge.T_change_to_W) then     ! ADDED 2007 07 19
       call change_to_w_SW                 ! small modif MP 22 July 2022
       T_change_to_W = 1.e10
   endif
endif
     
!      ------  Check for BCs in NS

if  (SPH_ProblemType.eq.2) then
!     call BCs_NS
endif 


END SUBROUTINE time_integration_sph



!-------------------------------------------------------------------

       SUBROUTINE adaptive_dt_sph

!-------------------------------------------------------------------

 
!      ------  Obtains dt_sph, the minimum dt for sph stability

implicit none

if     (SPH_ProblemType.eq.1) then
        call Adaptive_dt_sph_SW
elseif (SPH_ProblemType.eq.10) then
!        call Adaptive_dt_sph_TGSolid
endif  



END SUBROUTINE adaptive_dt_sph

!-------------------------------------------------------------------

       SUBROUTINE Out_plot_sph

!-------------------------------------------------------------------

implicit none

if (SPH_ProblemType.eq.1) then
    call OutputRes_SW  
elseif (SPH_ProblemType.eq.2) then
!    call OutputRes_NS
elseif (SPH_ProblemType.eq.3) then
!    call OutputRes_DF
elseif (SPH_ProblemType.eq.10) then     ! TGSolid
!     call OutputRes_TGSolid
!     call OutputRes_GMESH
else
   write(*,*) ' wrong type of problem= ', SPH_ProblemType
   pause
endif 


END SUBROUTINE Out_plot_sph


!-------------------------------------------------------------------

       Subroutine Out_Plot_Trigger_SPH  
       
!-------------------------------------------------------------------

    implicit none

!      ------  Output for TRIGGER

 Call Out_Plot_Trigger_SW

End Subroutine Out_Plot_Trigger_SPH  


!-------------------------------------------------------------------

       SUBROUTINE Out_save_sph

!-------------------------------------------------------------------

implicit none

if (SPH_ProblemType.eq.1) then
    call OutputSave_SW
endif 


END SUBROUTINE Out_save_sph

!-------------------------------------------------------------------

       SUBROUTINE Out_print_sph

!-------------------------------------------------------------------

implicit none

integer(ink) ipoin
real   (irk) total_mass, vtemp, vtempx, vtempy

total_mass=0.0
vtemp = 0.0
do ipoin = 1, npoin
   total_mass = total_mass + mass(ipoin)
   if (If_Out_Domain(ipoin).eq.1) cycle
   vtempx = vx(1,ipoin)
   vtempy = 0.0
   if (ndimn.eq.2) vtempy = vx(2,ipoin)
   vtemp  = vtemp + (vtempx*vtempx + vtempy*vtempy)**0.5
enddo
vtemp = vtemp/npoin

if (Ic_TRIGGER.EQ.0 ) then
   write(*,*)'   '
   write(*,*)'  sph output '
   write(*,*)' time step_sph is ', itimestep_sph, ' time_sph =', time_sph 
   write(*,*)'  dt_sph=', dt_sph , ' MASS = ', total_mass 
   write(*,*) '   vtemp =  ', vtemp      !  MP July 2022
   write(*,*)'.........................................'
endif
   
if (SPH_ProblemType.eq.1) then
   if (Ic_cases_win.eq.0) then
      call Moni_SW
      call output_points_SW
   else
      call Output_cases_win_SW
   endif
endif
         
101 format(1x,3(2x,a12))     
100 format(1x,4(2x,e13.6))  


END SUBROUTINE Out_print_sph     


!---------------------------------------------------------------------- 

      subroutine clean_up_sph
      
!----------------------------------------------------------------------  

!     Subroutine for deallocating variables and closing files

implicit none   

deallocate (     x0 ,       x ,      vx0,    vx                     )
deallocate (     dx ,     dvx ,       av                            )
deallocate ( indvxdt,  exdvxdt,  ardvxdt                            )
deallocate (    mass,      rho,        p,    u, hsml, c,     s, e   )
deallocate (      u0,     rho0,       du, drho, ds  , t, tdsdt      ) 
deallocate (  ahdudt,   avdudt,      eta                            ) 
deallocate (   itype, countiac                                      )

If (SPH_t_Integ_Alg.ge.4) then
   deallocate ( Dx_RK,  Dvx_RK, Du_RK,Drho_RK )
Endif 

if (SPH_ProblemType.eq.1) call Clean_SW  ! for variables dimensioned in SPH

close(dat_file)
close(chk_file)
close(gid_res)  
close(gid_msh)
close(gmesh)
close(res_file1)
close(res_file2)
close(res_file3)

END SUBROUTINE clean_up_sph 


!----------------------------------------------------------------------
!
       subroutine TG2
!
!----------------------------------------------------------------------
!
!*   RECOVER FROM TB AFTER CHECKING SW

END subroutine  TG2 


!----------------------------------------------------------------------

subroutine First_step_TG

!----------------------------------------------------------------------

!*   RECOVER FROM TB AFTER CHECKING SW

END subroutine  First_step_TG


!----------------------------------------------------------------------

subroutine Second_step_TG

!----------------------------------------------------------------------

!*   RECOVER FROM TB AFTER CHECKING SW

END Subroutine Second_step_TG

!----------------------------------------------------------------------

subroutine MLS_for_Boundaries(na,BNneighb_vct,BNneighb,BN_count,nmlspts,nfact,flux,diveF)

!----------------------------------------------------------------------


!*   RECOVER FROM TB AFTER CHECKING SW


END SUbroutine  MLS_for_Boundaries

!----------------------------------------------------------------------

subroutine MLS_for_fn(na,BNneighb_vct,BNneighb,BN_count,nmlspts,nfact,fi,unkno)

!----------------------------------------------------------------------



!*   RECOVER FROM TB AFTER CHECKING SW



END SUbroutine  MLS_for_fn


!-------------------------------------------------------------------------

        subroutine get_sources_RK4_Solid  (n1_RK)

!-------------------------------------------------------------------------



!*   RECOVER FROM TB AFTER CHECKING SW



END SUBROUTINE get_sources_RK4_Solid  

!-------------------------------------------------------------------------

        subroutine get_sources_RK3_Solid  (n1_RK)

!-------------------------------------------------------------------------



!*   RECOVER FROM TB AFTER CHECKING SW

END SUBROUTINE get_sources_RK3_Solid  


!---------------------------------------------------------------

     subroutine  art_visc

!---------------------------------------------------------------


!*   RECOVER FROM TB AFTER CHECKING SW


 
END SUBROUTINE art_visc


!----------------------------------------------------------------------

      subroutine Trigger00     

!----------------------------------------------------------------------
 
implicit none          !  MP March 2022 Will be called when ic_t_integ_Alg = 200 NEW NEW 

 
if (ic_Mcarlo .eq.6) then
    call Trigger_FOSM_SW  ! new FOSM computations  24 april 2022
else   
    call Trigger_FoS_SW  
endif

END subroutine  Trigger00 


!----------------------------------------------------------------------

      subroutine RK00 (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) 

!----------------------------------------------------------------------

!   Random Walk  

implicit none

integer (ink) icrnd
integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
real   (irk) Dx_RK(ndimn,npoin), Dvx_RK(ndimn,npoin), Du_RK(npoin), Drho_RK(npoin) ! These are local!
										                                           ! defined as allocat 										                                           ! in MAIN Vars
real   (irk) dummy1, vmod, xrnd, xrnd0 , pi2
real   (irk) dxt(ndimn), dvxt(ndimn), dxtmod, dvxtmod, nt(ndimn), nn(ndimn), dxmc, dymc, dvxmc, dvymc 
integer(ink) t_clock
integer(ink), dimension(1):: seed

data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

!icrnd = 0                     ! MP 7 Aug 2022
!if (rndm_fact.GT.1.e-3) then
!   icrnd = 1
! endif
!      ------  Initialize ( done also in subroutines, strictly unnecessary )

if (ndimn.ne.2.AND.SPH_t_Integ_alg.eq.100) then
   write (*,*) ' random only for 2D'
   read  (*,*) dummy1
   PAUSE
endif
 
exdvxdt    = 0.     
vmod = 0.0; nt=0.0; dxt=0.0; xrnd=0.0

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0 

t0 = time_sph; vx0 = vx; x0 = x

DO i = 1,4                                         ! ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph
   do ipoin = 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      x (:,ipoin)  = x0(:,ipoin) + f1rk(i)*(dt_sph)*vx(:,ipoin)
      vx(:,ipoin) = vx0(:,ipoin) + f1rk(i)*(dt_sph)*(exdvxdt(:,ipoin) )         !    +   indvxdt(:,ipoin)+ardvxdt(:,ipoin))
   enddo

   Call Check_Out_Domain 
   call ExtForces_SW (0)                !       External Forces SW (slope, friction,virt.parts)
   
   DO ipoin = 1, npoin      
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      Dvx_RK(:,ipoin) = Dvx_RK  (:,ipoin) + f2rk(i)*(exdvxdt(:,ipoin) )        !    +indvxdt(:,ipoin)+ardvxdt(:,ipoin)) MP Aug 2022
      Dx_RK (:,ipoin) = Dx_RK   (:,ipoin) + f2rk(i)*vx(:,ipoin)
   ENDDO
   ipoin = 1 !  just for debug
   
ENDDO

pi2 = 2*atan(1.0) 
dxt = 0.0;  dvxt = 0.0; nt = 0. ; nn = 0.

!if (SPH_t_Integ_Alg.EQ.100) then
if (icrnd.eq.1) then
   call system_clock (t_clock)   ! CHGD MP 4 April 2022
   seed = t_clock
   call random_seed (PUT=seed)
endif 

DO ipoin =1,npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   dxt (:) =  dt_sph/6. * Dx_RK (:,ipoin) 
   dvxt(:) =  dt_sph/6. * Dvx_RK (:,ipoin)
   dxtmod  =  dxt(1)*dxt(1) 
   dvxtmod =  dvxt(1)*dvxt(1)
   if (ndimn.eq.2) then
      dxtmod  = dxtmod  + dxt(2)*dxt(2) 
      dvxtmod = dvxtmod + dvxt(2)*dvxt(2) 
   endif 
   dxtmod  = dxtmod **0.5
   dvxtmod = dvxtmod**0.5
   if (dxtmod.ge.1.e-6) then
      nt = dxt/dxtmod
      if (ndimn.eq.2) then
         nn(1) = -nt(2)
         nn(2) =  nt(1)
         ! if (SPH_t_Integ_Alg.EQ.100) call random_number (xrnd0)  ! MP August 2022
         if (icrnd.EQ.1) call random_number (xrnd0)  ! MP August 2022
         xrnd     = (-1.+ 2*xrnd0)*pi2
         dxt      =  dxtmod * ( nt + rndm_fact *xrnd0*nn )
      else
         dxt      =  dxtmod 
      endif
!     dxt      = dxtmod* ( nt*cos(xrnd) + 0.5* nn*sin(xrnd) )   !   ellipse  
   endif
!   if (dvxtmod.ge.1.e-6) then   ! deactivated
!      nt = dvxt/dvxtmod
!      nn(1) = -nt(2)
!      nn(2) =  nt(1)
 !     call system_clock (t_clock)
 !     seed = t_clock
 !     call random_seed (PUT=seed)
 !     call random_number (xrnd0)
 !     xrnd   = (-1.+ 2*xrnd0)*pi2
 !     dvxt    = dvxtmod* ( nt*cos(xrnd) + nn*sin(xrnd) )
 !  endif

   x  (:,ipoin) =  x0 (:,ipoin) + dxt (:)      
   vx (:,ipoin) = vx0 (:,ipoin) + dvxt(:)
ENDDO

 Call Check_Out_Domain 

call ExtForces_SW (1)                !       External Forces SW (slope, friction,virt.parts)

1991 format (12(f10.4,2x))

!      ------   Clears linked list containing interaction data

!   call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  RK00 


!----------------------------------------------------------------------

      subroutine RK45TRI  (itimestep_sph ) 

!----------------------------------------------------------------------

 

implicit none

integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0  

if (nnps.eq.1) then 
    call direct_find(itimestep_sph)
else if (nnps.eq.2) then
    call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
elseif (nnps.eq.3) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.4) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.100) then                          ! CHGD MP april 2020
    call grid_find_NEW_Tools_TEST  (itimestep_sph) ! CHGD MP Aug 2022 test Increm Mb0
endif 

call Pint_Update_SW_TOOLS   !   Updates interactions to 0 (default) 

if (itimestep_sph.eq.1) then
   DO ipoin =1,npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE      
      vx (:,ipoin) =   x (:,ipoin) 
   ENDDO         
else
   DO ipoin =1,npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE      
      x  (:,ipoin) =  x  (:,ipoin) + dt_sph* vx (:,ipoin) 
   ENDDO        
endif 

call sum_density_TOOLS        ! CHGD MP April 2020 

if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML
                                        
call add_w_to_s_SW_TRI                  !    MP for mixing water to soil 21 July 2023    
 

END subroutine  RK45TRI 
!----------------------------------------------------------------------

      subroutine RK45 (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) 

!----------------------------------------------------------------------

!   Subroutine to integrate in time. 
!   First, (dt_sph/2) RK4 for sources  
!   Then   convective terms (pressure)
!   Second (dt_sph/2) RK4 for sources

!   *** integrating 4 and 5 1st June 2018

implicit none

integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
                                                                                   ! These are local!
real   (irk) Dx_RK(ndimn,npoin), Dvx_RK(ndimn,npoin), Du_RK(npoin), Drho_RK(npoin) ! defined as allocat										    
										   ! in MAIN Vars
data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

real   (irk) zz, dhh, xi	!      used to update pwp with dh										   
integer(ink) icUaux, iaux	!      icUaux=0 we do not use PWP + FD's
real   (irk) temp1
!      ------  Initialize ( done also in subroutines, strictly unnecessary )

avdudt     = 0.
ahdudt     = 0.
indvxdt    = 0.
ardvxdt    = 0.
exdvxdt    = 0.
av         = 0.     

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0
   
!      ------   Interaction parameters, calculating neighboring particles
!               and optimizing smoothing length

if (nnps.eq.1) then 
    call direct_find(itimestep_sph)
else if (nnps.eq.2) then
    call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
elseif (nnps.eq.3) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.4) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.100) then                          ! CHGD MP april 2020
    call grid_find_NEW_Tools_TEST  (itimestep_sph) ! CHGD MP Aug 2022 test Increm Mb0
endif 

if (SPH_ProblemType.eq.1) then  ! CHGD MP April 2020
 !  call Pint_Update_SW         !   Updates interactions to 0 (default) 
    call Pint_Update_SW_TOOLS   !   Updates interactions to 0 (default) 
endif    

t0 = time_sph; vx0 = vx; x0 = x; u0 = u
rho0 = rho				           ! **chnged for pwp

if (itimestep_sph.eq.1) then       !   we need normalized hs at the beginning 
   call Vn0_BCs_SW (0.0) 
endif

DO i = 1,4                                         !  ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph
   do ipoin = 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      x (:,ipoin)  = x0(:,ipoin) + f1rk(i)*(dt_sph)*vx(:,ipoin)
      vx(:,ipoin) = vx0(:,ipoin) + f1rk(i)*(dt_sph)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      
      if (.not.summation_density) then                  
          rho(ipoin) = rho0 (ipoin)   + f1rk(i)*(dt_sph)*drho(ipoin)
      endif
      if (SPH_t_Integ_Alg.eq.5.and.itype(ipoin).eq.2) then                  
          rho(ipoin) = rho0 (ipoin)   + f1rk(i)*(dt_sph)*drho(ipoin)
      endif
      
      u   (ipoin) = u0 (ipoin)   + f1rk(i)*(dt_sph)*du(ipoin)    !  modify only for icpwp=1
       
      if (SPH_ProblemType.eq.1.and.nuaux.eq.0) then             !   SW 
           if ( u(ipoin).gt.rho(ipoin) ) then 
               u(ipoin) = rho(ipoin)   ! ** CHK pwp
           elseif ( u(ipoin).lt.0.0 ) then      
               u(ipoin) = 0.0 
           endif
        endif 

    enddo
   
   Call Check_Out_Domain 
   
!      ------   Get alpha correction factor

   if (pa_sph.eq.3.or.pa_sph.eq.4) then                    !      Alpha correction factor JB                                                                
       call Get_AlphaJB
   endif   
   
   if (SPH_t_Integ_Alg.eq.4) then
       if (summation_density) then          
           call sum_density_TOOLS    ! CHGD MP April 2020         
       else             
           call con_density_TOOLS    ! CHGD MP April 2020       
       endif
   elseif (SPH_t_Integ_Alg.eq.5) then
       call get_wir_density_TOOLS    ! CHGD MP April 2020 
   endif    
 
 
!  test for divergence  and gradient 2020 october...

!call get_div_TOOLS  (x00,  1, drho)
!call get_grad_TOOLS (drho, 1, x00)

!STOP
   temp1 = 1               !  BUG MP & VZ 21 feb 2024
   ! chuan to detect NAN
   ! temp1 = isnan(rho(ipoin))
   if (temp1.lt.0)then
       write(*,*)'in RK45 rho(ipoin)=NAN when ipoin=',ipoin
       pause    
   endif   
   
   if (SPH_ProblemType.eq.1) then           !     SW
!          write(*,*) ' -----------  irkutta, time_sph ', i, time_sph
       call IntForces_SW                    !       Internal Forces SW.. Gets hw at s nodes
       
       if (i.eq.1.and.allocated(auxS1)) then!
!           call Get_s_at_w_DF_SW            !     called inside intforces
           auxS10 = auxS1                    !     transfer initial value of porosity to auxS10
       endif    
       
       call ExtForces_SW (0)                !       External Forces SW (slope, friction,virt.parts)
                                            !       March 2014 0 means standard, get dh/dt
   elseif (SPH_ProblemType.eq.2) then       !       NS
!       call IntForces_NS                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_NS                   !       
   elseif (SPH_ProblemType.eq.3) then       !       DF
!       call IntForces_DF                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_DF                   !
   else
       write(*,*) ' wrong type of problem= ', SPH_ProblemType
       pause
   endif
   
   DO ipoin = 1, npoin      
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      Dvx_RK(:,ipoin) = Dvx_RK  (:,ipoin) + f2rk(i)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      Dx_RK (:,ipoin) = Dx_RK   (:,ipoin) + f2rk(i)*vx(:,ipoin)
      if (SPH_t_Integ_Alg.eq.4) then
         if (.not.summation_density) then     !   ** wir
            Drho_RK (ipoin) = Drho_RK (ipoin)   + f2rk(i)*drho(ipoin) !       
         endif
      elseif (SPH_t_Integ_Alg.eq.5.and. itype(ipoin).eq.2) then
         Drho_RK (ipoin) = Drho_RK (ipoin)   + f2rk(i)*drho(ipoin) !
      endif
      Du_RK   (ipoin) = Du_RK   (ipoin)   + f2rk(i)*du  (ipoin)
   ENDDO
   
ENDDO

DO ipoin =1,npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   x  (:,ipoin) =  x0 (:,ipoin) + (dt_sph/6.)* Dx_RK (:,ipoin)       
   vx (:,ipoin) = vx0 (:,ipoin) + (dt_sph/6.)* Dvx_RK(:,ipoin)
   if (SPH_t_Integ_Alg.eq.4) then
      if (.not.summation_density) then
         rho(ipoin)  = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)        
      endif
   elseif (SPH_t_Integ_Alg.eq.5) then
      if (itype(ipoin).eq.2) then
         rho(ipoin) = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)
      endif      
   endif   
   u  (ipoin)   =  u0   (ipoin) + (dt_sph/6.)* Du_RK (ipoin)
ENDDO

Call Vn0_BCs_SW ( dt_sph)             ! change test MP

if (SPH_t_Integ_Alg.eq.4) then
    if (summation_density) then       ! ------  MP march 2018 DBKs   
        call sum_density_TOOLS        ! CHGD MP April 2020 
        Drho_RK = rho - rho0          !     CH MP march 2019
    endif
endif  
        
! Call Vn0_BCs_SW ( dt_sph)

Call Check_Out_Domain                                            

if (SPH_ProblemType.eq.1) then            

   call Get_s_at_w_DF_SW                ! *** get hw at soil pics and viceversa   
             
   call ExtForces_SW (1)                !  External Forces SW (slope, friction,virt.parts)
                                        !    March 2014 1 means erosion, once we know v  
                                        !    MP for Stefan  erosion march 2021 
!  call add_w_to_s_SW                   !    MP for mixing water to soil 21 July 2023    
                                        !  commented and moved to the end MP Jan 2023

   if (nUaux.GE.1) then
      
      if (allocated(auxS1)) then
         dauxS1 = auxS1 - auxS10
         call Get_Pwp_SW (Drho_RK, dauxS1)   !    compute pwp DFs and pwp
      else 
         call Get_Pwp_SW (Drho_RK)		!    compute pwp
      endif
      
      do ipoin = 1,npoin
	  u(ipoin) = UAux(1,ipoin)
      enddo
            
   endif
             
!    call ExtForces_SW (1)                !       External Forces SW (slope, friction,virt.parts)
                                         !       March 2014 1 means erosion, once we know v  
                                         !       Replaces old Erosion_SW 
    
            
    Call Vn0_BCs_SW ( dt_sph )           !  **  Moved here MP 30.12.2018 changhd dt/6 to dt
    

    Call Abs_BCs_SW                      !    remember this should go also to other time integration routines
!      
!      ------   
!
   if(nuaux.eq.0) then
     DO ipoin=1,npoin
        if ( if_Out_Domain(ipoin).eq.1) CYCLE
        if ( u(ipoin).gt.rho(ipoin) ) u(ipoin) = rho(ipoin)  ! ** CHK pwp
        if ( u(ipoin).lt.0.0 )        u(ipoin) = 0.0
     Enddo
   endif

   
endif      
   
if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML
                                        
if (average_velocity) then      !     March 2014 Moved here, more logical for erosion and BCs
   call av_vel		        !     Average v of particles to avoid penetration
   do ipoin = npoin0 + 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      vx (:,ipoin) = vx(:,ipoin) + av(:,ipoin)
   enddo   
endif

   call add_w_to_s_SW                   !    MP for mixing water to soil 21 July 2023    
                                        !  commented and moved to the end MP Jan 2023

!      ------   Clears linked list containing interaction data

!   call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  RK45 

!----------------------------------------------------------------------

      subroutine RK10_slow (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) 

!----------------------------------------------------------------------

!   Subroutine to integrate in time. 
!   First, (dt_sph/2) RK4 for sources  
!   Then   convective terms (pressure)
!   Second (dt_sph/2) RK4 for sources

!   *** integrating 4 and 5 1st June 2018

implicit none

integer(ink) itimestep_sph, i, ipoin , in_step
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
                                                                                   ! These are local!
real   (irk) Dx_RK(ndimn,npoin), Dvx_RK(ndimn,npoin), Du_RK(npoin), Drho_RK(npoin) ! defined as allocat										    
										   ! in MAIN Vars
data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

real   (irk) zz, dhh, xi	!      used to update pwp with dh										   
integer(ink) icUaux, iaux	!      icUaux=0 we do not use PWP + FD's
real   (irk) temp1
!      ------  Initialize ( done also in subroutines, strictly unnecessary )

avdudt     = 0.
ahdudt     = 0.
indvxdt    = 0.
ardvxdt    = 0.
exdvxdt    = 0.
av         = 0.     

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0
   
!      ------   Interaction parameters, calculating neighboring particles
!               and optimizing smoothing length

in_step = itimestep_sph
If (ic_search_step.eq.1) then 
   in_step = itimestep_sph/nsearch_step   
endif
                           !  ---------------------  !  MP Jan 2023  skip search when v is very small
! If ((in_step*nsearch_step.eq.itimestep_sph).OR.itimestep_sph.EQ.1)  then    

   if (nnps.eq.1) then 
       call direct_find(itimestep_sph)
   else if (nnps.eq.2) then
       call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
   elseif (nnps.eq.3) then                            ! CHGD MP april 2020
       call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
   elseif (nnps.eq.4) then                            ! CHGD MP april 2020
       call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
   elseif (nnps.eq.100) then                          ! CHGD MP april 2020
       call grid_find_NEW_Tools_TEST  (itimestep_sph) ! CHGD MP Aug 2022 test Increm Mb0
   endif 

   if (SPH_ProblemType.eq.1) then  ! CHGD MP April 2020 
       call Pint_Update_SW_TOOLS   !   Updates interactions to 0 (default) 
   endif    

!  endif

t0 = time_sph; vx0 = vx; x0 = x; u0 = u
rho0 = rho				           ! **chnged for pwp

if (itimestep_sph.eq.1) then       !   we need normalized hs at the beginning 
   call Vn0_BCs_SW (0.0) 
   call Get_Rain_SW                ! to get value at time 0 reactivated JMS 14th March 2024
endif

DO i = 1,4                                         !  ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph
   do ipoin = 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      x (:,ipoin)  = x0(:,ipoin) + f1rk(i)*(dt_sph)*vx(:,ipoin) 
      
      if (.not.summation_density) then                  
          rho(ipoin) = rho0 (ipoin)   + f1rk(i)*(dt_sph)*drho(ipoin)
      endif
!      if (SPH_t_Integ_Alg.eq.5.and.itype(ipoin).eq.2) then                  
!          rho(ipoin) = rho0 (ipoin)   + f1rk(i)*(dt_sph)*drho(ipoin)
!      endif
      
      u   (ipoin) = u0 (ipoin)   + f1rk(i)*(dt_sph)*du(ipoin)    !  modify only for icpwp=1
       
      if (SPH_ProblemType.eq.1.and.nuaux.eq.0) then             !   SW 
           if ( u(ipoin).gt.rho(ipoin) ) then 
               u(ipoin) = rho(ipoin)   ! ** CHK pwp
           elseif ( u(ipoin).lt.0.0 ) then      
               u(ipoin) = 0.0 
           endif
        endif 

    enddo
   
   Call Check_Out_Domain 
 
   ! if (SPH_t_Integ_Alg.eq.4) then
       if (summation_density) then          
           call sum_density_TOOLS    ! CHGD MP April 2020         rain12
       else             
           call con_density_TOOLS    ! CHGD MP April 2020       
       endif
   !elseif (SPH_t_Integ_Alg.eq.5) then
   !    call get_wir_density_TOOLS    ! CHGD MP April 2020 
   !endif    
 
 
!  test for divergence  and gradient 2020 october...

!call get_div_TOOLS  (x00,  1, drho)
!call get_grad_TOOLS (drho, 1, x00)

!STOP
   
   ! chuan to detect NAN
   ! temp1 = isnan(rho(ipoin))
   if (temp1.lt.0)then
       write(*,*)'in RK10 rho(ipoin)=NAN when ipoin=',ipoin
       pause    
   endif   
   
   if (SPH_ProblemType.eq.1) then           !     SW
!          write(*,*) ' -----------  irkutta, time_sph ', i, time_sph
       call IntForces_SW                    !       Internal Forces SW.. Gets hw at s nodes
       
       if (i.eq.1.and.allocated(auxS1)) then!
!           call Get_s_at_w_DF_SW            !     called inside intforces
           auxS10 = auxS1                    !     transfer initial value of porosity to auxS10
       endif    
                                            !       External Forces SW (slope, friction,virt.parts)
       call ExtForces_SW (0)                !       March 2014 0 means standard, get dh/dt
                                            !       Slow: obtains VX
   else
       write(*,*) ' wrong type of problem= ', SPH_ProblemType
       pause
   endif
   
   DO ipoin = 1, npoin      
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      Dx_RK (:,ipoin) = Dx_RK   (:,ipoin) + f2rk(i)*vx(:,ipoin)
 !     if (SPH_t_Integ_Alg.eq.4) then
         if (.not.summation_density) then     !   ** wir
            Drho_RK (ipoin) = Drho_RK (ipoin)   + f2rk(i)*drho(ipoin) !       
         endif
 !     elseif (SPH_t_Integ_Alg.eq.5.and. itype(ipoin).eq.2) then
         Drho_RK (ipoin) = Drho_RK (ipoin)   + f2rk(i)*drho(ipoin) !
 !     endif
      Du_RK   (ipoin) = Du_RK   (ipoin)   + f2rk(i)*du  (ipoin)
   ENDDO
   
ENDDO

DO ipoin =1,npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   x  (:,ipoin) =  x0 (:,ipoin) + (dt_sph/6.)* Dx_RK (:,ipoin) 
   if (.not.summation_density) then
      rho(ipoin)  = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)        
   endif 
!   if (SPH_t_Integ_Alg.eq.4) then
!      if (.not.summation_density) then
!         rho(ipoin)  = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)        
!      endif
 !  elseif (SPH_t_Integ_Alg.eq.5) then
 !     if (itype(ipoin).eq.2) then
 !        rho(ipoin) = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)
 !     endif      
 !  endif   
   u  (ipoin)   =  u0   (ipoin) + (dt_sph/6.)* Du_RK (ipoin)
ENDDO

Call Vn0_BCs_SW ( dt_sph)             ! change test MP

! if (SPH_t_Integ_Alg.eq.4) then
    if (summation_density) then       ! ------  MP march 2018 DBKs   
        call sum_density_TOOLS        ! CHGD MP April 2020 
        Drho_RK = rho - rho0          !     CH MP march 2019
    endif
! endif  
        
! Call Vn0_BCs_SW ( dt_sph)

Call Check_Out_Domain                                            

if (SPH_ProblemType.eq.1) then            

   call Get_s_at_w_DF_SW                ! *** get hw at soil pics and viceversa   
             
   call ExtForces_SW (1)                ! --- check slow case MP MARCH 2023
                                        !    March 2014 1 means erosion, once we know v  
                                        !    MP for Stefan erosion march 2021 
   call Get_Rain_SW        !  Updates phreatic level, PWp
   
   if (nUaux.GE.1) then
      
      if (allocated(auxS1)) then
         dauxS1 = auxS1 - auxS10
         call Get_Pwp_SW (Drho_RK, dauxS1)   !    compute pwp DFs and pwp
      else 
         call Get_Pwp_SW (Drho_RK)		!    compute pwp
      endif
      
      do ipoin = 1,npoin
	  u(ipoin) = UAux(1,ipoin)
      enddo
            
   endif
            
    Call Vn0_BCs_SW ( dt_sph )           !  **  Moved here MP 30.12.2018 changhd dt/6 to dt
    

    Call Abs_BCs_SW                      !    remember this should go also to other time integration routines
!      
!      ------   
!
   if(nuaux.eq.0) then
     DO ipoin=1,npoin
        if ( if_Out_Domain(ipoin).eq.1) CYCLE
        if ( u(ipoin).gt.rho(ipoin) ) u(ipoin) = rho(ipoin)  ! ** CHK pwp
        if ( u(ipoin).lt.0.0 )        u(ipoin) = 0.0
     Enddo
   endif

   
endif      
   
if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML
                                        
if (average_velocity) then      !     March 2014 Moved here, more logical for erosion and BCs
   call av_vel		        !     Average v of particles to avoid penetration
   do ipoin = npoin0 + 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      vx (:,ipoin) = vx(:,ipoin) + av(:,ipoin)
   enddo   
endif

!      ------   Clears linked list containing interaction data

!   call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  RK10_slow 


!----------------------------------------------------------------------

      subroutine RK4 (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) 

!----------------------------------------------------------------------

!   Subroutine to integrate in time. 
!   First, (dt_sph/2) RK4 for sources  
!   Then   convective terms (pressure)
!   Second (dt_sph/2) RK4 for sources

!   *** changed step by step 15 June *****

implicit none

integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
real   (irk) Dx_RK(ndimn,npoin), Dvx_RK(ndimn,npoin), Du_RK(npoin), Drho_RK(npoin) ! These are local!
										   ! defined as allocat 
										   ! in MAIN Vars
real   (irk) zz, dhh, xi	!      used to update pwp with dh										   
integer(ink) icUaux, iaux	!      icUaux=0 we do not use PWP + FD's
real   (irk) beta  ! **** **** **** **** **** ****

data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

!      ------  Initialize ( done also in subroutines, strictly unnecessary )

avdudt     = 0.
ahdudt     = 0.
indvxdt    = 0.
ardvxdt    = 0.
exdvxdt    = 0.
av         = 0.     

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0
   
!      ------   Interaction parameters, calculating neighboring particles
!               and optimizing smoothing length

if (nnps.eq.1) then 
    call direct_find(itimestep_sph)
else if (nnps.eq.2) then
    call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
elseif (nnps.eq.3) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.4) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
endif

if (SPH_ProblemType.eq.1) then 

  ! call Pint_Update_SW           !      Updates interactions to 0 (default) 
    call Pint_Update_SW_TOOLS     !      CHGD MP April 2028

!    Call Vn0_BCs_SW  ( 0.0)      !      SW necessary for normalizing density  *** deactivated 7 march
    
!    if (summation_density) then  !     -TMP delete after debug -----  If continuity is used, then get drho                      
!       call sum_density       
!    endif
    
endif    

t0 = time_sph; vx0 = vx; x0 = x; u0 = u
!if (.not.summation_density) rho0 = rho
rho0 = rho				           ! **chnged for pwp

if (itimestep_sph.eq.1) then       ! we need normalized hs at the beginning MP chgd 6-june-17
   call Vn0_BCs_SW (0.0) 
endif

DO i = 1,4                                         ! ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph
!   Call Vn0_BCs_SW ( f1rk(i)*dt_sph )             !      ***** CH added 6 March. *** deactivated 7 march
   do ipoin = 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      x (:,ipoin)  = x0(:,ipoin) + f1rk(i)*(dt_sph)*vx(:,ipoin)
      vx(:,ipoin) = vx0(:,ipoin) + f1rk(i)*(dt_sph)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      if (.not.summation_density) then                  
          rho(ipoin) = rho0 (ipoin)   + f1rk(i)*(dt_sph)*drho(ipoin)
      endif
      u   (ipoin) = u0 (ipoin)   + f1rk(i)*(dt_sph)*du(ipoin)                   !  modify only for icpwp=1
!      
!      ------   
!
        if (SPH_ProblemType.eq.1.and.nuaux.eq.0) then                           !   SW 
           if     ( u(ipoin).gt.rho(ipoin) ) then 
               u(ipoin) = rho(ipoin)   ! ** CHK pwp
           elseif ( u(ipoin).lt.0.0 ) then      
               u(ipoin) = 0.0 
           endif
        endif 

    enddo
   
   Call Check_Out_Domain 
   
!      ------   Get alpha correction factor

   if (pa_sph.eq.3.or.pa_sph.eq.4) then                    !      Alpha correction factor JB                                                                
       call Get_AlphaJB
   endif   

   if (summation_density) then       ! ------  If continuity is used, then get drho    
       call sum_density_TOOLS    ! CHGD MP April 2020         
     else             
       call con_density_TOOLS    ! CHGD MP April 2020       
   endif
   
   if (SPH_ProblemType.eq.1) then           !     SW
       call IntForces_SW                    !       Internal Forces SW.. Gets hw at s nodes
       
       if (i.eq.1.and.allocated(auxS1)) then!
           call Get_s_at_w_DF_SW            !     transfer initial value of porosity to auxS10
           auxS10 = auxS1
       endif    
       
       call ExtForces_SW (0)                !       External Forces SW (slope, friction,virt.parts)
                                            !       March 2014 0 means standard, get dh/dt
   elseif (SPH_ProblemType.eq.2) then       !     NS
!       call IntForces_NS                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_NS                   !       
   elseif (SPH_ProblemType.eq.3) then       !     DF
!       call IntForces_DF                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_DF                   !
   else
       write(*,*) ' wrong type of problem= ', SPH_ProblemType
       pause
   endif
   
   DO ipoin = 1, npoin      
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      Dvx_RK(:,ipoin) = Dvx_RK  (:,ipoin) + f2rk(i)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      Dx_RK (:,ipoin) = Dx_RK   (:,ipoin) + f2rk(i)*vx(:,ipoin)
      if (.not.summation_density) then
         Drho_RK (ipoin) = Drho_RK (ipoin)   + f2rk(i)*drho(ipoin) !       
      endif
      Du_RK   (ipoin) = Du_RK   (ipoin)   + f2rk(i)*du  (ipoin)
   ENDDO
   ipoin = 1 !  just for debug
   
ENDDO

DO ipoin =1,npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   x  (:,ipoin) =  x0 (:,ipoin) + (dt_sph/6.)* Dx_RK (:,ipoin)       
   vx (:,ipoin) = vx0 (:,ipoin) + (dt_sph/6.)* Dvx_RK(:,ipoin)
   if (.not.summation_density) then
      rho(ipoin)  = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)        
   endif
   u  (ipoin)   =  u0   (ipoin) + (dt_sph/6.)* Du_RK (ipoin)
ENDDO

if (summation_density) then       ! ------  MP march 2018 DBKs   
       call sum_density_TOOLS     ! CHGD MP April 2020 
endif
        
!     Call Vn0_BCs_SW ( dt_sph/6. )    deactivated and moved down MP 30.12.2018  

Call Check_Out_Domain                                            

if (SPH_ProblemType.eq.1) then            

   call Get_s_at_w_DF_SW                ! *** get hw at soil pics and viceversa
   
   Drho_RK = rho - rho0			!    we will use drho_rk here to avoid increasing memory
   
   if (nUaux.GE.1) then
      
      if (allocated(auxS1)) then
         dauxS1 = auxS1 - auxS10
         call Get_Pwp_SW (Drho_RK, dauxS1)   !    compute pwp DFs and pwp
      else 
         call Get_Pwp_SW (Drho_RK)		!    compute pwp
      endif
      
      do ipoin = 1,npoin
	  u(ipoin) = UAux(1,ipoin)
      enddo
      
!      u (:) = Uaux (1,:)
   endif
             
    call ExtForces_SW (1)                !       External Forces SW (slope, friction,virt.parts)
                                         !       March 2014 1 means erosion, once we know v  
                                         !       Replaces old Erosion_SW 
        
    Call Vn0_BCs_SW ( dt_sph )           !  **  Moved here MP 30.12.2018 changhd dt/6 to dt

    Call Abs_BCs_SW                      !    remember this should go also to other time integration routines

      
!      ------   
!
   if(nuaux.eq.0) then
     DO ipoin=1,npoin
        if ( if_Out_Domain(ipoin).eq.1) CYCLE
        if ( u(ipoin).gt.rho(ipoin) ) u(ipoin) = rho(ipoin)  ! ** CHK pwp
        if ( u(ipoin).lt.0.0 )        u(ipoin) = 0.0
     Enddo
   endif

   
endif      
                                        
if (average_velocity) then      !     March 2014 Moved here, more logical for erosion and BCs
   call av_vel		            !     Average v of particles to avoid penetration
   do ipoin = npoin0 + 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      vx (:,ipoin) = vx(:,ipoin) + av(:,ipoin)
   enddo   
endif
   
if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML

!if (SPH_ProblemType.eq.1) then            
!   call Get_s_at_w_DF_SW               ! *** get hw at soil pics and viceversa
!endif                                  ! now done in internal forces ***

ipoin = 1        ! *** to debug

!write(res_file3,1991) time_sph+dt_sph, (Uaux(i,1),i=1,nUaux-1), exp(-3.141593*3.141593*0.40528473*(time_sph+dt_sph)/4)
1991 format (12(f10.4,2x))


!      ------   Clears linked list containing interaction data

!   call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  RK4 


!----------------------------------------------------------------------

      subroutine RK4_WIR (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) 

!----------------------------------------------------------------------

!   Subroutine to integrate in time. 
!   First, (dt/2) RK4 for sources
!   Then   convective terms (pressure)
!   Second (dt/2) RK4 for sources

implicit none
integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
real   (irk) Dx_RK(ndimn,ntotal), Dvx_RK(ndimn,ntotal), Du_RK(ntotal), Drho_RK(ntotal)

data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

!      ------  Initialize ( done also in subroutines, strictly unnecessary )

avdudt     = 0.
ahdudt     = 0.
indvxdt    = 0.
ardvxdt    = 0.
exdvxdt    = 0.
av         = 0.

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0
   
!      ------   Interaction parameters, calculating neighboring particles
!               and optimizing smoothing length

if (nnps.eq.1) then 
    call direct_find(itimestep_sph)
else if (nnps.eq.2) then
    call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
elseif (nnps.eq.3) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.4) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
endif

! if (SPH_ProblemType.eq.1) call Pint_Update_SW         !      Updates interactions to 0 (default) 
if (SPH_ProblemType.eq.1) call Pint_Update_SW_TOOLS     !      CHGD MP April 2020

t0 = time_sph; vx0 = vx; x0 = x; u0 = u
where (itype.eq.2) rho0 = rho                           !  **** changed 30 Juin  water

DO i = 1,4                           ! ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph
   do ipoin = 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      x (:,ipoin)  = x0(:,ipoin) + f1rk(i)*(dt_sph)*vx(:,ipoin)
      vx(:,ipoin) = vx0(:,ipoin) + f1rk(i)*(dt_sph)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      if (itype(ipoin).eq.2) then
          rho(ipoin) = rho0(ipoin) + f1rk(i)*(dt_sph) * drho(ipoin)                     !  **** changed 30 Juin  water
      endif          
      u (ipoin) = u0 (ipoin)   + f1rk(i)*(dt_sph)*du(ipoin)
   enddo
   Call Check_Out_Domain 
   if (SPH_ProblemType.eq.1) then          !     SW
      DO ipoin=1,npoin
         if ( u(ipoin).gt.rho(ipoin) ) u(ipoin) = rho(ipoin)
         if ( u(ipoin).lt.0.0 )        u(ipoin) = 0.0
      Enddo
    endif      
   
!      ------   Get alpha correction factor

   if (pa_sph.eq.3.or.pa_sph.eq.4) then                 !      Alpha correction factor JB                                                                
       call Get_AlphaJB
   endif   

   call get_wir_density
   
   if (SPH_ProblemType.eq.1) then          !     SW
       call IntForces_SW                   !       External Forces SW (slope, friction,virt.parts) 
       call ExtForces_SW  (0)              !       Internal forces SW (pressure gradient)
   elseif (SPH_ProblemType.eq.2) then      !     NS
!       call IntForces_NS                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_NS                   !       
   elseif (SPH_ProblemType.eq.3) then      !     DF
!       call IntForces_DF                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_DF                   !
   else
       write(*,*) ' wrong type of problem= ', SPH_ProblemType
       pause
   endif
   
   DO ipoin = 1, npoin      
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      Dvx_RK(:,ipoin) = Dvx_RK  (:,ipoin) + f2rk(i)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      Dx_RK (:,ipoin) = Dx_RK   (:,ipoin) + f2rk(i)*vx(:,ipoin)
      Du_RK   (ipoin) = Du_RK   (ipoin)   + f2rk(i)*du  (ipoin)
      if (itype(ipoin).eq.2) then
         Drho_RK(ipoin) = Drho_RK(ipoin) + f2rk(i)*drho(ipoin) 
      endif                     
   ENDDO

   
ENDDO

DO ipoin=1,npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   x  (:,ipoin) =  x0 (:,ipoin) + (dt_sph/6.)* Dx_RK (:,ipoin)
   vx (:,ipoin) = vx0 (:,ipoin) + (dt_sph/6.)* Dvx_RK(:,ipoin)
!    if (.not.summation_density) then                            !    MP 2013
!      rho(ipoin)  = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)    !  water done twice      
!    endif                                                       !  if sum=false
   u  (ipoin)   =  u0   (ipoin) + (dt_sph/6.)* Du_RK (ipoin)
   if (itype(ipoin).eq.2) then
      rho(ipoin) = rho0(ipoin) + (dt_sph/6.)*Drho_RK(ipoin)
   endif      
ENDDO

Call Check_Out_Domain                                           !      ******* MP



if (SPH_ProblemType.eq.1) then          !     SW

   Call Abs_BCs_SW                      !    remember this should go also to other time integration routines

   DO ipoin=1,npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      if ( u(ipoin).gt.rho(ipoin) ) u(ipoin) = rho(ipoin)
      if ( u(ipoin).lt.0.0 )        u(ipoin) = 0.0
   Enddo
   
endif      
   
if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML

if (average_velocity) then              !     Average v of particles to avoid penetration
   call av_vel
   vx (:,npoin0 + 1:npoin) = vx(:,npoin0 + 1:npoin) + av(:,npoin0 + 1:npoin)
endif

!      ------   Clears linked list containing interaction data

!    call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  RK4_WIR 


!----------------------------------------------------------------------

      subroutine chk_rho (itimestep_sph, Dx_RK, Dvx_RK, Du_RK, Drho_RK ) 

!----------------------------------------------------------------------


implicit none
integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
real   (irk) Dx_RK(ndimn,npoin), Dvx_RK(ndimn,npoin), Du_RK(npoin), Drho_RK(npoin)

data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

!      ------  Initialize ( done also in subroutines, strictly unnecessary )

av         = 0.

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0
   
!      ------   Interaction parameters, calculating neighboring particles
!               and optimizing smoothing length

if (nnps.eq.1) then 
    call direct_find(itimestep_sph)
else if (nnps.eq.2) then
    call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
elseif (nnps.eq.3) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.4) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
endif

! if (SPH_ProblemType.eq.1) call Pint_Update_SW           !      Updates interactions to 0 (default)
if (SPH_ProblemType.eq.1) call Pint_Update_SW_TOOLS       !      CHGD MP april 2020 

t0 = time_sph; vx0 = vx; x0 = x; u0 = u
if (.not.summation_density) rho0 = rho

DO i = 1,4                              ! ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph
   x (:,1:npoin)  = x0(:,1:npoin) + f1rk(i)*(dt_sph)*vx(:,1:npoin)
   
   call get_vx                          ! ------    get velocity depending on position 
   
   if (.not.summation_density) then
       rho(1:npoin) = rho0 (1:npoin)   + f1rk(i)*(dt_sph)*drho(1:npoin)
   endif 
   Call Check_Out_Domain 
   
!      ------   Get alpha correction factor

if (pa_sph.eq.3.or.pa_sph.eq.4) then                    !      Alpha correction factor JB                                                                
    call Get_AlphaJB
endif   

   if (summation_density) then       ! ------  If continuity is used, then get drho    
       call sum_density_TOOLS    ! CHGD MP April 2020        
     else             
       call con_density_TOOLS    ! CHGD MP April 2020       
   endif 
   
   Dvx_RK(:,1:npoin) = 0.0
   Dx_RK (:,1:npoin) = Dx_RK   (:,1:npoin) + f2rk(i)*vx(:,1:npoin)
   if (.not.summation_density) then
      Drho_RK (1:npoin) = Drho_RK (1:npoin)   + f2rk(i)*drho(1:npoin)   
   endif 
   
ENDDO

x  (:,1:npoin) =  x0 (:,1:npoin) + (dt_sph/6.)* Dx_RK (:,1:npoin)
Call Check_Out_Domain                                           ! ******* MP 
if (.not.summation_density) then
   rho(1:npoin)  = rho0(1:npoin) + (dt_sph/6.)*Drho_RK(1:npoin)
endif  
   
if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML

if (average_velocity) then              !     Average v of particles to avoid penetration
   call av_vel
   vx (:,1:npoin) = vx(:,1:npoin) + av(:,1:npoin)
endif

!      ------   Clears linked list containing interaction data

!    call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  chk_rho 


!----------------------------------------------------------------------

      subroutine get_vx  

!----------------------------------------------------------------------


implicit none
integer(ink) ipoin 

DO ipoin = 1, npoin
!   vx(1,ipoin) = 0.1
   vx(1,ipoin) = 1.0*(x(1,ipoin)+10.00)/20.
ENDDO


END subroutine get_vx  


! -----------------------------------------------------------------
 
       SUBROUTINE Check_Out_Domain  
 
! -----------------------------------------------------------------

!   Subroutine to check wheter a given point is outside the computational domain  
 
implicit none

integer(ink)  ipoin, idimn
real   (irk)  dxx 

!      ------  Initialization 

! If_Out_Domain = 0               CHANGED 19 July 2007

!      ------  

DO ipoin = 1,npoin
   Do idimn = 1, ndimn
      dxx   = (x(idimn,ipoin)-Xmin_Domain(idimn))*(x(idimn,ipoin)-Xmax_Domain(idimn))
      if ( dxx.gt.0.0 ) If_Out_Domain(ipoin) = 1
   Enddo
Enddo

END  SUBROUTINE Check_Out_Domain


! -----------------------------------------------------------------
 
       SUBROUTINE Check_Out_Domain_SAVE  
 
! -----------------------------------------------------------------

!   Subroutine to check wheter a given point is outside the computational domain  
 
implicit none

integer(ink)  ipoin, idimn
real   (irk)  dxx, rr

!      ------  Initialization 

!If_Out_Domain = 0    CHANGED 19 July 2007

!      ------  

DO ipoin = 1,npoin
   rr = 0.0
   Do idimn = 1, ndimn
      dxx   = (x(idimn,ipoin)-Xmin_Domain(idimn))*(x(idimn,ipoin)-Xmax_Domain(idimn))
      rr    = rr + x(idimn,ipoin)*x(idimn,ipoin)
      if ( dxx.ge.0.0 ) If_Out_Domain(ipoin) = 1
   Enddo
   if (rr.gt.1.1) then
      If_Out_Domain(ipoin) = 1
   endif
Enddo

END  SUBROUTINE Check_Out_Domain_SAVE


!----------------------------------------------------------------------

      subroutine direct_find (itimestep_sph ) 

!----------------------------------------------------------------------

implicit none

integer(ink)  i,j, idimn, niac, itimestep_sph 
integer(ink) sumiac, maxiac, miniac, noiac 
integer(ink) maxp, minp, scale_k
real   (irk) driac, r, mhsml, w 
real   (irk), allocatable::  dxiac(:), tdwdx(:)

if (.NOT.allocated (dxiac) ) then 
   allocate (dxiac(ndimn),tdwdx(ndimn) )
endif

nullify  (last) 

niac = 0
countiac(:) = 0

if (skf.eq.1.or.skf.eq.4) then  ! Saeidmt 24Oct2020 Wendland kernel (skf=4)
   scale_k = 2 
else if (skf.eq.2) then 
   scale_k = 3 
else if (skf.eq.3) then 
   scale_k = 3 
endif 
   
DO i = 1,ntotal-1
   if ( If_Out_Domain(i).eq.1 ) CYCLE
   DO j = i+1,ntotal
      if ( If_Out_Domain(j).eq.1 ) CYCLE
      dxiac(1) = x(1,i)-x(1,j)
      driac    = dxiac(1)*dxiac(1)
      do idimn = 2, ndimn
         dxiac(idimn) = x(idimn,i)-x(idimn,j)
         driac = driac + dxiac(idimn)*dxiac(idimn)
      enddo
      mhsml =  (hsml(i)+hsml(j))/2.
      r = sqrt(driac)
      if (r.lt.scale_k*mhsml) then
         niac = niac + 1
         countiac(i) = countiac(i) + 1
         countiac(j) = countiac(j) + 1
         call kernel(r, dxiac, mhsml, w, tdwdx)  !      using r, dxiac and mshml obtains w and its gradient
         allocate (current)                      !      create next pair structure
         allocate (current%dwdx(ndimn))          !      allocate dwdx
         current%niac      = niac                !      this is the pair number
         current%pair_i    = i                   !      first particle
         current%pair_j    = j                   !      and second particle
         current%Pint_type = 0                   !      Sets interaction type to default (same material)
         current%w         = w                   !      Weigth
         current%dwdx(:)   = tdwdx(:)            !      Gradient of weight
         current%next      => last               !      
         last              => current            !
      endif
   enddo
enddo

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
 
if (mod(itimestep_sph,print_step).eq.0) then
       print *,' >> Statistics: interactions per particle:'
       print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
       print *,'**** Particle:',minp, ' minimal interactions:',miniac
       print *,'**** Average :',real(sumiac)/real(ntotal)
       print *,'**** Total pairs : ',niac
       print *,'**** Particles with no interactions:',noiac
endif   

deallocate (dxiac, tdwdx )


END SUBROUTINE  direct_find 


!----------------------------------------------------------------------
 
       subroutine grid_find_NEW_OLD  (itimestep_sph)  
 
!----------------------------------------------------------------------

!      ------  This is the new routine
 
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
 
if (skf.eq.1.or.skf.eq.4) then  ! Saeidmt 24Oct2020 Wendland kernel (skf=4)
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

if (itimestep_sph.eq.1) then
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
               
               call kernel(r, dxiac, mhsml, w, tdwdx)  !      using r, dxiac and mshml obtains w and its gradient
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
endif    


END SUBROUTINE  grid_find_NEW_OLD 

 
 
!----------------------------------------------------------------------
 
       subroutine grid_find_OLD (itimestep_sph)  
 
!----------------------------------------------------------------------
 
 implicit none
 
 integer(ink)  idimn, itotal, jtotal
 integer(ink)  idx, idy, idz, jdx, jdy, jdz, idx0, idy0, idz0, idx1, idy1, idz1
 integer(ink)    i,   j,   k
 integer(ink)  niac, itimestep_sph 
 integer(ink)  sumiac, maxiac, miniac, noiac 
 integer(ink)  maxp, minp, scale_k 
 
 real   (irk)  length, length_new 
 real   (irk)  driac, r, mhsml, w 
 real   (irk), allocatable::  dxiac(:), tdwdx(:)
 real   (irk), allocatable::  xmin(:), xmax(:), deltx(:)
 
 integer(ink)  ndivx(3),idivx(3),jdivx(3)
 integer(ink)  ndivt,   idivt,   jdivt
 
 type cell
      integer(ink) ipic
      type   (cell), pointer:: next
 end type cell
 
 type list_of_cells
      integer(ink) npics
      type   (cell), pointer:: current, last, keepit
 endtype list_of_cells
 
 type (list_of_cells), allocatable:: grid(:)
 
 allocate (xmin(ndimn),  xmax(ndimn),  deltx(ndimn))
 
 xmin = 1.e+10; xmax = -xmin; deltx = 0.0; ndivx = 1; idivx = 1; jdivx = 1
 
 if (skf.eq.1) then 
    scale_k = 2 
 else if (skf.eq.2) then 
    scale_k = 3 
 else if (skf.eq.3) then 
    scale_k = 3 
 endif 
 
 !      ------  Task 1: Setup Grid parameters
 
  
 Do idimn = 1,ndimn
    Do itotal = 1, ntotal
       if ( If_Out_Domain(itotal).eq.1 ) CYCLE
       if (x(idimn,itotal).lt.xmin(idimn)) xmin(idimn)  = x(idimn,itotal)
       if (x(idimn,itotal).gt.xmax(idimn)) xmax(idimn)  = x(idimn,itotal)
       if (deltx(idimn).lt.hsml(itotal))   deltx(idimn) = hsml(itotal)
    enddo
    deltx(idimn) = deltx(idimn)*2    !   was multiplied by scale_k*
    length       = xmax(idimn) - xmin(idimn)
    ndivx(idimn) = (length/deltx(idimn)) + 1
    length_new   = ndivx(idimn)*deltx(idimn)
    xmin(idimn)  = xmin(idimn) - (length_new-length)/2 - 0.001*length
    xmax(idimn)  = xmax(idimn) + (length_new-length)/2 + 0.001*length
 Enddo
 
 !      ------  Task 2: Fill Grid structure with particles
 
 ndivt = ndivx(1)*ndivx(2)*ndivx(3)    
 allocate ( grid(ndivt) )
 
 Do idivt = 1, ndivt
    nullify (grid(idivt)%last)
 enddo
 
 DO itotal = 1, ntotal
    if ( If_Out_Domain(itotal).eq.1 ) CYCLE
    DO idimn = 1, ndimn
       idivx(idimn) = (x(idimn,itotal)-xmin(idimn))/deltx(idimn) + 1
          if (idivx(idimn).gt.ndivx(idimn)) idivx(idimn)=ndivx(idimn)
    enddo
    idivt =  ndivx(1)*ndivx(2)*(idivx(3)-1) + ndivx(1)*(idivx(2)-1) + idivx(1)
    allocate (grid(idivt)%current)
    grid(idivt)%npics = grid(idivt)%npics +1
    grid(idivt)%current%ipic = itotal
    grid(idivt)%current%next => grid(idivt)%last
    grid(idivt)%last => grid(idivt)%current
 Enddo
 
 !      ------  Task 3: Get list of interactions
 
 allocate (dxiac(ndimn), tdwdx(ndimn) )
 nullify  (last) 
 niac = 0
 countiac(:) = 0
 
 DO idx = 1, ndivx(1)
 DO idy = 1, ndivx(2)
 DO idz = 1, ndivx(3)
    idivt =  ndivx(1)*ndivx(2)*(idz-1) + ndivx(1)*(idy-1) + idx 
    grid(idivt)%current => grid(idivt)%last
    DO while ( associated (grid(idivt)%current))
       itotal = grid(idivt)%current%ipic
       idx0   = max (       1,idx-1); idy0   = max (       1,idy-1); idz0   = max (       1,idz-1)
       idx1   = min (ndivx(1),idx+1); idy1   = min (ndivx(2),idy+1); idz1   = min (ndivx(3),idz+1)
       grid(idivt)%keepit=>grid(idivt)%current
       DO jdx = idx0, idx1  
       DO jdy = idy0, idy1
       DO jdz = idz0, idz1
          jdivt = ndivx(1)*ndivx(2)*(jdz-1) + ndivx(1)*(jdy-1) + jdx
          if (jdivt.lt.idivt) CYCLE
          grid(jdivt)%current => grid(jdivt)%last
          if (jdivt.eq.idivt) THEN                                !  ------  If we are at same cell, search starts at 
              grid(jdivt)%current => grid(idivt)%keepit           !          position next to idivt-current
              grid(jdivt)%current => grid(idivt)%current%next
          endif    
          DO while ( associated (grid(jdivt)%current))
             jtotal = grid(jdivt)%current%ipic
             dxiac(1) = x(1,itotal)-x(1,jtotal)                     !          note jdivt is always greater than idivt
             driac    = dxiac(1)*dxiac(1)
             do idimn = 2, ndimn
               dxiac(idimn) = x(idimn,itotal)-x(idimn,jtotal)
               driac = driac + dxiac(idimn)*dxiac(idimn)
             enddo
             mhsml =  (hsml(itotal)+hsml(jtotal))/2.
             r = sqrt(driac)
             if (r.lt.scale_k*mhsml) then
               niac = niac + 1
               countiac(itotal) = countiac(itotal) + 1
               countiac(jtotal) = countiac(jtotal) + 1
               call kernel(r, dxiac, mhsml, w, tdwdx)  !      using r, dxiac and mshml obtains w and its gradient
               allocate (current)                      !      create next pair structure
               allocate (current%dwdx(ndimn))          !      allocate dwdx
               current%niac      = niac                !      this is the pair number
               current%pair_i    = itotal              !      first particle
               current%pair_j    = jtotal              !      and second particle
               current%Pint_type = 0                   !      Sets interaction type to default (same material)         
               current%w         = w                   !      Weigth
               current%dwdx(:)   = tdwdx(:)            !      Gradient of weight
               current%next      => last               !      
               last              => current            !
             endif
             grid(jdivt)%current => grid(jdivt)%current%next      
          ENDDO
          grid(idivt)%current=>grid(idivt)%keepit           !      This is done to avoid spoiling j loop 
       ENDDO
       ENDDO
       ENDDO
       
       grid(idivt)%current => grid(idivt)%current%next
       
    ENDDO
    
 ENDDO
 ENDDO
 ENDDO
 
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
 
if (mod(itimestep_sph,print_step).eq.0) then
       print *,' >> Statistics: interactions per particle:', time, dt_sph
       print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
       print *,'**** Particle:',minp, ' minimal interactions:',miniac
       print *,'**** Average :',real(sumiac)/real(ntotal)
       print *,'**** Total pairs : ',niac
       print *,'**** Particles with no interactions:',noiac
endif    
 
 !      ------  Clear GRID, because there is no longer needed
 
 DO idivt = 1,ndivt
    grid(idivt)%current => grid(idivt)%last
    DO while (associated(grid(idivt)%last))
       grid(idivt)%last => grid(idivt)%current%next
       deallocate (grid(idivt)%current)
       grid(idivt)%current=>grid(idivt)%last
    ENDDO
 ENDDO
 
 
 deallocate (dxiac, tdwdx )
 deallocate (xmin, xmax, deltx, grid )
 
 
END SUBROUTINE  grid_find_OLD 
 

!--------------------------------------------
      subroutine cross_link
!--------------------------------------------

!      ------  Crossing the list

integer(ink) i,j

!  open(1,file="link.res")
!  Write(*,*) ' ----  Now cross the list and print  '
Write(res_file2,*) ' ----  Now cross the list and print  '
current=> last
Do while (associated(current))
   i =  min(current%pair_i,current%pair_j)
   j =  max(current%pair_i,current%pair_j)
   write(res_file2,100) i,j,current%niac
!  write(1,100) current%niac, current%pair_i, current%pair_j, current%w, current%dwdx(1:ndimn)
!  write(*,100) current%niac, current%pair_i, current%pair_j, current%w, current%dwdx(1:ndimn) 
   current=>current%next
enddo
100 format(3(2x,I6),5(2x,e11.4))

end subroutine cross_link

! -------------------------------------------------------------

     subroutine kernel(r,dx,hsml,w,dwdx)   

! -------------------------------------------------------------

!   Subroutine to calculate the smoothing kernel wij and its 
!              derivatives dwdxij.
!   if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!          = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!          = 3, Quintic kernel (Morris 1997)

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
endif 
        
END SUBROUTINE kernel

!---------------------------------------------------------------

     subroutine Get_AlphaJB

!---------------------------------------------------------------
 

!     Subroutine to calculate alpha factor
!
!     alpha(I) = - SUM J ( Mj * Rij * Grad Wij)  with Rij = Xj-Xi
!                then divide it by ndimn*rho(i) 

implicit none
  

integer(ink) i ,j ,idimn    
real   (irk) xcc  
real   (irk), allocatable:: dtx(:) 

allocate ( dtx(ndimn)  )

alphaJB = 0.      

current=> last
Do while (associated(current))
   i = current%pair_i 
   j = current%pair_j
   if (current%Pint_type.eq.0) then 
      do idimn = 1,ndimn
         dtx(idimn) =  x(idimn,j) -  x(idimn,i) 
      enddo        
      xcc = dtx(1)*current%dwdx(1)        
      do idimn = 2,ndimn
         xcc = xcc + dtx(idimn)*current%dwdx(idimn)
      enddo 
   
      if (pa_sph.eq.3) then
         alphaJB(i) = alphaJB(i) + mass(j)*xcc
         alphaJB(j) = alphaJB(j) + mass(i)*xcc
      elseif (pa_sph.eq.4) then
         alphaJB(i) = alphaJB(i) + mass(j)*xcc/rho(j)
         alphaJB(j) = alphaJB(j) + mass(i)*xcc/rho(i)
      endif
   endif
   
   current=>current%next
enddo

if (pa_sph.eq.3) then
   Do i = 1,ntotal
      alphaJB(i) =  alphaJB(i)/(ndimn*rho(i))
   enddo
endif

deallocate ( dtx  ) 

END SUBROUTINE Get_AlphaJB


!---------------------------------------------------------------

     subroutine get_wir_density

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

END    subroutine get_wir_density



!---------------------------------------------------------------

     subroutine sum_density

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

if (.NOT.allocated (hv) ) then 
   allocate ( hv(ndimn), wi(ntotal) )
endif
hv = 0.0; wi = 0.0
 
hv = 0.e0
r  = 0.
      
!      ------  First  calculate the integral of the kernel over the space

do i = 1,ntotal
   if (if_out_domain(i).eq.1) cycle
   call kernel(r,hv,hsml(i),selfdens,hv)
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
   call kernel(r,hv,hsml(i),selfdens,hv)
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

 END SUBROUTINE sum_density
 

!---------------------------------------------------------------

     subroutine con_density

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
       call kernel(r,hv,hsml(i),selfdens,hv)
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
       call kernel(r,hv,hsml(i),selfdens,hv)
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

 END SUBROUTINE con_density


!---------------------------------------------------------------

     subroutine  h_upgrade 

!---------------------------------------------------------------
!     Subroutine to evolve smoothing length  Type 3 by Lingang 11th Feb 2022

!     dt     : time step                                            [in]
!     ntotal : Number of particles                                  [in]
!     mass   : Particle masses                                      [in]
!     vx     : Velocities of all particles                          [in]
!     rho    : Density                                              [in]
!     niac   : Number of interaction pairs                          [in]
!     pair_i : List of first partner of interaction pair            [in]
!     pair_j : List of second partner of interaction pair           [in]
!     dwdx   : Derivative of kernel with respect to x, y and z      [in]
!     hsml   : Smoothing Length                                 [in/out]
        
implicit none

integer(ink) i,j,idimn,B
integer(ink), allocatable:: N_neighb(:)
real   (irk) fac, hvcc, Ls 
real   (irk), allocatable:: dvx(:), vcc(:), dhsml(:), h_by_m(:)

if (.NOT.allocated (dvx) ) then
   allocate ( dvx(ndimn), vcc(npoin), dhsml(npoin)  )
endif
dvx = 0.0; vcc = 0.0; dhsml = 0.0

if (sle.eq.0 ) then         !      Keep smoothing length unchanged    
    return
else if (sle.eq.2) then     !      dh/dt = (-1/dim)*(h/rho)*(drho/dt).    
    vcc(:) = 0.0 
    current=> last          !    Initializes linked list to last element 
    DO WHILE (associated(current))
       i = current%pair_i 
       j = current%pair_j 
       if (current%Pint_type.eq.0) then
           if ( (itype(i).gt.0).and.(itype(j).gt.0) ) then 
              do idimn =1,ndimn
                 dvx(idimn) = vx(idimn,j) - vx(idimn,i) 
              enddo
              hvcc = dvx(1)*current%dwdx(1)
              do idimn = 2,ndimn
                 hvcc = hvcc + dvx(idimn)*current%dwdx(idimn)     
              enddo    
              vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
              vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
           endif
       endif
       current=>current%next
    ENDDO  
        
    do i = 1, npoin
   if (if_out_domain(i).eq.1) cycle
       dhsml(i) = (hsml(i)/ndimn)*vcc(i)
       hsml(i) = hsml(i) + dt_sph*dhsml(i)      
       if (hsml(i).le.0) hsml(i) = hsml(i) - dt_sph*dhsml(i) 
    enddo
    
else if(sle.eq.1) then
            
     fac = 2.0
     do i = 1, npoin
        if (if_out_domain(i).eq.1) cycle          
        hsml(i) = fac * (mass(i)/rho(i))**(1./ndimn)
     enddo

else if(sle.eq.3) then                       ! Lingang 11Feb2022   **** McDougall & Hungr 2004
    if (.NOT.allocated (N_neighb) ) then
        allocate(N_neighb(npoin), h_by_m(npoin))
    endif
    N_neighb(:)=0
    h_by_m(:)=0
    B=4
    current=> last
    do while (associated(current))
       i = current%pair_i
       j = current%pair_j
       if ( current%Pint_type.eq.0.OR.current%Pint_type.eq. 1 ) THEN
          N_neighb(i)=N_neighb(i)+1
          N_neighb(j)=N_neighb(j)+1               !  MP feb 
          h_by_m(i)=h_by_m(i)+rho(i)/mass(i)
          h_by_m(j)=h_by_m(j)+rho(j)/mass(j)
       endif
       current=>current%next
    enddo
    
    do i=1,npoin
        if (if_out_domain(i).eq.1) cycle
        Ls=(h_by_m(i)/N_neighb(i))**0.5
        hsml(i)=B/Ls
    enddo
       
endif 

END SUBROUTINE h_upgrade





!---------------------------------------------------------------

     subroutine  h_upgrade_OLD 

!---------------------------------------------------------------
!     Subroutine to evolve smoothing length

!     dt     : time step                                            [in]
!     ntotal : Number of particles                                  [in]
!     mass   : Particle masses                                      [in]
!     vx     : Velocities of all particles                          [in]
!     rho    : Density                                              [in]
!     niac   : Number of interaction pairs                          [in]
!     pair_i : List of first partner of interaction pair            [in]
!     pair_j : List of second partner of interaction pair           [in]
!     dwdx   : Derivative of kernel with respect to x, y and z      [in]
!     hsml   : Smoothing Length                                 [in/out]
        
implicit none

integer(ink) i,j,idimn
real   (irk) fac, hvcc 
real   (irk), allocatable:: dvx(:), vcc(:), dhsml(:)

if (.NOT.allocated (dvx) ) then
   allocate ( dvx(ndimn), vcc(npoin), dhsml(npoin)  )
endif
dvx = 0.0; vcc = 0.0; dhsml = 0.0

if (sle.eq.0 ) then         !      Keep smoothing length unchanged    
    return
else if (sle.eq.2) then     !      dh/dt = (-1/dim)*(h/rho)*(drho/dt).    
    vcc(:) = 0.0 
    current=> last          !    Initializes linked list to last element 
    DO WHILE (associated(current))
       i = current%pair_i 
       j = current%pair_j 
       if (current%Pint_type.eq.0) then
           if ( (itype(i).gt.0).and.(itype(j).gt.0) ) then 
              do idimn =1,ndimn
                 dvx(idimn) = vx(idimn,j) - vx(idimn,i) 
              enddo
              hvcc = dvx(1)*current%dwdx(1)
              do idimn = 2,ndimn
                 hvcc = hvcc + dvx(idimn)*current%dwdx(idimn)     
              enddo    
              vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
              vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
           endif
       endif
       current=>current%next
    ENDDO  
        
    do i = 1, npoin
   if (if_out_domain(i).eq.1) cycle
       dhsml(i) = (hsml(i)/ndimn)*vcc(i)
       hsml(i) = hsml(i) + dt_sph*dhsml(i)      
       if (hsml(i).le.0) hsml(i) = hsml(i) - dt_sph*dhsml(i) 
    enddo
    
else if(sle.eq.1) then
            
     fac = 2.0
     do i = 1, npoin
        if (if_out_domain(i).eq.1) cycle          
        hsml(i) = fac * (mass(i)/rho(i))**(1./ndimn)
     enddo
       
endif 

END SUBROUTINE h_upgrade_OLD


!---------------------------------------------------------------

     subroutine  av_vel 

!---------------------------------------------------------------


!     Subroutine to calculate the average velocity to correct velocity
!     for preventing.penetration (Monaghan, 1992)

!     ntotal : Number of particles                                  [in]
!     mass   : Particle masses                                      [in]
!     niac   : Number of interaction pairs                          [in]
!     pair_i : List of first partner of interaction pair            [in]
!     pair_j : List of second partner of interaction pair           [in]
!     w      : Kernel for all interaction pairs                     [in]
!     vx     : Velocity of each particle                            [in]
!     rho    : Density of each particle                             [in]
!     av     : Average velocityof each particle                    [out]

implicit none
      
integer(ink)  i,j, idimn
real   (irk)  epsilon
real   (irk), allocatable :: dvx(:)

if (.NOT.allocated (dvx) )  allocate ( dvx(ndimn) )
dvx = 0.0

!     epsilon --- a small constants chosen by experence, may lead to instability.
!             for example, for the 1 dimensional shock tube problem, the E <= 0.3

epsilon = 0.3    

av(:,:) = 0.0

    current=> last          !    Initializes linked list to last element 
    DO WHILE (associated(current))
       i = current%pair_i 
       j = current%pair_j 
        if (current%Pint_type.eq.0.and.i.gt.npoin0.and.j.gt.npoin0) then 
          do idimn = 1,ndimn
             dvx(idimn) = vx(idimn,i) - vx(idimn,j)            
             av(idimn, i) = av(idimn,i) - 2*mass(j)*dvx(idimn)/(rho(i)+rho(j))*current%w 
             av(idimn, j) = av(idimn,j) + 2*mass(i)*dvx(idimn)/(rho(i)+rho(j))*current%w                      
          enddo
       endif
       current=>current%next
    ENDDO 
 
    do i = 1, npoin  				! *** TB uses npoin1
       if (if_out_domain(i).eq.1) cycle
       do idimn = 1, ndimn
          av(idimn,i) = epsilon * av(idimn,i)
       enddo 
    enddo   

deallocate ( dvx )
    
END SUBROUTINE av_vel


!---------------------------------------------------------------------- 

      subroutine clear_interactions
      
!----------------------------------------------------------------------  

!     Subroutine for clearing linked list containing interaction data

implicit none      

!     ------  Deleting

current=> last
DO while (associated(last))
   last=> current%next
   deallocate (current%dwdx)
   deallocate (current)
   current=>last
enddo                           
     
END SUBROUTINE clear_interactions 

!----------------------------------------------------------------------

      subroutine RK4_FS  (itimestep_sph, dt_sph_factor, if_IntForces_SW, if_density_comput)

!----------------------------------------------------------------------

implicit none

integer(ink) if_IntForces_SW, if_density_comput  ! ** MP if=0 do not compute int forces or density
real   (irk) dt_sph_factor                       ! multiplies dt_factor. Used if 2 TKs in FracTStep

integer(ink) itimestep_sph, i, ipoin 
real   (irk) t, t0, f1rk(4), f2rk(4), dtby2
real   (irk) Dx_RK(ndimn,npoin), Dvx_RK(ndimn,npoin), Du_RK(npoin), Drho_RK(npoin) ! These are local!
						    ! defined as allocat 
						    ! in MAIN Vars
real   (irk) zz, dhh, xi	!      used to update pwp with dh										   
integer(ink) icUaux, iaux	!      icUaux=0 we do not use PWP + FD's
real   (irk) beta           ! **** **** **** **** **** ****

data f1rk/0., 0.5, 0.5, 1.0/
data f2rk/1., 2., 2., 1.0/

!      ------  Initialize ( done also in subroutines, strictly unnecessary )

avdudt     = 0.
ahdudt     = 0.
indvxdt    = 0.
ardvxdt    = 0.
exdvxdt    = 0.
av         = 0.     

Dx_RK = 0.0; Dvx_RK= 0.0 ;Du_RK= 0.0; Drho_RK = 0.0
   
!      ------   Interaction parameters, calculating neighboring particles
!               and optimizing smoothing length

if (nnps.eq.1) then 
    call direct_find(itimestep_sph)
else if (nnps.eq.2) then
    call grid_find_new_TOOLS  (itimestep_sph)      ! CHGD MP april 2020
elseif (nnps.eq.3) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS2  (itimestep_sph)     ! CHGD MP april 2020
elseif (nnps.eq.4) then                            ! CHGD MP april 2020
    call grid_find_new_TOOLS3  (itimestep_sph)     ! CHGD MP april 2020
endif

if (SPH_ProblemType.eq.1) then 
    ! call Pint_Update_SW        !      CHGD MP april 2020
    call Pint_Update_SW_TOOLS    !      Updates interactions to 0 (default) 
                                 !      Gets list of interactions and nr of interactions if FS
                                 !      MP ** FS 29 april 2017  
    if (itimestep_sph.eq.1) then       !   we need normalized hs at the beginning ONLY MP 2021 29 jan
        call Vn0_BCs_SW (0.0) 
    endif
!    Call Vn0_BCs_SW  ( 0.0)      !      SW necessary for normalizing density  *** deactivated 7 march ***reactivated for FS Nov 7 2017
    
!    if (summation_density) then  !     -TMP delete after debug -----  If continuity is used, then get drho                      
!       call sum_density       
!    endif
    
endif    

t0 = time_sph; vx0 = vx; x0 = x; u0 = u
!if (.not.summation_density) rho0 = rho
rho0 = rho				           ! **chnged for pwp

DO i = 1,4                                         ! ------    RK loop 

   t  = t0  + f1rk(i)*dt_sph*dt_sph_factor         ! ** MP semiimplicit
   do ipoin = 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      x (:,ipoin)  = x0(:,ipoin) &
                     + f1rk(i)*(dt_sph*dt_sph_factor)*vx(:,ipoin)
      vx(:,ipoin) = vx0(:,ipoin) & 
                     + f1rk(i)*(dt_sph*dt_sph_factor)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      if (if_density_comput.eq.1) then    !  ** MP semiimplicit not for FStep
         if (.not.summation_density) then                  
            rho(ipoin) = rho0 (ipoin)   + f1rk(i)*(dt_sph*dt_sph_factor)*drho(ipoin)
         endif
      endif
      u   (ipoin) = u0 (ipoin)   + f1rk(i)*(dt_sph*dt_sph_factor)*du(ipoin)     ! ** MP FS
!      
!      ------   
!
        if (SPH_ProblemType.eq.1.and.nuaux.eq.0) then                           !   SW 
           if     ( u(ipoin).gt.rho(ipoin) ) then 
               u(ipoin) = rho(ipoin)   ! ** CHK pwp
           elseif ( u(ipoin).lt.0.0 ) then      
               u(ipoin) = 0.0 
           endif
        endif 

    enddo
   
   Call Check_Out_Domain 
   
!      ------   Get alpha correction factor

   if (pa_sph.eq.3.or.pa_sph.eq.4) then                    !      Alpha correction factor JB                                                                
       call Get_AlphaJB
   endif   

   if (if_density_comput.eq.1) then     !  ** MP FStep
      if (summation_density) then       ! ------  If continuity is used, then get drho    
         call sum_density         
      else             
         call con_density       
      endif
   endif
   
   if (SPH_ProblemType.eq.1) then           !     SW

       if (if_IntForces_SW.eq.1) then          !       ** MP not computed in FS
          call IntForces_SW                    !       Internal Forces SW.. Gets hw at s nodes      
          if (i.eq.1.and.allocated(auxS1)) then!
              call Get_s_at_w_DF_SW            !     transfer initial value of porosity to auxS10
              auxS10 = auxS1
          endif 
       endif   
       
       call ExtForces_SW (0)                !       External Forces SW (slope, friction,virt.parts)
                                            !       March 2014 0 means standard, get dh/dt
   elseif (SPH_ProblemType.eq.2) then       !     NS
!       call IntForces_NS                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_NS                   !       
   elseif (SPH_ProblemType.eq.3) then       !     DF
!       call IntForces_DF                   !
!       if (visc_artificial) call art_visc  !       Artificial viscosity 
!       call ExtForces_DF                   !
   else
       write(*,*) ' wrong type of problem= ', SPH_ProblemType
       pause
   endif
   
   DO ipoin = 1, npoin      
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      Dvx_RK(:,ipoin) = Dvx_RK  (:,ipoin) + f2rk(i)*(exdvxdt(:,ipoin)+indvxdt(:,ipoin)+ardvxdt(:,ipoin))
      Dx_RK (:,ipoin) = Dx_RK   (:,ipoin) + f2rk(i)*vx(:,ipoin)
      if (if_density_comput) then 
          if (.not.summation_density) then
             Drho_RK (ipoin) = Drho_RK (ipoin)   + f2rk(i)*drho(ipoin) !       
          endif
      endif
      Du_RK   (ipoin) = Du_RK   (ipoin)   + f2rk(i)*du  (ipoin)
   ENDDO
   ipoin = 1 !  just for debug
   
ENDDO

DO ipoin =1,npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   x  (:,ipoin) =  x0 (:,ipoin) + (dt_sph*dt_sph_factor/6.)* Dx_RK (:,ipoin)       
   vx (:,ipoin) = vx0 (:,ipoin) + (dt_sph*dt_sph_factor/6.)* Dvx_RK(:,ipoin)
   if (if_density_comput) then   !  ** MP fract Step
      if (.not.summation_density) then
         rho(ipoin)  = rho0(ipoin) + (dt_sph*dt_sph_factor/6.)*Drho_RK(ipoin)        
      endif
   endif
   u  (ipoin)   =  u0   (ipoin) + (dt_sph*dt_sph_factor/6.)* Du_RK (ipoin)
ENDDO

Call Vn0_BCs_SW ( dt_sph*dt_sph_factor/6. )

Call Check_Out_Domain                                            

if (SPH_ProblemType.eq.1) then            

   call Get_s_at_w_DF_SW        !    get hw at soil pics and viceversa
   
   Drho_RK = rho - rho0			!    we will use drho_rk here to avoid increasing memory
   
   if (nUaux.GE.1) then
      
      if (allocated(auxS1)) then
         dauxS1 = auxS1 - auxS10
         call Get_Pwp_SW (Drho_RK, dauxS1)   !    compute pwp DFs and pwp
      else 
         call Get_Pwp_SW (Drho_RK)		     !    compute pwp
      endif
      
      do ipoin = 1,npoin
	     u(ipoin) = UAux(1,ipoin)
      enddo
      
!      u (:) = Uaux (1,:)
   endif
             
    call ExtForces_SW (1)                !       External Forces SW (slope, friction,virt.parts)
                                         !       March 2014 1 means erosion, once we know v  
                                         !       Replaces old Erosion_SW 

    Call Abs_BCs_SW                      !    remember this should go also to other time integration routines

      
!      ------   
!
   if(nuaux.eq.0) then
     DO ipoin=1,npoin
        if ( if_Out_Domain(ipoin).eq.1) CYCLE
        if ( u(ipoin).gt.rho(ipoin) ) u(ipoin) = rho(ipoin)  ! ** CHK pwp
        if ( u(ipoin).lt.0.0 )        u(ipoin) = 0.0
     Enddo
   endif

   
endif      
                                        
if (average_velocity) then      !     March 2014 Moved here, more logical for erosion and BCs
   call av_vel		            !     Average v of particles to avoid penetration
   do ipoin = npoin0 + 1, npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      vx (:,ipoin) = vx(:,ipoin) + av(:,ipoin)
   enddo   
endif
   
if (sle.ne.0) call h_upgrade            !     neighboring particles & update HSML

!if (SPH_ProblemType.eq.1) then            
!   call Get_s_at_w_DF_SW               ! *** get hw at soil pics and viceversa
!endif                                  ! now done in internal forces ***

ipoin = 1        ! *** to debug

!write(res_file3,1991) time_sph+dt_sph, (Uaux(i,1),i=1,nUaux-1), exp(-3.141593*3.141593*0.40528473*(time_sph+dt_sph)/4)
1991 format (12(f10.4,2x))


!      ------   Clears linked list containing interaction data

!   call clear_interactions

101   format(1x,4(2x,a12))      
100   format(1x,9(2x,e11.4)) 

END subroutine  RK4_FS



!----------------------------------------------------------------------

      subroutine FS_semi_implicit_FEM (itimestep_sph)

!----------------------------------------------------------------------

implicit none

integer(ink) itimestep_sph, ipoin 
integer(ink) type_problem      ! 1 steagy heat  2 transient 3..SW FS

if (itimestep_sph.eq.1) then
   write(*,*) ' type problem type'
   read (*,*) type_FS_problem
endif

!type_FS_problem = 4               !  CHANGE if necessary...21 jan 2021

if (type_FS_problem.eq.1) then     ! ----------------  Steady state heat
    if (itimestep_sph.eq.1) then
       dt_sph = 0.01                  
    endif
    call Second_Step_FS_SW_FEM    !  (i) vstar  (ii) div vstar  (iii) rrhsp
                                  !  (iv) Get Laplacian  (v) Get BCs  (vi) Solve
                                  !  output is p at n+1  
     rho = p                      !  
     call OutputRes_SW            !    
     STOP                         !  
elseif (type_FS_problem.eq.2) then   ! ----------------  transient heat    
     if (itimestep_sph.eq.1) then
         p   = 10.                   !   
         p(1) = 0. ; p(npoin)=0.     !   
         dt_sph = 0.01 
     endif
     call Second_Step_FS_SW_FEM   !  (i) vstar  (ii) div vstar  (iii) rrhsp
                                  !  (iv) Get Laplacian  (v) Get BCs  (vi) Solve
                                  !  output is p at n+1  
     rho = p                      ! 
elseif (type_FS_problem.eq.3) then   ! ----------------  sw
                       
     call RK4_FS (itimestep_sph, 0.5, 0, 1)
     call Second_Step_FS_SW_FEM   !  (i) vstar  (ii) div vstar  (iii) rrhsp
                                  !  (iv) Get Laplacian  (v) Get BCs  (vi) Solve
                                  !  output is p at n+1 
                                  ! do not call grad pn1, we do it inside 2nd step
     call RK4_FS (itimestep_sph, 0.5, 0, 1)   
     
elseif (type_FS_problem.eq.4) then   ! ----------------  sw
                       
     call RK4_FS (itimestep_sph, 1.0, 0, 1)  ! ADDED to check(itimestep_sph, 0.5, 0, 1)
     !x = x0; vx = vx0    !  ADDED to check
     call Second_Step_FS_SW_FEM   !  (i) vstar  (ii) div vstar  (iii) rrhsp
                                  !  (iv) Get Laplacian  (v) Get BCs  (vi) Solve
                                  !  output is p at n+1 
                                  ! do not call grad pn1, we do it inside 2nd step
    !  call RK4_FS (itimestep_sph, 0.5, 0, 1)  ! ADDED to check(itimestep_sph, 0.5, 0, 1)  

elseif (type_FS_problem.eq.5) then   ! ----------------  shao sw
                       
     call RK4_FS (itimestep_sph, 1.0, 0, 1)
     call Second_Step_FS_SW_FEM   !  (i) vstar  (ii) div vstar  (iii) rrhsp
                                  !  (iv) Get Laplacian  (v) Get BCs  (vi) Solve
                                  !  output is p at n+1 
endif 
  
END subroutine FS_semi_implicit_FEM 


!----------------------------------------------------------------------

      subroutine FS_semi_implicit_SFM (itimestep_sph)

!----------------------------------------------------------------------

implicit none

integer(ink) itimestep_sph, ipoin 
integer(ink) type_problem      ! 1 basic

if (itimestep_sph.eq.1) then
   write(*,*) ' type problem type'
   read (*,*) type_FS_problem
endif

if (type_FS_problem.eq.1) then   ! ----------------  shao sw
                       
     call RK4_FS (itimestep_sph, 1.0, 0, 1)
     call Second_Step_FS_SW_SFM    
                                                                     
endif 
  
END subroutine FS_semi_implicit_SFM 

!-------------------------------------------------------------------

       Subroutine Update_Const_SPH (icode, index_mat, Cval, Cval0_mat) 
       
!-------------------------------------------------------------------


implicit none
integer(ink) icode, index_mat
real   (irk) Cval, Cval0_mat
call Update_Const_SW (icode, index_mat, Cval, Cval0_mat) 

End Subroutine Update_Const_SPH 



!-------------------------------------------------------------------

       Subroutine  Stop_case_SPH
       
!-------------------------------------------------------------------


implicit none 

call Stop_case_SW  

End Subroutine Stop_case_SPH 

!-------------------------------------------------------------------

       Subroutine  Get_Prob_Propag_SPH
       
!-------------------------------------------------------------------


implicit none 

call Get_Prob_Propag_SW  

End Subroutine Get_Prob_Propag_SPH 


!-------------------------------------------------------------------

       Subroutine Manage_pixel_Hazard_SPH (iopt)
       
!-------------------------------------------------------------------

implicit none 

integer (ink) iopt

if (iopt.EQ.1) then
   call Get_pixel_hVars_SW
elseif (iopt.eq.2) then
   call Get_pixel_FOSM_SW 
elseif (iopt.eq.3) then
   call Plot_pixel_RES_SW 
endif

!  ---  calls for mcarlo...to be used wiTh 1 pic

if (iopt.EQ.11) then
   call Get_MC_pixel_Vars_SW
elseif (iopt.eq.12) then
   call Norm_MC_pixel_Vars_SW (1)   !  changed to 1 always print every problem
   call Plot_MC_pixel_Vars_SW (1)
elseif (iopt.eq.13) then
   call Stats_MC_pixel_Vars_SW
elseif (iopt.eq.14) then
   call Plot_MC_pixel_Vars_SW
endif

!if (iopt.eq.999) then  !  ********************  DEBUG
 9999 format (i4, 9(2x, f10.4))
!    write(34,9999) nproblems_index , dt_sph, time_sph, t_plot_reset, time_plot,   &
!                        x(1,9)        ,vx(1,9)         , rho(9)    
!                       x(1,9), x(2,9),vx(1,9), vx(2,9), rho(9)    
! endif

! if (iopt.eq.998) then  !  ********************  DEBUG
!    write(34,*) ' ******************************' 
!    write(34,9999) nproblems_index , dt_sph, time_sph, t_plot_reset, time_plot,   &
!                          x(1,9), x(2,9),vx(1,9), vx(2,9), rho(9)    
!  endif

End Subroutine Manage_pixel_Hazard_SPH

!-------------------------------------------------------------------

       Subroutine Close_files_SPH
       
!-------------------------------------------------------------------

!      ------    

implicit none

call Close_files_SW

End Subroutine Close_files_SPH 

!-------------------------------------------------------------------

       Subroutine deallocate_list_SPH
       
!-------------------------------------------------------------------

!      ------    

implicit none

call deallocate_list_Tools            !  MP Aug 2022 to avoid increase of size in multicases

End Subroutine deallocate_list_SPH


!-------------------------------
END MODULE SPH_MAIN_2019


 