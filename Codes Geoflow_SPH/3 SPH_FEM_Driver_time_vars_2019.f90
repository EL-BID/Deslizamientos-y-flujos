!--------------------------------------------------------


                 MODULE SPH_FEM_Driver_time_vars_2019  
                 

!--------------------------------------------------------

!                version 2006 02 22...time stepping
!                included if out domain last 27 09 2007

                       
USE SPH_FEM_variable_types_2019

implicit none

integer (ink) MASTER_dat, MASTER_chk  	! Files with global information
integer (ink) MASTER_cases, MASTER_res  ! input cases loop info and output

character(16) SPH_SW_problem_name       ! name of SW problem being solved.  CHG MP 27 March 2020
character(16) GFL_SW_problem_name       ! name of SW problem being solved. Moved to Driver vars MP 1 dic 2021 
character(16) problem_name_Global       ! name of problem being solved Brought from time vars MP sept 2022

integer(ink) if_sph, if_gfl, if_tgf     !  1 if we will use sph and geoflow, zero otherwise

real   (irk) time, dt, initial_dt		        !  Saeidmt 19Oct2023    ! time increment, time, 
real   (irk) time_sph, dt_sph, initial_dt_sph   !  Saeidmt 19Oct2023	! local time and dt for sph

real   (irk) time_gfl, dt_gfl	! local time and dt for gfl
real   (irk) time_tgf, dt_tgf	! local time and dt for tgf

real   (irk) time_sph_restart   ! time at which we have restarted SPH computation

integer(ink) current_ts		! absolute number of time step. Accounts several restarts
integer(ink) nstart		    ! number of timesteps in previous cycle
integer(ink) itimestep		! local time step in present restart cycle
integer(ink) itimestep_sph	! same for sph (there is an inner loop)
integer(ink) itimestep_gfl	! same for gfl (there is an inner loop)
integer(ink) itimestep_tgf	! same for gfl (there is an inner loop)

real   (irk) dt_limit_inf	! minimum allowed dt, set to 1.e-7*dt
real   (irk) time_end		! max time for computation
integer(ink) maxtimestep	! maximum number of timesteps in computation

integer(ink) ic_adapt_dt_sph		! controls SPH adapt. t stepping (0...no     1...yes)
integer(ink) ic_adapt_dt_gfl		! controls GFL adapt. t stepping (0...no     1...yes)
integer(ink) ic_adapt_dt_elem_gfl	! dt local at elements in TG 2 steps
integer(ink) ic_adapt_dt_tgf		! controls tgf adapt. t stepping (0...no     1...yes)
integer(ink) ic_adapt_dt_elem_tgf	! dt local at elements in TG 2 steps

integer(ink) print_step     	!  Print every print_step    (On Screen)
integer(ink) save_step      	!  Save  every  save_step    (To Disk File)
integer(ink) plot_step      	!  Plot  every  plot_step    (For GID v.7)
integer(ink) moni_part      	!  Number of control particle (only one)
integer(ink) isearch_step       !  auxiliar


real   (irk) time_plot		 ! we plot every time_plot time
real   (irk) t_plot_reset	 ! time counter reset to zero after plotting
real   (irk) time_save		 ! we save every time_plot time
real   (irk) t_save_reset	 ! time counter reset to zero after saving
real   (irk) time_print		 ! we print every time_print time
real   (irk) t_print_reset	 ! time counter reset to zero after printing
real   (irk) t_casewin       ! same when considering cases output from windows
real   (irk) t_casewin_reset !           
!   --------------------  Time curves info for geoflow. Sph uses hist.
 						
integer(ink) ntcurves 				        !     number of time curves (0 if none)
integer(ink) mptstcurves		            !     max number of time pts in all t curves
integer(ink),allocatable:: nptstcurves(:)	!     nr. of time pts in each curve
real   (irk),allocatable::   ttcurves(:,:)	!     time values for all time curves
real   (irk),allocatable::   ftcurves(:,:)	!     factors at all times. 
						                    !     Multiply a0 in prer variable
						
real   (irk)    T_change_to_W	!  in WIR computations, we changed algorithm after this
						        !  used in sph_main module, read in SW module shared 
						        !           trough this driver-time_vars module

integer(ink) SPH_t_Integ_Alg ! Type of integration for sources 
                             ! Source_Integration = 0   Central, original version
                             !                    = 4   Runge Kutta 4t
                             !                    = 11  TG 2 steps
                             !                    = 6   semi implicit RK  ** MP s_implicit 17.03.18
                             !                          calls RK_FS_semi_implicit
                             !                    = 7   semi implicit RK  ** FE
                             !                    = 10  4 DM, stops intforces until t > T_end_solid
                             !                    = 6   semi implicit RK  ** MP s_implicit 17.03.18
                             !                    = 100  random walk 
real   (irk) T_end_solid
                             !
                             ! **** comes from source_integration *** TAKE CARE ThB codes

integer(ink) SPH_Problemtype !  1 for SW, 2 for NS, 3 for DF, > 10 for TGSolid 						!       trough this driver-time_vars module

integer(ink) GFL_t_Integ_Alg ! Type of integration for sources 
                             ! Source_Integration = 0   Central, original version
                             !                    = 2   Runge Kutta 4th 
integer(ink) GFL_Problemtype !  1 for SW  (not used yet) 
  
                                                     !     if or not cases loop, nr of sets of parameters,                   
integer(ink) ic_cases_mat , ncases_mat, ninfo_mat    !     nr of info data: index, initial val, nvals, min, max   
real   (irk),allocatable::  Ccases_mat (:,:)         !     (ncases_mat, ninfo_mat)

integer(ink) nproblems           !  total number of problems to solve: 
integer(ink) nproblems_index     !  active case in list of cases * values
integer(ink) ndat_out            !  how many data columns in csv 
integer(ink) Mdat_out            !  max data columns in csv 
integer(ink) ic_stop_cases       !  when 1 code will stop and proceed to next case
real   (irk) vstop_cases         !  the velocity value to stop computation 
real   (irk) timestop_cases      !  the velocity value to stop computation
real   (irk) dtimestop_cases     !  the velocity value to stop computation
real   (irk) xmax_cases          ! to chk 1d cases
real (irk), allocatable::  OC_mat(:,:)   ! for all cases info for csv file
                                         ! (nproblems, ndat_out: 3+ncases_mat+2 (vectors)+ nwin*(6)+ 1dist to target))
integer(ink), allocatable::  OC_index(:,:) ! for all problems info indexes for variables
                                         ! (nproblems, ncases_mat)  made integer MP 15 March 2022

real (irk), allocatable::  List_Vals (:,:) ! list of values for each case   2022 03 07
real (irk), allocatable::  List_Probs(:,:) ! probabilities for each case
                                           ! (ncase_mat, Nvals_max)
integer(ink) Nvals_max                     ! above: max nr of vals in all cases

integer(ink) ic_FOSM                       ! to 1 when FOSm  MP 25th July
integer(ink) ic_MCarlo                     ! set to 1 when ic_cases_mat=3. Then Ic_cases_mat -> 2
                                           !  MP 15th March 2022
integer(ink) ic_cases_win                  !  1 when we are writting vars at windows (Qpeak, hpeak, etc)

integer(ink) ic_Trigger                    !  set to 1 ehen SPH_t_integ_alg = 200  MP 18th 03 2022

character(16) TOPO_problem_name        ! name of problem being solved in Topo. CHGD MP March 2020.
integer (ink) if_Read_TopoDat          ! set to 1 after reading topo.dat file

integer (ink) nprob_props                         !  MP 24 03 2022 nr of prob parameters in a law
integer (ink), allocatable:: probab_law       (:) ! (ncases_mat)  for each magnitude, law
real    (irk), allocatable:: Param_prob_law (:,:) ! (ncases_mat, nprob_props)  table of props
real    (irk), allocatable:: Tmp_List         (:) ! (ncases_mat, nprob_props)  table of props

integer (ink) np_h, ne_h                 ! histograms 1D points and elems. Values at elems
real    (irk), allocatable:: x_h (:)     ! (np_h)  magnitude values at nodes
real    (irk), allocatable:: p_h (:)     ! (ne_h)  prob values at elems = n/Ntotal,or n
real    (irk), allocatable:: F_h (:)     ! (np_h)  Accum values at nodes, integral p_h

integer (ink) n_GZones                   ! MP April 2022  geo zones for trigger ic_cases = 5
integer (ink) n_MCarlo_problems          !  nproblems = nMCarlo_pr * n_GZones   ic MCarlo = 5

real    (irk), allocatable:: FOSMms(:,:) ! (1+2*N_GZones,ncases_mat)  table of props
                                         ! columns are cases
                                         ! row 1 index, (row 2 mean, row 3 sigma) for all zones

real    (irk), allocatable:: FOSM_Pr2   (:,:) !  Probability propagation (2,4)  MP 30th May 2022                                         
                                              !  files mu sigma, columns dist, hpk, vpk, Qpk 
real    (irk), allocatable:: MCarlo_Pr2 (:,:) !  Probability propagation (2,4)  MP J09 une 2022  
                                                                                       
integer (ink) ic_Prob_Propag               !  MP 4 June 2022 if 1 propag + Fosm
integer (ink) ic_xtp                     ! 1 if prob(x,t) activated else 0
integer (ink) ma_xtp, mh_xtp             ! number of values of PGA and alphas to compute time probability
integer (ink) Nyra_xtp, Nyrh_xtp         ! nr of yearsfor which Prob [Ai] is computed
real    (irk), allocatable:: hrs_xtp (:), Acc_xtp (:)    ! values of Acc and alphas
real    (irk), allocatable:: Trh_xtp (:), Tra_xtp (:)    ! values of return periods
real    (irk), allocatable:: PAi_xtp (:), Phi_xtp(:)     ! probabilities P[Ai]

integer (ink) ic_pixel_Hazard            ! MP 21 July 2022  for FOSm and Propag using pixeles if set=1
                                         ! code does it when ic_cases_mat=21. Line 358-366
                                         ! Then back to 20 and ic_pixel hazard set to 1    
integer (ink) npixel_hVars               ! nr of vars to characterize hazard at pixels MP 23 July 2022

TYPE Aux_H_pixel                                         ! stores info for hazard at ipoigs
   real    (irk), dimension (:,:), pointer:: OCMM        ! same info than OC_mat columns 3+ncases+2 on  
   real    (irk), dimension (:,:), pointer:: FOSM_Pr     ! same than FOSM_Pr2 (2, npixel_hVars) h,v,q 
   integer (ink), dimension (:),   pointer:: nr_pts_pix  ! different for each problem            
   integer (ink)                          :: if_alloc    ! the only way to allocate ointer dims w/o asking 
   integer (ink)                          :: if_all_FOSM ! 
END TYPE Aux_H_pixel

TYPE (Aux_H_pixel),  allocatable:: FOSM_pixel (:)
                                                             ! 3rd August 2022 MP 
TYPE Aux_MC_pixel                                            ! stores info for hazard at ipoigs
   real    (irk), dimension (:)  , pointer:: OCMC            ! same info than OC_mat columns 3+ncases+2 on  
   real    (irk), dimension (:,:), pointer:: MC_Pr           ! (2, npixel_hVars) stores mean and sigma 
   integer (ink)                          :: nr_pts_pix      ! different for each problem 
   integer (ink)                          :: nr_problem_pix  ! different for each problem 
   integer (ink)                          :: if_OCMC_alloc   ! the only way to allocate ointer dims w/o asking 
   integer (ink)                          :: if_MC_Pr_alloc  ! 
END TYPE Aux_MC_pixel 

TYPE (Aux_MC_pixel),  allocatable:: MC_pixel (:)

integer(ink) ic_vuln                       ! set to 1 when performing vulnerability analysis  MP Sept 2022 

                                           !  MP Jan 2023

integer(ink) ic_search_step, nsearch_step  ! 1 if we search only every nsearch_steps... vels are very small
integer(ink) ic_time_series              ! 1 if a cuve for rain intensity is provided  
real   (ink) tscale_code, tscale_series  !  13 th march Jw MP SM 

TYPE ts_info                                     
   character(60)                         :: description    ! 60 char describing time series
   integer (ink)                         :: np_ts          ! nr of points in time series   
   real    (ink)                         :: unit_ts        ! seconds in time unit (60, 3600, 86400) 
   real    (ink)                         :: type_Val       ! 1 Rain 2 basal Pwp abs (Pa) Pwp rel (0-1) 
   real    (irk),dimension (:), pointer  :: time_ts        ! list of times in days      
   real    (irk),dimension (:), pointer  :: val_ts         ! list of values
   integer (ink)                         :: if_pwp         !  1 means pwp and rain
   integer (ink)                         :: if_drain       !  1 means drainage of upper
                                                           !  2   drainage lower layer 
   real    (irk)                         :: pwp_presc_basal!  prescribed pwp at base [L] type_val=2 MP Jw SM 13th March      
   real    (irk)                         :: h_basal        !   h of basal layer   
   real    (irk)                         :: poros_basal    !   porosity of basal layer  
   real    (irk)                         :: tanTh_basal    !  avg slope for drainage    
   real    (irk)                         :: L_basal        !   avg length of basal layer  
   real    (irk)                         :: kw_basal       !   horz perm code units1 of draining layer  

   real    (irk)                         :: h_upper        !   h of upper layer   
   real    (irk)                         :: poros_upper    !   porosity of upperl layer  
   real    (irk)                         :: tanTh_upper    !  avg slope for drainage    
   real    (irk)                         :: L_upper        !   avg length of upper layer  
   real    (irk)                         :: kw_upper       !   horz perm ms-1 of draining layer                                   ! 
END TYPE ts_info 

TYPE (ts_info) :: ts_Rain
                                         ! MP 5 March 2023
integer (ink) ic_slow                    ! for slow formulation  
integer (ink) ic_FS,   ic_DF             ! indicate wheter or not (1,9) is a FS or a DF                                            

END MODULE SPH_FEM_Driver_time_vars_2019
 
