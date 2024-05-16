!-------------------------------------------------------


                MODULE SPH_SW_Vars_2019  


!-------------------------------------------------------

!         This module includes all variables used by module MODULE SPH_2008_SW
!         August 2016 Pwp_Avgd and dPwp_Avgd  averaged pwps to add in DF to hydrostatic
!          ...jan 2018 crushing

USE SPH_FEM_variable_types_2019

character(72) text

integer(ink) :: icpwp                           ! -- Control of PWP 
real   (irk) :: pwp0rel

integer(ink)    ic_K0                           ! --  If 1 we will use active and pasives K0 in internal forces
real   (irk)    K0_active, K0_pasive
!      ------

real   (irk), allocatable:: Z00   (:)           ! --  Initial Z coordinate of particles. Used in plotting routines
real   (irk), allocatable:: rho00 (:)           ! --  Initial h  of particles. Used in plotting routines

real   (ink),allocatable:: const(:)             !   Mat constants. In versions up to 2014, 15. Now 45
integer(ink):: nconst                           !   dimensioned in SW module

integer(ink):: nfric, if_frictional, if_V2byR   !   MP 7th december 2017 to include new frict las mu(I)
                                                !   nfric is the basal tau law, if_frictional = 1 if frictional,
                                                !   if_V2byR if centripetal forces added
integer(ink):: NslMasses                        ! -- Sliding masses info
integer(ink), allocatable::Pin_slm(:)           !    Number of Points in sliding masses
real   (irk), allocatable::T0_slm(:)            !    Triggering instant

!      =======  For soil-water interactions: sph waves in reservoirs and 2 phase DFs =============
                                
integer(ink)                ic_SWAlg            ! --  0 is h, 1 is (h2-Z2) classical for SW computations
                                                !     2 tmp value for Drum
                                                !     3 tmp for DF 2 phases

integer(ink)                ic_semi_implicit    !     1 if semi implicit SW algorithm ** MP s_implicit 17.03.18 
                                                !      SPH_t_Integ_Alg=6 (MAIN SPH Vars) makes ic_semi_implicit =1  
                                                
integer(ink)                ic_ws_Interact      !     0 no interaction (when npoin_s and npoin_w >0= 1 interact
                                                !     1 waves in reservoirs (when npoin_s and npoin_w >0= 1)
                                                !     2 2 phase DFs

integer(ink) :: ic_WIR                          ! --  Always 1 Obsolete

real   (irk), allocatable:: hs_at_W      (:)    !     height of soil at water nodes
real   (irk), allocatable:: gradhs_at_W(:,:)    !            and its gradient

real   (irk), allocatable:: hs_plus_Z_topo(:)   !     height of soil at topo mesh
real   (irk), allocatable:: hs_plus_Z_w   (:)   !     hs  at water points Dim is npoin!!
real   (irk), allocatable:: gradhs_plus_Z_w(:,:)!         and its gradient

real   (irk), allocatable:: h_DF (:)            !     hs + hw in DF hs=(1-n)h_DF hw=n.h_DF
real   (irk), allocatable:: n_DF (:)            !     depth averaged porosity
real   (irk), allocatable:: h_sw (:)            !        1st npoin_s posn hw at s pts, 
                                                !        2nd npoin_w posn hs at w pts
real   (irk), allocatable:: v_sw (:,:)          !        1st npoin_s posn hw at s pts, 
                                                !        2nd npoin_w posn hs at w pts

real   (irk), allocatable:: Pwp_avgd (:)        !     averaged pwp to add to hydr.component in DFs

!      ===================================================================================================                                                 

!      ------  for ABS BCs
                                                                
integer(ink)  ic_abs, nabs                              ! --  Info for absorbing BCs  Just for water, sub critical (Froude < 1.)                                                                
integer(ink), allocatable:: labs(:)                     !     Ic_abs=0 no Abs BCs, Ic_abs=1 yes,  nabs number of abs nodes                              
real   (irk), allocatable:: Zext_abs(:),  vn_abs(:,:)   !     labs(nabs)  List of abs nodes             
real   (irk), allocatable:: x_abs(:,:), hsml_abs(:)     !     Absolute Z of smwl at exterior of computational domain

!      ------  for Vn=0 BCs
                                                                
integer(ink)  ic_vn0, nvn0, max_nvn0                    ! --  Info for vn=0 BCs***Chgd MP 28th Nov 2022 set max pts in walls               
integer(ink), allocatable:: lvn0(:)                     !     Ic_vn0=0 no  BCs, Ic_vn0=1 yes,  nvn0 number of vn0 nodes                         
real   (irk), allocatable:: vn_vn0(:,:)                 !     lvn0(nvn0)  List of vn0 nodes             
real   (irk), allocatable:: x_vn0(:,:), hsml_vn0(:)     !              
integer(ink), allocatable:: type_vn0(:)                 !  type_vn0(nvn0)  1 ..not perm  2 perm to water ***Chgd MP 28th Nov 2022   


!      ------  If we use PWP FD meshes for points, here are BC's same for all points

integer(ink)  BC_pwp_Type(2)
real   (irk)  BC_pwp_Val (2)

!      ------  If we use PWP FD meshes with erosion (Stefan problem)

integer(ink)  BC_Adv_Type(2)
real   (irk)  BC_Adv_Val (2)
real   (irk)  safetpwp         !  safety factor multiplying pwp should be >1

!      ------  for wir points in water and in soil

integer(ink)  ic_water, npoin_w 
real   (irk), allocatable:: x_w  (:,:), rho_w(:), mass_w(:), hsml_w(:)
real   (irk), allocatable:: vx_w (:,:), p_w  (:), u_w   (:)
integer(ink), allocatable:: itype_w(:)

integer(ink)  ic_soil, npoin_s
real   (irk), allocatable:: x_s  (:,:), rho_s(:), mass_s(:), hsml_s(:)
real   (irk), allocatable:: vx_s (:,:), p_s  (:), u_s   (:)
integer(ink), allocatable:: itype_s(:)

!      ------  for Virtual Particles  

integer(ink)  ic_vps, npoin_vps 
real   (irk), allocatable:: x_vps  (:,:), rho_vps(:), mass_vps(:), hsml_vps(:)
real   (irk), allocatable:: vx_vps (:,:), p_vps  (:), u_vps   (:)
integer(ink), allocatable:: itype_vps(:)

real   (irk)  vp_r0, vp_D, vp_n1, vp_n2                 ! --  Info for virtual particles (dx/2, 5gh, 12 , 4)

!      ------  for histogram input

integer(ink) ic_histogram               !   1 if there are histograms
integer(ink) n_hist                     !   nr of histograms
integer(ink) type_of_histogram          !   type of histog  1...h   2...v or q
integer(ink) mp_hist                    !   max nr of points in histograms
integer(ink), allocatable:: np_hist(:)  !   nr of pts in each histogram
real   (irk), allocatable:: t_his(:,:)  !   time series for each histogram
real   (irk), allocatable:: h_his(:,:)  !   h series
real   (irk), allocatable:: v_his(:,:)  !   v or q series

!      ------  for gates definition. Each gate is associated to a time curve

integer(ink) npoin_in_model             !   number of nodes active in the model
integer(ink) np_bakG_Total              !   number of nodes behind all gates
integer(ink), allocatable:: np_bakG(:)  !   number of nodes behind each gate

integer(ink)                ic_sgates                          !   1 if there are gates
integer(ink)                n_sgates                           !   number of gates
integer(ink)                mp_gate, ml_gate                   !   max nr of (i) points in gate (ii) lines to be output
integer(ink), allocatable:: np_gate(:), nl_gate(:) 
real   (irk), allocatable:: xl_gate(:)
real   (irk), allocatable:: x_gate(:,:,:), nv_gate(:,:)         !   --- the gate (line) has points where particles are entered
integer(ink), allocatable:: np_track(:), ndof_track(:,:)        !   --- nr of points exited and list with them
real   (irk) hsml_gfact  

!   New coarse mesh for plotting  -------------------------------------
!
!     if_coarse_plt = 1                 1 means we activate it
!     fc_coarse_plt                     is a factor which multiplies deltxg and deltyg in topo mesh
!     npoi_plt                          points in new mesh
!     xmin_plt...ymax_plt               inscribing rectangle
!     deltx_plt, delty_plt              new spacings
!     npoix_plt, npoiy_plt              number of grid lines in x and y
!     coor_plt (3,npoi_plt)             stores x, y and Z
!     poi_in_model_plt(npoi_plt) = 1    point in computational domain

integer(ink) if_coarse_plt
real   (irk) fc_coarse_plt
integer(ink) npoi_plt          !  when it works, take out all PUBLIC here
integer(ink) npoix_plt, npoiy_plt
real   (irk) xmin_plt, xmax_plt, ymin_plt, ymax_plt
real   (irk) deltx_plt, delty_plt
integer(ink), allocatable:: poi_in_model_plt(:)
real   (irk), allocatable:: coor_plt(:,:)

!      ------------------------------------------

!    output at selected points

  
integer(ink) if_chk_pts, nchk_pts                !  if_chk_pts = 1 if there are points to output
                                                 !  nchk_pts   number of points to be output
real   (irk), allocatable :: coor_chk_pts(:,:)   !  coordinates stored here: ndimn x nchk_pts

!      ------------------------------------------

!      Filter for GID

integer(ink)                nGidMaskSw          !   We will use a mask to plot variables, nGidMask = 20
integer(ink), allocatable:: Gid_Mask_SW(:)      !   Gid_MAsk(nGidMask) can be zero (no plot) or one (plot)
                                                !   SW case:
                                                !          
                                                !   1 ... h soil       
                                                !   2 ... displacement          
                                                !   3 ... velocity          
                                                !   4 ... Pwp          
                                                !   5 ... erosion          
                                                !   6 ... Z         
                                                !   7 ... h soil relative  (used in Lausanne)         
                                                !   8 ... h water          
                                                !   9 ... eta          
                                                !   10 ... hs + hw
                                                !   11 ... hsml (if -1 Froude)
                                                !   12 ... PWP full
                                                
!      ----------  In SW module we limit inf(h)

real   (irk) h_inf_SW
                                                
!      ----------  For drum

real   (irk) R_drum, w_drum_radps               !    Radius and rads per second

!      ------  Inclined planes: gravity components

real   (irk) gravx, gravy, gravz

!      ------   crushing

integer(ink) ic_crush                     !  default 0 , 1 if activated  
                                                 
real   (irk), allocatable :: r_crush(:)   ! current representative radius of grains 
real   (irk), allocatable :: w_crush(:)   ! current crush work 

!      ------------------------------------------

!    solid behaviour via no int forces
 
integer(ink) ic_end_solid, law_end_solid, aG_end_only
real   (irk) t1_end_solid, fct_exp_end_solid

!      ------------------------------------------

!    Solid water interfaces in WIR and submarine landslide problems
 
integer(ink) ic_ws_erosion                   !  set to 1 if there is interface erosion 
real   (irk), allocatable :: n_wir   (:)     !  porosity of soil in WIR. Used for erosion interface
real   (irk), allocatable :: dens_wir(:)     !  density dependent on porosity of soil in WIR. Used for erosion interface
 
integer(ink) ic_ws_friction                  !  set to 1 if there is interface friction 

!      ---------------------------------------

!      erosion control, including Stefan pwp FDs ALE

integer(ink) ic_erosion
integer(ink) ic_Stef_eros                               ! will be set to 1 when pwp+FD and erosion
integer(ink) ic_Activ_Pwp_Stef_eros                     ! 1 active, 0 deactivated. Read in Initialize PWP
integer(ink), allocatable :: Activ_Pwp_Stef_Eros_Pt (:) ! same for all pics
real   (irk), allocatable :: erosion_rate_V (:)         ! array with erosion rates at ipoin

!      ---------------------------------------

!     restart 

integer(ink) ic_restart

!      ------------------------------------

!    DF with variable saturation, Upper layer(1-a)h Sr lower a.h sat

integer(ink) ic_hrelSat
real   (irk), allocatable :: hrelSat_DF(:)

real   (irk), allocatable:: fact_dens_DF(:)   
real   (irk), allocatable:: dens_DF (:), densw_DF(:), dens1_DF(:), dens2_DF(:)  
real   (irk), allocatable:: dens_bar_DF(:), dens_bar_eff_DF(:)
real   (irk), allocatable:: fact_pwp_max_DF (:), fact_pwp_min_DF (:)
real   (irk), allocatable:: fpwp_max_DF     (:), fpwp_min_DF     (:)

real   (irk), allocatable:: basal_flux (:)   !  MP CHGD Aug 2021 Stores fluxes at basal grids

!      ------   Shiomi with or w/0 pwp control   MP 30th June

integer(ink) ic_gradn      !  if 0 gradn computations are deactivated.
                           !   Porosity changes not taken into account

!      --------   Info for windows when considering cases

integer (ink) nwin_vars   !  nr of variables in windows. 20th Luly MP **

TYPE window_case_info
     real   (irk) x0 (2)
     real   (irk) ab (2)
     real   (irk) Theta
     real   (irk) nva (2)
     real   (irk) nvb (2)
     integer(ink) is_target !   when 1, the window is of target type, to assess vulnerab
END TYPE

TYPE (window_case_info), allocatable :: info_cw(:) 
integer (ink) ncase_win 

integer, allocatable::    ic_T0_cw (:)   ! used to determine start and end of aval. in window
real (irk), allocatable::       T0_cw (:)   ! initial time of passing at window
real (irk), allocatable::       Tf_cw (:)   !   final time of passing at window
real (irk), allocatable::    qpeak_cw (:)   !  peak discharge
real (irk), allocatable::   vnpeak_cw (:)   !  peak normal v to section
real (irk), allocatable::    hpeak_cw (:)   !  peak h
real (irk), allocatable::     Qtot_cw (:)   ! total volume trough section
real (irk), allocatable:: dstminpeak_cw (:) ! min distance to avalanche

!      ------  parameters for TRIGGER    MP 18 March 2022

integer (ink) N_Slopes                            !  number of slopes to be considered in a topo model
integer (ink) N_Stats                             !  nr of props of the distributio to keep
real    (irk) Slope_Max_cut, Slope_Min_cut        !  cut values which will limit values from DEM
real    (irk), allocatable::    Slope (:)         !  list of N_solpes slopes
real    (irk), allocatable::    FoS_Slope (:,:)   !  (nproblems, N_Slopes) FoS for each problem and slope
real    (irk), allocatable::    FoS_Stats (:,:)   !  (N_Slopes,  N_Stats ) FoS prob properties for each slope
real    (irk), allocatable::    FoS_DEM   (:,:)   !  (Npoig,     N_Stats ) prob properties at each DTM point

!      ------  for pts generation type 8  (same than 7) read h_expans, mass_expans, mass_shift 

real    (irk)  h_expans, mass_expans  
real    (irk), allocatable:: mass_shift (:)   
integer (ink)  ic_geomm_pts     ! set to 1 when activated

!      ------  water injection into soil in Flow Slides  MP 21st July 2023

integer (ink)  ic_w_inflow     ! set to 1 when activated
real    (irk), allocatable::  ns_dep_props(:,:)  ! soem 2 props at every point dep in (1-n)

TYPE w_inflow_info
     integer(ink)  flow_type
     real   (irk)  Zmin, Zmax  
     integer(ink)  flow_law
     real   (irk)  flow_q 
END TYPE

TYPE (w_inflow_info):: w_inflow

END MODULE SPH_SW_Vars_2019
