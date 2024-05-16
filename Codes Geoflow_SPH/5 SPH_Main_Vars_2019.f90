!--------------------------------------------------------


                 MODULE SPH_Main_Vars_2019    


!--------------------------------------------------------


!                version 2008 08 22...includes TGSolid
!                from tb 14-01-2011 large def
!                included auxS1, auuxS10, dauxS1(:) scalar 

USE SPH_FEM_variable_types_2019
implicit none

TYPE pairs
     integer (ink) niac                             ! Used in a linked list describing part interactions
     integer (ink) Pint_type			            ! Water in reservoirs: 0 same material 1 soil and water 
     integer (ink) pair_i, pair_j                   ! two linked particles
     real    (irk) w                                ! Weight function and
     real    (irk), dimension(:), pointer :: dwdx   !         derivatives
     type    (pairs),             pointer :: next   ! Pointer to next interacting pair of particles 
END TYPE pairs

type (pairs), pointer :: current, last              ! List described by pointers to current and last of list
type (pairs), pointer :: keepit00, keepit11         ! housekeeping 26 Aug 2022 MP  only for find test

real   (irk), allocatable:: x0       (:,:)   !  initial coordinates of particles. Updated at every tn 
real   (irk), allocatable:: x00      (:,:)   !  initial coordinates of particles. Used for relative displacements 
real   (irk), allocatable:: x        (:,:)   !  coordinates of particles 
real   (irk), allocatable:: dx       (:,:)   !  dx = vx = dx/dt  
real   (irk), allocatable:: vx0      (:,:)   !  velocities at time n 
real   (irk), allocatable:: vx       (:,:)   !  velocities of particles 
real   (irk), allocatable:: dvx      (:,:)   !  dvx = dvx/dt, force per unit mass 
real   (irk), allocatable:: av       (:,:)   !  Monaghan average velocity
real   (irk), allocatable:: indvxdt  (:,:)   !  Internal forces at particles
real   (irk), allocatable:: exdvxdt  (:,:)   !  External forces at particles
real   (irk), allocatable:: ardvxdt  (:,:)   !  ... forces at particles
real   (irk), allocatable:: mass     (:)     !  mass of particles 
real   (irk), allocatable:: mass0    (:)     !  mass of particles at time 0 (before any erosion)
real   (irk), allocatable:: rho0     (:)     !  density at time tn   
real   (irk), allocatable:: rho      (:)     !  densities of particles
real   (irk), allocatable:: drho     (:)     !  drho =  drho/dt 
real   (irk), allocatable:: p        (:)     !  pressure  of particles
real   (irk), allocatable:: u0       (:)     !  internal energy at time tn
real   (irk), allocatable:: u        (:)     !  internal energy of particles
real   (irk), allocatable:: du       (:)     !  du  = du/dt 
real   (irk), allocatable:: ahdudt   (:)     !  artificial heat du/dt (energy)
real   (irk), allocatable:: avdudt   (:)     !  artificial visc du/dt (energy) 
real   (irk), allocatable:: hsml     (:)     !  smoothing lengths of particles
real   (irk), allocatable:: c        (:)     !  sound velocity of particles
real   (irk), allocatable:: s        (:)     !  entropy of particles
real   (irk), allocatable:: e        (:)     !  total energy of particles 
real   (irk), allocatable:: t        (:)     !  Temperature
real   (irk), allocatable:: ds       (:)     !  ds  = ds/dtds  = ds/dt
real   (irk), allocatable:: tdsdt    (:)     !  Production of viscous entropy t*ds/dt
real   (irk), allocatable:: eta      (:)     !  kinematic viscosity at particles 
real   (irk), allocatable:: alphaJB  (:)     !  kinematic viscosity at particles 		 

integer(ink)                nUaux            !  we will use a new var Uaux(nUaux,npoin) (Pwp)
real   (irk), allocatable:: Uaux    (:,:)    !      
real   (irk), allocatable:: Uaux0   (:,:)    ! 
real   (irk), allocatable:: dUaux   (:,:)    ! 
real   (irk), allocatable:: dUaux_RK(:,:)    !  we will allocate it if necessary in RK4/WIR

real   (irk), allocatable:: Dx_RK (:,:)	     ! These variables only for RK4 integrator
real   (irk), allocatable:: Dvx_RK(:,:) 
real   (irk), allocatable:: Du_RK   (:)
real   (irk), allocatable:: Drho_RK (:) 

real   (irk), allocatable:: auxS1   (:)      ! Auxiliar variables for all modules     
real   (irk), allocatable:: auxS10  (:)      ! in Sw DFs: n_DF, n0_DF and dn_DF
real   (irk), allocatable:: dauxS1  (:)      !  

real   (irk), allocatable:: other1   (:)     ! variable to plot for testing the model
real   (irk), allocatable:: other2   (:)     ! variable to plot for testing the model


integer(ink), allocatable:: itype    (:)     !  type of particles 
integer(ink), allocatable:: countiac (:)     !  List of particles connected to each one.  
    
integer(ink) if_coupled_pw                   ! 1 = Coupled fractional step  ; 0 = Drained TG

integer(ink) pa_sph          !  SPH algorithm for particle approximation (pa_sph)
                             !      pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
                             !               2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
                             !               3 : uses 2 + alpha correction factor 	 	
integer(ink) nnps            !  Nearest neighbor particle searching (nnps) method
                             !  nnps = 1 : Simplest and direct searching
                             !         2 : Sorting grid linked list
                             !         3 : same than 2 but rotating the cells IXX IYY PXY
integer(ink) sle             !  Smoothing length evolution (sle) algorithm
                             !  sle = 0 : Keep unchanged,
                             !        1 : h = fac * (m/rho)^(1/dim)
                             !        2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
                             !        3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) )
integer(ink) skf             !  Smoothing kernel function 
                             !  skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
                             !        2, Gauss kernel   (Gingold and Monaghan 1981) 
                             !        3, Quintic kernel (Morris 1997)
                             !        4, Wendland kernel     Saeidmt 24Oct2020
integer(ink) nsym            !  Symmetry of the problem 0 : no symmetry, 1 axis, 2 center

logical summation_density    ! .TRUE. : Use density summation model in the code,
                             ! .FALSE.: Use continuity equation
logical average_velocity     ! .TRUE. : Monaghan treatment on average velocity,
                             ! .FALSE.: No average treatment.
logical visc                 ! .TRUE. : Consider viscosity,
                             ! .FALSE.:  
logical virtual_part         ! .TRUE. : Use virtual particle,
                             ! .FALSE.:  
logical ex_force             ! .TRUE. : Consider external force,
                             ! .FALSE.:  
logical visc_artificial      ! .TRUE. : Use artificial viscosity
                             ! .FALSE.: 
logical Heat_artificial      ! .TRUE. : Use artificial heat
                             ! .FALSE.:  
logical self_gravity         ! .TRUE. : Self_gravity active
                             ! .FALSE.:   
logical nor_density          ! .TRUE. : normalized density active
                             ! .FALSE.: 
                                                   
integer(ink) ntotal, npoin, nvirt         ! ntotal  : real + virtual particles
                                          ! npoin   : number of real particles
                                          ! nvirt   : virtual particles for bcs
integer(ink) ndimn                        ! dimension of the problem 1,2 3
integer(ink) niac                         ! Maximum number of interacting pairs. 
 
integer(ink) dat_file, chk_file, gid_msh, gid_res, gmesh    ! files for input/output
integer(ink) gis_csv1, gis_csv2                             ! Adding .CSV for GiS    Saeidmt 30March2022
integer(ink) QGIS_asset, QGIS_avg, QGIS_peaks               ! MP 2023 11 09 for QGIS Vuln
integer(ink) res_file1,res_file2, res_file3                 ! files for input/output
integer(ink) res_file4,res_file5                            ! files for input/output
integer(ink) res_file6, res_file7                           ! files for input/output  ! Saeidmt 12Nov2023
integer(ink) res_file_cw                                    ! files for input/output window cases
integer(ink) init_dat, init_chk                             ! files for read sigma0
integer(ink) restart_file                                   ! file  for restart
integer(ink) restart_Top_file                               ! file for topo 1,14 and 15 restart if erosion exists

integer(ink) ic_random_walk                                 !  1 if random walk  activated
real   (irk) rndm_fact                                      !  factor for normal v , see RK00

integer(ink), allocatable:: If_Out_Domain(:)                ! value is 1 if particle is outside control domain 
integer(ink), allocatable:: If_Out_Domain_LAST(:)           ! las known value, to check wheteher part is just exited              
real   (irk), allocatable:: Xmin_Domain(:), Xmax_Domain(:)  ! Xmin_Domain(ndimn) and Xmax_Domain(ndimn)

real(irk) xmin_geom, xmax_geom        ! max and min. X coordinates of domain
real(irk) ymin_geom, ymax_geom        ! max and min. Y coordinates of domain
real(irk) zmin_geom, zmax_geom        ! max and min. Z coordinates of domain

real(irk) pi

!character(16) SPH_SW_problem_name     ! name of SW problem being solved.  MOVED TO Driver_time_vars
! character(60) text                   ! general purpose char string

integer(ink) n_avail_nodes                       !  we store here the list of nodes which we can inject. 
integer(ink), allocatable::list_avail_nodes(:)   !  set in set global vars, updated in injection
integer(ink) npoin0					             !  sometimes nodes 1 to npoin 0 are special (like behind gates)

!      --------  new vars for fractional step ** MP 2017/03/24

 
integer(ink)  ss_niac, ss_niac0                             ! Maximum number of ss interacting pairs.    ** MP FS Chuan 
integer(ink), allocatable:: list_ss_interacts(:,:)          ! Maximum number of ss interacting pairs.    ** MP FS Chuan 

real    (irk), allocatable:: vstar    (:,:)  !  ** MP s_implicit 24.03.17 V after explicit step (ndimn, npoin)
real    (irk), allocatable:: dvstar   (:,:)  !  ** MP s_implicit 05.02.21 dV after explicit step (ndimn, npoin)
real    (irk), allocatable:: vn_plus1 (:,:)  !  ** MP s_implicit 24.03.17  V after laplacian    (ndimn, npoin)
real    (irk), allocatable:: grad_p   (:,:)  !  ** MP s_implicit 11.04.17  gradient of pressure (ndimn, npoin)

real    (irk), allocatable:: gstif    (:,:)  !  ** MP s_implicit 24.03.17 K_laplacian (npoin,npoin)
real    (irk), allocatable:: gstif2 (:,:,:)  !  ** MP s_implicit 24.03.17 K_laplacian (2,2,total_ss_niac)
real    (irk), allocatable:: rhs_p      (:)  !  ** MP s_implicit 24.03.17  (npoin)...but in SFM version is npoin_FS
integer (irk), allocatable:: iffix      (:)  !  ** MP s_implicit 24.03.17  (npoin)
real    (irk), allocatable:: fixed      (:)  !  ** MP s_implicit 24.03.17  (npoin)
                                             
!      ------  back ground mesh covering active and not active elems and nodes

integer (ink)  nelem_BG, npoin_BG, nnode_BG    !  nodes per element and elements
integer (ink)  melem_BG, mpoin_BG              ! maximun dimension reached so far
integer (ink), allocatable :: intmat_BG (:,:)  !  (nnode_FS, nelem_FS)
real    (irk), allocatable :: coord_BG  (:,:)  !  (ndimn, npoin_FS)
integer (ink), allocatable :: act_nodes_BG (:) !  (npoin_FS)
integer (ink), allocatable :: act_elems_BG (:) !  (nelem_FS)
integer (ink)  ngeom_BG                        !  = nnode*ndimn +1 (3 in 1d, 7 in 2d)
real    (irk), allocatable :: geome_BG (:,:)   !  (ngeom_FS,nelem_FS) N1x N2x N3x N1y N2y N3y Area 2d 
integer (ink)  nvars_BG                        ! nr of vars transferred      
real    (ink), allocatable :: vars_BG  (:,:)   !                   
! real    (irk)  xminbg, xmaxbg, deltxbg,  Lbg

integer (ink)  nelem_FS, npoin_FS, nnode_FS    !  nodes per element and elements
integer (ink)  melem_FS, mpoin_FS              ! maximun dimension reached so far
integer (ink), allocatable :: intmat_FS (:,:)  !  (nnode_FS, nelem_FS)
real    (irk), allocatable :: coord_FS  (:,:)  !  (ndimn, npoin_FS)
integer (ink), allocatable :: act_nodes_FS (:) !  (npoin_FS)
integer (ink), allocatable :: act_elems_FS (:) !  (nelem_FS)                     
integer (ink)  ngeom_FS                        !  N1x N2x area 1d
real    (irk), allocatable :: geome_FS (:,:)   !  (ngeom_FS,nelem_FS) N1x N2x N3x N1y N2y N3y Area 2d
real    (ink), allocatable :: vars_FS  (:,:)   !                        
integer (ink)  nvars_FS                        ! nr of vars transferred                        

real    (irk), allocatable :: rho_FS      (:)  
real    (irk), allocatable :: p_FS        (:)
real    (irk), allocatable :: vstar_FS  (:,:)  
real    (irk), allocatable :: grad_p_FS (:,:)    

integer(ink)  if_left_wall, if_right_wall      !  if waalls at left or right
real   (irk)  X_left_wall,  X_right_wall       !  position

    
integer(ink)  nor_div_FS      ! 1  Normalize  div_v*        in s/implicit  ** MP s_implicit 24.03.18
integer(ink)  nor_gradP_FS    ! 1  Normalize  grad p t(n+1) in s/implicit  ** MP s_implicit
real   (irk)  tolcg_MAIN      !  ** MP s_implicit 24.03.18. read together with nor_dens_SI
integer(ink)  type_FS_problem 

!      --------  new vars added for TG2  by TB -------

integer(ink)  ntype_solid !  Type of solid problem 
						  !       0 ... 1D
						  !       1 ... plane stress
						  !       2 ... plane strain
						  !       3 ... axial symm.
						  !       4 ... 3D

integer(ink)  nprop, nmats                      !  number of mat props and number of mats
real   (irk), allocatable:: props       (:,:)   !  Material properties nmats * nprop
integer(ink), allocatable:: matno         (:)   !  Material types of all nodes


integer(ink)  nint_Vars                         !  number of internal vars
real   (irk), allocatable:: Internal_Vars (:,:) !  Internal Vars at nodes               nint_Vars * npoin1
real   (irk), allocatable:: Ddev_strn      (:), Dvol_strn(:),  DPc(:)

integer(ink) IC_norm_TG				            !  0 nothing 1 fi normalized 2 Fi and gradient normalized

integer(ink)  npoin1, npoin2, npoin2out         !  points in TG nodes and elements
integer(ink)  namaTG, nstre, nstr1              !  components of Fi in TG and nr. stress components

real   (irk), allocatable:: unknoS      (:,:)   !  vector of unks at main nodes						namaTG * npoin1
real   (irk), allocatable:: unknoS0     (:,:)   !  vector of unks at main nodes at the beginging of each time step
real   (irk), allocatable:: unkneS      (:,:)   !  vector of unks at intermediate nodes				namaTG * npoin2
real   (irk), allocatable:: unkneS0     (:,:)   !  vector of unks at intermediate nodes	at the beginging of each time step
real   (irk), allocatable:: fluxN     (:,:,:)   !  fluxes at main nodes							    ndimn * namaTG * npoin1
real   (irk), allocatable:: fluxe     (:,:,:)   !  fluxes at intr nodes							    ndimn * namaTG * npoin2
real   (irk), allocatable:: divefp      (:,:)   !  div of fluxes at main nodes					    namaTG * npoin1
real   (irk), allocatable:: divefe      (:,:)   !  div of fluxes a intr nodes					    namaTG * npoin2
real   (irk), allocatable:: SourceP     (:,:)   !  Sources at main nodes						    namaTG * npoin1
real   (irk), allocatable:: SourceE     (:,:)   !  Sources at intr nodes						    namaTG * npoin2
real   (irk), allocatable:: Sigma_p     (:,:)   !  vector of principal stress at main nodes         ndimn  * npoin1
real   (irk), allocatable:: unkno_12   (:,:)    !  used in 2nd step to store Fi (n+1/2)    

integer(ink)  ic_MLS
integer(ink), allocatable:: BNneighb_nod(:,:), BN_count_nod(:)
integer(ink), allocatable:: BNnei_vct_nod(:,:)
integer(ink)  nmlspts_nod
integer(ink), allocatable:: BNneighb_ele(:,:), BN_count_ele(:)
integer(ink), allocatable:: BNnei_vct_ele(:,:)
integer(ink)  nmlspts_ele

integer(ink) ic_intsource

integer(ink) ic_update_x

real   (irk), allocatable:: displ(:,:)

integer(ink) ic_stress_integ                ! we have a new stress integration algorithm
real   (irk) theta_stress_integ             ! Consists of obtaining v(n+theta), frem here gradsv(n+theta)
                                            ! and s(n+1)= s(n) + dt*th+De.grad V(nth)- De. evp(nth) ( No RK here)
                                            
!      ----------  vn = coditions -> normalize rho Var if_correct (npoin)=1 for these points

integer(ink), allocatable :: if_correct (:)
integer(ink) ic2_vn0                                             

END MODULE SPH_Main_Vars_2019
