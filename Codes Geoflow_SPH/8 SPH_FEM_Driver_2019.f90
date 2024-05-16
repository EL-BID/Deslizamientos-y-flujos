       
      
!      -------------------------------------------------------
!       
!
                      PROGRAM SPH_FEM_Driver_2019
                      
!
!      -------------------------------------------------------


USE SPH_FEM_variable_types_2019     !  Variable precision defined here.

USE SPH_FEM_Driver_time_vars_2019   !  This is a module shared by SPH and Geoflow
                                    !  Main sph driver...connects to SW,..
                                
USE SPH_MAIN_2019, ONLY:   Init_SPH,              	&
                           time_integration_sph,  	& 
                           adaptive_dt_sph,       	&
                           Out_plot_sph,          	&
                           Out_print_sph,         	& 
                           Out_save_sph,          	& 
                           clean_up_sph,            &
                           stop_case_SPH,           &
                           close_files_SPH,         &
                           Out_Plot_Trigger_SPH,    & 
                           Get_Prob_Propag_SPH,     &                                    
                           Update_Const_SPH,        & 
                           Manage_pixel_Hazard_SPH, &
                           deallocate_list_SPH  

USE GFL_MAIN_2019, ONLY: Init_GFL,                    &
                         Time_Integration_GFL,        &
                         Adaptive_dt_GFL,             &
                         Out_Plot_GFL ,               &
                         Out_print_GFL 

USE SPH_FEM_topo_2019    ! changed MP 01 10 2022

!* USE SPH_GFL_SW_Interactions

USE SPH_Vulnerability_2022              !  access to vulnerability module

implicit none

real   (irk) dt_max                     !  max of dt_sph and gfl. Interval to be subdivided 

integer(ink) icase, ivals_mat, iproblem ! loopcontrol cases MP CHGD 19th March 2020
real   (irk) Cval0_mat, minval_mat, maxval_mat, dval_mat, Cval, Cval0

integer(ink) index_mat, nvals_mat 

character(72) text 

1005 format(a60)
1006 format(a16)

if_Read_TopoDat = 0   !  we will read topo only the 1st time which is called

Call Driver_input

IF (ic_cases_mat.eq.1) THEN

   nproblems_index = 0
   DO icase = 1, ncases_mat
      index_mat  = Ccases_mat (icase,1)
      Cval0_mat  = Ccases_mat (icase,2) 
      nvals_mat  = Ccases_mat (icase,3)
      DO ivals_mat = 1, nvals_mat 
         nproblems_index = nproblems_index +1
         if (if_sph.eq.1) call Init_sph   !   we read here all input
         if (if_gfl.eq.1) call Init_Gfl   !  open files should be closed before reusing
         Cval = OC_mat(nproblems_index, 3 + icase)
         call Driver_Flow 
         call update_const_SPH (2, index_mat, Cval, Cval0_mat) !   reset constant(index)
         OC_mat (nproblems_index,3 + icase) = Cval             ! if not modified by update is the same
                                                               ! if density, stores porosity
          if (ic_stop_cases.eq.1) then 
             call Write_CSV
             call Out_PLOT     !  MP May 2020 this is to plot ONLY at the very end of computation                  
             call Close_files_SPH
             CYCLE
          endif
         !   close files after stop Call stop_flow          
      ENDDO   
   ENDDO

ELSEIF (ic_cases_mat .GE.2) then     !  MP 24 03 2022

   DO iproblem = 1, nproblems ! --------------------------  problems loop
      write (*,*) ' -------------  iproblem = ' , iproblem
      nproblems_index = iproblem
      call init_sph   
                                  ! changed MP 01 10 2022
      topol = topol0              ! for cases with erosion, resets topol to the original value

                                  !   check if vuln analysis required...file.vuln exists -> ic_vuln = 1 
      if (iproblem.eq.1.AND.ic_pixel_hazard.eq.1) then                         
          call init_vuln     !  MP sept 2022
      endif
      
      dt = initial_dt                !  Saeidmt 19Oct2023             
      dt_SPH = initial_dt_SPH        !  Saeidmt 19Oct2023    
      time_print    = print_step*dt; time_save  = save_step*dt; time_plot  = plot_step*dt
      t_print_reset = 0.; t_save_reset = 0.; t_plot_reset = 0. !  MP21Oct2023
      dt_max        = dt

      call  Driver_Flow_23   !  Driver_Flow_23   
      if (ic_stop_cases.eq.1) then  !  consider the case of finishing w/o stop=1 
         if (ic_TRIGGER.NE.1) then 
            if (ic_cases_win.eq.1)    call Write_CSV
            if (ic_pixel_hazard.eq.0) call Out_PLOT  !  MP May 2020 this is to plot ONLY at the very end of computation                    
         endif                                       !  do not plot when hazad propag!
         call Close_files_SPH          ! close files for all MP March 2002                
      endif
      
        if (nproblems.GT.1) then 
!!!         call deallocate_list_SPH         ! August 2022 MP avoids increase of memory un multi-problems
      endif

   ENDDO    ! -------------------------- end of problems loop
   
   IF (ic_TRIGGER.eq.1) then 
      call Out_Plot_Trigger_SPH
   endif
   
   IF (ic_Prob_Propag.eq.1) then         ! block changed MP 24th July 2022
      if (ic_pixel_hazard.eq.0) then
         call Get_Prob_Propag_SPH
      elseif (ic_pixel_hazard.eq.1) then  ! when pixel hazard analysis is done
         if (ic_FOSM.eq.1) then
            call Manage_pixel_Hazard_SPH (2) ! we will get FOSm stats at all pixels (if active)
            call Manage_pixel_Hazard_SPH (3) ! we will plot FOSM hazard vars in GID
         else
!             call Manage_pixel_Hazard_SPH (13) ! we will get  MCarlo stats at all pixels (if active)
!             call Manage_pixel_Hazard_SPH (14) ! we will plot MCarlo hazard vars in GID
         endif
      endif                               ! code 2 implies calling Get_pixel_FOSM_SW
   endif                                  ! 3 is plotting FOSM_pixel(ipoig)%FOSM_Pr(npixel_HVars)

ELSEIF (ic_cases_mat .Eq.0) then
   if (if_sph.eq.1) call Init_sph   !   we read here all input
   if (if_gfl.eq.1) call Init_Gfl   !  open files should be closed before reusing
   nproblems       = 1              ! MP 26th Aug 2022...will be used in grid_fid test
   nproblems_index = 1
   call Driver_Flow 
ENDIF

close( MASTER_dat)
close( MASTER_chk)
if (ic_cases_mat.NE.0) then
   close (MASTER_cases)
   close (MASTER_res)
endif

if (allocated(nptstcurves)) deallocate (nptstcurves)
if (allocated   (ttcurves)) deallocate    (ttcurves)
if (allocated   (ftcurves)) deallocate    (ftcurves)

CONTAINS


!-------------------------------------------------------------------

       Subroutine Driver_Input 
       
!-------------------------------------------------------------------


implicit none

character(72) text	
! character(16) problem_name_Global       ! moved to time vars by MP Sept 2022

integer (ink) len1, i, j, my_iostat				       	                			       	
integer (ink) nline, iline
integer (ink) nval, iproblem, indx, icase, jcase
integer (ink) izone, ipr, idest, jdest
integer (ink) ival, ilaw
integer (ink) iRnd                     !  MP 15th march 2022 Nr of problems to be generated<nproblems
integer (ink) n_cases, mx_cases        ! MP oct 2023
real    (irk) rnd_indx
real    (irk) minval, maxval, dval, Cval0, Cval
real    (irk) Mu_i, sigma_i

integer (ink), allocatable:: AuxT(:,:)
integer (ink), allocatable:: OC_index_tmp(:,:)   !  temporal aux var for M-CArlo MP 15 March 2022

Type record_prob_gral   !  list of stochastic variables, some zeroed,                                 
   character(60)                         :: description        ! 60 char describing Asset
   integer (ink)                         :: if_active          ! active stoch var
   integer (ink)                         :: prob_var           ! index of the variable
   integer (ink)                         :: prob_law           ! which law 
   real    (irk),dimension (:), pointer  :: prob_props         ! list of properties of law
END TYPE record_prob_gral

TYPE (record_prob_gral), allocatable :: record_prob(:)

1005 format(a60)
1006 format(a16)

MASTER_dat   = 1 ; MASTER_chk = 2		! control data read here
MASTER_cases = 3 ; MASTER_res = 4       ! cases input and output

write(*,*) ' '
write(*,*) ' '
write(*,*) ' ------------------------------------------ '
write(*,*) ' '   
write(*,*) '        MASTER INPUT...reading MASTER.dat file name '
write(*,*) ' '
write(*,*) '        Input GLOBAL problem name?'
write(*,*) ' '
write(*,*) ' ------------------------------------------ '

read (*,*)  problem_name_Global
len1 = len_trim(problem_name_Global)

open( MASTER_dat,  file=problem_name_Global (1:len1)//'.master.dat'     )
open( MASTER_chk,  file=problem_name_Global (1:len1)//'.master.chk'     )
open( MASTER_res,  file=problem_name_Global (1:len1)//'.master.res'     )
if (ic_cases_mat.EQ.1) then
   open( MASTER_cases,  file=problem_name_Global (1:len1)//'.master.cases' )
endif

write(*,*) '        MASTER file name for data is    ', problem_name_Global
write(*,*) ' '
write(*,*) '   '

read (MASTER_dat,*)  nline
write(MASTER_chk,*)  nline
Do iline = 1, nline
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
enddo

read (MASTER_dat,1005) text                         
write(MASTER_chk,1005) text                         
read (MASTER_dat,*)  if_sph, if_gfl, ic_cases_mat   
write(MASTER_chk,*)  if_sph, if_gfl, ic_cases_mat    
           ! changed if_tgf (unused) by ic_cases_mat                                                     
           ! ic_cases_mat = 1 not used (sensib analysis)                                                      
           !  2  Hypercube including reading values                                                     
           !  3  Hyp.c. MCarlo based on Hypercube.                                                  
           !  4   Full MCarlo  ic_MCarlo = 2 
           !  5 same, with GeoZones MP April 2022
           !  6 FOSM with zones    
           ! 20 FOSMfor prob_propag
           
           ! 20 FOSMfor prob_propag
           ! 31 MCarlo with Vuln at pixels and assets  prob_propag
           ! 32same, user frindly input  MP Oct 2023
 
                                                
!* type_of_fem_sph_interaction = 0                  
!* if (if_sph.eq.1.and.if_gfl.eq.1) then
!*    write(*,*) 'input type of fem sph interaction '
!*    write(*,*) ' 0 none   1 WIR   2 gate histogram '
!*    read (*,*) type_of_fem_sph_interaction
!* endif

if (if_sph.eq.1) then           
   !  Input Type of Problem....SPH_ProblemType 
   !        1....SW'
   !        2....NS'
   !        3....DF'
   !
   !  .AND  Source Integration....SPH_t_Integ_Alg 
   !        0....Central, original scheme'
   !        1....Central, Modified'
   !        4....4th order Runge Kutta 
   !        5....4th order Runge Kutta using sum for soil and con for water'
   !        6,7,8,9  .... semiimplicit and heat cases (work in progress)
   !       10.... 
   !       100...Random walk
   !       200...TRIGGER *******
         
  read  (MASTER_dat,1005) text
  write (MASTER_chk,1005) text
  read  (MASTER_dat,*)  SPH_ProblemType, SPH_t_Integ_Alg 
  write (MASTER_chk,*)  SPH_ProblemType, SPH_t_Integ_Alg

  write (*,*) ' ------------------------------------------ '
  write (*,*) ' '   
  write (*,*) '        SPH INPUT...reading file.dat name '
  write (*,*) ' '   
  read  (MASTER_dat,1005) text
  write (MASTER_chk,1005) text
  read  (MASTER_dat,1006) SPH_SW_problem_name
  write (MASTER_chk,1005) SPH_SW_problem_name 
  write (*,*) ' sph file name for data is    ', SPH_SW_problem_name
  write (*,*) ' '
  write (*,*) '   '
  write (*,*) ' ------------------------------------------ '
  write (*,*) ' '
  write (*,*) '   '
    
endif

if (if_gfl.eq.1) THEN                                           ! added MP 2021 12 01
							                                    !Select source integr'
    read (MASTER_dat,1005) text				                    !     0-> standard '
    write(MASTER_chk,1005) text				                    !     1-> Obsolete'
    read (MASTER_dat,*)  GFL_ProblemType, GFL_t_Integ_Alg       !     2-> RK4 split'
    write(MASTER_chk,*)  GFL_ProblemType, GFL_t_Integ_Alg  

    write(*,*) ' ------------------------------------------ '
    write(*,*) ' '   
    write(*,*) '        GFL INPUT...reading file.dat name '
    write(*,*) ' '

    read (MASTER_dat,1005) text              
    write(MASTER_chk,1005) text
    read (MASTER_dat,1006) GFL_SW_problem_name
    write(MASTER_chk,1005) GFL_SW_problem_name 

    write(*,*) ' GFL file name for data is    ', GFL_SW_problem_name
    write(*,*) ' '
    write(*,*) '   '
    write(*,*) ' ------------------------------------------ '
    write(*,*) ' '
    write(*,*) '   '

endif

write(MASTER_chk,*) '  .....INPUT FROM DRIVER  ............ ' 
write(MASTER_chk,*) '  .......................................... '
write(MASTER_chk,*) '     Time is = ', time
write(MASTER_chk,*) '            input (1) dt  for reference in output'
write(MASTER_chk,*) '                  (2) Time for end of computation and '
write(MASTER_chk,*) '                  (3) Maximum number of time steps maxtimestep    '
write(MASTER_chk,*) '     type a value of dt < 0 to stop         '
write(MASTER_chk,*) '  .......................................... '

read (MASTER_dat,1005) text
write(MASTER_chk,1005) text
read (MASTER_dat,*)  dt, time_end, maxtimestep
write(MASTER_chk,*)  dt, time_end, maxtimestep

initial_dt = dt           !  Saeidmt 19Oct2023

write(MASTER_chk,*) '      '  
write(MASTER_chk,*) '.  nr of steps for print_step, save_step, plot_step  '   
read (MASTER_dat,1005) text
write(MASTER_chk,1005) text
read (MASTER_dat,*)  print_step, save_step, plot_step 
write(MASTER_chk,*)  print_step, save_step, plot_step 

nsearch_step = 1		!  search for neighbours every nsearch_steps
isearch_step = 0
if (print_step.lt.0) then   ! to keep compatibility and old input
    read (MASTER_dat,1005) text
    write(MASTER_chk,1005) text
    read (MASTER_dat,*)  nsearch_step
    write(MASTER_chk,*)  nsearch_step
    print_step = - print_step
endif
                
dt_sph = -999.	    ! --- SPH ---
if (if_sph.eq.1) then
   write(MASTER_chk,*) '  input (1) dt  for sph'
   write(MASTER_chk,*) '        (2) adaptive control   (1 yes 0 n0) for sph'      
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)  dt_sph, ic_adapt_dt_sph 
   write(MASTER_chk,*)  dt_sph, ic_adapt_dt_sph
   
   initial_dt_sph = dt_sph     !  Saeidmt 19Oct2023  

endif 
      
dt_gfl = -999.	    ! --- GFL ---
if (if_gfl.eq.1) then
   write(*,*) '          Next read    (1) dt  for tgf'
   write(*,*) '                       (2) adaptive control (1 yes 0 n0) for tgf'
   write(*,*) '                       (3) ELEMENT dt control (1 yes 0 n0) for tgf'
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)  dt_gfl, ic_adapt_dt_gfl,  ic_adapt_dt_elem_gfl
   write(MASTER_chk,*)  dt_gfl, ic_adapt_dt_gfl,  ic_adapt_dt_elem_gfl
   write         (*,*)  dt_gfl, ic_adapt_dt_gfl,  ic_adapt_dt_elem_gfl      
endif   

read (MASTER_dat,1005) text
write(MASTER_chk,1005) text
read (MASTER_dat,*)  ntcurves, mptstcurves
write(MASTER_chk,*)  ntcurves, mptstcurves
   
if (ntcurves.gt.0) then

   allocate ( nptstcurves(ntcurves) )               ;nptstcurves = 0        
   allocate ( ttcurves(ntcurves,mptstcurves) )      ;ttcurves    = 0.0      
   allocate ( ftcurves(ntcurves,mptstcurves) )      ;ftcurves    = 0.0      
   
   do i = 1, ntcurves
      read (MASTER_dat,1005) text
      write(MASTER_chk,1005) text
      read (MASTER_dat,*) nptstcurves(i)     !    nr of pts in each time curve
      write(MASTER_chk,*) nptstcurves(i)
      read (MASTER_dat,1005) text
      write(MASTER_chk,1005) text
      read (MASTER_dat,*) (ttcurves(i,j), j=1,nptstcurves(i))        !  t values
      write(MASTER_chk,*) (ttcurves(i,j), j=1,nptstcurves(i))  
      read (MASTER_dat,*) (ftcurves(i,j), j=1,nptstcurves(i))        !  f values
      write(MASTER_chk,*) (ftcurves(i,j), j=1,nptstcurves(i)) 
   enddo
      
endif
    
time_print    = print_step*dt; time_save  = save_step*dt; time_plot  = plot_step*dt
t_print_reset = 0.; t_save_reset = 0.; t_plot_reset = 0. 
dt_max        = dt

dt_limit_inf  = 1.e-7*dt  

! defaults for triggering and probabilistic propagation 

  ic_Trigger     = 0 
  ic_Prob_Propag = 0               !  changed 9th June 2022
  ic_pixel_Hazard= 0               !  MP 21 July 2022
  ic_MCarlo = 0
  ic_FOSM   = 0
  N_Gzones       = 1                                  
  npixel_HVars   = 10              !  MP 21 Sept 2022  h, v, Intens, Suscep, Vuln & their peaks                              
 ! npixel_HVars   = 6              !  MP 23 July 2022  h_av, hpeak v, vpeak, q, qpeak

  if (SPH_t_integ_Alg.EQ.200) ic_Trigger = 1                !  all block MP August 2022

  if (ic_cases_mat.eq.20.AND.SPH_t_integ_Alg.EQ.4)  then
     ic_Prob_Propag = 1 ! FOSM and prob_propagation 
     ic_FOSM        = 1
  endif
  if (ic_cases_mat.eq.21.AND.SPH_t_integ_Alg.EQ.4) then
     ic_Prob_Propag  = 1
     ic_pixel_Hazard = 1
     ic_cases_mat    = 20
     ic_FOSM         = 1
  endif 
  if (ic_cases_mat.eq.21.AND.SPH_t_integ_Alg.EQ.110) then
     ic_Prob_Propag  = 1
     ic_pixel_Hazard = 1
     ic_cases_mat    = 20
     ic_FOSM         = 1
  endif 
  if (ic_cases_mat.eq.31.AND.SPH_t_integ_Alg.EQ.4) then
     ic_Prob_Propag  = 1
     ic_pixel_Hazard = 1
     ic_MCarlo       = 31
  endif
  if (ic_cases_mat.eq.31.AND.SPH_t_integ_Alg.EQ.110) then
     ic_Prob_Propag  = 1
     ic_pixel_Hazard = 1
     ic_MCarlo       = 31
  endif 
  if (ic_cases_mat.eq.32.AND.SPH_t_integ_Alg.EQ.4) then
     ic_Prob_Propag  = 1
     ic_pixel_Hazard = 1
     ic_MCarlo       = 32
  endif
  if (ic_cases_mat.eq.32.AND.SPH_t_integ_Alg.EQ.110) then
     ic_Prob_Propag  = 1
     ic_pixel_Hazard = 1
     ic_MCarlo       = 32
  endif 
  if (ic_cases_mat.eq.4 .AND.SPH_t_integ_Alg.EQ.4)   ic_Prob_Propag = 1 ! MCarlo and prob_propagation 
  if (ic_cases_mat.eq.5 .AND.SPH_t_integ_Alg.EQ.4)   ic_Prob_Propag = 1 ! MCarlo and prob_propagation
  ! if (ic_cases_mat.eq.6 .AND.SPH_t_integ_Alg.EQ.4) ic_Prob_Propag = 1 ! FOSM  use 20,21 for propagation  
  if (ic_cases_mat.eq.4 .AND.SPH_t_integ_Alg.EQ.100) ic_Prob_Propag = 1 ! MCarlo and RWalk 
  if (ic_cases_mat.eq.5 .AND.SPH_t_integ_Alg.EQ.100) ic_Prob_Propag = 1 ! MCarlo and RWalk

  if (ic_cases_mat.EQ.3) then        !   MCarlo analysis run when ic_cases set to 3, then back to 1
     ic_MCarlo = 1                   !   the analysis continues with type 2
     ic_cases_mat = 2
  endif  
    
IF ( ic_cases_mat.eq.1) THEN       !   ------   first option explore along lines
  
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max, law 
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
   enddo

   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     !  

   nproblems = 0
   do icase = 1, ncases_mat
      nval = Ccases_mat(icase,3)
      nproblems = nproblems + nval
   enddo

   ! ndat_out = 3 + ncases_mat + 1  ! nr of columns in output csv file
   Mdat_out = 30   ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out))
   OC_mat = 0

   iproblem = 0
   do icase = 1, ncases_mat
      indx   = Ccases_mat (icase,1)
      Cval0  = Ccases_mat (icase,2)
      nval   = Ccases_mat (icase,3)
      minval = Ccases_mat (icase,4)
      maxval = Ccases_mat (icase,5)
      if (indx.EQ.13.OR.indx.EQ.27) then    !  13 is Bfact and 27 is kw permeability, Vt
         minval = log (minval)               !  6th April to include log scales
         maxval = log (maxval)
      endif
      if (nval.eq.1) then
         dval = 0.0
      else
         dval = ( maxval- minval )/(nval-1)
      endif 
      do ival = 1, nval
         iproblem = iproblem +1
         OC_mat (iproblem,1) = icase
         OC_mat (iproblem,2) = ival
         OC_mat (iproblem,3) = indx 
         Cval  = minval + dval*(ival-1)
         if (indx.EQ.13.OR.indx.EQ.27)  Cval = exp(Cval)   !  log scales 
         do jcase = 1, ncases_mat
            if (jcase.eq.icase) then
               OC_mat (iproblem, 3+jcase) = Cval
            else
               OC_mat (iproblem, 3+jcase) = Ccases_mat (jcase,2)
            endif
         enddo
      enddo
   enddo

ELSEIF (ic_cases_mat.eq.2) then   !   --- analyze N1 * N2 *...Nncases 
  
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   Nvals_max = 0                            ! max nr of stats values read
   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max  
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      if (Ccases_mat(icase,6).GE.100) then 
          ival      = Ccases_mat(icase,3)
          Nvals_max = max (Nvals_max, ival)  ! MP 2022 03 07
      endif
   enddo

   if (Nvals_max.GT.0) then                  ! allocate matrixes with stat values and their probs
      read (MASTER_dat,1005) text               ! reads heading
      write(MASTER_chk,1005) text
      allocate ( List_Vals (ncases_mat,Nvals_max+1) ) ! we use 1st position to store index
      allocate ( List_Probs(ncases_mat,Nvals_max+1) )
      List_Vals = 0.0; List_Probs = 0.0
      Do icase = 1, ncases_mat
         if (Ccases_mat(icase,6).GE.100 ) then
            nval   = Ccases_mat(icase,3)    ! number of values to read, 1st var index
            read (MASTER_dat,*)    (List_Vals (icase, ival), ival=1, nval+1)
            read (MASTER_dat,*)    (List_Probs(icase, ival), ival=1, nval+1)
            write(MASTER_chk,1005) (List_Vals (icase, ival), ival=1, nval+1)
            write(MASTER_chk,1005) (List_Probs(icase, ival), ival=1, nval+1)
         endif
      enddo
   endif

   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     ! 

   nproblems = 1             !   --- compute nproblems as product of vars in each case
   Do icase= 1, ncases_mat
      nproblems = nproblems * Ccases_mat(icase,3)
   enddo 

!      ------  get indexes matrix OC_index. For each problem the set of ival indexes

   if (ic_MCarlo.eq.1) then
      read (MASTER_dat,1005) text
      write(MASTER_chk,1005) text
      read (MASTER_dat,*)    n_MCarlo_problems    !  read the nr of simulations, MP 15 March 2022
      write(MASTER_chk,*)    n_MCarlo_problems    !  should be smaller than nproblems
      if (n_Mcarlo_problems.GT.nproblems) then
          n_Mcarlo_problems = nproblems
          !  ic_MCarlo = 0
          !  exit loop
      endif 
      allocate (OC_index_tmp(n_Mcarlo_problems,ncases_mat))
      OC_index_tmp = 0
   endif
   
   allocate (OC_index(nproblems, ncases_mat))
   allocate (AuxT    (4, ncases_mat))
   OC_index = 0 ; AuxT = 0
   do icase = 1, ncases_mat
      auxT(1,icase) = Ccases_mat (icase,3)   !  nvals of variable icase
   enddo
   auxT(2,1) = 1
   do icase = 2, ncases_mat
      auxT(2,icase)= auxT(2,icase-1)*auxT(1,icase-1)
   enddo
   DO iproblem = 0, nproblems-1
      auxT (3,ncases_mat) = iproblem
      auxT (4,ncases_mat) = auxT (3,ncases_mat)/auxT (2,ncases_mat)
      DO icase = ncases_mat -1, 1,-1
         auxT (3,icase) = auxT (3,icase+1) - auxT (2,icase+1) *auxT (4,icase+1) 
         auxT (4,icase) = auxT (3,icase) / auxT (2,icase) 
      enddo
      OC_index (iproblem+1, :) = auxT(4,:) +1
   enddo
   
   if(ic_MCarlo.eq.1) then
      do iproblem = 1, n_MCarlo_problems
         call random_number(rnd_indx)
         iRnd = nint (rnd_indx * nproblems)
         if (iRnd.eq.0) iRnd = 1
         if (iRnd.Gt.nproblems) irnd = nproblems
         OC_index_tmp(iproblem,:) = OC_index(iRnd,:)
      enddo
      deallocate (OC_index)
      allocate   (OC_index(n_MCarlo_problems,ncases_mat))
      OC_index = OC_index_tmp
      deallocate (OC_index_tmp)
      nproblems = n_MCarlo_problems
   endif

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)

   Mdat_out = 30                           ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0                              ! 30 is OK for 3 win and 13 variables

   DO iproblem = 1, nproblems              ! ---   fill the props
      DO icase = 1, ncases_mat
         ilaw   = Ccases_mat(icase,6)
         ival   = OC_index (iproblem, icase)   !  this is the i-th value for var icase in iproblem
         
         IF (ilaw.GE.100) THEN                 ! statistic list 
            OC_mat (iproblem, 3 + icase) = list_Vals(icase, ival+1)
         
         ELSE
         
            indx   = Ccases_mat (icase,1)   ! ------  this is the index saying which var is it
            Cval0  = Ccases_mat(icase,2)
            nval   = Ccases_mat(icase,3)
            MinVal = Ccases_mat(icase,4)
            Maxval = Ccases_mat(icase,5)
            if (ilaw.eq.2) then             !   ---  logaritmic spacing (permeability)
               MinVal = log (MinVal)
               MaxVal = log (MaxVal)
            endif
            if (nval.eq.1) then
               dval = 0.0
            else
               dval = (MaxVal - MinVal) / (nval-1)
            endif
            cval = MinVal + dval*(ival-1)
            if (ilaw.eq.2) then
               Cval = exp (CVal)
            endif
            OC_mat (iproblem, 3 + icase) = Cval

            if (indx.eq.521) then   !     ------  if var is a vector (dx, dy) 
               MinVal = Ccases_mat(icase,7)
               Maxval = Ccases_mat(icase,8)
               if (nval.eq.1) then
                  dval = 0.0
               else
                  dval = (MaxVal - MinVal) / (nval-1)
               endif
               cval = MinVal + dval*(ival-1)
               OC_mat (iproblem, 3 + ncases_mat +1 ) = Cval
            endif
         ENDIF   !  ends list or direct computation if
        
      enddo      !  icase loop
   enddo         !  iproblem loop
   
   icase= 1 ;   !   for debugging STOP
     

ELSEIF (ic_cases_mat.eq.4) then   !   --- we will generate sets of nproblems 
                                  ! stochastic sets for all ncases_mat props
   ic_MCarlo = 2
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max 
   Ccases_mat = 9.9
    
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
   enddo
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     ! 

!      ------  get indexes matrix OC_index. For each problem the set of ival indexes


   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    n_MCarlo_problems    !  read the nr of simulations, MP 15 March 2022
   write(MASTER_chk,*)    n_MCarlo_problems    !  should be smaller than nproblems
   nproblems = n_MCarlo_problems 
   allocate (OC_index(nproblems, ncases_mat))
   OC_index = 0 

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)

   Mdat_out = 30                           ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0.0                            ! 30 is OK for 3 win and 13 variables

!      ------  generate the sets of stochastic values for each case magnitude

   nprob_props = 10
   allocate (probab_law (ncases_mat))                  ; probab_law = 0
   allocate (Param_prob_law (ncases_mat, nprob_props)) ; Param_prob_law = 0.0
   allocate (Tmp_List(nproblems))                       ; Tmp_List = 0.0

   DO icase = 1, ncases_mat                       ! generate a set of nproblems rnd values for icase
      
      read (MASTER_dat,1005) text                 !   type of prob distribution law
      write(MASTER_chk,1005) text                 !   0...    constant value
      read (MASTER_dat,*)    probab_law (icase)   !   1...    triangle 
      write(MASTER_chk,*)    probab_law (icase)   !   100 ... histogram    

      if (probab_law (icase) .LT.100) then        !   continuous

         read (MASTER_dat,1005) text                 !   parameters of prob distribution law
         write(MASTER_chk,1005) text                 !     
         read (MASTER_dat,*)    (Param_prob_law (icase, i), i = 1,nprob_props )   
         write(MASTER_chk,*)    (Param_prob_law (icase, i), i = 1,nprob_props )

         Tmp_List = 0.0
        
         call GetMCarlo_Cont_values (icase)          ! continuum distributions 
          
      else                                           

         Tmp_List = 0.0
                                
         call GetMCarlo_Dscrt_values (icase)        ! discrte laws...histograms        
                      
      endif

      Do iproblem = 1, nproblems                          ! store TMp_list in OC_Mat
         OC_Mat (iproblem, 3+icase ) = Tmp_List(iproblem) ! fills column 3 + icase 
      enddo                                               ! with the nproblems values

   ENDDO   !   ends loop icase       

ELSEIF (ic_cases_mat.eq.5) then   !   --- we will generate sets of nproblems 
                                  ! stochastic sets for all ncases_mat props
   ic_MCarlo = 5                  ! includes Geo Zones  
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max 
   Ccases_mat = 9.9
    
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
   enddo
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     ! 

!      ------  get indexes matrix OC_index. For each problem the set of ival indexes

   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    n_MCarlo_problems, n_GZones    !  read the nr of simulations 
   write(MASTER_chk,*)    n_MCarlo_problems, n_GZones    !  and geo zones  
   nproblems = n_MCarlo_problems * n_GZones              !  MP April 2022
   allocate (OC_index(nproblems, ncases_mat))
   OC_index = 0 

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)

   Mdat_out = 30                           ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0.0                            ! 30 is OK for 3 win and 13 variables

!      ------  generate the sets of stochastic values for each case magnitude

   nprob_props = 10
   allocate (probab_law (ncases_mat))                  ; probab_law = 0
   allocate (Param_prob_law (ncases_mat, nprob_props)) ; Param_prob_law = 0.0
   allocate (Tmp_List (n_MCarlo_problems))             ; Tmp_List = 0.0

   DO icase = 1, ncases_mat                       ! generate a set of nproblems rnd values for icase
   DO izone = 1, n_GZones                         ! for each case , nGeo sets of vars MP April 2022
      
      read (MASTER_dat,1005) text                 !   type of prob distribution law
      write(MASTER_chk,1005) text                 !   0...    constant value
      read (MASTER_dat,*)    probab_law (icase)   !   1...    triangle 
      write(MASTER_chk,*)    probab_law (icase)   !   100 ... histogram    

      if (probab_law (icase) .LT.100) then        !   continuous

         read (MASTER_dat,1005) text              !   parameters of prob distribution law
         write(MASTER_chk,1005) text              !     
         read (MASTER_dat,*)    (Param_prob_law (icase, i), i = 1,nprob_props )   
         write(MASTER_chk,*)    (Param_prob_law (icase, i), i = 1,nprob_props )

         Tmp_List = 0.0
        
         call GetMCarlo_Cont_values (icase)        ! continuum distributions 
          
      else                                           

         Tmp_List = 0.0
                                
         call GetMCarlo_Dscrt_values (icase)        ! discrte laws...histograms        
                      
      endif

      Do ipr    = 1, n_MCarlo_problems               ! store TMp_list in OC_Mat
         idest  = (izone-1)* n_MCarlo_problems + ipr ! sets for zones                      
         OC_Mat (idest, 3+icase ) = Tmp_List(ipr)    ! fills column 3 + icase 
      enddo                                          ! with the npMCarlo values
   
   ENDDO   !   ends loop geo zones
   ENDDO   !   ends loop icase  
  

ELSEIF (ic_cases_mat.eq.6) then   !  FOSM method   24 April 2022 MP
                                  !  
   ic_MCarlo = 6                  ! includes Geo Zones  
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max 
   Ccases_mat = 9.9
    
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
   enddo
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     ! 

!      ------  get indexes matrix OC_index. For each problem the set of ival indexes

   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    n_MCarlo_problems, n_GZones    !  read the nr of simulations 
   write(MASTER_chk,*)    n_MCarlo_problems, n_GZones    !  and geo zones  
   nproblems = n_MCarlo_problems * n_GZones              !  MP April 2022
   allocate (OC_index(nproblems, ncases_mat))
   OC_index = 0 

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)

   Mdat_out = 30                           ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0.0                            ! 30 is OK for 3 win and 13 variables

!      ------  generate the sets of stochastic values for each case magnitude

   nprob_props = 10
   allocate (FOSMms( 2*N_GZones + 1, ncases_mat)) ! (1+2*N_GZones,ncases_mat)   
                                                  ! columns are cases row 1 index,
   FOSMms = 0.0                                   ! (row 2 mean, row 3 sigma) for all zones 

   DO icase = 1, ncases_mat                           !  
      FOSMms(1,icase) = CCases_mat(icase,1)  ! this is the index of the icase-th variable
   ENDDO   
      
   DO izone = 1, n_GZones                ! for each case , nGeo sets of vars MP April 2022
      
      read (MASTER_dat,1005) text                  
      write(MASTER_chk,1005) text                  
      idest = 1 + (izone-1)*2                              ! index to get position of files 
      iproblem =  1 + (izone-1)                            ! iproblem in OC_MAt
      read (MASTER_dat,*)  (FOSMms(idest + 1 , i), i= 1, ncases_mat)    ! Means  
      read (MASTER_dat,*)  (FOSMms(idest + 2 , i), i= 1, ncases_mat)    ! sigmas 
      write(MASTER_chk,*)  (FOSMms(idest + 1 , i), i= 1, ncases_mat)    ! Means  
      write(MASTER_chk,*)  (FOSMms(idest + 2 , i), i= 1, ncases_mat)    ! sigmas    

   ENDDO   !   ends loop geo zones

   DO  iproblem = 1, nproblems
   DO  icase    = 1, ncases_mat
       OC_mat (iproblem, 3+icase) = FOSMms (2*iproblem, icase)
   ENDDO
   ENDDO
 
   read (MASTER_dat,1005) text    ! if 1 prob(x,t) activated               
   write(MASTER_chk,1005) text    !  ONLY for FOSm , hence block moved here
   read (MASTER_dat,*)    ic_xtp  !   MP 11 May 2022              
   write(MASTER_chk,*)    ic_xtp

   if (ic_xtp.eq.1) then
      read (MASTER_dat,1005) text               
      write(MASTER_chk,1005) text   
      read (MASTER_dat,*)    ma_xtp, mh_xtp  !  nr of values of PGA and alphas used                
      write(MASTER_chk,*)    ma_xtp, mh_xtp
      if (ma_xtp.GT.0) then
         allocate ( Acc_xtp(mA_xtp), TrA_xtp(mA_xtp) , PAi_xtp  (mA_xtp) )
         read (MASTER_dat,1005) text               
         write(MASTER_chk,1005) text     
         read (MASTER_dat,*)    (Acc_xtp(i), i = 1, ma_xtp)   ! accelerations             
         write(MASTER_chk,*)    (Acc_xtp(i), i = 1, ma_xtp) 
         read (MASTER_dat,1005) text               
         write(MASTER_chk,1005) text  
         read (MASTER_dat,*)    (Tra_xtp(i), i = 1, ma_xtp)   !  Return periods             
         write(MASTER_chk,*)    (Tra_xtp(i), i = 1, ma_xtp)
         read (MASTER_dat,1005) text               
         write(MASTER_chk,1005) text  
         read (MASTER_dat,*)    Nyra_xtp                      ! period (years) for P[Ai]               
         write(MASTER_chk,*)    Nyra_xtp
!          call Get_Ai (ma_xtp, TrA_xtp, Nyra_xtp, PAi_xtp  )  
         call Get_Ai (1, ma_xtp, TrA_xtp, Nyra_xtp, PAi_xtp  )  ! MP 15th July use exp for accel
      endif
      if (mh_xtp.GT.0) then
         allocate ( Hrs_xtp (mh_xtp), Trh_xtp(mh_xtp) , PHi_xtp  (mh_xtp) )
         read (MASTER_dat,1005) text               
         write(MASTER_chk,1005) text     
         read (MASTER_dat,*)    (hrs_xtp(i), i = 1, ma_xtp)   ! hsat/h            
         write(MASTER_chk,*)    (hrs_xtp(i), i = 1, ma_xtp) 
         read (MASTER_dat,1005) text               
         write(MASTER_chk,1005) text  
         read (MASTER_dat,*)    (Trh_xtp(i), i = 1, ma_xtp)   !  Return periods             
         write(MASTER_chk,*)    (Trh_xtp(i), i = 1, ma_xtp)
         read (MASTER_dat,1005) text               
         write(MASTER_chk,1005) text  
         read (MASTER_dat,*)    Nyrh_xtp                      ! period (years) for P[Ai]               
         write(MASTER_chk,*)    Nyrh_xtp
         call Get_Ai (2,mh_xtp, Trh_xtp, Nyrh_xtp, Phi_xtp  ) ! MP 15th July use exp for accel
      endif
   endif

ELSEIF (ic_cases_mat.eq.20) then   !  FOSM for Propagation (1 pic RW aor full SPH)
                                   !  MP 30th May 2022
   ic_MCarlo      = 0              ! no Geo Zones, so is 1  
   ic_Prob_Propag = 1              
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max 
   Ccases_mat = 0.0
    
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
   enddo
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     ! 

!      ------  get indexes matrix OC_index. For each problem the set of ival indexes

   n_MCarlo_problems = 2*ncases_mat + 1
   n_GZones  = 1
   nproblems = n_MCarlo_problems * n_GZones              

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)

   Mdat_out = 30                           ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0.0                            ! 30 is OK for 3 win and 13 variables

!      ------  generate the sets of stochastic values for each case magnitude

   nprob_props = 10
   allocate (FOSMms( 3, ncases_mat)) ! (1+2*N_GZones,ncases_mat)   
                                     ! columns are cases row 1 index,
   FOSMms = 0.0                      ! (row 2 mean, row 3 sigma) for all zones
                                     ! rows mu and sigma, cols dist, hpk, vpk, Qpk 
   ! allocate (FOSM_Pr2 (2,7)) 
   ! FOSM_Pr2 = 0.0

   DO icase = 1, ncases_mat                           !  
      FOSMms(1,icase) = CCases_mat(icase,1)  ! this is the index of the icase-th variable
   ENDDO   
   
   read (MASTER_dat,1005) text                  
   write(MASTER_chk,1005) text  
   read (MASTER_dat,*)  (FOSMms (2 , i), i= 1, ncases_mat)    ! Means  
   read (MASTER_dat,*)  (FOSMms( 3 , i), i= 1, ncases_mat)    ! sigmas 
   write(MASTER_chk,*)  (FOSMms( 2 , i), i= 1, ncases_mat)    ! Means  
   write(MASTER_chk,*)  (FOSMms( 3 , i), i= 1, ncases_mat)    ! sigmas    

   DO icase = 1, ncases_mat          ! 1st line contains mu1...mu_ncase (means)
      OC_mat (1 , 3 + icase) = FOSMms (2 , icase)
   ENDDO
   
   DO icase = 1, ncases_mat  
      idest = 2*icase
      DO jcase = 1, ncases_mat
         OC_mat (idest  , jcase + 3 ) = FOSMms (2, icase)
         OC_mat (idest+1, jcase + 3 ) = FOSMms (2, icase)
      ENDDO
   ENDDO

   DO icase = 1, ncases_mat
      idest = 2*icase
      jdest = 3 + icase
      mu_i    = FOSMms (2 , icase)
      sigma_i = FOSMms (3 , icase)
      OC_mat (idest     , jdest) = mu_i + sigma_i/10.
      OC_mat (idest + 1 , jdest) = mu_i - sigma_i/10.
   ENDDO

ELSEIF (ic_cases_mat.eq.31) then   !  MCarlo for Propagation (1 pic RW aor full SPH)
                                   !  MP 4th August 2022
   ic_MCarlo      = 31               
   ic_Prob_Propag = 1              
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   read (MASTER_dat,*)    ncases_mat, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat, ninfo_mat    ! matrix defining changes

   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max 
   Ccases_mat = 0.0
    
   read (MASTER_dat,1005) text              ! posn t for dymin dymax. Overrides ninfo_mat
   write(MASTER_chk,1005) text
   do icase = 1, ncases_mat
      read (MASTER_dat,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
      write(MASTER_chk,*) (Ccases_mat(icase, ivals_mat), ivals_mat=1, ninfo_mat)
   enddo
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     ! 

!      ------  get indexes matrix OC_index. For each problem the set of ival indexes

   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    n_MCarlo_problems, n_GZones    !  read the nr of simulations 
   write(MASTER_chk,*)    n_MCarlo_problems, n_GZones    !  and geo zones  
   nproblems = n_MCarlo_problems * n_GZones              !  MP April 2022
   allocate (OC_index(nproblems, ncases_mat))
   OC_index = 0 
                        ! 30 is OK for 3 win and 13 variables

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)

   Mdat_out = 30                           ! 3 + ncases_mat + 2 + 6*ncase_win 
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0.0                            ! 30 is OK for 3 win and 13 variables

!      ------  generate the sets of stochastic values for each case magnitude

   nprob_props = 10
   allocate (probab_law (ncases_mat))                  ; probab_law = 0
   allocate (Param_prob_law (ncases_mat, nprob_props)) ; Param_prob_law = 0.0
   allocate (Tmp_List (n_MCarlo_problems))             ; Tmp_List = 0.0

   DO icase = 1, ncases_mat                       ! generate a set of nproblems rnd values for icase
   DO izone = 1, n_GZones                         ! for each case , nGeo sets of vars MP April 2022
      
      read (MASTER_dat,1005) text                 !   type of prob distribution law
      write(MASTER_chk,1005) text                 !   0...    constant value
      read (MASTER_dat,*)    probab_law (icase)   !   1...    triangle 
      write(MASTER_chk,*)    probab_law (icase)   !   100 ... histogram    

      if (probab_law (icase) .LT.100) then        !   continuous

         read (MASTER_dat,1005) text              !   parameters of prob distribution law
         write(MASTER_chk,1005) text              !     
         read (MASTER_dat,*)    (Param_prob_law (icase, i), i = 1,nprob_props )   
         write(MASTER_chk,*)    (Param_prob_law (icase, i), i = 1,nprob_props )

         Tmp_List = 0.0
        
         call GetMCarlo_Cont_values (icase)        ! continuum distributions 
          
      else                                           

         Tmp_List = 0.0
                                
         call GetMCarlo_Dscrt_values (icase)        ! discrte laws...histograms        
                      
      endif

      Do ipr    = 1, n_MCarlo_problems               ! store TMp_list in OC_Mat
         idest  = (izone-1)* n_MCarlo_problems + ipr ! sets for zones                      
         OC_Mat (idest, 3+icase ) = Tmp_List(ipr)    ! fills column 3 + icase 
      enddo                                          ! with the npMCarlo values
   
   ENDDO   !   ends loop geo zones
   ENDDO   !   ends loop icase  
  

ELSEIF (ic_cases_mat.eq.32) then   !  MCarlo for Propagation (1 pic RW aor full SPH)
                                   !  MP 4th August 2022
   ic_MCarlo      = 32               
   ic_Prob_Propag = 1  
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    vstop_cases, timestop_cases, dtimestop_cases     !  stop problem conditions
   write(MASTER_chk,*)    vstop_cases, timestop_cases, dtimestop_cases     !  

   nprob_props    = 10
   Mdat_out       = 30  
            
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    ncases_mat,mx_cases, ninfo_mat    !  CHGD MP 19th March
   write(MASTER_chk,*)    ncases_mat,mx_cases, ninfo_mat    ! matrix defining changes
  
   allocate (Ccases_mat(ncases_mat, 10))    ! Const index, initvalue, nr vals, min, max 
   Ccases_mat = 0.0
   
!      ------  get  n_MCarlo_problems, n_GZones

   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text
   read (MASTER_dat,*)    n_MCarlo_problems, n_GZones                                !  read  nr of simulations 
   write(MASTER_chk,*)    n_MCarlo_problems, n_GZones    !  and geo zones = 1 here
   n_GZones  = 1                                         ! just in case user fails 
   nproblems = n_MCarlo_problems * n_GZones              !  MP April 2022                   

!      ------  Obtain the nproblems sets of values (note , 521 has dx and dy)
                          ! 3 + ncases_mat + 2 + 6*ncase_win 
   
   allocate ( OC_mat (nproblems,Mdat_out)) !  Still dont know ncase_win read in INPUT
   OC_mat = 0.0                            ! 30 is OK for 3 win and 13 variables

   allocate (record_prob(mx_cases))
   DO icase = 1, mx_cases
      allocate (record_prob(icase)%prob_props(nprob_props)  )
      record_prob(icase)%prob_var   = 0
      record_prob(icase)%prob_law   = 0
      record_prob(icase)%prob_props = 0.0
   enddo
!      ------  generate the sets of stochastic values for each case magnitude

   allocate (probab_law (ncases_mat))                  ; probab_law = 0
   allocate (Param_prob_law (ncases_mat, nprob_props)) ; Param_prob_law = 0.0
   allocate (Tmp_List (n_MCarlo_problems))             ; Tmp_List = 0.0

   n_cases = 0
   
   read (MASTER_dat,1005) text
   write(MASTER_chk,1005) text

   DO icase = 1, mx_cases
      read (MASTER_dat,*) record_prob(icase)%prob_var,& 
                          record_prob(icase)%prob_law,&   
                          (record_prob(icase)%prob_props(ival), ival = 1, nprob_props)
      write(MASTER_chk,*) record_prob(icase)%prob_var,& 
                          record_prob(icase)%prob_law,&  
                          (record_prob(icase)%prob_props(ival), ival = 1, nprob_props)
      if (record_prob(icase)%prob_law.GT.0 ) then
          n_cases = n_cases + 1
          Ccases_mat(n_cases, 1) = record_prob(icase)%prob_var 
          probab_law (n_cases)   = record_prob(icase)%prob_law
          do ival = 1, nprob_props
             Param_prob_law(n_cases,ival)= record_prob(icase)%prob_props(ival) 
          enddo       
      endif 
   enddo

   DO icase = 1, ncases_mat                       ! generate a set of nproblems rnd values for icase
   DO izone = 1, n_GZones                         ! for each case , nGeo sets of vars MP April 2022
      
      if (probab_law (icase) .LT.100) then        !   continuous 

         Tmp_List = 0.0
        
         call GetMCarlo_Cont_values (icase)        ! continuum distributions 
          
      else                                           

         write (*,*) " code stopped --- no discrete vals for MCarlo 32"
         PAUSe
         STOP
                                
         call GetMCarlo_Dscrt_values (icase)        ! discrte laws...histograms        
                      
      endif

      Do ipr    = 1, n_MCarlo_problems               ! store TMp_list in OC_Mat
         idest  = (izone-1)* n_MCarlo_problems + ipr ! sets for zones                      
         OC_Mat (idest, 3+icase ) = Tmp_List(ipr)    ! fills column 3 + icase 
      enddo                                          ! with the npMCarlo values
   
   ENDDO   !   ends loop geo zones
   ENDDO   !   ends loop icase  
      
ENDIF            !  ends (ic_cases_mat.eq.32) case 

read (MASTER_dat,1005, IOSTAT = my_iostat) text
if (my_iostat.LT.0) then
    ic_cases_win = 0 
else
    write(MASTER_chk,1005) text
    read (MASTER_dat,*)  ic_cases_win
    write(MASTER_chk,*)  ic_cases_win
endif 

END SUBROUTINE  Driver_Input  ! ============================================END DRIVER INPUT


!-------------------------------------------------------------------

       SUBROUTINE Get_Ai (accel_or_rain, m, TRi, N, PAi )
       
!-------------------------------------------------------------------

implicit none

integer (ink) accel_or_rain, m, N            ! accel or rain 1 accel 2 rain MP 15th Juky 2022
real    (irk), dimension (:) :: TRi, PAi

integer (ink) i
real    (irk) Pi (m)

Pi (1) = 1.0
do i = 2, m
   if (accel_or_rain.EQ.1) then 
       Pi (i) = 1. - (1.-1./TRi(i))**N
   elseif (accel_or_rain.EQ.2) then
       Pi (i) =  1. - exp(-N/TRi(i))    !  1. - (1.-1./TRi(i))**N 
   endif
enddo

do i = 1, m-1
   PAi (i) = Pi(i) - Pi(i+1)
enddo
PAi (m) = Pi(m)

END SUBROUTINE Get_Ai  

!-------------------------------------------------------------------

       SUBROUTINE GetMCarlo_Cont_values  (icc)
       
!-------------------------------------------------------------------

implicit none

integer (ink) p_law, icc
real    (irk) yrnd, xrnd, trnd
real    (irk) a, b, c, ba, ca, cb, baca, cbca, yb, aa, aainv
real    (irk) G_mean, G_sigma, LgN_mean, LgN_sigma 

integer (ink) Seed(1), Count

p_law = probab_law(icc)

if (p_law.eq.0) then         !   --- Uniform

   Tmp_List = Param_prob_law (icc,1)

elseif (p_law.eq.1) then     !   --- Triangular
   a = Param_prob_law (icc,1)
   b = Param_prob_law (icc,2)
   c = Param_prob_law (icc,3)
   ba = b-a ; ca = c-a ; cb = c-b
   baca = ba*ca ; cbca = cb*ca
   yb   = ba/ca
   Call System_clock (Count)
   Seed(1) = Count
   call RANDOM_SEED (PUT=Seed)
   do iproblem = 1, n_MCarlo_problems           !  nproblems
      call random_number(yrnd)
      if (yrnd.LE.yb) then
         xrnd = a + (baca * yrnd)**0.5
      else
         xrnd = c - (cbca * (1.-yrnd))**0.5
      endif   
      Tmp_List (iproblem) = xrnd   
   enddo  

elseif (p_law.eq.2) then               !   Normal distribution 

   G_mean  =  Param_prob_law (icc,1)   !  using Edous & Eidous 2018 approx
   G_sigma =  Param_prob_law (icc,2)   !  aa alternatives: 
   aa      = 5./8.  ; aainv = 1/aa     ! ( (2/pi), (5/8); 81/130; 0.647-0.0021x )
   Call System_clock (Count)
   Seed(1) = Count
   call RANDOM_SEED (PUT=Seed)
   do iproblem = 1, n_MCarlo_problems  !  nproblems !  FI = 0.5* ( 1 + sqrt(1-exp(-aa*x2))
      call random_number(yrnd)
      if (yrnd.LT.0.00001) yrnd = 0.00001
      if (yrnd.GT.0.99999) yrnd = 0.99999 
      b = log (4*yrnd*(1.-yrnd))
      if (yrnd.GE.0.5) then
         trnd =  sqrt (-aainv * b )
      else
         trnd = -sqrt (-aainv * b )
      endif   
      xrnd =  G_mean + G_sigma * trnd 
      Tmp_List (iproblem) = xrnd   
   enddo   

elseif (p_law.eq.3) then               !   Log Normal distribution 
   G_mean  =  Param_prob_law (icc,1)   !  using Edous & Eidous 2018 approx
   G_sigma =  Param_prob_law (icc,2)   !  aa alternatives:
                                       ! ( (2/pi), (5/8); ** 81/130; 0.647-0.0021x )
   LgN_sigma = (log (1. + (G_sigma/G_mean)**2) )**0.5     ! See Fenton and Griffiths book
   LgN_mean  = log (G_mean) - 0.5*LgN_sigma*LgN_sigma
   aa      = 5./8.  ; aainv = 1/aa 
   Call System_clock (Count)
   Seed(1) = Count
   call RANDOM_SEED (PUT=Seed)    
   do iproblem = 1, n_MCarlo_problems  !  nproblems !  FI = 0.5* ( 1 + sqrt(1-exp(-aa*x2))
      call random_number(yrnd)
      if (yrnd.LT.0.00001) yrnd = 0.00001
      if (yrnd.GT.0.99999) yrnd = 0.99999 
      b = log (4*yrnd*(1.-yrnd))
      if (yrnd.GE.0.5) then
         trnd =  sqrt (-aainv * b )
      else
         trnd = -sqrt (-aainv * b )
      endif  
      xrnd = exp (LgN_mean + trnd*LgN_sigma)
      Tmp_List (iproblem) = xrnd   
   enddo  
endif


END SUBROUTINE GetMCarlo_Cont_values

!-------------------------------------------------------------------

    SUBROUTINE GetMCarlo_Dscrt_values  (icc) 
       
!-------------------------------------------------------------------

implicit none
integer (ink) p_law, icc, i, ip, ie, iproblem, iexit
real    (irk) yrnd, xrnd, trnd
real    (irk) discrim, dd 

integer (ink) Seed(1), Count

1005 format(a60)
1006 format(a16)

p_law = probab_law(icc)
 
IF (p_law.eq.100) then                          !  read the histogram Xs and Ns
                                               
   read (MASTER_dat,1005) text                  !   histogram
   write(MASTER_chk,1005) text                  !   
   read (MASTER_dat,*)    np_h                  ! number of points
   write(MASTER_chk,*)    np_h                  ! ne_h is nr of elems
   ne_h = np_h-1
   allocate ( x_h (np_h) ) ; x_h  = 0.0         ! list of values of magnitudes
   allocate ( p_h (ne_h) ) ; p_h  = 0.0         ! p nr of observations (not normalized)
   allocate ( F_h (np_h) ) ; F_h  = 0.0         ! F  accumulated, norm to 1
   read (MASTER_dat,1005) text                  !   histogram nodes
   write(MASTER_chk,1005) text                  !   
   read (MASTER_dat,*)    (x_h (i), i = 1, np_h)
   write(MASTER_chk,*)    (x_h (i), i = 1, np_h) 
   read (MASTER_dat,1005) text                  !   histogram elems vslues
   write(MASTER_chk,1005) text                  !   
   read (MASTER_dat,*)    (p_h (i), i = 1, ne_h)
   write(MASTER_chk,*)    (p_h (i), i = 1, ne_h)   

!   ---  TASK 1: Generate the accumulated prob law at nodes
   
   Do ip = 2, np_h 
      ie = ip -1 
      F_h (ip) = F_h (ip-1) + p_h (ie)
   enddo
   F_h = F_h /F_h (np_h)   ! normalize, last value will be 1

!   ---  TASK 2:  generate rnd values of y, and obtain xrnd
   
   Call System_clock (Count)
   Seed(1) = Count
   call RANDOM_SEED (PUT=Seed) 
   Do iproblem = 1, n_MCarlo_problems           !  nproblems
      call random_number(yrnd)
      iexit = 0
      Do ip = 1, np_h - 1 
         discrim = (F_h(ip)-yrnd)*(F_h(ip+1)-yrnd)
         if (discrim.LT.0.0) then
            dd   = (yrnd - F_h(ip) ) / (F_h(ip+1)-F_h(ip))
            xrnd = x_h (ip) + dd * (x_h (ip+1) - x_h (ip))
            iexit = 1
         endif
         if (iexit.eq.1) EXIT
      enddo
      Tmp_List (iproblem) = xrnd   
   enddo  
ELSEIF (p_law.eq.200) then ! we will use it to generate histograms from 
                           ! continuous p(x) not invertible        
ENDIF

END SUBROUTINE GetMCarlo_Dscrt_values


!-------------------------------------------------------------------

       Subroutine Driver_Flow_23  
       
!-------------------------------------------------------------------

implicit none

real   (irk) tpr_try , dt_keep   !  MP 28 Oct 2023  search for exact timeplot
integer(ink) ntpr_try            !  MP 28 Oct 2023   

1005 format(a60)
1006 format(a16)

nstart   = 0.0; current_ts = 0          !  Initialize for timestepping
time     = 0.0;   
time_sph = 0.0  

if (time_sph_restart.gt.0.0) then
   time     = time_sph_restart
   time_sph = time_sph_restart
endif

if (ic_Trigger.eq.1) maxtimestep = 1     !   MP 19 March 2022
      
DTGT: DO itimestep = nstart+1, nstart + maxtimestep   
           
      current_ts    = current_ts + 1   
      
      if (dt.lt.0.0.or.dt.lt.dt_limit_inf) EXIT
            
      if (if_sph.eq.1.AND.ic_TRIGGER.eq.0) THEN	! MP MArch 2022 avoids unnecessary computations
         call adaptive_dt_sph                   !  Obtains dt_sph. MAIN calls to SW, NS,.. 
      endif 
      
      dt_keep  = dt_sph
      ntpr_try = 1
      tpr_try  = t_plot_reset + dt_sph
      if (tpr_try.GE.Time_plot) then
         ntpr_try = INT(tpr_try/time_plot)
         dt_sph   = ntpr_try*time_plot -  t_plot_reset
      endif

      t_plot_reset  = t_plot_reset  + dt_sph 

      t_print_reset = t_print_reset + dt_sph  ! new MP 27Oct 2023
      t_save_reset  = t_save_reset  + dt_sph 

      itimestep_sph = itimestep_sph + 1 

      call time_integration_sph      
                        
      time_sph = time_sph +  dt_sph   
      time     = time     +  dt_sph   
      
      !   debug case T07wMCarlo follow point 9  along problems ----------
      !  call Manage_pixel_Hazard_SPH (999)  !  ********************  DEBUG
    
      if (t_plot_reset.ge.time_plot) then
         if (ic_cases_mat.EQ.0.OR.nproblems.eq.1.AND.ic_pixel_hazard.EQ.0 ) then 
                                            !  MP May 2020...if cases, only at the end
            call Out_PLOT                   !  MP Oct 2023 
         elseif (ic_pixel_hazard.GT.0) then
             if (ic_FOSM.eq.1) then
                call Manage_pixel_Hazard_SPH (1)
             else
                call Manage_pixel_Hazard_SPH (11)        ! MP August 2022 Mcarlo 
             endif
         endif
         dt_sph = dt_keep      !  MP 29th Oct 2023 avoid dt decreasing 
         t_plot_reset = 0.0 
       endif 
   
      if (ic_trigger.NE.1) then                    ! 25 April 2022 MP
         if (t_save_reset.ge.time_save) then
            call Out_Save                          ! Routine in sph plus gfl. Router
            t_save_reset = 0.0 
          endif 
   
         if (t_print_reset.ge.time_print) then
             Call Out_PRN                         !  Routine in sph plus gfl. Router
             t_print_reset = 0.0 
         endif
         call stop_case_SPH
         if (ic_stop_cases.eq.1) then              ! skipp when ictrigger.ne.1 April 2022 MP 
             EXIT
         endif
         if (time.ge.time_end)   EXIT 
      endif

   101 format(1x,3(2x,a12))     
   100 format(1x,4(2x,e14.6)) 

   ENDDO DTGT  

nstart = current_ts

if (ic_pixel_hazard.GE.1.AND.ic_MCarlo.GE.1) then   ! normalizes accum values to obtain mean
    call Manage_pixel_Hazard_SPH (12)               ! MP August 2022 Mcarlo 
!     call Manage_pixel_Hazard_SPH (998)            !  ********************  DEBUG 
endif

if (ic_TRIGGER.eq.1) ic_stop_cases = 1       !  only 1 iter for Trigger problems MP 03 2022

! if (if_sph.eq.1) call clean_up_sph

End Subroutine Driver_Flow_23



!-------------------------------------------------------------------

       Subroutine Driver_Flow  
       
!-------------------------------------------------------------------

implicit none

integer(ink) nls_sph, nls_gfl, nls_tgf  !  number of substeps
integer(ink) ils_sph, ils_gfl, ils_tgf

1005 format(a60)
1006 format(a16)

nstart   = 0.0; current_ts = 0          !  Initialize for timestepping
time     = 0.0;                         !  MP 22Oct2023  dt = 0.001
time_gfl = 0.0; itimestep_gfl = 0       ! counts total iters in sph and gfl
time_sph = 0.0; itimestep_sph = 0 

dt_max   = dt                      !  Saeidmt 19Oct2023

if (time_sph_restart.gt.0.0) then
   time     = time_sph_restart
   time_sph = time_sph_restart
endif

if (ic_Trigger.eq.1) maxtimestep = 1     !   MP 19 March 2022
      
DTGT: DO itimestep = nstart+1, nstart + maxtimestep   
           
      current_ts    = current_ts + 1      
      t_print_reset = t_print_reset + dt_max
      t_plot_reset  = t_plot_reset  + dt_max
      t_save_reset  = t_save_reset  + dt_max
      
      if (dt.lt.0.0.or.dt.lt.dt_limit_inf) EXIT
            
      if (if_sph.eq.1.AND.ic_TRIGGER.eq.0) THEN	! MP MArch 2022 avoids unnecessary computations
         call adaptive_dt_sph                   !  Obtains dt_sph. MAIN calls to SW, NS,.. 
      endif
      
       if (if_gfl.eq.1) THEN	       ! --- GFLF ---
            call adaptive_dt_gfl       !  Obtains dt_gfl
       endif
      
      dt_max = max (dt_sph, dt_gfl)

      nls_sph = 0                              !  Obtain number of subdivisions and dt for sph
      if (if_sph.eq.1) then
          nls_sph = max(1.,(dt_max/dt_sph))
!         nls_sph = int( (0.99*dt_max/dt_sph)+1 )
          dt_sph  = dt_max/nls_sph
      endif

      nls_gfl = 0                               !  Obtain number of subdivisions and dt for gfl
      if (if_gfl.eq.1) then
          nls_gfl = max(1.,(dt_max/dt_gfl))
          dt_gfl  = dt_max/nls_gfl
      endif

!*      IF (type_of_fem_sph_interaction.eq.0) THEN! The standard loop when no interactions wished 

      if (if_sph.eq.1) then                  !  if SPH present, perform nls local steps
         DO ils_sph = 1, nls_sph
            itimestep_sph = itimestep_sph +1                    
                call time_integration_sph
            time_sph = time_sph + dt_sph     !  update local time for sph
         ENDDO
      endif
       
      if (if_gfl.eq.1) then                  !  if GFL present, perform nls local steps
          DO ils_gfl = 1, nls_gfl
             itimestep_gfl = itimestep_gfl +1                  
             call Time_Integration_GFL
             time_gfl = time_gfl + dt_gfl     !  update local time for gfl            
          ENDDO
      endif
!*    
!*      ELSEIF (type_of_fem_sph_interaction.eq.1) THEN    ! WIR interactions 
!*    
!*          call get_pics_in_element               ! T1 obtains which_element 
!*          
!*          DO ils_gfl = 1, nls_gfl
!*             itimestep_gfl = itimestep_gfl +1
!*             if (ils_gfl.eq.1.and.WIR_Interact) then 
!*                 call get_gradhs_w               ! T2 obtains hs_w and gradhs_w
!*             endif
!*             call time_integration_gfl           !  Timestep
!*             time_gfl = time_gfl + dt_gfl        !  update local time for sph
!*          ENDDO
!*          
!*          DO ils_sph = 1, nls_sph
!*             itimestep_sph = itimestep_sph +1
!*             if (WIR_Interact) call get_gradhw_s ! T3 obtains hs_w and gradhs_w
!*             call time_integration_sph           ! performs timestep
!*             call get_pics_in_element            ! T1 obtains which_element 
!*             time_sph = time_sph + dt_sph        !  update local time for sph
!*          ENDDO  
!*       
!*    ELSE
!*       write(*,*) ' type_of_fem_sph_interaction', type_of_fem_sph_interaction, ' NOT VALID'
!*       pause
!*       
!*    ENDIF
      
      time = time + dt_max
      
      !   debug case T07wMCarlo follow point 9  along problems ----------
 
!     call Manage_pixel_Hazard_SPH (999)  !  ********************  DEBUG
    
      if (t_plot_reset.ge.time_plot) then
         if (ic_cases_mat.EQ.0.OR.nproblems.eq.1.AND.ic_pixel_hazard.EQ.0 ) then 
                                            !  MP May 2020...if cases, only at the end
            call Out_PLOT                   !  MP Oct 2023 
         elseif (ic_pixel_hazard.GT.0) then
             if (ic_FOSM.eq.1) then
                call Manage_pixel_Hazard_SPH (1)
             else
                call Manage_pixel_Hazard_SPH (11)        ! MP August 2022 Mcarlo 
             endif
         endif
         t_plot_reset = 0.0 
       endif 
   
      if (ic_trigger.NE.1) then                    ! 25 April 2022 MP
         if (t_save_reset.ge.time_save) then
            call Out_Save                          ! Routine in sph plus gfl. Router
            t_save_reset = 0.0 
          endif 
   
         if (t_print_reset.ge.time_print) then
             Call Out_PRN                         !  Routine in sph plus gfl. Router
             t_print_reset = 0.0 
         endif
         call stop_case_SPH
         if (ic_stop_cases.eq.1) then              ! skipp when ictrigger.ne.1 April 2022 MP 
             EXIT
         endif
         if (time.ge.time_end)   EXIT 
      endif

   101 format(1x,3(2x,a12))     
   100 format(1x,4(2x,e14.6)) 

   ENDDO DTGT  

nstart = current_ts

if (ic_pixel_hazard.GE.1.AND.ic_MCarlo.GE.1) then   ! normalizes accum values to obtain mean
    call Manage_pixel_Hazard_SPH (12)               ! MP August 2022 Mcarlo 
!     call Manage_pixel_Hazard_SPH (998)              !  ********************  DEBUG 
endif

if (ic_TRIGGER.eq.1) ic_stop_cases = 1       !  only 1 iter for Trigger problems MP 03 2022

! if (if_sph.eq.1) call clean_up_sph

End Subroutine Driver_Flow



!-------------------------------------------------------------------

       Subroutine Driver_Flow_changed  
       
!-------------------------------------------------------------------

implicit none

integer(ink) nls_sph, nls_gfl, nls_tgf  !  number of substeps
integer(ink) ils_sph, ils_gfl, ils_tgf

integer(ink) if_reset_plot, ndt_target  !  MP 25 Oct 2023 Changes for exact timeplot
real   (irk) t_plot_target, t_plot_reset_target  !  MP 25 Oct 2023  exact timeplot

1005 format(a60)
1006 format(a16)

nstart   = 0.0; current_ts = 0          !  Initialize for timestepping
time     = 0.0;                         !  MP 22Oct2023  dt = 0.001
time_gfl = 0.0; itimestep_gfl = 0       ! counts total iters in sph and gfl
time_sph = 0.0; itimestep_sph = 0 

dt_max   = dt                      !  Saeidmt 19Oct2023

if (time_sph_restart.gt.0.0) then
   time     = time_sph_restart
   time_sph = time_sph_restart
endif

if (ic_Trigger.eq.1) maxtimestep = 1     !   MP 19 March 2022
      
DTGT: DO itimestep = nstart+1, nstart + maxtimestep   
           
      current_ts    = current_ts + 1  
          
!      t_print_reset = t_print_reset + dt_max   27thOct2023
!      t_plot_reset  = t_plot_reset  + dt_max
!      t_save_reset  = t_save_reset  + dt_max

      if_reset_plot = 0     !  MP 25 Oct 2023 Changes for exact timeplot
      
      if (dt.lt.0.0.or.dt.lt.dt_limit_inf) EXIT
            
      if (if_sph.eq.1.AND.ic_TRIGGER.eq.0) THEN	! MP MArch 2022 avoids unnecessary computations
         call adaptive_dt_sph                   !  Obtains dt_sph. MAIN calls to SW, NS,.. 
      endif
      
       if (if_gfl.eq.1) THEN	       ! --- GFLF ---
            call adaptive_dt_gfl       !  Obtains dt_gfl
       endif
      
      dt_max = max (dt_sph, dt_gfl)   
      t_print_reset = t_print_reset + dt_max  ! new MP 27Oct 2023
      t_plot_reset  = t_plot_reset  + dt_max
      t_save_reset  = t_save_reset  + dt_max

      nls_sph = 0                              !  Obtain number of subdivisions and dt for sph
      if (if_sph.eq.1) then
          nls_sph = max(1.,(dt_max/dt_sph))
!         nls_sph = int( (0.99*dt_max/dt_sph)+1 )
          dt_sph  = dt_max/nls_sph
      endif

      nls_gfl = 0                               !  Obtain number of subdivisions and dt for gfl
      if (if_gfl.eq.1) then
          nls_gfl = max(1.,(dt_max/dt_gfl))
          dt_gfl  = dt_max/nls_gfl
      endif

!*      IF (type_of_fem_sph_interaction.eq.0) THEN! The standard loop when no interactions wished 

      if (if_sph.eq.1) then                  !  if SPH present, perform nls local steps    
         
         !   -----------------     !  MP 25 Oct 2023 Changes for exact timeplot
         if_reset_plot       = 0           
         t_plot_target       = time_sph +  dt_sph
         t_plot_reset_target = t_plot_reset +  dt_sph 
         if ( t_plot_reset_target.GT.time_plot) then
            ndt_target          = (t_plot_target/time_plot)
            t_plot_target       =  ndt_target * time_plot
!!!!!!!            dt_sph              =  t_plot_target - time_sph
            if_reset_plot       =  1
         endif
     !  ----------------------  !  MP 25 Oct 2023 Changes for exact timeplot

         DO ils_sph = 1, nls_sph
            itimestep_sph = itimestep_sph +1                    
                call time_integration_sph
            time_sph = time_sph + dt_sph     !  update local time for sph
         ENDDO

      endif                 !  END of if SPH present, perform nls local steps 

      if (if_gfl.eq.1) then                  !  if GFL present, perform nls local steps
          DO ils_gfl = 1, nls_gfl
             itimestep_gfl = itimestep_gfl +1                  
             call Time_Integration_GFL
             time_gfl = time_gfl + dt_gfl     !  update local time for gfl            
          ENDDO
      endif
!*    
!*      ELSEIF (type_of_fem_sph_interaction.eq.1) THEN    ! WIR interactions 
!*    
!*          call get_pics_in_element               ! T1 obtains which_element 
!*          
!*          DO ils_gfl = 1, nls_gfl
!*             itimestep_gfl = itimestep_gfl +1
!*             if (ils_gfl.eq.1.and.WIR_Interact) then 
!*                 call get_gradhs_w               ! T2 obtains hs_w and gradhs_w
!*             endif
!*             call time_integration_gfl           !  Timestep
!*             time_gfl = time_gfl + dt_gfl        !  update local time for sph
!*          ENDDO
!*          
!*          DO ils_sph = 1, nls_sph
!*             itimestep_sph = itimestep_sph +1
!*             if (WIR_Interact) call get_gradhw_s ! T3 obtains hs_w and gradhs_w
!*             call time_integration_sph           ! performs timestep
!*             call get_pics_in_element            ! T1 obtains which_element 
!*             time_sph = time_sph + dt_sph        !  update local time for sph
!*          ENDDO  
!*       
!*    ELSE
!*       write(*,*) ' type_of_fem_sph_interaction', type_of_fem_sph_interaction, ' NOT VALID'
!*       pause
!*       
!*    ENDIF
      
      time = time +  dt_max
      
      !   debug case T07wMCarlo follow point 9  along problems ----------
 
!       call Manage_pixel_Hazard_SPH (999)  !  ********************  DEBUG
    
      if (t_plot_reset.ge.time_plot) then
         if (ic_cases_mat.EQ.0.OR.nproblems.eq.1.AND.ic_pixel_hazard.EQ.0 ) then 
                                            !  MP May 2020...if cases, only at the end
            call Out_PLOT                   !  MP Oct 2023 
         elseif (ic_pixel_hazard.GT.0) then
             if (ic_FOSM.eq.1) then
                call Manage_pixel_Hazard_SPH (1)
             else
                call Manage_pixel_Hazard_SPH (11)        ! MP August 2022 Mcarlo 
             endif
         endif
         t_plot_reset = 0.0 
       endif 
   
      if (ic_trigger.NE.1) then                    ! 25 April 2022 MP
         if (t_save_reset.ge.time_save) then
            call Out_Save                          ! Routine in sph plus gfl. Router
            t_save_reset = 0.0 
          endif 
   
         if (t_print_reset.ge.time_print) then
             Call Out_PRN                         !  Routine in sph plus gfl. Router
             t_print_reset = 0.0 
         endif
         call stop_case_SPH
         if (ic_stop_cases.eq.1) then              ! skipp when ictrigger.ne.1 April 2022 MP 
             EXIT
         endif
         if (time.ge.time_end)   EXIT 
      endif

   101 format(1x,3(2x,a12))     
   100 format(1x,4(2x,e14.6)) 

   ENDDO DTGT  

nstart = current_ts

if (ic_pixel_hazard.GE.1.AND.ic_MCarlo.GE.1) then   ! normalizes accum values to obtain mean
    call Manage_pixel_Hazard_SPH (12)               ! MP August 2022 Mcarlo 
!      call Manage_pixel_Hazard_SPH (998)              !  ********************  DEBUG 
endif

if (ic_TRIGGER.eq.1) ic_stop_cases = 1       !  only 1 iter for Trigger problems MP 03 2022

! if (if_sph.eq.1) call clean_up_sph

End Subroutine Driver_Flow_changed


!-------------------------------------------------------------------

       Subroutine Out_Plot  
       
!-------------------------------------------------------------------


    implicit none

!      ------  Output for geoflow and sph

!* if (if_sph.eq.0.and.if_gfl.eq.0.and.if_tgf.eq.1) then
!*    call Out_Plot_tgf
!* endif

if (if_sph.eq.1.and.if_gfl.eq.0) then
   call Out_Plot_SPH
endif   

 if (if_sph.eq.0.and.if_gfl.eq.1) then
    if (ic_pixel_hazard.eq.0) call Out_Plot_GFL
 endif

!* if (if_sph.eq.1.and.if_gfl.eq.1) then
!*    if (type_of_fem_sph_interaction.eq.0) then
!*       call Out_Plot_SPH
!*       call Out_Plot_GFL
!*    elseif (type_of_fem_sph_interaction.eq.1) then
!*       call plot_wir_res
!*    endif
!* endif

End Subroutine Out_Plot 


!-------------------------------------------------------------------

       Subroutine OUT_PRN   
       
!-------------------------------------------------------------------


    implicit none

!      ------  Output for geoflow and sph
 

!* if (if_sph.eq.0.and.if_gfl.eq.0.and.if_tgf.eq.1) then
!*    call Out_PRINT_tgf
!* endif 

if (if_sph.eq.1.and.if_gfl.eq.0) then
   call Out_PRINT_SPH
endif 

!* if (if_sph.eq.0.and.if_gfl.eq.1) then
!*    call Out_PRINT_GFL
!* endif 
  
!* if (if_sph.eq.1.and.if_gfl.eq.1) then
!*    call Out_PRINT_SPH
!*    call Out_PRINT_GFL
!* endif


End Subroutine OUT_PRN


!-------------------------------------------------------------------

       Subroutine Out_Save  
       
!-------------------------------------------------------------------


    implicit none

!      ------  Restart for   sph

if (if_sph.eq.1.and.if_gfl.eq.0) then
   call Out_Save_SPH
endif  

End Subroutine Out_Save 


!-------------------------------------------------------------------

       Subroutine Write_CSV  
       
!-------------------------------------------------------------------


implicit none
integer(ink) iproblem,jdat, ishift, i0,jf
real   (irk) ndout
character(32) my_fmt
character(3) spac

ishift = 0
if (ic_cases_mat.eq.2) ishift = 3

spac = ' , '
! my_fmt = ' 3(I6, " , ", 1X), 20 (G8.4," , ", 1X) '
my_fmt = '  (20 (G8.4," , ", 1X)  )'

1001 format (20 (G10.4,A3))

if (ic_cases_win.eq.0) then
   ndout = 3 + ncases_mat + 2 + 1
   write(MASTER_res, 1001) (OC_mat(nproblems_index,jdat), spac, jdat = 1 +ishift, ndout)
elseif  (ic_cases_win.eq.1) then
   write(MASTER_res, 1001) (OC_mat(nproblems_index,jdat), spac, jdat = 1 +ishift, ndat_out)
endif  


End Subroutine Write_CSV



END PROGRAM SPH_FEM_Driver_2019

