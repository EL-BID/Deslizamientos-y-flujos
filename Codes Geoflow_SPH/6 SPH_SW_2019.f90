 MODULE  SPH_SW_2019  

!----------------------------
!
!   2016: added d pwp for DFs
!
!   Note:  implemented erosion rate dh/dt = Es.h.v  Es estimated from ln(Volfinal/Vol 0) /distance
!          new constants: nrs 4 where coriolis was
!          21 June erosion goes also to linear momentum
!          22 June ellipsoidal for Stromboli icunk_s = 12
!          02 July Voellmy uses xi->cmanning tau = frict + rho*g*V*V/xi
!          24 July injection via if_out_domain
!          19 September abs and vn0 BCs
!          ...
!          16  March 2008  fixed Voellmy bug and new CF model
!          September 2008  new coarse mesh for output
!                          if_coarse_plt = 1 means we activate it
!                           f_coarse_plt   is a factor which multiplies deltxg and deltyg in topo mesh
!                          npoi_plt        points in new mesh
!                          xmin_plt...ymax_plt  inscribing rectangle
!                          deltx_plt, delty_plt new spacings
!                          npoix_plt, npoiy_plt number of grid lines in x and y
!                          coor_plt (3,npoi_plt) stores x, y and Z
!                          poi_plt_in_model(npoi_plt) = 1 point in computational domain
!
! EROSION ------  Read properties associated to terrain, such as basal friction and max erodible depth
!
!     ntopo_props = 6
!               1 ... Zmin
!               2 ... Zmax
!               3 ... type 1 for friction 2 for erosion 3 for both
!               4 ... Hungr's constant which modifies c4
!               5 ... max erosion depth
!               6 ... free

!
! See also changes in log file
!
! Changed IntForces_SW for ic_interact_SW = 2 (2 phase DFs)
!
!---------------------------------

USE SPH_FEM_variable_types_2019
USE SPH_FEM_Driver_time_vars_2019
USE SPH_Main_Vars_2019 
USE SPH_SW_Vars_2019
USE SPH_FEM_topo_2019
USE GCsph_solvers_2019                   !  CG for semi.implicit SW from TB ** MP 2017/03/25
USE quicksort_2019                       ! sort for 1D fem FS Chuan
USE My_SPH_Tools_2019                    !  New experiment
USE SPH_Vulnerability_2022               !  Vulnerability module Sept 2022

implicit none

PRIVATE       !      ------  No variable will be accesible from outside the module
              !              unless declared public (none indeed)

public:: Get_Init_SPH_SW
public:: Input_SW     ,  Pint_Update_SW 
public:: IntForces_SW ,  ExtForces_SW 
public:: OutputMesh_SW,  OutputRes_SW, OutputSave_SW, output_points_SW, Moni_SW
public:: Adaptive_dt_sph_SW  
public:: Abs_BCs_SW, Vn0_BCs_SW    
public:: change_to_w_SW
public:: Get_Pwp_SW
public:: Get_Rain_SW              ! MP Feb 2023 rain, pwp and drainage slow landslides
public:: Drum_drive_SW            ! for drum cases

public:: clean_SW                 !      ------  for histogram input
public:: Inject_histogram_SW

public:: Get_s_at_w_DF_SW         ! get hw at soil pics and viceversa  DFs
public:: Get_s_at_w_Wir_SW        ! get hw at soil pics and viceversa  WIRs

public:: Second_Step_FS_SW_FEM
public:: Second_Step_FS_SW_SFM

public:: Update_Const_SW
public:: Stop_case_SW
public:: Output_cases_win_SW
public:: Close_files_SW

public:: Trigger_FoS_SW             ! triggering March 16 2022
public:: Trigger_FOSM_SW 
public:: Out_Plot_Trigger_SW
public:: Get_Prob_Propag_SW 

public:: Get_pixel_hVars_SW         !  MP 23 July 2022
public:: Plot_pixel_RES_SW          !  MP 25th July 2022
public:: Get_pixel_FOSM_SW          !  MP 24 July 2022

public:: Get_MC_pixel_Vars_SW       !  MP 3rd August 2022
public:: Norm_MC_pixel_Vars_SW
public:: Stats_MC_pixel_Vars_SW
public:: Plot_MC_pixel_Vars_SW

public:: add_w_to_s_SW              !  MP 21st Jully 2023 to ad water. Called by RK45
public:: add_w_to_s_SW_TRI 

CONTAINS


!-------------------------------------------------------------------

       Subroutine add_w_to_s_SW_TRI
       
!-------------------------------------------------------------------


!      ------   adds water to soil Only FS  MP July 2021
  
implicit none 
                                            !  
real   (irk) dhw
real   (irk)  dMassT , dumm, TMass          !  MP for debugg only 2023 01 13
integer(ink) ipoin

if (.not.allocated(ds)) then  !  MP for debugg only 2023 01 13
  allocate ( ds(npoin))       ! ds is dmass
  ds  = 0.0
endif
 
IF (ic_w_inflow.EQ.1) THEN   !  only if inflow of waater activated 
      
      ds = 0.0;  dMassT= 0.0
      dhw    = w_inflow%flow_q * dt_sph
      Tmass  = 0.0
      do ipoin = 1, npoin_s
         ds (ipoin)  = dhw*mass(ipoin)/rho(ipoin)
         rho(ipoin)  = rho(ipoin) + dhw
         dMassT      = dMassT + ds (ipoin) 
         mass(ipoin) = mass(ipoin) + ds (ipoin)   !  dhw*mass(ipoin)/rho(ipoin)  !  mp Nov 2023
         TMass       = Tmass + mass(ipoin)
      enddo       
             
ENDIF    ! flow slide model ends

end subroutine  add_w_to_s_SW_TRI 



!-------------------------------------------------------------------

       Subroutine add_w_to_s_SW 
       
!-------------------------------------------------------------------


!      ------   adds water to soil Only FS  MP July 2021
  
implicit none 
                                            !  
real   (irk)  poros, hw, hw0, dhw, dn, xg, yg, Zp, Zmin, Zmax
real   (irk)  n_try, n_max                                  !  MP 5 Jan 2024
real   (irk)  dMassT , dumm          !  MP for debugg only 2023 01 13
integer(ink) ipoin

if (.not.allocated(ds)) then  !  MP for debugg only 2023 01 13
  allocate ( ds(npoin))       ! ds is dmass
  ds  = 0.0
endif
 
IF (ic_w_inflow.EQ.1) THEN   !  only if inflow of waater activated 

   IF (ic_FS.EQ.1) THEN      !  Flowslide model
      ds = 0.0;  dMassT= 0.0
      do ipoin = 1, npoin_s
         xg = x(1,ipoin)
         if (ndimn.eq.2) then 
            yg = x(2,ipoin)
         else
         yg = 0.5*(yming + ymaxg)
         endif
         Call Get_Z_Topo (ipoin,xg, yg, Zp )
         Zmin = w_inflow % Zmin; Zmax = w_inflow % Zmax
         if (Zp.GE.Zmin.AND.ZP.LE.Zmax) THEN
             poros  = n_DF(ipoin)
             hw0    = rho(ipoin)*poros  
             dhw    = w_inflow%flow_q * dt_sph
             dn     = (1-poros)*dhw/rho(ipoin)
             n_try  = poros + dn              ! MP 5th Jan 2024 limits n max
             if ( const(49).GT.1.e-6) then 
                n_max = const(49)
             else
                n_max = 0.8
             endif
             if (n_try.GT.n_max) then
                n_try = n_max
                dn = n_try -  poros
             endif
                                             
             n_DF(ipoin) = n_DF(ipoin) + dn  
             ds (ipoin)  = dhw*mass(ipoin)/rho(ipoin)
             rho(ipoin)  = rho(ipoin) + dhw  ! MP for debugg only 2023 01 13
             dMassT      = dMassT + ds (ipoin) 
             mass(ipoin) = mass(ipoin) +  ds (ipoin) !!!   MP JAN 2024 dhw*mass(ipoin)/rho(ipoin)  
         endif
         dMassT  = dMassT     
      enddo
     Dumm = 0.0  
   ENDIF    ! flow slide model ends

   IF (ic_DF.eq.1) THEN      !  MP Dec 2022 DFs only add h and mass to water nodes
      do ipoin = npoin_s + 1, npoin_s +  npoin_w
         xg = x(1,ipoin)
         if (ndimn.eq.2) then 
            yg = x(2,ipoin)
         else
         yg = 0.5*(yming + ymaxg)
         endif
         Call Get_Z_Topo (ipoin,xg, yg, Zp )
         Zmin = w_inflow % Zmin; Zmax = w_inflow % Zmax
         if (Zp.GE.Zmin.AND.ZP.LE.Zmax) THEN
             dhw   = w_inflow%flow_q * dt_sph 
             rho(ipoin)  = rho(ipoin) + dhw
             mass(ipoin) = mass(ipoin) + dhw*mass(ipoin)/rho(ipoin)  !  mp Nov 2023
         endif     
      enddo
      
      call Get_s_at_w_DF_SW         ! kmowing hw new update at soil ptes and N-DF etc

   ENDIF  ! DFs model ends

ENDIF  !  endif of only if inflow of waater activated 

end subroutine  add_w_to_s_SW 


!-------------------------------------------------------------------

       Subroutine add_w_to_s_SW_OLD
       
!-------------------------------------------------------------------


!      ------   adds water to soil Only FS  MP July 2021
  
implicit none 
                                            !  
real   (irk)  poros, hw, hw0, dhw, dn, xg, yg, Zp, Zmin, Zmax

integer(ink) ipoin
 
if (ic_w_inflow.EQ.1) then 

if (ic_FS.NE.1) then
   write(*,*) '  -----*************************************** '
   write(*,*) '  -----*************************************** '
   write(*,*) '  to use add water to soil, problem should be of FS type '
   write(*,*) '  code will STOP '
   write(*,*) '  -----*************************************** '
   write(*,*) '  -----*************************************** '
   PAUSE
   STOP
endif

do ipoin = 1, npoin_s
   xg = x(1,ipoin)
   if (ndimn.eq.2) then 
      yg = x(2,ipoin)
   else
   yg = 0.5*(yming + ymaxg)
   endif
   Call Get_Z_Topo (ipoin,xg, yg, Zp )
   Zmin = w_inflow % Zmin; Zmax = w_inflow % Zmax
   if (Zp.GE.Zmin.AND.ZP.LE.Zmax) THEN
      poros = n_DF(ipoin)
      hw0   = rho(ipoin)*poros 
      dhw   = w_inflow%flow_q * dt_sph
      dn    = (1-poros)*dhw/rho(ipoin)
      n_DF(ipoin) = n_DF(ipoin) + dn  
      rho(ipoin)  = rho(ipoin) + dhw
      mass(ipoin) = mass(ipoin) + dhw*mass(ipoin)/rho(ipoin)  !  mp Nov 2023
   endif     
enddo

endif

end subroutine  add_w_to_s_SW_OLD 


!-------------------------------------------------------------------

       Subroutine Get_Init_SPH_SW
       
!-------------------------------------------------------------------

!      ------  Just input and call to topo. Avoids link sph-topo

call Topo_main            !      ------  Get Coorg, Topol and Xming,Xmaxg,..npoig, deltas... 
if (ic_TRIGGER.EQ.1) then !  MP 22 March new input for trigger 
   call Input_SW_TRIGGER
else
   call Input_SW    
endif 
if ( ic_cases_win.GT.0) then
   call Input_cases_WIN_SW
endif

END SUBROUTINE Get_Init_SPH_SW


!-------------------------------------------------------------------

       Subroutine Output_cases_win_SW_tst
       
!-------------------------------------------------------------------

!      ------  writes window info 

implicit none 
                                            !   to file RES_file_cw every time_print steps
real   (irk), allocatable:: h_cw(:)         !   h_DF or rho (averaged of points in window)
real   (irk), allocatable:: v_cw(:,:)       !   velocity
real   (irk), allocatable:: q_cw(:)         !   outflow trough region
real   (irk), allocatable:: vn_cw(:)        !   normal velocity to side a 
real   (irk), allocatable:: dist_min_cw(:)  !   normal velocity to side a 
integer(ink), allocatable:: in_win(:)       !   point in window

integer(ink) icw, ipoin, ip_in, np_win
character(3) spac                               ! is " , "
real   (irk) xi, eta, h_cw_MAX, dstmp           ! max in windows
real   (irk) dstmpx, dstmpy                     ! MP chgd 15 June 2022 for -ve dist
real   (irk) x0_cw(2), a, b, na(2), nb(2)       ! simpler than  direct structure
real   (irk) xtg (2)
integer(ink) is_target

if (.NOT.allocated(h_cw)) then
   allocate (h_cw(ncase_win), v_cw(ndimn,ncase_win), q_cw(ncase_win), vn_cw(ncase_win) )
   allocate (in_win(ncase_win), dist_min_cw(ncase_win))
endif

spac = ' , '

DO icw = 1, ncase_win   !   loop in active windows 
   np_win = 0; h_cw(icw) = 0.; v_cw(:,icw) = 0.; q_cw(icw) = 0.; vn_cw(icw) = 0. ; dist_min_cw(icw) = 1.e18
   x0_cw     = info_cw(icw) % x0    !  recover info from structure, case icw of ncase_win
   a         = info_cw(icw) % ab(1)
   b         = info_cw(icw) % ab(2)
   na        = info_cw(icw) % nva 
   nb        = info_cw(icw) % nvb 
   is_target = info_cw(icw) % is_target
   Do ipoin = 1, npoin
      ip_in = 0   ! default, point is not in
      xi  = (x(1,ipoin)-x0_cw(1)) * na(1)   !  determine if ipoint is in a X b
      eta = (x(1,ipoin)-x0_cw(1)) * nb(1)
      if (ndimn.eq.2) then
         xi  = xi  + (x(2,ipoin)-x0_cw(2)) * na(2)
         eta = eta + (x(2,ipoin)-x0_cw(2)) * nb(2)
      endif
      if (ndimn.eq.2) then
          if ( (xi.ge.0.0).and.(xi.le.a ).and.(eta.ge.0.0).and.(eta.le.b ) ) ip_in = 1
      elseif (ndimn.eq.1) then
          if ( (eta.ge.0.0).and.(eta.le.b ) ) ip_in = 1
      endif   
      if (ip_in.eq.1)   then     !  point is inside
         np_win = np_win +1  ! increment number of points in rectangle a X b
         if (allocated(h_DF) ) then
            h_cw (icw) = h_cw (icw) + h_DF(ipoin)
         else
            h_cw (icw) = h_cw (icw) + rho (ipoin)
         endif
         v_cw (:,icw) = v_cw (:, icw) + vx(:,ipoin)
      endif
      if (is_target.eq.1) then     !  dist new CHGD MP 15th June 2022
          dstmp = 0.0; dstmpx = 0.0; dstmpy = 0.0  
          xtg (1) = x0_cw (1) + a*na(1)/2.
          dstmp = (x(1,ipoin)-xtg(1))*nb(1)
          if (ndimn.eq.2) then
               xtg (2) = x0_cw (2) + b*nb(2)/2.
               dstmp = (x(2,ipoin)-xtg(2))*nb(2)
          endif
          if (dstmp.lt. dist_min_cw(icw)) then
              dist_min_cw(icw) = dstmp
          endif  
      endif
   enddo       !    ------- ipoin loop
   if (np_win.GT.0 ) then
      h_cw (icw)   = h_cw (icw) / np_win
      v_cw (:,icw) = v_cw (:,icw) / np_win
      vn_cw (icw)  = - v_cw(1,icw)* nb(1)
      if (ndimn.eq.2) then 
         vn_cw (icw) = vn_cw (icw) - v_cw(2,icw)* nb(2)
      endif 
      q_cw (icw) = vn_cw (icw) * h_cw (icw) * a 
   endif
ENDDO          !    ------- icw loop     

!       We have obtained h_sw, v_cw, vn_cw and q_cw at all ncase_win windows at t = time
!       time to write to file filename.CaseWin.RES (unit is res_file_cw) if h diff from zero

is_target = info_cw(ncase_win) % is_target
if (is_target.eq.1) then          !   target window print always dist 
!    write(RES_file_cw,1001) (  (time, h_cw(icw),vn_cw(icw),q_cw(icw), dist_min_cw(icw)),icw= 1, ncase_win )
    write(RES_file_cw,1001)     time,  (( h_cw(icw),vn_cw(icw),q_cw(icw)),icw= 1, ncase_win) , dist_min_cw(ncase_win)  
else
   h_cw_MAX = 0.0
   DO icw = 1, ncase_win
      if (h_cw(icw).GT. h_cw_MAX) h_cw_MAX = h_cw(icw) 
   enddo
   if (h_cw_MAX.GE.const(15)) then   !   write only when max h is GT threshold C15 *****
      write(RES_file_cw,1001) (  (time, h_cw(icw),vn_cw(icw),q_cw(icw)),icw= 1, ncase_win )
   endif 
endif
1001 format (20 (G10.4,2x))

!  from current variables, update 
!  T0_cw, Tf_cw, hpeak_cw, vnpeak_cw, qpeak_cw, QTotal_cw (ncase_win)

DO icw = 1, ncase_win
   if (h_cw(icw).GE.const(15)) then  !   active value, there is DF at the window    
      if (ic_T0_cw(icw).eq.0) then   !   but before there was not, so update
         ic_T0_cw(icw) = 1           !   ic_T0(icw) to 1, i.e. active, and continue
         T0_cw   (icw) = time
      endif
      if (ic_T0_cw(icw).eq.1) then
         if (h_cw (icw).GT.hpeak_cw (icw)) hpeak_cw (icw) = h_cw (icw) ! select maxima
         if (q_cw (icw).GT.qpeak_cw (icw)) qpeak_cw (icw) = q_cw (icw)
         if (vn_cw(icw).GT.vnpeak_cw(icw)) vnpeak_cw(icw) = vn_cw(icw)
         Qtot_cw  (icw) = Qtot_cw   (icw) + t_print_reset *q_cw(icw)
      endif
   else
      if (ic_T0_cw(icw).EQ.1) then        ! DF just finished, set time control to -999
         ic_T0_cw(icw) = -999
         Tf_cw (icw) = time
      endif
   endif
   if (is_target.eq.1) then
      if (dist_min_cw (icw).LT.dstminpeak_cw(icw)) then 
         dstminpeak_cw(icw) = dist_min_cw(icw)
      endif
   endif
ENDDO !  loop ncase_win

END SUBROUTINE Output_cases_win_SW_tst

!-------------------------------------------------------------------

       Subroutine Output_cases_win_SW
       
!-------------------------------------------------------------------

!      ------  writes window info 

implicit none 
                                            !   to file RES_file_cw every time_print steps
real   (irk), allocatable:: h_cw(:)         !   h_DF or rho (averaged of points in window)
real   (irk), allocatable:: v_cw(:,:)       !   velocity
real   (irk), allocatable:: q_cw(:)         !   outflow trough region
real   (irk), allocatable:: vn_cw(:)        !   normal velocity to side a 
real   (irk), allocatable:: dist_min_cw(:)  !   normal velocity to side a 
integer(ink), allocatable:: in_win(:)       !   point in window

integer(ink) icw, ipoin, ip_in, np_win
character(3) spac                               ! is " , "
real   (irk) xi, eta, h_cw_MAX, dstmp           ! max in windows
real   (irk) dstmpx, dstmpy                     ! MP chgd 15 June 2022 for -ve dist
real   (irk) x0_cw(2), a, b, na(2), nb(2)       ! simpler than  direct structure
real   (irk) xtg (2)
integer(ink) is_target

if (.NOT.allocated(h_cw)) then
   allocate (h_cw(ncase_win), v_cw(ndimn,ncase_win), q_cw(ncase_win), vn_cw(ncase_win) )
   allocate (in_win(ncase_win), dist_min_cw(ncase_win))
endif

spac = ' , '

DO icw = 1, ncase_win   !   loop in active windows 
   np_win = 0; h_cw(icw) = 0.; v_cw(:,icw) = 0.; q_cw(icw) = 0.; vn_cw(icw) = 0. ; dist_min_cw(icw) = 1.e18
   x0_cw     = info_cw(icw) % x0    !  recover info from structure, case icw of ncase_win
   a         = info_cw(icw) % ab(1)
   b         = info_cw(icw) % ab(2)
   na        = info_cw(icw) % nva 
   nb        = info_cw(icw) % nvb 
   is_target = info_cw(icw) % is_target
   Do ipoin = 1, npoin
      ip_in = 0   ! default, point is not in
      xi  = (x(1,ipoin)-x0_cw(1)) * na(1)   !  determine if ipoint is in a X b
      eta = (x(1,ipoin)-x0_cw(1)) * nb(1)
      if (ndimn.eq.2) then
         xi  = xi  + (x(2,ipoin)-x0_cw(2)) * na(2)
         eta = eta + (x(2,ipoin)-x0_cw(2)) * nb(2)
      endif
      if (ndimn.eq.2) then
          if ( (xi.ge.0.0).and.(xi.le.a ).and.(eta.ge.0.0).and.(eta.le.b ) ) ip_in = 1
      elseif (ndimn.eq.1) then
          if ( (eta.ge.0.0).and.(eta.le.b ) ) ip_in = 1
      endif   
      if (ip_in.eq.1)   then     !  point is inside
         np_win = np_win +1  ! increment number of points in rectangle a X b
         if (allocated(h_DF) ) then
            h_cw (icw) = h_cw (icw) + h_DF(ipoin)
         else
            h_cw (icw) = h_cw (icw) + rho (ipoin)
         endif
         v_cw (:,icw) = v_cw (:, icw) + vx(:,ipoin)
      endif
      if (is_target.eq.1) then
          dstmp = 0.0  ; dstmpx = 0.0 ; dstmpy=0.0
!          xtg (1) = x0_cw (1) + a*na(1)/2.                 MP chg 15th June 2022
          xtg (1) = x0_cw (1) + a*na(1)/2. + b*nb(1)/2.   ! MP chg 15th June 2022
          dstmpx = (x(1,ipoin)-xtg(1))
!          dstmp = dstmpx**2.
          if (ndimn.eq.2) then
!               xtg (2) = x0_cw (2) + b*nb(2)/2.
               xtg (2) = x0_cw (2) + a*na(2)/2. + b*nb(2)/2.
               dstmpy   = (x(2,ipoin)-xtg(2))
!               dstmp = dstmp + dstmpy**2.   ! MP chg 15th June 2022
          endif
!          dstmp =dstmp**0.5                       !    MP chg 15th June 2022
          dstmp =  dstmpx*nb(1)  + dstmpy*nb(2)    !    MP chg 15th June 2022
!          if (ip_in.eq.1) dstmp = 0.111111        !   OUT when points reach target, d=0.1111
          if (dstmp.lt. dist_min_cw(icw)) then
              dist_min_cw(icw) = dstmp
          endif  
      endif
   enddo       !    ------- ipoin loop
   if (np_win.GT.0 ) then
      h_cw (icw)   = h_cw (icw) / np_win
      v_cw (:,icw) = v_cw (:,icw) / np_win
      vn_cw (icw)  = - v_cw(1,icw)* nb(1)
      if (ndimn.eq.2) then 
         vn_cw (icw) = vn_cw (icw) - v_cw(2,icw)* nb(2)
      endif 
      q_cw (icw) = vn_cw (icw) * h_cw (icw) * a 
   endif
ENDDO          !    ------- icw loop     

!       We have obtained h_sw, v_cw, vn_cw and q_cw at all ncase_win windows at t = time
!       time to write to file filename.CaseWin.RES (unit is res_file_cw) if h diff from zero

is_target = info_cw(ncase_win) % is_target
if (is_target.eq.1) then          !   target window print always dist 
!    write(RES_file_cw,1001) (  (time, h_cw(icw),vn_cw(icw),q_cw(icw), dist_min_cw(icw)),icw= 1, ncase_win )
    write(RES_file_cw,1001)     time,  (( h_cw(icw),vn_cw(icw),q_cw(icw)),icw= 1, ncase_win) , dist_min_cw(ncase_win)  
else
   h_cw_MAX = 0.0
   DO icw = 1, ncase_win
      if (h_cw(icw).GT. h_cw_MAX) h_cw_MAX = h_cw(icw) 
   enddo
   if (h_cw_MAX.GE.const(15)) then   !   write only when max h is GT threshold C15 *****
      write(RES_file_cw,1001) (  (time, h_cw(icw),vn_cw(icw),q_cw(icw)),icw= 1, ncase_win )
   endif 
endif
1001 format (20 (G10.4,2x))

!  from current variables, update 
!  T0_cw, Tf_cw, hpeak_cw, vnpeak_cw, qpeak_cw, QTotal_cw (ncase_win)

DO icw = 1, ncase_win
   if (h_cw(icw).GE.const(15)) then  !   active value, there is DF at the window    
      if (ic_T0_cw(icw).eq.0) then   !   but before there was not, so update
         ic_T0_cw(icw) = 1           !   ic_T0(icw) to 1, i.e. active, and continue
         T0_cw   (icw) = time
      endif
      if (ic_T0_cw(icw).eq.1) then
         if (h_cw (icw).GT.hpeak_cw (icw)) hpeak_cw (icw) = h_cw (icw) ! select maxima
         if (q_cw (icw).GT.qpeak_cw (icw)) qpeak_cw (icw) = q_cw (icw)
         if (vn_cw(icw).GT.vnpeak_cw(icw)) vnpeak_cw(icw) = vn_cw(icw)
         Qtot_cw  (icw) = Qtot_cw   (icw) + t_print_reset *q_cw(icw)
      endif
   else
      if (ic_T0_cw(icw).EQ.1) then        ! DF just finished, set time control to -999
         ic_T0_cw(icw) = -999
         Tf_cw (icw) = time
      endif
   endif
   if (is_target.eq.1) then
      if (dist_min_cw (icw).LT.dstminpeak_cw(icw)) dstminpeak_cw(icw) = dist_min_cw(icw)
   endif
ENDDO !  loop ncase_win

END SUBROUTINE Output_cases_win_SW

!-------------------------------------------------------------------

       Subroutine Input_cases_WIN_SW
       
!-------------------------------------------------------------------

!      ------  Reads basic control data for Depth Integrated problems

implicit none

integer (ink) len1
integer (ink) icase_win, idimn
real    (irk) diffab0(2), cosTh, sinTh
1005 format (A60)

len1 = len_trim(SPH_SW_problem_name)
open(res_file_cw,  file=SPH_SW_problem_name(1:len1)//'.CaseWin.res'   )

read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*)    ncase_win
write(chk_file,*)    ncase_win  

nwin_vars = 7      !   ------  to dimension all windows as target  MP 20 July 2022

if (.NOT.allocated (ic_T0_cw) ) then

   allocate (ic_T0_cw(ncase_win) ) ! we will use to obtain t0 for outflow trough section
                                   !  0 still not arrived 1 active 999 finished
   allocate ( T0_cw     (ncase_win) ) 
   allocate ( Tf_cw     (ncase_win) ) 
   allocate ( qpeak_cw  (ncase_win) )  
   allocate ( vnpeak_cw (ncase_win) ) ;  
   allocate ( hpeak_cw  (ncase_win) )  
   allocate ( Qtot_cw   (ncase_win) )  
   allocate ( dstminpeak_cw (ncase_win) )    

endif

! ndat_out = 3 + ncases_mat + 2 + 6*(ncase_win-1) + 7    ! CHGD MP 20 Julty 2022
ndat_out = 3 + ncases_mat + 2 + ncase_win*nwin_vars      ! nr of columns in OC_mat 
                                                         ! last arget window has 7 vars
! ndat_out  = 3 + ncases_mat + 2 + 6*ncase_win        
 
ic_T0_cw  = 0   ;  T0_cw     = 0.0  ; Tf_cw     = 0.0
qpeak_cw  = 0.0 ;  vnpeak_cw = 0.0  ; hpeak_cw  = 0.0 ; Qtot_cw = 0.0 ; dstminpeak_cw = 1.e18

if (ndat_out.GT.Mdat_out) then
   write (*,*) ' --------------------------------'
   write (*,*) '  max number ofcases and windows too large ', ncases_mat, ncase_win
   write (*,*) ' recommended 6 cases and 2 windows max ...or change code line 327 of driver '
   write (*,*) ' enter any character to stp '
   read  (*,*) text
   STOP
endif

if (.NOT.allocated (info_cw ) ) then
   allocate (info_cw(ncase_win))   ! allocate the window structure for ncase_win windows
endif

DO icase_win = 1, ncase_win
   read (dat_file,1005) text
   write(chk_file,1005) text 
   read (dat_file,*) (info_cw(icase_win)% x0 (idimn), idimn = 1,2 )  ! (x0,y0) origin
   write(chk_file,*) (info_cw(icase_win)% x0 (idimn), idimn = 1,2 ) 
   read (dat_file,1005) text
   write(chk_file,1005) text 
   read (dat_file,*) (info_cw(icase_win)% ab (idimn), idimn = 1,2 )  !  sides a, b
   write(chk_file,*) (info_cw(icase_win)% ab (idimn), idimn = 1,2 )
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)  info_cw(icase_win)% Theta                          !  angle a with OX
   write(chk_file,*)  info_cw(icase_win)% Theta  
   cosTh = cos (info_cw(icase_win)% Theta)
   sinTh = sin (info_cw(icase_win)% Theta)
   info_cw(icase_win)% nva (1) =  cosTh
   info_cw(icase_win)% nvb (1) = -sinTh
   !if (ndimn.eq.2) then
      info_cw(icase_win)% nva (2) =  sinTh
      info_cw(icase_win)% nvb (2) =  cosTh
   !endif
   info_cw(icase_win)% is_target = 0
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)  info_cw(icase_win)% is_target                 !  1 if of target type
   write(chk_file,*)  info_cw(icase_win)% is_target 
enddo

END SUBROUTINE Input_cases_WIN_SW

!-------------------------------------------------------------------

       Subroutine Input_SW
       
!-------------------------------------------------------------------

!      ------  Reads basic control data for Depth Integrated problems

implicit none

logical conected

integer(ink) len1,  nline, icunkno   !  , nfric made global SW var 7th dec 2017 MP 
integer(ink) iline, i , j, itotv, ipoin, idimn, iiabs, iivn0, ic_abss
integer(ink) ipoin_w, my_iostat, ip_ts

integer(ink) ic_Sat             ! if FS activated and rel sat h (x,t) MP march 2023
integer(ink) ic_fct_dens        !  march 2023 MP  used to set ic_FS = 1 C20

integer(ink) BCt1,BCt2          ! type of BCs at bottom and top 1 Dirich 2 N
real   (irk) BCv1,BCv2          ! values of BCs 
integer(ink) pwp0_type          ! type of law for init pwp(z) 0 ct 1 lin 2 cos
integer(ink) iaux               !         value at bottom is hrelpw*h
real   (irk) xi, hh             ! z/h xi=0 at bottom 1 at top of st layer

real   (irk) w_drum_rpm         ! auxiliar, revs per min
real   (irk) tmp_avail_eros     ! will be used to initialize topol(15

real   (irk) dens, denss, densw, nporos, hhs, hhw, porostmp
                                    !  CHGD MP 8 April 2020
integer (ink) icase_mat, index_mat  !  dummmy indexes for multiple case computation
real    (irk) Cval_mat , Cval0_mat  !  see after readind const line 414 approx      

real(irk), allocatable:: coord_cg(:)     ! center of mass coords for expansion of mass
real(irk), allocatable:: shift_case(:)   ! shift dx dy for cases analysis  MP 12th May 2020

real    (irk) tmp_val

character(72) text
!              Assign units to files, ask for a name and open them

dat_file     = 10
chk_file     = 21
restart_file     = 22   ! defined in main vars module
restart_TOP_file = 24   ! defined in main vars module  MP 06 06 19
res_file1    = 31
res_file2    = 32
res_file3    = 33
res_file4    = 34
res_file5    = 35
res_file6    = 36   ! MP 17 March 2020 CHGD
res_file_cw  = 37   ! MP for window output cases 15th May 2020. See imput_win_cases SUB
res_file7    = 38   ! Saeidmt 12Nov2023
gid_msh      = 41
gid_res      = 42 
QGIS_asset   = 45   !  mp 2023 11 09 QGIS files
QGIS_avg     = 46
QGIS_peaks   = 47


!write(*,*) ' ------------------------------------------ '
!write(*,*) ' '   
!write(*,*) '        SPH INPUT...reading file.dat name '
!write(*,*) ' '   
!read (MASTER_dat,1005) text
!write(MASTER_chk,1005) text
!read (MASTER_dat,1006) SPH_SW_problem_name
!write(MASTER_chk,1005) SPH_SW_problem_name 
!write(*,*) ' sph file name for data is    ', SPH_SW_problem_name
!write(*,*) ' '
!write(*,*) '   '
!write(*,*) ' ------------------------------------------ '
!write(*,*) ' '
!write(*,*) '   '

1005 format (a60)
1006 format (a16)

!len1 = len_trim(SPH_SW_problem_name)
len1 = len_trim (problem_name_Global)    !  BUG in platform excel file output MP 7 Nov 2023  ************
SPH_SW_problem_name = problem_name_Global

open( dat_file,  file=SPH_SW_problem_name(1:len1)//'.dat'     )
open( chk_file,  file=SPH_SW_problem_name(1:len1)//'.chk'     )
open(res_file1,  file=SPH_SW_problem_name(1:len1)//'.1.res'   )
open(res_file2,  file=SPH_SW_problem_name(1:len1)//'.2.res'   )
open(res_file3,  file=SPH_SW_problem_name(1:len1)//'.3.res'   )
open(res_file4,  file=SPH_SW_problem_name(1:len1)//'.4.res'   )     !  ********************  DEBUG
!    INQUIRE ( res_file4, OPENED= conected)               
!    if (.NOT.conected) then                                         !  ********************  DEBUG
!       open(res_file4,  file=SPH_SW_problem_name(1:len1)//'.4.res' )!  ********************  DEBUG
!   endif
open(res_file5,  file=SPH_SW_problem_name(1:len1)//'.5.res'   )
open(res_file6,  file=SPH_SW_problem_name(1:len1)//'.6.res'   )
open(res_file7,  file=SPH_SW_problem_name(1:len1)//'.7.res'   ) ! Saeidmt 12Nov2023

open(  gid_msh,  file=SPH_SW_problem_name(1:len1)//'.post.msh')
if (ic_Cases_mat.GT.0) then              ! MP May 2020 Cases. We will append at the end of the file each final situation for the new case
   INQUIRE ( gid_res, OPENED= conected)
   if (.NOT.conected) then
      open( gid_res,  file=SPH_SW_problem_name(1:len1)//'.post.res')
      write (gid_res,*) 'GiD Post Results File 1.0'
   endif
   INQUIRE ( QGIS_asset, OPENED= conected)
   if (.NOT.conected) then
      open  (QGIS_asset,  file=SPH_SW_problem_name(1:len1)//'.QGIS_assets')  ! MP 2023 11 09 QGIS
      open  (QGIS_avg,    file=SPH_SW_problem_name(1:len1)//'.QGIS_avg') 
      open  (QGIS_peaks,  file=SPH_SW_problem_name(1:len1)//'.QGIS_peaks') 
   endif
else
   open( gid_res,  file=SPH_SW_problem_name(1:len1)//'.post.res')
   write (gid_res,*) 'GiD Post Results File 1.0'
   open  (QGIS_asset,  file=SPH_SW_problem_name(1:len1)//'.QGIS_assets') 
   open  (QGIS_avg,    file=SPH_SW_problem_name(1:len1)//'.QGIS_avg') 
   open  (QGIS_peaks,  file=SPH_SW_problem_name(1:len1)//'.QGIS_peaks') 
endif


write(chk_file,*) ' Attention*** '
write(chk_file,*) 'If you are using frictional, a negative Manning is new Mu_CF=50 Pas2?'
write(chk_file,*)  '     '
write(chk_file,*)  '     '

 
!      ------ Initialize some vars 

pi        = 4*atan(1.0)
npoin_w   = 0
npoin_s   = 0
npoin_vps = 0
nabs      = 0
nvn0      = 0

ic_sgates      = 0

n_hist         = 0
ic_histogram   = 0  

!      ------  Ask whether we are solving an standard or a general problem

ic_WIR = 1
write(chk_file,*)  ' ic_WIR = ', ic_WIR

!      ------  read lines of text describing the problem 

nline = 0
read (dat_file,1005, IOSTAT = my_iostat) text
if (my_iostat.LT.0) then
    close(dat_file)
    open( dat_file,  file=SPH_SW_problem_name(1:len1)//'.dat'     )
else
     rewind(dat_file)
endif 

read (dat_file,*)  nline
write(chk_file,*)  nline    
do iline = 1,nline
   read (dat_file,1005) text
   write(chk_file,1005) text
enddo
!
! ic_SWAlg algorithm for SW computations. There are several:
!       0 for landslide algorithm p=0.5 g h2  
!       1 for coastal   SW        p=0.5 g (h2-H2) 
!       2 for rotating drum
!       3,4 for DFs w/o and with pwp
!       5,6 for Shiomi model, grad n is not computed

read  (dat_file,1005) text
write (chk_file,1005) text  
read  (dat_file,*)    ic_SWAlg    
write (chk_file,*)    ic_SWAlg 

!      ------ In case ic_SWAlg=2 we have a 1D rotating drum

if (ic_SWAlg.eq.2) then

   read  (dat_file,1005) text
   write (chk_file,1005) text  
   read  (dat_file,*)    R_drum, w_drum_rpm    
   write (chk_file,*)    R_drum, w_drum_rpm 
   w_drum_radps = w_drum_rpm*pi/30.
   ic_SWAlg     = 0

endif 

!      ------  Ask whether we have histogram type input

ic_histogram   = 0
rndm_fact      = 0.0  ! MP August 2022
ic_random_walk = 0
ic_w_inflow      = 0  ! MP July 2023 read water for mix with soil

write(chk_file,*) '  Input nr of gates with histograms, 0 otherwise AND random wak fct'
read  (dat_file,1005) text
write (chk_file,1005) text  
if (SPH_t_Integ_Alg.EQ.100.or.SPH_t_Integ_Alg.EQ.110 ) then
   ic_random_walk = 1
    read  (dat_file,*)    n_hist, rndm_fact
    write (chk_file,*)    n_hist, rndm_fact  
else
   read  (dat_file,*)    n_hist  ! chuan added if_inject_new
   write (chk_file,*)    n_hist       
endif

if (n_hist.ne.0)  ic_histogram = 1

npoin_in_model = 0  
                      !  -----------------------------------------------------
ic_slow = 0           ! default     

if (SPH_t_Integ_Alg.EQ.10 ) then     !   read info of time curves for rain intensity and 
   ic_slow = 1                       !   skipping search routine calls
   ic_time_series = 0; ic_search_step = 0
   read  (dat_file,1005) text
   write (chk_file,1005) text  
   read  (dat_file,*)    ic_time_series, tscale_code    ! Changed J+M+S march 2024
   write (chk_file,*)    ic_time_series, tscale_code 
   
   nsearch_step = 0
   If (ic_search_step.eq.1) then
      read  (dat_file,1005) text
      write (chk_file,1005) text  
      read  (dat_file,*)    nsearch_step 
      write (chk_file,*)    nsearch_step    
   endif    

   IF (ic_time_series.eq.1) then
      read  (dat_file,1005) text                !  description of time curve
      write (chk_file,1005) text 
      read  (dat_file,1005) ts_Rain % description    
      write (chk_file,1005) ts_Rain % description 
       
      read  (dat_file,1005) text     ! number of values in t curve  ! units: secs in units of time series                       
      write (chk_file,1005) text     ! (1, 60, 3600, 86400) Type_Val 1 Rain Intens mm/h 2 Pwp abs(Pa) 3 rel(0,1)
      read  (dat_file,*)    ts_Rain % np_ts, ts_Rain % unit_ts, ts_Rain % type_Val             
      write (chk_file,*)    ts_Rain % np_ts, ts_Rain % unit_ts, ts_Rain % type_Val

      !if (.NOT.allocated (ts_Rain % time_ts)) then
          allocate ( ts_Rain % time_ts(ts_Rain % np_ts) ) 
          allocate ( ts_Rain % val_ts (ts_Rain % np_ts) ) 
      !endif   
          
      read  (dat_file,1005) text                 !  series of time arain intensity values
      write (chk_file,1005) text  
      read  (dat_file,*)    (ts_Rain % time_ts(ip_ts), ip_ts = 1, ts_Rain % np_ts) 
      write (chk_file,*)    (ts_Rain % time_ts(ip_ts), ip_ts = 1, ts_Rain % np_ts) 
      read  (dat_file,*)    (ts_Rain % val_ts (ip_ts), ip_ts = 1, ts_Rain % np_ts) 
      write (chk_file,*)    (ts_Rain % val_ts (ip_ts), ip_ts = 1, ts_Rain % np_ts)

      read  (dat_file,1005) text                                   !  if_pwp 0 (no pwp induced) 1 yes
      write (chk_file,1005) text                                   !  if_drain 0 (no pwp induced) 1 yes
      read  (dat_file,*)    ts_Rain % if_pwp , ts_Rain % if_drain  !  1 if drainage is hrz trough top layer   
      write (chk_file,*)    ts_Rain % if_pwp , ts_Rain % if_drain  !  2 if drainage is hrz trough basal layer
      
      read  (dat_file,1005) text                ! h of basal drainage layer(only ifdrain=2), L length   
      write (chk_file,1005) text                !  tanth slope (avg), kw horiz permeability 
      read  (dat_file,*)    ts_Rain %h_basal, ts_Rain % L_basal , ts_Rain %tanTh_basal      !  
      write (chk_file,*)    ts_Rain %h_basal, ts_Rain % L_basal , ts_Rain %tanTh_basal
      
      read  (dat_file,1005) text                ! horiz perm and porosity   
      write (chk_file,1005) text                ! horiz perm and porosity
      read  (dat_file,*)    ts_Rain % kw_basal, ts_Rain %poros_basal        !  
      write (chk_file,*)    ts_Rain % kw_basal, ts_Rain %poros_basal


   endif  

ENDIF
 
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) ndimn   
write(chk_file,*) ndimn  
    
!      ------  Control parameters

read (dat_file,1005) text
write(chk_file,1005) text  
read (dat_file,*)  ic_soil, ic_water, ic_vps, ic_abss 
write(chk_file,*)  ic_soil, ic_water, ic_vps, ic_abss 

ic_ws_Interact = 0                          !  1 means s and w interact 
ic_DF = 0 ; ic_FS = 0                       ! defaults   
if (ic_soil.eq.1.AND.ic_water.eq.1) then 
   ic_ws_Interact = 1           ! Waves in rservoirs 
   ic_gradn = 0  
   if     (ic_SWAlg.eq.3) then
       ic_gradn       = 1
       ic_ws_Interact = 2; ic_DF = 1
       ic_SWAlg       = 0
   elseif (ic_SWAlg.eq.4) then   !  If ic_SWAlg=4 -> DFs + pwp Set pwp, read later icpwp
       ic_gradn       = 1
       ic_ws_Interact = 2 ; ic_DF = 1        
       icpwp          = 3
       ic_SWAlg       = 0
   elseif (ic_SWAlg.eq.5) then   !   
       ic_gradn       = 0
       ic_ws_Interact = 2 ; ic_DF = 1
       ic_SWAlg       = 0
   elseif (ic_SWAlg.eq.6) then   !  If ic_SWAlg=6 -> DFs + pwp Set pwp, read later icpwp
       ic_gradn       = 0
       ic_ws_Interact = 2 ; ic_DF = 1        
       icpwp          = 3
       ic_SWAlg       = 0
   endif     
endif   

if     (ic_abss.eq.0) then
           ic_abs = 0 ; ic_vn0 = 0
elseif (ic_abss.eq.1) then
           ic_abs = 1 ; ic_vn0 = 0
elseif (ic_abss.eq.2) then
           ic_abs = 1 ; ic_vn0 = 1
elseif (ic_abss.eq.3) then
           ic_abs = 0 ; ic_vn0 = 1
elseif (ic_abss.eq.4) then            !   new water permeable walls MP 20th dec 2018
           ic_abs = 0 ; ic_vn0 = 2
elseif (ic_abss.eq.10) then           !   new input for walls MP 28 Nov 2020
           ic_abs = 0 ; ic_vn0 = 10
endif  
    
if (ic_soil.eq.1)  Call Get_soil_points      !      ------  Read points in avalanching soil mass
    
if (ic_water.eq.1) Call Get_water_points     !      ------  Read points in reservoir or fjord
    
if (ic_vps.eq.1)   Call Get_vps_points       !      ------  Read virtual particles
    
if (ic_abs.eq.1)   call Read_AbsBcs          !      ------  Read Abs BCs    
    
if (ic_vn0.GE.1)   call Read_vn0Bcs          !      ------  MP 20th dec 18   
                                             !              Read vn0 BCs 1 for wall-solid &w  2 water passes 
npoin  = npoin_w + npoin_s 
nvirt  = npoin_vps
ntotal = npoin + nvirt + nabs + nvn0     
                                         !  Allocate global arrays,    
Call Setup_Global_Arrays_SW              !  pass local to global, deallocate, initialize
                                         !  deal with h_DF, n_DF, h_sw, v_SW                                          
!      ------  we start here case analysis for h (code 501), expansion factor (506)

if (.not.allocated (coord_cg))   allocate (coord_cg(ndimn))  !      ! for case studies
if (.not.allocated (shift_case)) allocate (shift_case(ndimn))

if (ic_geomm_pts.eq.1) then               !  MP CHGD June 2023 expand pts
   if (h_expans.GT.1.e-4) then
      Do ipoin = 1, npoin
         rho (ipoin) = rho(ipoin) * h_expans
         mass(ipoin) = mass(ipoin)* h_expans
      enddo 
      rho0 = rho ; rho00 = rho; mass0 = mass       
   endif 
   if (mass_expans.GT.1.e-4) then
       coord_cg = 0.0
       do ipoin = 1, npoin
          coord_cg(:) = coord_cg(:) + x(:,ipoin)
          mass(ipoin) = mass(ipoin)*mass_expans
       enddo
       coord_cg = coord_cg/npoin
       do ipoin = 1, npoin
          x(:,ipoin) = coord_cg(:) + mass_expans * (x(:,ipoin)-coord_cg(:))
          mass(ipoin) = mass(ipoin) * mass_expans  
       enddo  
       x0 = x; x00 = x ; mass0 = mass
   endif 
   tmp_val = tmp_val + abs(mass_shift(1))
   if (ndimn.eq.2) tmp_val = tmp_val + abs(mass_shift(2))
   if (tmp_val.GT.1.e-4) then  !   ------  shift of the whole mass case MP June 2023
      do ipoin = 1, npoin
         x(:,ipoin) = x(:,ipoin) + mass_shift(:)
      enddo
   endif 
endif

IF (ic_cases_mat.eq.1) THEN    ! new h cases
    icase_mat = OC_mat (nproblems_index, 1)
    index_mat = Ccases_mat ( icase_mat, 1) 
    if (index_mat.eq.501) then
        Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
        Do ipoin = 1, npoin
           rho (ipoin) = rho(ipoin) *Cval_mat
           mass(ipoin) = mass(ipoin)*Cval_mat
        enddo 
        rho0 = rho ; rho00 = rho; mass0 = mass       
    endif
    if (index_mat.eq.511) then
        Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)  ! expansion fact
        coord_cg = 0.0
        do ipoin = 1, npoin
           coord_cg(:) = coord_cg(:) + x(:,ipoin)
           mass(ipoin) = mass(ipoin)*Cval_mat
        enddo
        coord_cg = coord_cg/npoin
        do ipoin = 1, npoin
           x(:,ipoin) = coord_cg(:) + Cval_mat* (x(:,ipoin)-coord_cg(:))
           mass(ipoin) = mass(ipoin)*Cval_mat  
        enddo  
        x0 = x; x00 = x ; mass0 = mass
    endif
    if (index_mat.eq.521) then      !      ------  shift of the whole mass case
       shift_case(1) = OC_mat (nproblems_index, 3 +icase_mat)           !  this is DX
       if (ndimn.eq.2) then
          shift_case(2) = OC_mat (nproblems_index, 3 +ncases_mat + 1)   !  this is DY
       endif
        do ipoin = 1, npoin
           x(:,ipoin) = x(:,ipoin) + shift_case
        enddo
    endif

ELSEIF (ic_cases_mat.GE.2) THEN       !   change mat props according to cases MP 24-03 2022
   Do icase_mat = 1, ncases_mat
      index_mat = Ccases_mat ( icase_mat, 1)
      if (index_mat.eq.501) then
          Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
          Do ipoin = 1, npoin
             rho (ipoin) = rho (ipoin)*Cval_mat
             mass(ipoin) = mass(ipoin)*Cval_mat
          enddo   
          rho0 = rho ; rho00 = rho; mass0 = mass
      endif  
      if (index_mat.eq.511) then
          Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)  ! expansion fact
          coord_cg = 0.0
          do ipoin = 1, npoin
             coord_cg(:) = coord_cg(:) + x(:,ipoin)
          enddo
          coord_cg = coord_cg/npoin
          do ipoin = 1, npoin
             x(:,ipoin) = coord_cg(:) + Cval_mat* (x(:,ipoin)-coord_cg(:))
             mass(ipoin) = mass(ipoin)*Cval_mat    
          enddo  
          x0 = x; x00 = x ; mass0 = mass
      endif 

      if (index_mat.eq.521) then      !      ------  shift of the whole mass case
          shift_case(1) = OC_mat (nproblems_index, 3 +icase_mat)           !  this is DX
          if (ndimn.eq.2) then
             shift_case(2) = OC_mat (nproblems_index, 3 +ncases_mat + 1)   !  this is DY
          endif
          do ipoin = 1, npoin
            x(:,ipoin) = x(:,ipoin) + shift_case
          enddo  
          x0 = x; x00 = x ; mass0 = mass
      endif      
            
   enddo !  -------------  end loop icases         
ENDIF    !  -------------- end cases ic_case_mat

countiac = 0        
    
Xmin_Domain(1) = Xming
Xmax_Domain(1) = Xmaxg   
If (ndimn.eq.2) then
    Xmin_Domain(2) = Yming
    Xmax_Domain(2) = Ymaxg  
endif  
!      ------  read control parameters

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)  pa_sph, nnps, sle, skf 
write(chk_file,*)  pa_sph, nnps, sle, skf 
if ((ndimn.ne.2).and.(nnps.ge.3)) then ! rotation only for 2D
   nnps = 2                            ! CHGD MP 9th April 2020                       
endif
!      ------  read constants
        
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) summation_density, average_velocity, virtual_part , nor_density
write(chk_file,*) summation_density, average_velocity, virtual_part , nor_density

ic_semi_implicit = 0
if (SPH_t_Integ_Alg.eq.7.or.SPH_t_Integ_Alg.eq.8) then
    ic_semi_implicit = 1     
endif

if (ic_semi_implicit.EQ.1) then
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)  nor_div_FS, nor_gradP_FS, tolcg_MAIN
   write(chk_file,*)  nor_div_FS, nor_gradP_FS, tolcg_MAIN
endif

!      ------ If SPH_t_Integ_Alg.eq.6 -> semi_implicit  !  ** MP s_implicit 17.03.18
    
!      ------  Allocate and Read material constants
!               1..cgra       2..dens      3..cmanning    4..Erosion     5..nfrict   
!               6..Tauy0      7.. constK   8..visco       9..tanfi8      10..hfrict0  
!              11..          12..tanfi0   13..Bfact      14..h_pwp_rel   15.. Comp
!
!              We will read the 2nd group when C11 < 0 i.e. -999
!
!              16  iccosTh   17.. denss   18..densw      19  ----        20 ic_fct_dens
!              21  eros law  22   eros1   23  eros2      24  eros3       25 eros4
!              26  perm law  27   c1      28  c2         29              30  if -999 read another block
!
!              We will read the 2nd group when C30 < 0 i.e. -999
!
!              31  Kv skeleton vol stiff 32      33     34 ic_hrelSat.....45 min_porosity   
!              36..law_crush   37..r0_crush 38..rf_crush   39..B_crush  45 if -99 read
!
!              46  n0_eros        47  alpha_slip factor   48       49    
!              50   interface erosion law   51  Hungr extension   52 -> 60 free
!               
nconst = 110               ! before it was 15 and 30 and 60 C102 cohes, C103 expf, c104 mu. c105 muespf  
!
if (.NOT.allocated(const)) then 
        allocate (const(nconst))
        const = 0.0
endif

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) (const(i), i = 1,15)
write(chk_file,*) (const(i), i = 1,15)
if (const(4).gt.1.e-4) then   ! keeps compatibility with old versions
     const(21) = 1            ! allows us to use current data with 
     const(22) = const(4)     ! old files
endif

if (const(11).lt.0) then
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) (const(i), i = 16,30)
    write(chk_file,*) (const(i), i = 16,30)
endif    
if (const(30).lt.0) then
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) (const(i), i = 31,45)
    write(chk_file,*) (const(i), i = 31,45)
endif
   
if (const(45).lt.0) then          ! Changed MP 3rd April 2018
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) (const(i), i = 46,60)
    write(chk_file,*) (const(i), i = 46,60)
endif

ic_fct_dens = const(20)             ! we set unsat FS when fact dens and ic_hrelSat.Eq.1 
ic_hrelSat  = const (34)            !  MP march 2023
if (ic_fct_dens.eq.1.AND.ic_hrelSat.eq.1) ic_FS = 1

!      ------   Modify constants if we are running multiple case data CHGD MP 8April 2020

IF (ic_cases_mat.eq.1) THEN    ! CHGD MP 8April 2020
    icase_mat = OC_mat (nproblems_index, 1)
    Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat) 
    Cval0_mat = Ccases_mat ( icase_mat, 2)
    index_mat = Ccases_mat ( icase_mat, 1)
    if (icase_mat.LE.nconst) then               ! prepare options 501 (h fact), etc 
       Call update_const_SW ( 1, index_mat, Cval_mat, Cval0_mat)  
    endif
ELSEIF (ic_cases_mat.GE.2) THEN       !   change mat props according to cases MP 24-03 2022
   Do icase_mat = 1, ncases_mat
      index_mat = Ccases_mat ( icase_mat, 1)
      
      Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
      if (index_mat.LE.nconst) then               ! prepare options 501 (h fact), etc
         IF (index_mat.eq.2) THEN                 !   change porosity to density of the mixture
            if (Cval_mat.LE.1.01) then            !   we give poros in table. Need to get dens
               if (const(17).ge.1.e-4) then
                  const(2) = (1.-Cval_mat)*const(17) +Cval_mat*Const(18)
               else
                  write(*,*) ' *** provide dens solid and water. Type any number '
                  read  (*,*) Cval0_mat           ! just to input a keystroke
                  stop
               endif
            else                         ! we are providing density in table
               Const (index_mat) = Cval_mat      
            endif
         ELSE                             ! general material property
            Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
            Const (index_mat) = Cval_mat
         ENDIF
      endif              ! end prepare options 501 (h fact), etc                 
   enddo          
ENDIF   
!      ------   Obtain nfric and
!                      if_frictional  (= 1 for frictional laws: 5 to 10, 26 and 27
!                      if_CV2byR      ( =1  for 7,8,10, 26 and 27)
    
nfric = const(5) + 0.01
if_frictional = 0
if_V2byR      = 0
if ( ((nfric.ge.5).AND.(nfric.le.10)).OR.(nfric.eq.26).OR.(nfric.eq.27) ) then
    if_frictional = 1
     if ( (nfric.eq.7).OR.(nfric.eq.8).OR.(nfric.eq.10) &
          .OR.(nfric.eq.26).OR. (nfric.eq.27)) then
          if_V2byR = 1
     endif
endif  
                         !  Model is Julien with variable (1-n) C102 cohes_0   MP 21st July 2023
if (nfric.eq.201) THEN   !  C102 cohes_0  C103 C_Expf  C104 Mu0  C105 Mu_expf   
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) (const(i), i = 102,105)
    write(chk_file,*) (const(i), i = 102,105)   
    ic_w_inflow      =  1            
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) w_inflow%flow_q , w_inflow%Zmin, w_inflow%Zmax
    write(chk_file,*) w_inflow%flow_q , w_inflow%Zmin, w_inflow%Zmax  
ELSEIF(nfric.eq.202.OR.nfric.eq.203) THEN   !  C101 tanFi C102 xi_0  C103 xi_expf  MP Nov 19 2023
    if (.not.allocated(ns_dep_props)) then   !  C101 tanFi C102 xi_0  C103 xCv max 0.615
       allocate (ns_dep_props(1,npoin))
       ns_dep_props = 0.0
    endif
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) (const(i), i = 101,103)
    write(chk_file,*) (const(i), i = 101,103)  
    ic_w_inflow      =  1  
    if_frictional    =  1  
    if_V2byR         =  1          
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) w_inflow%flow_q , w_inflow%Zmin, w_inflow%Zmax
    write(chk_file,*) w_inflow%flow_q , w_inflow%Zmin, w_inflow%Zmax                                             
ENDIF
!
ic_erosion = 0                                               ! erosion control  MP 06 06 19
If ((const(21).gt.1.e-9).or.(ic_basal_erosion.eq.1)) then    ! 1st condition body erosion
    ic_erosion = 1                                           ! 2nd topo dependent erosion 
    if (ic_restart.eq.1) then
         call Get_Topo_restart (ic_erosion)                 ! in case erosion + restart read topo
    endif
    tmp_avail_eros = 999.
    if (const(21).gt.1.e-9) then             !  C21>0 means surface eros law present
       if (const(48).gt.1.e-6) then          !      C48>0 means we give max eros. depth
          tmp_avail_eros = const(48)         !          SET TEMP = c48
       else                                  !      C48 = 0 means we have surf erosion unlimited
          tmp_avail_eros = 999               !           Set temp = 999          
       endif                                 !
       topol(15,:) = tmp_avail_eros          ! ** CH Saeid January 2020 *****
    else                                     !  C21 = 0 no surf erosion law                                     
       tmp_avail_eros = 0.0                  !      set temp = 0
    endif                                    !
                                             !  if basal erosion, we have obtained T15 in some points 
    if (ic_basal_erosion.eq.1) then          !                    all others are zero
       where (topol(15,:).lt.1.e-6) topol(15,:) = tmp_avail_eros
    endif
endif              
!
ic_crush = 0                  !   Crushing law vars
if (const(36).gt.0.01) then
   ic_crush = 1
   if (.NOT.allocated (r_crush) ) allocate ( r_crush(npoin), w_crush(npoin) )
   r_crush = 0.0; w_crush = 0.0
endif  
!
ic_end_solid = 0              !  -------  solid behaviour laws
if (const(43).gt.0 ) then
   ic_end_solid       = 1  
   Law_end_solid      = const(42)   
   t1_end_solid       = const(43)     
   fct_exp_end_solid  = const(44)      
endif

!    ---------  case of interfacial erosion C50 law C51 Hungr constant

ic_ws_erosion = 0
if (ic_ws_interact.EQ.1.AND.const(50).gt.0) then           !  CH MP Feb 19
   ic_ws_erosion =  1 
   if (.NOT.allocated(h_DF)) then
      allocate (h_sw(ntotal), v_sw(ndimn,ntotal), h_DF(ntotal), n_wir(ntotal), dens_wir(ntotal))
   endif
   h_sw = 0.0; v_sw = 0.0; h_DF = 0.0; n_wir = 0.0; dens_wir = 0.0
   dens = const(2) ; denss = const(17); densw = const(18)
   n_wir =  (denss - dens)/(denss - densw)                 !  bug MP 25 May 2018
   dens_wir = dens  
endif

!    ---------  case of interfacial friction C55 law C56 mu  C57 Li mix length

ic_ws_friction = 0
if (const(55).gt.0) then
   ic_ws_friction =  1 
endif   
!
!      ------  gravity at inclined planes. OX along plane, downhill

gravx = 0.0                !      ------  default values
gravy = 0.0
gravz = const(1) 
if ( const(1).lt.0) then   !      ------ inclined plane set ic_costh = 2
    const (16) = 2
    read (dat_file,1005) text
    write(chk_file,1005) text
    read (dat_file,*) gravx, gravy, gravz
    write(chk_file,*) gravx, gravy, gravz
    const(1) = -const(1) 
endif     
            
!     ** we will limit h<comp as in ext forces. comp = const(15)

if (const(15).le.1.e-8 ) then
    write(*,*) ' const15 is comp should not be zero!!! '
    PAUSE
endif

if (const(10).le.1.e-8 ) then
    write(*,*) ' const10 is hfrict0 should not be zero!!! '
    PAUSE
endif  

ic_k0 = 0

if (if_frictional.eq.1) then       !  december 7th 2017 MP

   read  (dat_file,1005) text
   write (chk_file,1005) text 
   write (chk_file,*)  '   To use K active K pasive input 1,  otherwise 0 ' 
   read  (dat_file,*)    ic_K0    
   write (chk_file,*)    ic_K0  

   if (ic_k0.eq.1) then
      read  (dat_file,1005) text
      write (chk_file,1005) text  
      write (chk_file,*) '  Input active  << and pasive coefficients '
      read  (dat_file,*)    K0_active, K0_pasive    
      write (chk_file,*)    K0_active, K0_pasive 
   endif
   
endif

!      ------ if DFs normalize rho  DFs correction.  CH MP March 2019

if (summation_density.AND.ic_ws_interact.EQ.2) then  ! we normalize density here 
                                                     ! if ic_DF.eq.1 subst above
    if (nnps.eq.2) then
        call grid_find_new_TOOLS  (0)      ! CHGD MP april 2020
    elseif (nnps.eq.3) then                            ! CHGD MP april 2020
        call grid_find_new_TOOLS2 (0)      ! CHGD MP april 2020
    elseif (nnps.eq.4) then                            ! CHGD MP april 2020
        call grid_find_new_TOOLS3 (0)      ! CHGD MP april 2020
    endif
    call Pint_Update_SW_Tools                        ! ensures consistency with h_DF h_sw
    call Vn0_BCs_SW (0.0)           
    call sum_density_Tools    
    call DFs_initial_correction
    rho0  = rho                                     ! MP CH 24th may 2019
    rho00 = rho 
endif

ic_sat  = const(34) + 0.001
if (ic_FS.eq.1.and. ic_Sat.eq.1 ) then
   call FSs_initial
endif
                               
nUaux = 0                                            !  pwp info and correction
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) icpwp
write(chk_file,*) icpwp 

ic_Stef_eros = 0           !  default
if (icpwp.eq.20) then      !  this indicates ALE will be activated  June 2021
   icpwp = 2
   ic_stef_eros = 1
    if (.NOT.allocated(erosion_rate_V)) then
       allocate ( erosion_rate_V (npoin)); erosion_rate_V = 0.0
    endif
endif     

if (icpwp.ne.0) then 
                    
    if (ic_ws_interact.NE.2 ) then        ! we normalize hsbefore defining pwp
        if (summation_density) then
           if (nnps.eq.2) then
               call grid_find_new_TOOLS  (0)      ! ! needed for sum_density
           elseif (nnps.eq.3) then                ! CHGD MP april 2020
               call grid_find_new_TOOLS2 (0)      ! CHGD MP april 2020
           elseif (nnps.eq.4) then                ! CHGD MP april 2020
               call grid_find_new_TOOLS3 (0)      ! CHGD MP april 2020
           endif
           call Pint_Update_SW_Tools 
           call Vn0_BCs_SW (0.0)          ! changed 25th Feb 2020 MP         
           call sum_density_Tools                    
           rho0  = rho                    !  MP CH 23th may 19
           rho00 = rho
        endif
    endif

   call Initialize_Pwp                    !  March 2019

   IF (ic_cases_mat.eq.1) THEN    ! new pwp cases code is 401
       icase_mat = OC_mat (nproblems_index, 1)
       index_mat = Ccases_mat ( icase_mat, 1) 
       if (index_mat.eq.401) then
           Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
           Do ipoin = 1, npoin
              u (ipoin) = u(ipoin) * Cval_mat
              if (allocated (Uaux)) then
                 Uaux(:,ipoin) = Uaux(:,ipoin)* Cval_mat
              endif
           enddo 
           Uaux0 = Uaux     
       endif
   ELSEIF (ic_cases_mat.GE.2) THEN       !   change mat props according to cases MP 24 03 2022
      Do icase_mat = 1, ncases_mat
         index_mat = Ccases_mat ( icase_mat, 1)
         if (index_mat.eq.401) then
            Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
            Do ipoin = 1, npoin
               u (ipoin) = u(ipoin) * Cval_mat
               if (allocated (Uaux)) then
                  Uaux(:,ipoin) = Uaux(:,ipoin)* Cval_mat
               endif
            enddo 
            Uaux0 = Uaux     
         endif      
      enddo          
   ENDIF 

ENDIF     !   pwp     

write (chk_file,*) ' --- if you want to use a coarse mesh for output, type 1 '
read  (dat_file,1005) text
write (chk_file,1005) text  
read  (dat_file,*)    if_coarse_plt    
write (chk_file,*)    if_coarse_plt  

if (if_coarse_plt.eq.1) then
    fc_coarse_plt = 0.0
    write (chk_file,*) ' give factor to multiply topog deltxg  '
    read  (dat_file,1005) text
    write (chk_file,1005) text  
    read  (dat_file,*)    fc_coarse_plt    
    write (chk_file,*)    fc_coarse_plt
    DO while (fc_coarse_plt.lt.1)  
        write(*,*) ' wrong value, factor should be >= 1 code stopped '
        write(*,*) ' reenter it trough keyboard '
        read (*,*) fc_coarse_plt
    enddo      
    call setup_coarse_plt
endif  

!      ------  what follows is for check points

read  (dat_file,1005) text
write (chk_file,1005) text  
read  (dat_file,*)    if_chk_pts 
write (chk_file,*)    if_chk_pts  
write(chk_file,*) ' control for output points at which to check solution '

if (if_chk_pts.eq.1) then

   if (ndimn.eq.1) then
       write(*,*) ' -------  STOP: control points only 2D '
       STOP
   endif

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    nchk_pts
   write(chk_file,*)    nchk_pts
   if (.NOT.allocated (coor_chk_pts) ) allocate (coor_chk_pts(ndimn,nchk_pts))   !   REMEMBER to DEALLOCATE 
   
   read (dat_file,1005) text
   write(chk_file,1005) text
   do i = 1, nchk_pts
      read (dat_file,*) (coor_chk_pts(idimn,i),idimn=1,ndimn)
      write(chk_file,*) (coor_chk_pts(idimn,i),idimn=1,ndimn)
   enddo
   
endif
!
!      Filter for GID
!
!   We will use a mask to plot variables, nGidMaskSw = 12
!   Gid_Mask_SW (nGidMask) can be zero (no plot) or one (plot)
!   SW case:
!          
!       1 ... h soil       
!       2 ... displacement          
!       3 ... velocity          
!       4 ... Pwp          
!       5 ... erosion          
!       6 ... Z         
!       7 ... h soil relative  (used in Lausanne)         
!       8 ... h water          
!       9 ... eta          
!      10 ... hs + hw
!      11 ... hsml
!      12 ... full FD pwp?

nGidMaskSw = 12
if (.NOT.allocated (Gid_Mask_SW) ) allocate (Gid_Mask_SW(nGidMaskSw) )  !  ** remember to deallocate
read  (dat_file,1005) text
write (chk_file,1005) text 
read  (dat_file,*) ( Gid_Mask_SW(i),i=1,nGidMaskSw)
write (chk_file,*) ( Gid_Mask_SW(i),i=1,nGidMaskSw)


!     finally, for SW problems, obtain T_change_to_w  
!     stored in driver_time_vars

read  (dat_file,1005) text
write (chk_file,1005) text  
read  (dat_file,*)    T_change_to_W    
write (chk_file,*)    T_change_to_W 
  
Call Check_Out_Domain_SW

deallocate (coord_cg )

CONTAINS

!-------------------------------------------------------------------

       Subroutine FSs_initial 
       
!-------------------------------------------------------------------

!      ------  Flow slides with vriable degree of sat height MP March 2023 

implicit none

integer(ink) ipoin 
real   (irk) Sr, alphas, dens, denss, densw, dens_bar , dens_bar_eff 
real   (irk) fact_dens       
!      ------  FSs values.

dens       = const   (2) 
denss      = const  (17)
densw      = const  (18)
Sr         = const  (19)
alphas     = const  (14)
ic_hrelSat = const  (34)

dens_bar     = dens
dens_bar_eff = dens_bar - alphas*densw
fact_dens   = dens_bar_eff/dens_bar

if (.NOT.allocated (hrelSat_DF) )  then
   allocate (hrelSat_DF     (npoin))
   allocate (n_DF           (npoin))
   allocate (dens_bar_DF    (npoin)) 
   allocate (dens_bar_eff_DF(npoin)) 
   allocate (fact_dens_DF    (npoin)) 
   allocate (fact_pwp_max_DF (npoin) , fact_pwp_min_DF (npoin)) 
   allocate (fpwp_max_DF     (npoin) , fpwp_min_DF   (npoin))  
   fact_pwp_max_DF (:) = 0.0 
   fact_pwp_min_DF (:) = 0.0  
   fpwp_max_DF     (:) = 0.0 
   fpwp_min_DF     (:) = 0.0                                          
endif

hrelSat_DF      (:) = alphas
n_DF            (:) = (denss-dens)/(denss-alphas*densw)
fact_dens_DF    (:) = fact_dens
dens_bar_DF     (:) = dens_bar
dens_bar_eff_DF (:) = dens_bar_eff 

END Subroutine FSs_initial


!-------------------------------------------------------------------

       Subroutine DFs_initial_correction
       
!-------------------------------------------------------------------

!      ------   

implicit none

integer(ink) ipoin, ipoinw
real   (irk) dens, denss, densw, Sr, alphas, alph 
       
!      ------  DFs correction.

dens       = const   (2) 
denss      = const  (17)
densw      = const  (18)
Sr         = const  (19)
alphas     = const (14)
alph       = alphas
ic_hrelSat = const (34)

if (ic_hrelSat.eq.1) then                               !  allocate if necessary relative h sat
   if (.NOT.allocated (hrelSat_DF) )  then
      allocate (hrelSat_DF(npoin))                         ! CH MP bug 2nd April 19
   endif
   hrelSat_DF = alphas
endif

IF (npoin_s.eq.npoin_w) THEN                             !  we have given same nr of points
   IF( abs(rho(1)-rho(npoin_s +1)).LT.1.e-6) THEN        !  with same values hs = hw -> is h_DF
      do ipoin  = 1, npoin_s                             !  loop about soil poins...water is the same
         ipoinw = ipoin + npoin_s                        !  for ech s ipoin, w is ipoinw
         if (ic_hrelSat.eq.1) alph = hrelSat_DF (ipoin)  !  when desaturation, alpha depends on poin
         
         n_DF (ipoin)  = (denss-dens) / (denss - ((1.-alph)*Sr + alph)*densw)
         n_DF (ipoinw) = n_DF (ipoin)

         h_DF (ipoin) =  rho(ipoin)                      !   the values we have read for hs and hw are h_DFs
         h_DF (ipoinw) = h_DF(ipoin)

         rho  (ipoin)  = rho  (ipoin)  * (1.- n_DF(ipoin))  ! h for soil
         mass (ipoin)  = mass (ipoin)  * (1.- n_DF(ipoin))  ! and mass
         rho  (ipoinw) = rho  (ipoinw) * n_DF(ipoin) * ((1. - alph)*Sr + alph)
         mass (ipoinw) = mass (ipoinw) * n_DF(ipoin) * ((1. - alph)*Sr + alph)
      enddo
   ELSE     
      do ipoin  = 1, npoin_s
         ipoinw = ipoin + npoin_s
         if (ic_hrelSat.eq.1) alph = hrelSat_DF (ipoin)
         n_DF (ipoin)  = (denss-dens) / (denss - ((1.-alph)*Sr + alph)*densw)
         n_DF (ipoinw) = n_DF (ipoin)
         h_DF (ipoin) = (rho(ipoin) + rho(ipoinw)) / ( 1.-n_DF(ipoin)*(1.-Sr)*(1.-alph) )
         h_DF (ipoinw) = h_DF(ipoin)
      enddo
   ENDIF
!   call get_s_at_w_DF_SW (0)    !  CH MP 5th April 19 arg means do not normalize hs close to walls
   call get_s_at_w_DF_SW (0)    
!  call Get_s_at_w_DF_SW_BAK  !  ***** check


ELSE      !   npoin_s .NE. npoin_w 

      call get_s_at_w_DF_SW     !  CH MP 5th April 19 arg means do not normalize hs close to walls
!      DO ipoin = 1, npoin      !  CH MP 14th april n_DF and h_df computed by get_s_at_w
!         if (ic_hrelSat.eq.1) alph = hrelSat_DF (ipoin)
!         n_DF (ipoin) =  (denss-dens) / (denss - ((1.-alph)*Sr + alph)*densw)   
!         h_DF (ipoin) = ( rho(ipoin) + h_sw(ipoin) ) / &
!                        ( 1.-n_DF(ipoin)*(1.-Sr)*(1.-alph) )
!      ENDDO
ENDIF


END Subroutine DFs_initial_correction



!-------------------------------------------------------------------

       Subroutine Initialize_Pwp
       
!-------------------------------------------------------------------

!      ------   

implicit none

integer(ink) BCt1, BCt2, iaux, ipoin

real   (irk) BCv1,   BCv2, pwp0_type, xi, hh
real   (irk) alphas, alfs, xirel

character(60) text 
1005 format (a60)

alphas = const (14)    
nUaux = 0

if (icpwp.eq.1) then                             ! Pwp at bottom cos dissipation Use u(npoin) 

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) pwp0rel
   write(chk_file,*) pwp0rel
   
   Do ipoin = 1, npoin
      alfs  = alphas
      if (ic_hrelSat.eq.1)       alfs = hrelSat_DF(ipoin)
      hh    = rho(ipoin) * alfs
      if (ic_ws_interact.EQ.2)   hh = h_DF(ipoin) * alfs
      u(ipoin) = pwp0rel * hh
   enddo
   Bc_pwp_type(1) = 2    !  MP Sept 2021 initialize for grids, avoid all zeros
   Bc_pwp_type(2) = 1
   Bc_pwp_val (1) = 0.
   Bc_pwp_val (2) = 0.
   
elseif (icpwp.eq.2) then    ! Pwp via FD mesh using Uaux(nuaux,npoin)

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)   nUaux
   write(chk_file,*)   nUaux
   if (.NOT.allocated (Uaux) ) then
      allocate ( Uaux  (nUaux,npoin) )
      allocate ( Uaux0 (nUaux,npoin) )
      allocate ( dUaux (nUaux,npoin) )
   endif
   Uaux = 0.0; Uaux0=0.0; dUaux = 0.0
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    BcT1,BcV1,BcT2,BcV2    ! Bc type at bottom and top
   write(chk_file,*)    BcT1,BcV1,BcT2,BcV2    !  type 0 ct 1 lin 2 cos 
   Bc_pwp_type(1) = BcT1
   Bc_pwp_type(2) = BcT2
   Bc_pwp_val (1) = BcV1
   Bc_pwp_val (2) = BcV2
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) pwp0_type, pwp0rel     ! Value is pwprel*h sat  
   write(chk_file,*) pwp0_type, pwp0rel
   if (pwp0_type.eq.4) then        !   new MP July 2021 
         read (dat_file,1005) text
         write(chk_file,1005) text
         read (dat_file,*)    xirel     ! Value is pwprel*h sat  
         write(chk_file,*)    xirel
   endif
   do ipoin = 1,npoin                       !  this block CH MP 22 Feb 19
      alfs  = alphas
      if (ic_hrelSat.eq.1)       alfs = hrelSat_DF(ipoin)
      hh    = rho(ipoin)                     !! * alfs
      if (ic_ws_interact.EQ.2)   hh = h_DF(ipoin) !!  * alfs
      if (pwp0_type.eq.0) then
         Uaux(:,ipoin) = pwp0rel*hh  ! pwp0rel*rho(1:npoin)  *** changed to pwp=pwp*rho*g*h
      elseif (pwp0_type.eq.1) then
         do iaux = 1, nUaux
            xi = 1.*(iaux-1)/(nUaux-1)
            Uaux(iaux,ipoin) = (1-xi)*pwp0rel*hh                
         enddo
      elseif (pwp0_type.eq.2) then
         do iaux = 1, nUaux
            xi = (iaux-1.)/(nUaux-1.)   ! Patetico
            pi = 4.*atan(1.0)
            Uaux(iaux,ipoin) = cos(pi*xi/2)*pwp0rel*hh       
         enddo
      elseif (pwp0_type.eq.3) then   !   new MP July 2021 
            Uaux(: ,ipoin) = 0.0
            Uaux(1 ,ipoin) = pwp0rel*hh
      elseif (pwp0_type.eq.4) then   !   new MP July 2021 
            Uaux(: ,ipoin) = 0.0
            Do iaux = 1, nuaux
               xi = (iaux-1.)/(nUaux-1.)
               if (xi.lt.xirel) then
                  Uaux(iaux ,ipoin) = hh*pwp0rel*(1-xi) 
               endif    
            enddo
      endif
   enddo
   
   Uaux0 = Uaux                  ! MP 11 May 2020
   dUaux = 0.0  
   u(1:npoin) = Uaux(1, 1:npoin)  ! for DFs, do it again after correcting pwps

   !      ------ dealing with stefan erosion + pwp (ic_erosion=1 nUaux>1)

   safetpwp = 1.0; ic_Activ_Pwp_Stef_eros= 0        !   1st may 2021 
   if (ic_erosion.eq.1.AND.ic_Stef_eros.Eq.1) then  !   CHGD 31st May 2021
       read (dat_file,1005) text
       write(chk_file,1005) text
       read (dat_file,*)    BC_Adv_Type(1), BC_Adv_Val(1)     ! Bc type at bottom  and Value
       write(chk_file,*)    BC_Adv_Type(1), BC_Adv_Val(1)     !  
       read (dat_file,1005) text
       write(chk_file,1005) text
       read (dat_file,*)    safetpwp, ic_Activ_Pwp_Stef_eros  ! safety fact * dt and 
       write(chk_file,*)    safetpwp, ic_Activ_Pwp_Stef_eros  !  activate 1 or not 0  Stefan ALE  
       if (ic_Activ_Pwp_Stef_eros.EQ.1) then
          if (.NOT.allocated(Activ_Pwp_Stef_Eros_Pt)) then
              allocate (Activ_Pwp_Stef_Eros_Pt(npoin))
          endif
          Activ_Pwp_Stef_Eros_Pt = 1
       endif
   endif
    
elseif (icpwp.eq.3) then    ! Pwp at bottom cos dissipation Uaux with nUaux=1 

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) pwp0rel
   write(chk_file,*) pwp0rel
   nUaux = 1
   if (.NOT.allocated (Uaux) ) then
      allocate ( Uaux  (nUaux,npoin) )
      allocate ( Uaux0 (nUaux,npoin) )
      allocate ( dUaux (nUaux,npoin) )
   endif
   Uaux = 0.0; Uaux0=0.0; dUaux = 0.0
   do ipoin = 1, npoin
      alfs  = alphas
      if (ic_hrelSat.eq.1)       alfs = hrelSat_DF(ipoin)
      hh    = rho(ipoin) * alfs
      if (ic_ws_interact.EQ.2)   hh = h_DF(ipoin) * alfs
      Uaux(1,ipoin) = pwp0rel*hh
   enddo

   Uaux0 = Uaux                   ! MP 11 May 2020
   dUaux = 0.0     
   u(1:npoin) = Uaux(1, 1:npoin)  ! for DFs, do it again after correcting pwps            

endif         

if ((ic_basal_flux.EQ.1).AND.(.NOT.allocated(basal_flux))) then   !  MP Aug 2021 CHGD
     allocate (basal_flux(npoin))
endif

END Subroutine Initialize_Pwp

END SUBROUTINE Input_SW



!-------------------------------------------------------------------

       Subroutine Input_SW_TRIGGER
       
!-------------------------------------------------------------------

!      ------  Reads basic control data for Depth Integrated problems

implicit none
logical conected

integer(ink) len1,  nline, icunkno   !  , nfric made global SW var 7th dec 2017 MP 
integer(ink) iline, i , j, itotv, ipoin, idimn, iiabs, iivn0, ic_abss
integer(ink) ipoin_w, nGidMaskSw

integer(ink) BCt1,BCt2          ! type of BCs at bottom and top 1 Dirich 2 N
real   (irk) BCv1,BCv2          ! values of BCs 
integer(ink) pwp0_type          ! type of law for init pwp(z) 0 ct 1 lin 2 cos
integer(ink) iaux               !         value at bottom is hrelpw*h
real   (irk) xi, hh             ! z/h xi=0 at bottom 1 at top of st layer
real   (irk) hrel
real   (irk) dens, denss, densw, nporos, hhs, hhw, porostmp
                                     
integer (ink) icase_mat, index_mat  !  dummmy indexes for multiple case computation
real    (irk) Cval_mat , Cval0_mat  !  see after readind const line 414 approx      

real(irk), allocatable:: coord_cg(:)     ! center of mass coords for expansion of mass
real(irk), allocatable:: shift_case(:)   ! shift dx dy for cases analysis  MP 12th May 2020

character(72) text
!              Assign units to files, ask for a name and open them

dat_file     = 10
chk_file     = 21
gid_msh      = 41
gid_res      = 42
gis_csv1     = 43          ! Saeidmt 30March2022   Adding .CSV for GiS
gis_csv2     = 44          ! Saeidmt 30March2022   Adding .CSV for GiS 

1005 format (a60)
1006 format (a16)

len1 = len_trim(SPH_SW_problem_name)

open( dat_file,  file=SPH_SW_problem_name(1:len1)//'.dat'     )
open( chk_file,  file=SPH_SW_problem_name(1:len1)//'.chk'     )

open(  gid_msh,  file=SPH_SW_problem_name(1:len1)//'.post.msh')
if (ic_Cases_mat.GT.0) then              ! MP May 2020 Cases. We will append at the end of the file each final situation for the new case
   INQUIRE ( gid_res, OPENED= conected)
   if (.NOT.conected) then
      open( gid_res,  file=SPH_SW_problem_name(1:len1)//'.post.res')
      write (gid_res,*) 'GiD Post Results File 1.0'
   endif
else
   open( gid_res,  file=SPH_SW_problem_name(1:len1)//'.post.res')
   write (gid_res,*) 'GiD Post Results File 1.0'
endif

if (ic_TRIGGER.EQ.1) then                                              ! Saeidmt 30March2022   Adding .CSV for GiS
    INQUIRE ( gis_csv1, OPENED= conected)
    if (.NOT.conected) then
        open(  gis_csv1,  file=SPH_SW_problem_name(1:len1)//'.1.csv')
        write (gis_csv1,*) 'Node No.,X,Y,Z,PoF,FoS mean,FoS sigma,FoS variab,Reliab index,PGA'
    end if
    INQUIRE ( gis_csv2, OPENED= conected)
    if (.NOT.conected) then
        open(  gis_csv2,  file=SPH_SW_problem_name(1:len1)//'.2.csv')
        write (gis_csv2,*) 'Slope,PoF,FoS mean,FoS sigma,FoS variab,Reliab index,PGA'
    end if
end if
 
!      ------ Initialize some vars

pi        = 4*atan(1.0)
npoin_w   = 0
npoin_s   = 0
ic_WIR    = 1

nGidMaskSw = 12
if (.NOT.allocated (Gid_Mask_SW) ) then
   allocate (Gid_Mask_SW(nGidMaskSw) ) 
   Gid_Mask_SW = 0.0
endif
!      ------  read lines of text describing the problem 

read (dat_file,*)  nline
write(chk_file,*)  nline    
do iline = 1,nline
   read (dat_file,1005) text
   write(chk_file,1005) text
enddo

ic_SWAlg = 0   
npoin_in_model = 0   
 
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) ndimn   
write(chk_file,*) ndimn  

npoin  =  1 
npoin_s = 1; ntotal = 1; 

! ---------------  global arrrays routine ------
  
    if (.NOT.allocated (x00) ) then
      allocate (     x00(ndimn,ntotal),    x0 (ndimn,ntotal),     x (ndimn,ntotal),   vx0 (ndimn,ntotal),   vx (ndimn,ntotal) )
      allocate (     Z00(ntotal), rho00(ntotal) )
      allocate (     dx (ndimn,ntotal),   dvx (ndimn,ntotal),    av (ndimn,ntotal) )
      allocate ( indvxdt(ndimn,ntotal),exdvxdt(ndimn,ntotal),ardvxdt(ndimn,ntotal) )
      allocate ( mass   (ntotal), rho(ntotal), p(ntotal), u(ntotal), hsml(ntotal), c(ntotal) )
      allocate ( mass0  (ntotal) )
      allocate ( u0     (ntotal), rho0(ntotal),  du(ntotal),drho(ntotal) )
      allocate ( itype(ntotal), countiac(ntotal) )
      allocate (If_Out_Domain(ntotal))
      allocate (Xmin_Domain(ndimn))
      allocate (Xmax_Domain(ndimn))
    endif 

        
!      ------  Initialize (changed, because we have to initialize x,vx,itype,rho,mass,hsml,u and p)   

        x0       = x  ;  x00      =  x0;  vx0     = vx ;  Z00      = 0.0 
        u0       = u  ;  rho0     = rho;  rho00   = rho;   
        
        dx       = 0.0;  dvx      = 0.0;  du      = 0.0;  drho     = 0.0                
        indvxdt  = 0.0;  exdvxdt  = 0.0;  ardvxdt = 0.0;  av       = 0.0
        c        = 0.0
        mass0 = mass   

IF (ic_cases_mat.GE.2) THEN       !   change mat props according to cases MP 24 03 2022
   Do icase_mat = 1, ncases_mat
      index_mat = Ccases_mat ( icase_mat, 1)
      if (index_mat.eq.501) then   !   Cval_mat will multiply rho and mass ..ctrl of diff hs
          Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
          Do ipoin = 1, npoin
             rho (ipoin) = rho (ipoin)*Cval_mat
             mass(ipoin) = mass(ipoin)*Cval_mat
          enddo   
          rho0 = rho ; rho00 = rho; mass0 = mass
      endif  
      if (index_mat.eq.511) then
          Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)  ! expansion fact
          coord_cg = 0.0
          do ipoin = 1, npoin
             coord_cg(:) = coord_cg(:) + x(:,ipoin)
          enddo
          coord_cg = coord_cg/npoin
          do ipoin = 1, npoin
             x(:,ipoin) = coord_cg(:) + Cval_mat* (x(:,ipoin)-coord_cg(:))
             mass(ipoin) = mass(ipoin)*Cval_mat    
          enddo  
          x0 = x; x00 = x ; mass0 = mass
      endif 

      if (index_mat.eq.521) then      !      ------  shift of the whole mass case
          shift_case(1) = OC_mat (nproblems_index, 3 +icase_mat)           !  this is DX
          if (ndimn.eq.2) then
             shift_case(2) = OC_mat (nproblems_index, 3 +ncases_mat + 1)   !  this is DY
          endif
          do ipoin = 1, npoin
            x(:,ipoin) = x(:,ipoin) + shift_case
          enddo  
          x0 = x; x00 = x ; mass0 = mass
      endif      
            
   enddo  !  -------------  end loop icases         
ENDIF                  !  -------------- end cases ic_case_mat

! pa_sph, nnps, sle, skf = defaults?
!               
nconst = 80               ! before it was 15 and 30  
!
if (.NOT.allocated(const)) then 
        allocate (const(nconst))
        const = 0.0
endif

read (dat_file,1005) text   !  N_slopes Slope_max_cut Slope_min_cut
write(chk_file,1005) text
read (dat_file,*)  const(71), const(72), const(73)
write(chk_file,*)  const(71), const(72), const(73) 

read (dat_file,1005) text   ! n_stats  
write(chk_file,1005) text
read (dat_file,*)  const(74) 
write(chk_file,*)  const(74) 

read (dat_file,1005) text   !  h  hrelSat   Sr  
write(chk_file,1005) text
read (dat_file,*)  const(75), const(14), const(19)
write(chk_file,*)  const(75), const(14), const(19)

read (dat_file,1005) text   ! ic_cosTheta, g, PGA  
write(chk_file,1005) text
read (dat_file,*)  const(16), const(1), const(76)
write(chk_file,*)  const(16), const(1), const(76) 

read (dat_file,1005) text   !  dens, denss densw ic_fct_dens 
write(chk_file,1005) text
read (dat_file,*)  const(2), const(17), const(18), const(20)
write(chk_file,*)  const(2), const(17), const(18), const(20)

read (dat_file,1005) text   !  nfric cohesion tanFi 
write(chk_file,1005) text
read (dat_file,*)  const(5), const(6), const(9)
write(chk_file,*)  const(5), const(6), const(9) 

nUaux = 0                   ! -------- pwp info 
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) icpwp
write(chk_file,*) icpwp 

if (icpwp.ne.0) then 
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) pwp0rel
   write(chk_file,*) pwp0rel
   const (77) = pwp0rel
endif 

IF (ic_cases_mat.GE.2) THEN       !   change mat props according to cases MP 24 03 2022
   Do icase_mat = 1, ncases_mat
      index_mat = Ccases_mat ( icase_mat, 1)
      Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
      if (index_mat.LE.nconst) then               ! prepare options 501 (h fact), etc
         IF (index_mat.eq.2) THEN                 !   change porosity to density of the mixture
            if (Cval_mat.LE.1.01) then            !   we give poros in table. Need to get dens
               if (const(17).ge.1.e-4) then
                  const(2) = (1.-Cval_mat)*const(17) +Cval_mat*Const(18)
               else
                  write(*,*) ' *** provide dens solid and water. Type any number '
                  read  (*,*) Cval0_mat           ! just to input a keystroke
                  stop
               endif
            else                         ! we are providing density in table
               Const (index_mat) = Cval_mat      
            endif
         ELSE                             ! general material property
            Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
            Const (index_mat) = Cval_mat
         ENDIF
      endif              ! end prepare options 501 (h fact), etc                 
   enddo          
ENDIF   

!      ------   Obtain nfric and
!                      if_frictional  (= 1 for frictional laws: 5 to 10, 26 and 27
!                      if_CV2byR      ( =1  for 7,8,10, 26 and 27)
    
nfric = const(5) + 0.01
if_frictional = 0
if_V2byR      = 0
if ( ((nfric.ge.5).AND.(nfric.le.10)).OR.(nfric.eq.26).OR.(nfric.eq.27) ) then
    if_frictional = 1
endif 

IF (icpwp.ne.0) then                  !  keep this old version using code 401 

   Do ipoin = 1, npoin
      hrel  = const(14)
      hh    = rho(ipoin) * hrel
      u(ipoin) = pwp0rel * hh
   enddo
   
   IF (ic_cases_mat.GE.2) THEN       !   change mat props according to cases MP 24 03 2022
      Do icase_mat = 1, ncases_mat
         index_mat = Ccases_mat ( icase_mat, 1)
         if (index_mat.eq.401) then
            Cval_mat  = OC_mat (nproblems_index, 3 +icase_mat)
            Do ipoin = 1, npoin
               u (ipoin) = u(ipoin) * Cval_mat
            enddo     
         endif      
      enddo          
   ENDIF 
ENDIF     !   pwp     

END Subroutine Input_SW_TRIGGER



! ----------------------------------------------------------------------   

      subroutine setup_coarse_plt

! ----------------------------------------------------------------------  

!   Used when ic_coarse_plt = 1 to build a new coarse mesh for plotting.

implicit none  
   
integer(ink) ipoi_plt, ipoix_plt, ipoiy_plt, ipoigx, ipoigy, ipoig, ic_imd, ic_tmp 
real   (irk) xg, yg, Zp, zming

xmin_plt   = xming;  xmax_plt = xmaxg
ymin_plt   = yming;  ymax_plt = ymaxg

deltx_plt  = deltxg*fc_coarse_plt
delty_plt  = deltyg*fc_coarse_plt
npoix_plt  = (xmaxg-xming)/deltx_plt + 1
npoiy_plt  = (ymaxg-yming)/delty_plt + 1
npoi_plt   = npoix_plt * npoiy_plt
deltx_plt  = (xmaxg-xming)/(npoix_plt-1)
delty_plt  = (ymaxg-yming)/(npoiy_plt-1)

if (.NOT.allocated (coor_plt) ) then
   allocate ( coor_plt(3,npoi_plt), poi_in_model_plt(npoi_plt) )  !  x,y,Z
endif
coor_plt = 0.0  ; poi_in_model_plt = 0
   
!  Get coordinates of structured mesh npoix_plt x npoiy_plt (divs)

zming = 1.e6
do ipoig = 1, npoig
   if (topol(1,ipoig).lt.zming) zming = topol(1,ipoig)
enddo

         
  ipoi_plt  = 0
  DO ipoiy_plt = 1,npoiy_plt
      ic_tmp = 0 
      yg     = yming + (ipoiy_plt-1)*delty_plt
      ipoigy = (yg-yming)/deltyg  + 1
      if (ipoigy.eq.npoigy) THEN       ! (a4) bug fixed
          ipoigy = npoigy-1
          yg     = 0.99*yg
      endif
      
      Do ipoix_plt = 1,npoix_plt 
         ipoi_plt  = ipoi_plt + 1
         xg     = xming + (ipoix_plt-1)*deltx_plt  
         coor_plt (1,ipoi_plt) = xg
         coor_plt (2,ipoi_plt) = yg
         coor_plt (3,ipoi_plt) = zming                             
         ipoigx = (xg-xming)/deltxg  + 1
         if (ipoigx.eq.npoigx) THEN       ! (a4) bug fixed
             ipoigx = npoigx -1 
             xg     = 0.99*xg
         endif
         ipoig   = (ipoigy-1)*npoigx + ipoigx
           if (ic_tmp.ne.1) then
            ic_imd  = poig_in_model(ipoig)*poig_in_model(ipoig+1)*poig_in_model(ipoig+npoigx+1)*poig_in_model(ipoig+npoigx)
            if (ic_imd.eq.1) then
                Call Get_Z_Topo (ipoi_plt,xg, yg, Zp )
                coor_plt      (3,ipoi_plt) = Zp
                poi_in_model_plt(ipoi_plt) = 1
            endif
         endif
      enddo
  ENDDO



END SUBROUTINE setup_coarse_plt

! ----------------------------------------------------------------------   

      subroutine Get_water_points

! ----------------------------------------------------------------------   

implicit none  

integer(ink) icunk_w, icc
character(72) text
!      ------   

if (npoin_s.gt.0) then
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    icunk_w    
   write(chk_file,*)    icunk_w  
else   
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    icunk_w, h_inf_SW   
   write(chk_file,*)    icunk_w, h_inf_SW
endif   

if     (icunk_w.eq.10.or.icunk_w.eq.11.or.icunk_w.eq.12)  then 
    call Read_WP_Topol(icunk_w)
elseif (icunk_w.eq.1) then 
    call Read_Cylinder_SW_w
elseif (icunk_w.eq.2) then 
    call Read_Multilinear_X_Water_Pts
elseif (icunk_w.eq.3) then 
    call Read_1Dbk_sw_W    
elseif (icunk_w.eq.7) then              ! --- MP June 2023 for modif pts generated
    call Read_Source_plus_New_WATER ( icunk_w, xshiftg, yshiftg)    
elseif (icunk_w.eq.8) then              ! --- MP June 2023 for modif pts generated
    call Read_Source_plus_New_WATER ( icunk_w, xshiftg, yshiftg)
    ic_geomm_pts = 1                    !  MP June 2023                             !                          
elseif (icunk_w.eq.0) then
    write(*,*) ' no water points in WIR code...want to continue? (0,1) '
    read (*,*) icc
    if (icc.eq.0) STOP                           
elseif (icunk_w.eq.200) then      !  ****  RESTART 9th aprilh 2106 
    call Get_Restart_Info_SW_w    !        new  
else
    write(*,*) ' incorrect icunk_w option ', icunk_w
endif     

1005     format(a60) 

END SUBROUTINE  Get_water_points

! ----------------------------------------------------------------------   

      subroutine Read_Cylinder_SW_w  

! ----------------------------------------------------------------------   

!      Reads cylider or ellipsoid circular base

implicit none  

character(72) text 
integer(ink) icdens, nlinep, ipoin, ilinep, ith   
integer(ink) icunk_w, icc, itp_ps 
   
real   (irk) xp0, yp0, rp0, hp0, cgra05, mss, drr, delth, rpp, hfact, theta
real   (irk) rp1, hp1, factp1     
real   (irk) rx ,  ry
real   (irk) facthsml, twopi 
 
!      ------   

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    xp0, yp0, rp0, hp0, icdens       
write(chk_file,*)    xp0, yp0, rp0, hp0, icdens 

factp1 = 1.0
if (icdens.eq.2) then      !      --- Inner cyl with a diff heigth...read radius and h
    read (dat_file,1005) text         
    write(chk_file,1005) text        
    read (dat_file,*)    rp1, hp1      
    write(chk_file,*)    rp1, hp1
    factp1 = hp1/hp0
endif

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    nlinep, facthsml       
write(chk_file,*)    nlinep, facthsml  

write(chk_file,*) ' --- we are reading from cylider/hemisphere point generation -----'
write(chk_file,*) ' give factor for smoothing length. Total is 2*facthsml'
write(chk_file,*)  facthsml

! cgra05  = const(1)/2.0
                                        ! nlinep  = (-3 + sqrt (9. + 12.*real(npoin-1)))/6
npoin_w = 1 + 3*nlinep*(nlinep+1)
mss     = pi*rp0*rp0*hp0/npoin_w
if (nlinep.ge.1) then  ! changed MP 28.12.2018 use 1 point only
    drr = (rp0/nlinep)  
else
    drr = rp0
endif

if (.NOT.allocated (x_w) ) then
   allocate ( x_w(ndimn,npoin_w),  rho_w(npoin_w), mass_w(npoin_w), hsml_w (npoin_w) )
   allocate (vx_w(ndimn,npoin_w),  u_w  (npoin_w), itype_w(npoin_w) )
endif                                       ! allocate p_w(npoin_w) deleted
ipoin        = 1
x_w (1,1)    = xp0
x_w (2,1)    = yp0
vx_w(:,1)    = 0.0
itype_w (1)  = 2
rho_w   (1)  = hp0*factp1
mass_w  (1)  = mss*factp1
hsml_w  (1)  = 2.*drr
u_w     (1)  = 0.0       

DO ilinep = 1,nlinep
   delth  = 2.*pi/(6.*ilinep)
   rpp    = (ilinep)*drr
   if ( (rpp/rp0).gt.1.00 ) rpp=0.9999999*rp0 ! avoids -ve root below
   factp1 = 1.0
   if (icdens.eq.2.and.rpp.le.rp1) factp1 = hp1/hp0
    
   hfact  = sqrt(1.0-(rpp/rp0)*(rpp/rp0))
   if (hfact.lt.0.1) hfact=0.1
   DO ith = 1, 6*ilinep
      ipoin = ipoin + 1
      if (ipoin.gt.npoin_w) exit !       ------  take care of roundoffs
      theta = ith*delth
      rx    = rpp*cos(theta)
      ry    = rpp*sin(theta)
      x_w   (1,ipoin)   = xp0 + rx
      x_w   (2,ipoin)   = yp0 + ry
      vx_w  (:,ipoin)   = 0.0
      itype_w (ipoin)   = 2
      if (icdens.eq.0) then
         rho_w  (ipoin)  = hp0   
         mass_w (ipoin)  = mss
      elseif (icdens.eq.1) then
         rho_w (ipoin)   = hp0*hfact  
         mass_w(ipoin)   = mss*hfact
      elseif (icdens.eq.2) then
         rho_w (ipoin)   = hp0*factp1  
         mass_w(ipoin)   = mss*factp1
      endif
     hsml_w  (ipoin) = 2.*drr*facthsml
     u_w     (ipoin) = 0.0
   ENDDO
ENDDO   

1005     format(a60)  

END SUBROUTINE  Read_Cylinder_SW_w
!------------------------------------------
 
   subroutine Get_Restart_Info_SW_w
 
!------------------------------------------

!      ------  Reads restart info; icunk_s option is 20

implicit none
 
integer(ink)  In_Source, Chk_Source

integer(ink)  len1,  i , ipoin,id , nliness

! time_sph_restart  is time at which computation is restarted 
! (driver time vars MODULE) 

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string
character( 2) textS

In_Source    = 86
Chk_Source   = 88

   1 format (a60)
1005 format (a60) 
1006 format (a16) 
1002 format (a2 )
  99 format  (i7,2x,i5,8(2x,g14.7))

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1006) problem_name 
  
len1 = len_trim(problem_name)
 
open  (  In_Source ,  file = 'restart.dat'   )
open  (  Chk_Source,  file = 'restart.chk'   )

if (npoin_s.gt.0) then       !  read but not use solid pts
    nliness = 3 + npoin_s
    do i = 1, nliness
       read  (In_Source, 1002) textS
       read  (Chk_Source,1002) textS
    enddo
endif    
       
read ( In_Source , 1005) text
write( Chk_Source, 1005) text
read ( In_Source , *)    npoin_w, time_sph_restart     
write( Chk_Source, *)    npoin_w, time_sph_restart   

time = time_sph_restart

!      ------ all Fi_s vars have been defined as allocatable in SW vars module

if (.NOT.allocated (x_w) )  then
   allocate ( x_w(ndimn,npoin_w),  rho_w(npoin_w), mass_w (npoin_w), hsml_w (npoin_w) )
   allocate (vx_w(ndimn,npoin_w),  u_w  (npoin_w), itype_w(npoin_w) )
endif

x_w = 0.0; vx_w=0.0; rho_w=0.0; mass_w=0.0;hsml_w=0.0; u_w = 0.0

read ( In_Source , 1005) text
write( Chk_Source, 1005) text

DO i = 1, npoin_w 

   read (In_Source, *) ipoin, itype_w(i), (x_w(id,i),id=1,ndimn), (vx_w(id,i),id=1,ndimn), &
                               rho_w (i), mass_w(i), hsml_w(i), u_w(i)
   write(Chk_Source,99)ipoin, itype_w(i), (x_w(id,i),id=1,ndimn), (vx_w(id,i),id=1,ndimn), &
                               rho_w (i), mass_w(i), hsml_w(i), u_w(i)
ENDDO

close ( In_Source)
close (Chk_Source)

END subroutine Get_Restart_Info_SW_w


! ----------------------------------------------------------------------   

      subroutine Read_1Dbk_SW_w
      
! ----------------------------------------------------------------------   

!     This subroutine is used to generate initial data for the 
!     1 d break dam problem (for DF option) Same than solid
!      
!     Reads:                Xl,Xc,Xr   (X0 is at the dam)
!                           hl and Hr  heights of water left and right of X0
!                           npoin      number of particles to be generated
! 
!      

implicit none  

integer(ink)  nl, nr, i , iomega  
real   (irk)  Xl,   Xc,   Xr,  hl,  hr, mss, pleft, pright
real   (irk)  xll, xrr, htot,  fl,  fr, dxl,   dxr, dxlr

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_w   
write(chk_file,*)    npoin_w 

if (.NOT.allocated (x_w) ) then 
   allocate ( x_w(ndimn,npoin_w),  rho_w(npoin_w), mass_w(npoin_w), hsml_w (npoin_w) )
   allocate (vx_w(ndimn,npoin_w),  u_w  (npoin_w), itype_w(npoin_w) )
endif

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    Xl, Xc, Xr, hl, hr, iomega   
write(chk_file,*)    Xl, Xc, Xr, hl, hr, iomega 

!      ------       Initialize

xll   = (Xc-Xl)
xrr   = (Xr-Xc)
htot  = xll*hl + xrr*hr
fl    = xll*hl/htot
fr    = 1.-fl
nl    = npoin_w * fl
nr    = npoin_w - nl
dxl   = xll/(nl-1)
dxr   = xrr/nr

!      ------       Generate real particles 

If (iomega.eq.0) then

   x_w     = 0.
   vx_w    = 0.

   mss  = hl*xll/nl
   do i = 1,nl
      mass_w (i) = mss
      hsml_w (i) = dxl*2.4   ! ****2
      itype_w(i) = 1
      x_w  (1,i) = xl + dxl*(i-1)
   enddo 

   mss  = hr*xrr/nr
   do i = 1,nr
      mass_w (nl+i) = mss
      hsml_w (nl+i) = dxr*2.5   ! ***** changed 7 June 2013 from 2.0
      itype_w(nl+i) = 1
      x_w  (1,nl+i) = xc + dxr*(i)
   enddo 

!   pleft = (const(1)/2)*hl*hl 
!   pright= (const(1)/2)*hr*hr 
   do i = 1,npoin_w
      if (x_w(1,i).le.(xc + 1.e-6)) then
          u_w  (i)  = 0.0
          rho_w(i)  = hl
 !         p_w  (i)  = pleft
      else  
          u_w  (i)  = 0.0
          rho_w(i)  = hr
 !         p_w  (i)  = pright
      endif        
   enddo   
   
elseif (iomega.eq.1) then

   if (hr.le.1.e-6) then
       dxlr = (xc-xl)/(npoin_w - 1)
   else
       dxlr = (xr-xl)/(npoin_w - 1)
   endif
   do i = 1,npoin_w
      itype_w  (i) = 1
      hsml_w   (i) = 2.*dxlr
      x_w    (1,i) = xl + (i-1)*dxlr
      if (x_w(1,i).le.xc) then
          rho_w(i) = hl
      else
          rho_w(i) = hr
      endif
      mass_w(i)  = dxlr*rho_w(i)
      u_w   (i)  = 0.0
   enddo

else
   write(*,*) ' iomega is incorrect only 0 or 1 ', iomega
   pause
endif 

vx_w = 0.0

1005     format(a60) 

END subroutine Read_1Dbk_SW_w


! ----------------------------------------------------------------------   

        subroutine Read_Multilinear_X_Water_Pts

! ---------------------------------------------------------------------- 

implicit none

integer(ink) nptsh , iptsh, ipoin_w, i
real   (irk) x0ptsh, x1ptsh, xxptsh, f0ptsh, f1ptsh, hhptsh, deltptsh
real   (irk) cgra05
real   (irk) facthsml
real   (irk), allocatable:: xptsh(:), hptsh(:)
       
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_w, facthsml ,nptsh     
write(chk_file,*)    npoin_w, facthsml ,nptsh 

if (.NOT.allocated (x_w) )  then
   allocate ( x_w(ndimn,npoin_w),  rho_w(npoin_w), mass_w(npoin_w), hsml_w (npoin_w) )
   allocate (vx_w(ndimn,npoin_w),  u_w  (npoin_w), itype_w(npoin_w) ) 
endif 
if (.NOT.allocated (xptsh) )  then
   allocate ( xptsh(nptsh), hptsh(nptsh) )  
endif                                                                    ! allocate p_w(npoin_w) deleted


read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*) (xptsh(i),i=1,nptsh)
write(chk_file,*) (xptsh(i),i=1,nptsh)
read (dat_file,*) (hptsh(i),i=1,nptsh)
write(chk_file,*) (hptsh(i),i=1,nptsh)
cgra05       = 9.8/2.0       ! ------ we do it here because properties have not been read!!

if (npoin_w.gt.1) then
   deltptsh     = (xptsh(nptsh)-xptsh(1))/(npoin_w - 1)!!
   do ipoin_w     = 1, npoin_w
      xxptsh    = xptsh(1) + (ipoin_w-1)*deltptsh
      do iptsh  = 1,nptsh-1
         x0ptsh = xptsh(iptsh)
         x1ptsh = xptsh(iptsh+1)
         if (xxptsh.ge.x0ptsh.and.xxptsh.le.x1ptsh) then
             f0ptsh = (xxptsh-x1ptsh)/(x0ptsh-x1ptsh)
             f1ptsh = 1.-f0ptsh
             hhptsh = f0ptsh*hptsh(iptsh) + f1ptsh*hptsh(iptsh+1)
             if (hhptsh.le.1.e-3) hhptsh = 1.e-3
         endif
      enddo
      x_w  (1,ipoin_w) = xxptsh
      vx_w (1,ipoin_w) = 0.0
      itype_w(ipoin_w) = 2
      rho_w  (ipoin_w) = hhptsh
      mass_w (ipoin_w) = hhptsh*deltptsh
      if ((ipoin_w.eq.1).or.(ipoin_w.eq.npoin_w)) mass_w(ipoin_w)=mass_w(ipoin_w)/2.0
      hsml_w (ipoin_w) = 2.*deltptsh*facthsml
      u_w    (ipoin_w) = 0.0
   enddo
elseif (npoin_w.eq.1) then
      x_w  (1,1) = xptsh(1)
      vx_w (1,1) = 0.0
      itype_w(1) = 2
      rho_w  (1) = hptsh(1)
      mass_w (1) = hptsh(1) 
      hsml_w (1) = 2.*facthsml
      u_w    (1) = 0.0
endif      
   
1005     format(a60) 

deallocate(xptsh)
deallocate(hptsh)

END subroutine Read_Multilinear_X_Water_Pts


! ----------------------------------------------------------------------   

        subroutine Read_Multilinear_X_Water_PtsOLD

! ---------------------------------------------------------------------- 

implicit none

integer(ink) nptsh , iptsh, ipoin_w, i
real   (irk) x0ptsh, x1ptsh, xxptsh, f0ptsh, f1ptsh, hhptsh, deltptsh
real   (irk) cgra05
real   (irk) facthsml
real   (irk), allocatable:: xptsh(:), hptsh(:)
       
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_w, facthsml ,nptsh     
write(chk_file,*)    npoin_w, facthsml ,nptsh 

if (.NOT.allocated (x_w) ) then
   allocate ( x_w(ndimn,npoin_w),  rho_w(npoin_w), mass_w(npoin_w), hsml_w (npoin_w) )
   allocate (vx_w(ndimn,npoin_w),  u_w  (npoin_w), itype_w(npoin_w) ) 
                                                                    ! allocate p_w(npoin_w) deleted
   allocate ( xptsh(nptsh), hptsh(nptsh) ) 
endif

read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*) (xptsh(i),i=1,nptsh)
write(chk_file,*) (xptsh(i),i=1,nptsh)
read (dat_file,*) (hptsh(i),i=1,nptsh)
write(chk_file,*) (hptsh(i),i=1,nptsh)

deltptsh     = (xptsh(nptsh)-xptsh(1))/(npoin_w - 1)
cgra05       = 9.8/2.0       !                  ------ we do it here because properties have not been read!!!!
do ipoin_w     = 1, npoin_w
   xxptsh    = xptsh(1) + (ipoin_w-1)*deltptsh
   do iptsh  = 1,nptsh-1
      x0ptsh = xptsh(iptsh)
      x1ptsh = xptsh(iptsh+1)
      if (xxptsh.ge.x0ptsh.and.xxptsh.le.x1ptsh) then
          f0ptsh = (xxptsh-x1ptsh)/(x0ptsh-x1ptsh)
          f1ptsh = 1.-f0ptsh
          hhptsh = f0ptsh*hptsh(iptsh) + f1ptsh*hptsh(iptsh+1)
          if (hhptsh.le.1.e-3) hhptsh = 1.e-3
      endif
   enddo
   x_w  (1,ipoin_w) = xxptsh
   vx_w (1,ipoin_w) = 0.0
   itype_w(ipoin_w) = 2
   rho_w  (ipoin_w) = hhptsh
   mass_w (ipoin_w) = hhptsh*deltptsh
   if ((ipoin_w.eq.1).or.(ipoin_w.eq.npoin_w)) mass_w(ipoin_w)=mass_w(ipoin_w)/2.0
   hsml_w (ipoin_w) = 2.*deltptsh*facthsml
   u_w    (ipoin_w) = 0.0
enddo

1005     format(a60) 

deallocate(xptsh)
deallocate(hptsh)

END subroutine Read_Multilinear_X_Water_PtsOLD


! ----------------------------------------------------------------------   

      subroutine Read_WP_Topol(icunk_w)

! ----------------------------------------------------------------------   

!     This routine gets water nodes info from topo routine...NEW CALL TO GET_Z_TOPO

implicit none  

integer(ink) ngx, ngy, igx, igy, icunk_w
integer(ink) ipoin_w

real   (irk) Z_swl, dgx, dgy, Rext, facthsml, lgx, lgy, xx, yy, RR 
real   (irk) Zp, dArea

real   (irk), allocatable:: xwtmp(:,:), rhowtmp(:)  ! tmp vars 

! to create water particls which are above soil particles  CHUAN
real   (irk), allocatable:: auxv(:,:)
integer(ink) ipoin, idimn
integer(ink) ipoigx, ipoigy, ipoig
integer(ink) inx0, iny0, inx1,iny1, inx, iny 
integer(ink) ieleg, n1, n2, n3, n4, kpoin
real   (irk) distx2, disty2, dist
real   (irk) xref, yref,xi, eta
real   (irk) sh1, sh2, sh3, sh4,  dZpdx, dZpdy
real   (irk) dsh1dx, dsh2dx, dsh3dx, dsh4dx
real   (irk) dsh1dy, dsh2dy, dsh3dy, dsh4dy
real   (irk) dxginv, dyginv
real   (irk) xwmin, xwmax, ywmin, ywmax
integer(ink) insid

if (.NOT.allocated (auxv) ) allocate (auxv(2,npoig))
auxv = 0.0

if (.not.allocated( hs_plus_Z_topo ) ) then
   allocate (    hs_plus_Z_topo    (npoig))
   hs_plus_Z_topo= 0.0 
endif

! rho_s at grid points

DO ipoin=1,npoin_s

   !if ( if_Out_Domain(ipoin).eq.1) CYCLE
   !if ( itype(ipoin).ne.1) CYCLE
   xx = x_s(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x_s(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
   !inx0      = max(1,ipoigx-1)
   !iny0      = max(1,ipoigy-1)
   !inx1      = min(npoigx,ipoigx+2)
   !iny1      = min(npoigy,ipoigy+2)
   inx0      = max(1,ipoigx)
   iny0      = max(1,ipoigy)
   inx1      = min(npoigx,ipoigx+1)
   iny1      = min(npoigy,ipoigy+1)   
   
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      distx2 = ( x_s(1,ipoin) - coorg(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x_s(2,ipoin) - coorg(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      auxv(1,ipoig) = auxv(1,ipoig) + rho_s(ipoin)/dist            
      auxv(2,ipoig) = auxv(2,ipoig) + 1.0/dist 
   enddo
   enddo
   
ENDDO      

! hs+Z at grid points

Do ipoig =1,npoig
   if (auxv(2,ipoig).gt.1.e-6) then
      hs_plus_Z_topo (ipoig) = topol(1,ipoig) +  auxv(1,ipoig)/auxv(2,ipoig)
   else
      hs_plus_Z_topo (ipoig) = topol(1,ipoig)
   endif
enddo


deallocate(auxv)   

! to create water particls which are above soil particles  CHUAN end***

read (dat_file,1005) text 
write(chk_file,1005) text
if (icunk_w.eq.10) then 
   read (dat_file,*)    Z_swl, dgx, dgy !      Z_swl is sea water level dgx is deltax for spacing nodes in X   
   write(chk_file,*)    Z_swl, dgx, dgy 
   Rext = 10. * xmaxg
   facthsml = 2.0    !  ** was 1.1
else if (icunk_w.eq.11) then 
   read (dat_file,*)    Z_swl, dgx, dgy, Rext, facthsml   
   write(chk_file,*)    Z_swl, dgx, dgy, Rext, facthsml
else if (icunk_w.eq.12) then             ! MP change July 2019
   read (dat_file,*)    Z_swl, dgx, dgy, xwmin, xwmax, ywmin, ywmax   
   write(chk_file,*)    Z_swl, dgx, dgy, xwmin, xwmax, ywmin, ywmax   
   Rext = 10. * xmaxg
   facthsml = 2.0    !  ** was 1.1
endif   

1005     format(a60) 

lgx   = (xmaxg-xming-dgx)
lgy   = 1.0
if (ndimn.eq.2) lgy   = (ymaxg-yming-dgy)

ngx   = (lgx/dgx) + 1
ngy   = 1
if (ndimn.eq.2) ngy   = (lgy/dgy) + 1

dgx   = lgx/(ngx-1)
dgy   = 1.0
if (ndimn.eq.2) dgy   = lgy/(ngy-1)

darea = dgx*dgy

npoin_w = ngx*ngy
if (.NOT.allocated (xwtmp) ) then
   allocate (xwtmp(ndimn,npoin_w), rhowtmp(npoin_w))
endif
rhowtmp = 0.0 ; xwtmp = 0.0

ipoin_w = 0

DO igy  = 1,ngy
DO igx  = 1,ngx
   xx     = xming + (igx-1)*dgx + dgx/2.        !      ------  Avoid points running out of the domain
   if (ndimn.eq.1) then
      yy  = ( yming + ymaxg )/2.0
   elseif (ndimn.eq.2) then
      yy  = yming + (igy-1)*dgy + dgy/2.
   endif
   RR = (xx*xx +yy*yy)**0.5
   call Get_Z_plus_hs (ipoin_w,xx,yy,Zp,hs_plus_Z_topo)
    ! avoid Zp=-10000   CHUAN
   if (abs(-1000000-zp).lt.1e-6)then
       Zp=Z_swl + 1000
   endif
    ! avoid Zp=-10000   end
   
   if (Zp.lt.Z_swl.and.RR.le.Rext) then
      insid = 1
      if (icunk_w.eq.12) then     ! MP to limit the reservoir domain
         insid = 0
         if ( xx.GE.xwmin.AND.xx.LE.xwmax.AND.yy.GE.ywmin.AND.yy.LE.ywmax) then
            insid = 1
         endif
      endif
      if (insid.eq.1) then
         ipoin_w           = ipoin_w + 1
         xwtmp (1,ipoin_w) = xx 
         if (ndimn.eq.2) then
            xwtmp (2,ipoin_w) = yy
         endif
         rhowtmp(ipoin_w) = Z_swl - Zp
      endif
   endif
   
ENDDO
ENDDO   

npoin_w = ipoin_w

if (.NOT.allocated (x_w) ) then
   allocate (x_w(ndimn,npoin_w), vx_w(ndimn,npoin_w), rho_w(npoin_w), mass_w(npoin_w) )
   allocate ( hsml_w  (npoin_w), itype_w   (npoin_w), p_w  (npoin_w), u_w   (npoin_w) ) 
endif

DO ipoin_w = 1, npoin_w
 
      rho_w  (ipoin_w)   = rhowtmp (ipoin_w) 
      x_w    (:,ipoin_w) = xwtmp(:,ipoin_w) 
      mass_w (ipoin_w)   = rho_w (ipoin_w)*darea
      hsml_w (ipoin_w)   = facthsml * amax1(dgx,dgy)
      itype_w(ipoin_w)   = 2
      vx_w (:,ipoin_w)   = 0.0
!      p_w    (ipoin_w)   = (const(1)/2.)*rho_w  (ipoin_w)*rho_w  (ipoin_w)
      u_w    (ipoin_w)   = 0.0
      
enddo

deallocate (rhowtmp, xwtmp)
 
END SUBROUTINE  Read_WP_Topol


! ----------------------------------------------------------------------   

      subroutine Read_WP_Topol_BACK (icunk_w)

! ----------------------------------------------------------------------   

!     This routine gets water nodes info from topo routine...NEW CALL TO GET_Z_TOPO

implicit none  

integer(ink) ngx, ngy, igx, igy, icunk_w
integer(ink) ipoin_w

real   (irk) Z_swl, dgx, dgy, Rext, facthsml, lgx, lgy, xx, yy, RR 
real   (irk) Zp, dArea

real   (irk), allocatable:: xwtmp(:,:), rhowtmp(:)  ! tmp vars 

! to create water particls which are above soil particles  CHUAN
real   (irk), allocatable:: auxv(:,:)
integer(ink) ipoin, idimn
integer(ink) ipoigx, ipoigy, ipoig
integer(ink) inx0, iny0, inx1,iny1, inx, iny 
integer(ink) ieleg, n1, n2, n3, n4, kpoin
real   (irk) distx2, disty2, dist
real   (irk) xref, yref,xi, eta
real   (irk) sh1, sh2, sh3, sh4,  dZpdx, dZpdy
real   (irk) dsh1dx, dsh2dx, dsh3dx, dsh4dx
real   (irk) dsh1dy, dsh2dy, dsh3dy, dsh4dy
real   (irk) dxginv, dyginv

if (.NOT.allocated (auxv) ) allocate (auxv(2,npoig))
auxv = 0.0

if (.not.allocated( hs_plus_Z_topo ) ) then
   allocate (    hs_plus_Z_topo    (npoig))
   hs_plus_Z_topo= 0.0 
endif

! rho_s at grid points

DO ipoin=1,npoin_s

   !if ( if_Out_Domain(ipoin).eq.1) CYCLE
   !if ( itype(ipoin).ne.1) CYCLE
   xx = x_s(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x_s(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
   !inx0      = max(1,ipoigx-1)
   !iny0      = max(1,ipoigy-1)
   !inx1      = min(npoigx,ipoigx+2)
   !iny1      = min(npoigy,ipoigy+2)
   inx0      = max(1,ipoigx)
   iny0      = max(1,ipoigy)
   inx1      = min(npoigx,ipoigx+1)
   iny1      = min(npoigy,ipoigy+1)   
   
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      distx2 = ( x_s(1,ipoin) - coorg(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x_s(2,ipoin) - coorg(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      auxv(1,ipoig) = auxv(1,ipoig) + rho_s(ipoin)/dist            
      auxv(2,ipoig) = auxv(2,ipoig) + 1.0/dist 
   enddo
   enddo
   
ENDDO      

! hs+Z at grid points

Do ipoig =1,npoig
   if (auxv(2,ipoig).gt.1.e-6) then
      hs_plus_Z_topo (ipoig) = topol(1,ipoig) +  auxv(1,ipoig)/auxv(2,ipoig)
   else
      hs_plus_Z_topo (ipoig) = topol(1,ipoig)
   endif
enddo


deallocate(auxv)   

! to create water particls which are above soil particles  CHUAN end***

read (dat_file,1005) text 
write(chk_file,1005) text
if (icunk_w.eq.10) then 
   read (dat_file,*)    Z_swl, dgx, dgy !      Z_swl is sea water level dgx is deltax for spacing nodes in X   
   write(chk_file,*)    Z_swl, dgx, dgy 
   Rext = 10. * xmaxg
   facthsml = 2.0    !  ** was 1.1
else if (icunk_w.eq.11) then 
   read (dat_file,*)    Z_swl, dgx, dgy, Rext, facthsml   
   write(chk_file,*)    Z_swl, dgx, dgy, Rext, facthsml  
endif   

1005     format(a60) 

lgx   = (xmaxg-xming-dgx)
lgy   = 1.0
if (ndimn.eq.2) lgy   = (ymaxg-yming-dgy)

ngx   = (lgx/dgx) + 1
ngy   = 1
if (ndimn.eq.2) ngy   = (lgy/dgy) + 1

dgx   = lgx/(ngx-1)
dgy   = 1.0
if (ndimn.eq.2) dgy   = lgy/(ngy-1)

darea = dgx*dgy

npoin_w = ngx*ngy
allocate (xwtmp(ndimn,npoin_w), rhowtmp(npoin_w))

ipoin_w = 0

DO igy  = 1,ngy
DO igx  = 1,ngx
   xx     = xming + (igx-1)*dgx + dgx/2.        !      ------  Avoid points running out of the domain
   if (ndimn.eq.1) then
      yy  = ( yming + ymaxg )/2.0
   elseif (ndimn.eq.2) then
      yy  = yming + (igy-1)*dgy + dgy/2.
   endif
   RR = (xx*xx +yy*yy)**0.5
   call Get_Z_plus_hs (ipoin_w,xx,yy,Zp,hs_plus_Z_topo)
    ! avoid Zp=-10000   CHUAN
   if (abs(-1000000-zp).lt.1e-6)then
       Zp=Z_swl + 1000
   endif
    ! avoid Zp=-10000   end
  
   
   
   if (Zp.lt.Z_swl.and.RR.le.Rext) then
      ipoin_w          = ipoin_w + 1
      xwtmp (1,ipoin_w) = xx 
      if (ndimn.eq.2) then
         xwtmp (2,ipoin_w) = yy
      endif
      rhowtmp(ipoin_w) = Z_swl - Zp
   endif
   
ENDDO
ENDDO   

npoin_w = ipoin_w

allocate (x_w(ndimn,npoin_w), vx_w(ndimn,npoin_w), rho_w(npoin_w), mass_w(npoin_w) )
allocate ( hsml_w  (npoin_w), itype_w   (npoin_w), p_w  (npoin_w), u_w   (npoin_w) ) 

DO ipoin_w = 1, npoin_w
 
      rho_w  (ipoin_w)   = rhowtmp (ipoin_w) 
      x_w    (:,ipoin_w) = xwtmp(:,ipoin_w) 
      mass_w (ipoin_w)   = rho_w (ipoin_w)*darea
      hsml_w (ipoin_w)   = facthsml * amax1(dgx,dgy)
      itype_w(ipoin_w)   = 2
      vx_w (:,ipoin_w)   = 0.0
!      p_w    (ipoin_w)   = (const(1)/2.)*rho_w  (ipoin_w)*rho_w  (ipoin_w)
      u_w    (ipoin_w)   = 0.0
      
enddo

deallocate (rhowtmp, xwtmp)
 
END SUBROUTINE  Read_WP_Topol_BACK



! ----------------------------------------------------------------------   

      subroutine Read_WP_Topol_OLD (icunk_w)

! ----------------------------------------------------------------------   

!     This routine gets water nodes info from topo routine...NEW CALL TO GET_Z_TOPO

implicit none  

integer(ink) ngx, ngy, igx, igy, icunk_w
integer(ink) ipoin_w

real   (irk) Z_swl, dgx, dgy, Rext, facthsml, lgx, lgy, xx, yy, RR 
real   (irk) Zp, dArea

real   (irk), allocatable:: xwtmp(:,:), rhowtmp(:)  ! tmp vars 

read (dat_file,1005) text 
write(chk_file,1005) text
if (icunk_w.eq.10) then 
   read (dat_file,*)    Z_swl, dgx, dgy !      Z_swl is sea water level dgx is deltax for spacing nodes in X   
   write(chk_file,*)    Z_swl, dgx, dgy 
   Rext = 10. * xmaxg
   facthsml = 2.0    !  ** was 1.1
else if (icunk_w.eq.11) then 
   read (dat_file,*)    Z_swl, dgx, dgy, Rext, facthsml   
   write(chk_file,*)    Z_swl, dgx, dgy, Rext, facthsml  
endif   

1005     format(a60) 

lgx   = (xmaxg-xming-dgx)
lgy   = 1.0
if (ndimn.eq.2) lgy   = (ymaxg-yming-dgy)

ngx   = (lgx/dgx) + 1
ngy   = 1
if (ndimn.eq.2) ngy   = (lgy/dgy) + 1

dgx   = lgx/(ngx-1)
dgy   = 1.0
if (ndimn.eq.2) dgy   = lgy/(ngy-1)

darea = dgx*dgy

npoin_w = ngx*ngy
allocate (xwtmp(ndimn,npoin_w), rhowtmp(npoin_w))

ipoin_w = 0

DO igy  = 1,ngy
DO igx  = 1,ngx
   xx     = xming + (igx-1)*dgx + dgx/2.        !      ------  Avoid points running out of the domain
   if (ndimn.eq.1) then
      yy  = ( yming + ymaxg )/2.0
   elseif (ndimn.eq.2) then
      yy  = yming + (igy-1)*dgy + dgy/2.
   endif
   RR = (xx*xx +yy*yy)**0.5
   call Get_Z_Topo (ipoin_w,xx,yy,Zp)
   
   if (Zp.lt.Z_swl.and.RR.le.Rext) then
      ipoin_w          = ipoin_w + 1
      xwtmp (1,ipoin_w) = xx 
      if (ndimn.eq.2) then
         xwtmp (2,ipoin_w) = yy
      endif
      rhowtmp(ipoin_w) = Z_swl - Zp
   endif
   
ENDDO
ENDDO   

npoin_w = ipoin_w

allocate (x_w(ndimn,npoin_w), vx_w(ndimn,npoin_w), rho_w(npoin_w), mass_w(npoin_w) )
allocate ( hsml_w  (npoin_w), itype_w   (npoin_w), p_w  (npoin_w), u_w   (npoin_w) ) 

DO ipoin_w = 1, npoin_w
 
      rho_w  (ipoin_w)   = rhowtmp (ipoin_w) 
      x_w    (:,ipoin_w) = xwtmp(:,ipoin_w) 
      mass_w (ipoin_w)   = rho_w (ipoin_w)*darea
      hsml_w (ipoin_w)   = facthsml * amax1(dgx,dgy)
      itype_w(ipoin_w)   = 2
      vx_w (:,ipoin_w)   = 0.0
!      p_w    (ipoin_w)   = (const(1)/2.)*rho_w  (ipoin_w)*rho_w  (ipoin_w)
      u_w    (ipoin_w)   = 0.0
      
enddo

deallocate (rhowtmp, xwtmp)
 
END SUBROUTINE  Read_WP_Topol_OLD 


! ----------------------------------------------------------------------   

      subroutine Get_soil_points

! ----------------------------------------------------------------------   

implicit none      

integer(ink) icunk_s, icc 
real   (irk) dx_sg

ic_geomm_pts = 0      ! default for modifying initial mass MPJune 2023

!      ------   

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    icunk_s,  h_inf_SW   
write(chk_file,*)    icunk_s,  h_inf_SW 

1005     format(a60)

if (icunk_s.eq.0)  then 
    write(*,*) ' no water points in WIR code...want to continue? (0,1) '
    read (*,*) icc
    if (icc.eq.0) STOP
elseif (icunk_s.eq.1) then
    call Read_Cylinder_SW 
elseif (icunk_s.eq.2) then
    call Read_Multilinear_X_SW  
elseif (icunk_s.eq.-2) then
    call Read_Multilinear_X_SW_Trams 
elseif (icunk_s.eq.3) then
    call Read_1Dbk_SW
elseif (icunk_s.eq.4) then
    call Read_3DCloud      (xshiftg,yshiftg)
elseif (icunk_s.eq.5) then
    call Read_Source       (xshiftg,yshiftg)
elseif (icunk_s.eq.6) then
    call Read_Source_plus66(xshiftg,yshiftg) 
elseif (icunk_s.eq.-6) then
    call Read_Source_plus  (xshiftg,yshiftg)      
elseif (icunk_s.eq.7) then                         ! Reads dx and dy    
                                                   ! uses a typical spacing of landslide points 
    call Read_Source_plus_New  (icunk_s, xshiftg, yshiftg)                                                          
                                                   !     avoids void areas         
                                                   !     use for points spacing  
                                                   !     larger than output grid 
    
elseif (icunk_s.eq.8) then         ! Reads dx and dy     MP June 2023 for modif pts generated     
    ic_geomm_pts = 1                                                         
    if (.not.allocated(mass_shift)) then
       allocate (mass_shift(ndimn)) ; mass_shift = 0.0     
    endif                                                   
    call Read_Source_plus_New  (icunk_s, xshiftg, yshiftg)  
                                                                                               
elseif (icunk_s.eq.10.or.icunk_s.eq.9) then      ! 9 is for Tuostolo
    call Read_SPH_Points_sliding_masses_SW (icunk_s, xshiftg,yshiftg)  ! added shift 26 april 2014
elseif (icunk_s.eq.11) then
    call Read_rect_SW (xshiftg,yshiftg) 
elseif (icunk_s.eq.12) then
    call Read_ellipsoidal_SW
elseif (icunk_s.eq.13) then
    call Read_rect_ellipsoid 
elseif (icunk_s.eq.14) then  
    call Read_histogram_SW
elseif (icunk_s.eq.20) then      !      This is used to fill the space of a reservoir with points
    call Read_reservoir_SW
elseif (icunk_s.eq.23) then      ! --- Valpola type. We give (i) cloud of points  
    call Read_Valpola_Type       !                           (ii) elipse + profile to model h 
elseif (icunk_s.eq.200) then     !  ****  RESTART 9th march 2106 
    call Get_Restart_Info_SW     !        new      
elseif (icunk_s.eq.202) then     !  ****  RESTART 1D to 2D 25th March 2016 
    call Get_Restart_1Dto2D      !        new      
elseif (icunk_s.eq.14) then  
    write(*,*) ' incorrect icunk_s option ', icunk_s
endif
  
         
END SUBROUTINE  Get_soil_points

!------------------------------------------
 
   subroutine Get_Restart_Info_SW
 
!------------------------------------------

!      ------  Reads restart info; icunk_s option is 20

implicit none
 
integer(ink)  In_Source, Chk_Source

integer(ink)  len1,  i , ipoin,id 

! time_sph_restart  is time at which computation is restarted 
! (driver time vars MODULE) 

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

ic_restart   = 1

In_Source    = 86
Chk_Source   = 88

   1 format (a60)
1005 format (a60) 
1006 format (a16)
  99 format  (i7,2x,i5,8(2x,g14.7))

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1006) problem_name 
  
len1 = len_trim(problem_name)
 
open (  In_Source ,  file = 'restart.dat'   )
open (  Chk_Source,  file = 'restart.chk'   )

read ( In_Source , 1005) text
write( Chk_Source, 1005) text
read ( In_Source , *)    npoin_s, time_sph_restart     
write( Chk_Source, *)    npoin_s, time_sph_restart   

time = time_sph_restart

!      ------ all Fi_s vars have been defined as allocatable in SW vars module

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

x_s = 0.0; vx_s=0.0; rho_s=0.0; mass_s=0.0;hsml_s=0.0; u_s = 0.0

read ( In_Source , 1005) text
write( Chk_Source, 1005) text

DO i = 1, npoin_s 

   read (In_Source, *) ipoin, itype_s(i), (x_s(id,i),id=1,ndimn), (vx_s(id,i),id=1,ndimn), &
                               rho_s (i), mass_s(i), hsml_s(i), u_s(i)

   write(Chk_Source,99) ipoin, itype_s(i), (x_s(id,i),id=1,ndimn), (vx_s(id,i),id=1,ndimn), &
                               rho_s (i), mass_s(i), hsml_s(i), u_s(i)
ENDDO

close ( In_Source)
close (Chk_Source)

END subroutine Get_Restart_Info_SW


!------------------------------------------
 
   subroutine Get_Restart_1Dto2D
 
!------------------------------------------

!      ------  Reads restart info; icunk_s option is 20

implicit none
 
integer(ink)  In_Source, Chk_Source
integer(ink)  npg, npyg
integer(ink)  len1,  i , ipoin,id, idb, ipoin_s 

real   (irk)  xd0, thetag, xgmin0, ygmin0, bg, db, bg2
real   (irk)  xgmin , ygmin, xi, eta, thetarad, costh, sinth 
real   (irk)  ax, bx
real, allocatable :: xg(:), vxg(:), rhog(:), massg(:), hsmlg(:), ug(:)
integer, allocatable:: itypeg(:)

! time_sph_restart  is time at which computation is restarted 
! (driver time vars MODULE) 

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

In_Source    = 86
Chk_Source   = 88

   1 format (a60)
1005 format (a60) 
1006 format (a16)
  99 format  (i7,2x,i5,8(2x,g14.7))

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1006) problem_name 
  
len1 = len_trim(problem_name)
 
open (  In_Source ,  file = 'restart.dat'   )
open (  Chk_Source,  file = 'restart.chk'   )

read ( In_Source , 1005) text
write( Chk_Source, 1005) text
read ( In_Source , *)    xd0, thetag, xgmin0, ygmin0, bg, npyg 
write( Chk_Source, *)    xd0, thetag, xgmin0, ygmin0, bg, npyg 
thetarad = thetag *4 * atan(1.0) /180.
costh = cos(thetarad)
sinth = sin(thetarad)
bg2 = bg/2.
xgmin = xgmin0 + bg2*sinth
ygmin = ygmin0 - bg2*costh

read ( In_Source , 1005) text
write( Chk_Source, 1005) text
read ( In_Source , *)    npg, time_sph_restart     
write( Chk_Source, *)    npg, time_sph_restart   

time = time_sph_restart

!      ------ allocate 1D Fi_s vars  

allocate ( xg(npg),  rhog(npg), massg (npg), hsmlg (npg) )
allocate (vxg(npg),  ug  (npg), itypeg(npg) )

xg = 0.0; rhog = 0.0; massg = 0.0; hsmlg = 0.0; vxg = 0.0; ug = 0.0; itypeg = 0

read ( In_Source , 1005) text
write( Chk_Source, 1005) text

DO i = 1, npg 

   read (In_Source, *)  ipoin, itypeg(i), xg(i), vxg(i), &
                               rhog (i), massg(i), hsmlg(i), ug(i)

   write(Chk_Source,99) ipoin, itypeg(i), xg(i), vxg(i), &
                               rhog (i), massg(i), hsmlg(i), ug(i)
ENDDO

npoin_s = npg*npyg
db      = bg/(npyg -1)

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate ( vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )
x_s = 0.0; rho_s = 0.0; mass_s=0.0; hsml_s = 0.0; vx_s = 0.0; u_s = 0.0; itype_s = 0

do i   = 1, npg
   xi  = xg(i) - xd0
    if (i.eq.1) then
     ax = (xg(i+1) - xg(i))/2
    else if (i.eq.npg) then
     ax = (xg(i) - xg(i-1))/2
    else
     ax = (xg(i+1) - xg(i-1))/2
    end if
   do idb = 1, npyg
      ipoin_s = (idb-1)*npg + i 
      eta = (idb-1)*db
      bx = db
       if ((idb.eq.1).or.(idb.eq.npyg)) then
        bx = db/2
       end if
      itype_s (ipoin_s) = itypeg(i)
      x_s   (1,ipoin_s) = xgmin + xi*costh - eta*sinth
      x_s   (2,ipoin_s) = ygmin + xi*sinth + eta*costh
      vx_s  (1,ipoin_s) = vxg(i) * costh
      vx_s  (2,ipoin_s) = vxg(i) * sinth 
      rho_s   (ipoin_s) = rhog  (i)
      mass_s  (ipoin_s) = ax * bx * rhog(i)  
      hsml_s  (ipoin_s) = (ax*ax + bx*bx)**0.5
      u_s     (ipoin_s) = ug    (i)
      write(Chk_Source,*) ' ----------------  chk generated pts ---------'
      write(Chk_Source,99) ipoin, itype_s (ipoin_s), &
                                  (x_s(id,ipoin_s),id=1,ndimn), (vx_s(id,ipoin_s),id=1,ndimn), &
                                   rho_s (ipoin_s), mass_s(ipoin_s), hsml_s(ipoin_s), u_s(ipoin_s)
   enddo
enddo      

close ( In_Source)
close (Chk_Source)
 
!      deallocate 1D vars 

deallocate ( xg,  rhog , massg , hsmlg  )
deallocate (vxg ,  ug  , itypeg  )

END subroutine Get_Restart_1Dto2D

!------------------------------------------
 
   subroutine Get_Restart_1Dto2D_v0
 
!------------------------------------------

!      ------  Reads restart info; icunk_s option is 20

implicit none
 
integer(ink)  In_Source, Chk_Source
integer(ink)  npg
integer(ink)  len1,  i , ipoin,id, idb, ipoin_s 

real   (irk)  xd0, thetag, xgmin0, ygmin0, bg, npyg, db, bg2
real   (irk)  xgmin , ygmin, xi, eta, thetarad, costh, sinth 
real, allocatable :: xg(:), vxg(:), rhog(:), massg(:), hsmlg(:), ug(:)
integer, allocatable:: itypeg(:)

! time_sph_restart  is time at which computation is restarted 
! (driver time vars MODULE) 

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

In_Source    = 86
Chk_Source   = 88

   1 format (a60)
1005 format (a60) 
1006 format (a16)
  99 format  (i7,2x,i5,8(2x,g14.7))

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1006) problem_name 
  
len1 = len_trim(problem_name)
 
open (  In_Source ,  file = 'restart.dat'   )
open (  Chk_Source,  file = 'restart.chk'   )

read ( In_Source , 1005) text
write( Chk_Source, 1005) text
read ( In_Source , *)    xd0, thetag, xgmin0, ygmin0, bg, npyg 
write( Chk_Source, *)    xd0, thetag, xgmin0, ygmin0, bg, npyg 
thetarad = thetag *4 * atan(1.0) /180.
costh = cos(thetarad)
sinth = sin(thetarad)
bg2 = bg/2.
xgmin = xgmin0 + bg2*sinth
ygmin = ygmin0 - bg2*costh

read ( In_Source , 1005) text
write( Chk_Source, 1005) text
read ( In_Source , *)    npg, time_sph_restart     
write( Chk_Source, *)    npg, time_sph_restart   

time = time_sph_restart

!      ------ allocate 1D Fi_s vars  

allocate ( xg(npg),  rhog(npg), massg (npg), hsmlg (npg) )
allocate (vxg(npg),  ug  (npg), itypeg(npg) )

xg = 0.0; vxg=0.0; rhog=0.0; massg=0.0;hsmlg=0.0; ug = 0.0

read ( In_Source , 1005) text
write( Chk_Source, 1005) text

DO i = 1, npg 

   read (In_Source, *)  ipoin, itypeg(i), xg(i), vxg(i), &
                               rhog (i), massg(i), hsmlg(i), ug(i)

   write(Chk_Source,99) ipoin, itypeg(i), xg(i), vxg(i), &
                               rhog (i), massg(i), hsmlg(i), ug(i)
ENDDO

npoin_s = npg*npyg
db      = bg/(npyg -1)

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate ( vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

do i   = 1, npg
   xi  = xg(i) - xd0
   do idb = 1, npyg
      ipoin_s = (idb-1)*npg + i 
      eta = (idb-1)*db
      itype_s (ipoin_s) = itypeg(i)
      x_s   (1,ipoin_s) = xgmin + xi*costh - eta*sinth
      x_s   (2,ipoin_s) = ygmin + xi*sinth + eta*costh
      vx_s  (1,ipoin_s) = vxg(i) * costh
      vx_s  (1,ipoin_s) = vxg(i) * sinth 
      rho_s   (ipoin_s) = rhog  (i)
      mass_s  (ipoin_s) = massg (i) 
      hsml_s  (ipoin_s) = hsmlg (i)
      u_s     (ipoin_s) = ug    (i)
      write(Chk_Source,*) ' ----------------  chk generated pts ---------'
      write(Chk_Source,99) ipoin, itype_s (ipoin_s), &
                                  (x_s(id,ipoin_s),id=1,ndimn), (vx_s(id,ipoin_s),id=1,ndimn), &
                                   rho_s (ipoin_s), mass_s(ipoin_s), hsml_s(ipoin_s), u_s(ipoin_s)
   enddo
enddo      

mass_s = mass_s/10.

close ( In_Source)
close (Chk_Source)
 
!      deallocate 1D vars 

deallocate ( xg,  rhog , massg , hsmlg  )
deallocate (vxg ,  ug  , itypeg  )

END subroutine Get_Restart_1Dto2D_v0

! ----------------------------------------------------------------------   

  subroutine Read_Valpola_Type

! ---------------------------------------------------------------------- 

integer(ink) nptsh
integer(ink) i, ipoin, iptsh, ick


real   (irk), allocatable:: xptsh(:), hptsh(:)
   
real   (irk)  deltpts, facthsml
real   (irk)  x0pts, y0pts,a0pts, b0pts, expmpts
real   (irk)  cgra05
real   (irk)  xpts , ypts, bpts, xi, eta
real   (irk)  fxpts, fypts, fhpts
real   (irk)  x0ptsh, x1ptsh, f0ptsh, f1ptsh 

read (dat_file,1005) text                               !      ------  Read spacing (used for mass) and hsml
write(chk_file,1005) text 
read (dat_file,*)    npoin_s, deltpts, facthsml       
write(chk_file,*)    npoin_s, deltpts, facthsml 

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s   (npoin_s), itype_s(npoin_s) )

                                                        ! allocate p_s(npoin_s) deleted

read (dat_file,1005) text                               !      ------  Read points (generated by GID)
write(chk_file,1005) text       

cgra05       = 9.81/2.0       !               

DO ipoin = 1,npoin_s
   read (dat_file,*) x_s(1,ipoin), x_s(2,ipoin)
   write(chk_file,*) x_s(1,ipoin), x_s(2,ipoin)
ENDDO
      
read (dat_file,1005) text                               !      ------  Properties of ellipse
write(chk_file,1005) text
read (dat_file,*)    x0pts, y0pts,a0pts, b0pts, expmpts     
write(chk_file,*)    x0pts, y0pts,a0pts, b0pts, expmpts
      
read (dat_file,1005) text                               !      ------  Profile along X
write(chk_file,1005) text
read (dat_file,*)    nptsh     
write(chk_file,*)    nptsh

allocate ( xptsh(nptsh) )
allocate ( hptsh(nptsh) ) 

read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*) (xptsh(i),i=1,nptsh)
write(chk_file,*) (xptsh(i),i=1,nptsh)
read (dat_file,*) (hptsh(i),i=1,nptsh)
write(chk_file,*) (hptsh(i),i=1,nptsh)

DO ipoin = 1, npoin_s                                      !      ------  for all points get h = profile * conoid

   xpts  = x_s(1,ipoin)
   ypts  = x_s(2,ipoin)
   xi    = (xpts - x0pts)/a0pts                         !      ------  Conoid starts -> get fypts
   bpts  = b0pts * ( 1.0 - xi*xi )**0.5
   eta   = (ypts - y0pts)/bpts
   fypts = ( 1.0 - eta**expmpts )
   ick   = 0 
   fxpts = 0.0
   iptsh = 0 
   DO WHILE (iptsh.lt.nptsh-1.AND.ick.eq.0)
      iptsh  = iptsh +1  
      x0ptsh = xptsh(iptsh)
      x1ptsh = xptsh(iptsh+1)
      if (xpts.ge.x0ptsh.and.xpts.le.x1ptsh) then
          f0ptsh = (xpts - x1ptsh)/(x0ptsh - x1ptsh)
          f1ptsh = 1.-f0ptsh
          fxpts  = f0ptsh*hptsh(iptsh) + f1ptsh*hptsh(iptsh+1)
                  ick    = 1
      endif
   enddo
   
   fhpts  = fxpts*fypts
   if (fhpts.lt.0.1) fhpts = 0.1

   vx_s (:,ipoin) = 0.0
   itype_s(ipoin) = 1
   rho_s  (ipoin) = fhpts
   mass_s (ipoin) = fhpts   * deltpts * deltpts
   hsml_s (ipoin) = deltpts * facthsml
   u_s    (ipoin) = 0.0  
   
   write(chk_file,1111) ipoin, xpts, ypts, xi, eta, fxpts, fypts, rho_s(ipoin)
1111 format(i5,2x,7(f10.4,2x))

ENDDO

1005     format(a60) 

deallocate(xptsh)
deallocate(hptsh)         

END subroutine Read_Valpola_Type  


! ----------------------------------------------------------------------   

      subroutine Read_SPH_Points_sliding_masses_SW (icunk_s, xshiftg, yshiftg)

! ---------------------------------------------------------------------- 

!     This subroutine is used to read sliding masses info. 

implicit none  
 
integer(ink)  NP_slm, In_source, Chk_source, len1 
integer(ink)  islm, ip_slm, ipoin0, ipoin, idimn, icunk_s 
real   (irk)  Area_slm, Height_slm, mass_slm,  twopi
real   (irk)  hsml_slm, facthsml2D
real   (irk)  xshiftg, yshiftg

character(60) problem_name            ! name of PFC3D file
character(60) text                    ! general purpose char string

twopi        = 4.0*atan(1.0)
Pin_slm(:)   = 0
ipoin0       = 0

In_Source    = 86
Chk_Source   = 88

! -----------------------  changed ---------------

1005 format (a60)
1006 format (a16)
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 

len1 = len_trim(problem_name)

open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'     )

! -----------------------  changed ---------------

!      ------  Read  

read (In_Source, 1005) text
write(Chk_Source,1005) text
read (In_Source, *)    NslMasses, npoin_s, facthsml2D
write(Chk_Source,*)    NslMasses, npoin_s, facthsml2D

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

                                                        ! allocate p_s(npoin_s) deleted
if (.not.allocated(t0_slm)) then                        ! MP 19 march 2022
   allocate( t0_slm(NslMasses), Pin_slm(npoin_s))
endif

If (icunk_s.eq.10) then

   DO islm = 1,Nslmasses
      read (In_Source, 1005) text
      write(Chk_Source,1005) text
      read (In_Source, *)    NP_slm, Area_slm, Height_slm, hsml_slm, T0_slm(islm) 
      write(Chk_Source,*)    NP_slm, Area_slm, Height_slm, hsml_slm, T0_slm(islm)
      mass_slm = Area_slm*Height_slm/NP_slm
      read (In_Source, 1005) text
      write(Chk_Source,1005) text
      DO ip_slm = 1, NP_slm
         ipoin  = ipoin0 + ip_slm
         read (In_Source, *) (x_s(idimn,ipoin),idimn=1,ndimn) 
         write(Chk_Source,*) (x_s(idimn,ipoin),idimn=1,ndimn) 
         rho_s  (ipoin) = Height_slm
         mass_s (ipoin) = mass_slm
         hsml_s (ipoin) = hsml_slm*facthsml2D
         vx_s (:,ipoin) = 0.0
         u_s    (ipoin) = 0.0
         itype_s(ipoin) = 1
         Pin_slm(ipoin) = ip_slm
      enddo
      ipoin0 = ipoin0 + NP_slm
   Enddo

ElseIf (icunk_s.eq.9) then  !  Assumes points equally distributed

   DO islm = 1,Nslmasses
      read (In_Source, 1005) text
      write(Chk_Source,1005) text
      read (In_Source, *)    NP_slm, Area_slm,  hsml_slm, T0_slm(islm) 
      write(Chk_Source,*)    NP_slm, Area_slm,  hsml_slm, T0_slm(islm)
      mass_slm = Area_slm / NP_slm
      read (In_Source, 1005) text
      write(Chk_Source,1005) text
      DO ip_slm = 1, NP_slm
         ipoin  = ipoin0 + ip_slm
         read (In_Source, *) (x_s(idimn,ipoin),idimn=1,ndimn), Height_slm 
         write(Chk_Source,*) (x_s(idimn,ipoin),idimn=1,ndimn), Height_slm 
         rho_s  (ipoin) = Height_slm
         mass_s (ipoin) = mass_slm*Height_slm
         hsml_s (ipoin) = hsml_slm*facthsml2D
         vx_s (:,ipoin) = 0.0
         u_s    (ipoin) = 0.0
         itype_s(ipoin) = 1
         Pin_slm(ipoin) = ip_slm
      enddo
      ipoin0 = ipoin0 + NP_slm
   Enddo

Endif

x_s(1,:) = x_s(1,:) - xshiftg
if (ndimn.eq.2) then
    x_s(2,:) = x_s(2,:) - yshiftg
endif    

close (  In_Source  )
close (  Chk_Source )
   

END subroutine Read_SPH_Points_sliding_masses_SW 



! ----------------------------------------------------------------------   

      subroutine Read_ellipsoidal_SW  

! ----------------------------------------------------------------------   

!      Reads ellipsoid 

implicit none  
 
integer(ink) ipoin, nxi, neta, ixi, ieta
   
real   (irk) xp0, yp0, a, b, theta0, hp0,  facthsml
real   (irk) cgra05
real   (irk) rpp, drr, bmin,xi, eta, dxi, deta
real   (irk) dArea, xx, yy, hfact 
 
!      ------   

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    xp0, yp0, a, b, theta0, hp0       
write(chk_file,*)    xp0, yp0, a, b, theta0, hp0 

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    drr, facthsml       
write(chk_file,*)    drr, facthsml  

!    cgra05  = const(1)/2.0
pi      = 4.*atan(1.0)

nxi     = 2*a/drr + 1
dxi     = 2*a/(nxi-1)
neta    = 2*b/drr + 1           !     -------  This is an upper value
npoin_s = nxi*neta              !              will be less, we are allocating for the rectangle
                                        
allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

                                                        ! allocate p_s(npoin_s) deleted

ipoin = 0
DO ixi   = 1, nxi
   xi    = -a + (ixi-1)*dxi
   bmin  = (1-xi/a)*(1+xi/a)
   if (bmin.gt.0.0) then
       bmin = b * sqrt(bmin)
   endif
   neta  = 2*bmin/drr + 1
   if (neta.gt.1) then
       deta = 2*bmin/(neta-1)
   else
       deta  = drr
   endif 
   dArea = dxi*deta
   DO ieta = 1, neta
      eta   = -bmin + (ieta-1)*deta
      rpp   =  (xi/a)*(xi/a) + (eta/b)*(eta/b) 
      IF ( rpp.le.1 ) THEN
         ipoin = ipoin +1 
         xx    = xp0 + xi*cos(theta0) - eta*sin(theta0)
         yy    = yp0 + xi*sin(theta0) + eta*cos(theta0)
         hfact = (1.-rpp) + 0.1*rpp             !   Varies between (0.1, 1)
         x_s   (1,ipoin)   = xx
         x_s   (2,ipoin)   = yy  
         vx_s  (:,ipoin)   = 0.0
         itype_s (ipoin)   = 1
         rho_s   (ipoin)   = hp0 * hfact
         mass_s  (ipoin)   = rho_s(ipoin) * dArea
         hsml_s  (ipoin)   = drr*facthsml
         u_s     (ipoin)   = 0.0
      ENDIF   
   ENDDO   
ENDDO 

npoin_s = ipoin

1005     format(a60)  


END SUBROUTINE  Read_ellipsoidal_SW



! ----------------------------------------------------------------------   

  subroutine Read_rect_ellipsoid

! ---------------------------------------------------------------------- 

implicit none  
 
integer(ink) npoinx, npoiny
integer(ink) ipoinx, ipoiny, ipoin
   
real   (irk) x0, y0, x1, y1,  b, d0, Zfact, facthsml
real   (irk) xs, ys    
real   (irk) xc, yc,  a, vnxi, vneta, vvxi, vveta 
real   (irk) xi, eta, area_Poin, cgra05  

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    npoinx, npoiny       
write(chk_file,*)    npoinx, npoiny

npoin_s = npoinx*npoiny

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )
rho_s  = 0.0
mass_s = 0.0            ! allocate p_s(npoin_s) deleted
                                                        
read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    x0,y0,x1,y1,b,d0,facthsml       
write(chk_file,*)    x0,y0,x1,y1,b,d0,facthsml

write(*,*) ' --- we are reading from cylider/hemisphere point generation -----'

! cgra05 = 9.81/2
 
xc    = (x0+x1)/2.
yc    = (y0+y1)/2.
a     = ( (x1-xc)**2.+(y1-yc)**2 )**0.5
vnxi  = (x1-xc)/a
vneta = (y1-yc)/a
vvxi  = - vneta
vveta =   vnxi

xs = x0 + b*vneta
ys = y0 - b*vnxi

area_Poin = 4*a*b/npoin_s

ipoin     = 0
do ipoiny = 1, npoiny
do ipoinx = 1, npoinx

   ipoin  = ipoin + 1
   
   x_s(1,ipoin) = xs + 2.*a*vnxi *(ipoinx-1)/(npoinx-1) + 2.*b*vvxi *(ipoiny-1)/(npoiny-1)
   x_s(2,ipoin) = ys + 2.*a*vneta*(ipoinx-1)/(npoinx-1) + 2.*b*vveta*(ipoiny-1)/(npoiny-1)
   
   xi    = ( x_s(1,ipoin)- xc )*vnxi + ( x_s(2,ipoin)-yc )*vneta
   eta   = ( x_s(1,ipoin)- xc )*vvxi + ( x_s(2,ipoin)-yc )*vveta
   if (abs(xi).le.a.and.abs(eta).le.b) then
      rho_s(ipoin) = d0*(1-xi/a)*(1+xi/a)*(1-eta/b)*(1+eta/b) + (d0*d0/npoin_s)**0.5
   else
      rho_s(ipoin) = (d0*d0/npoin_s)**0.5
   endif
   
   vx_s  (:,ipoin) = 0.0
   itype_s (ipoin) = 1
   mass_s  (ipoin) = rho_s(ipoin)* area_Poin
   hsml_s  (ipoin) = (area_Poin**0.5) * facthsml
   u_s     (ipoin) = 0.0     
   
enddo
enddo

1005     format(a60) 

END subroutine Read_rect_ellipsoid  

! ----------------------------------------------------------------------   

      subroutine Read_Cylinder_SW  

! ----------------------------------------------------------------------   

!      Reads cylider or ellipsoid circular base

implicit none  

 
integer(ink) icdens, nlinep, ipoin, ilinep, ith   
integer(ink) icunk_s, icc, itp_ps 
   
real   (irk) xp0, yp0, rp0, hp0, cgra05, mss, drr, delth, rpp, hfact, theta
real   (irk) rp1, hp1, factp1     
real   (irk) rx ,  ry
real   (irk) facthsml, twopi 
 
!      ------   

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    xp0, yp0, rp0, hp0, icdens       
write(chk_file,*)    xp0, yp0, rp0, hp0, icdens 

factp1 = 1.0
if (icdens.eq.2) then      !      --- Inner cyl with a diff heigth...read radius and h
    read (dat_file,1005) text         
    write(chk_file,1005) text        
    read (dat_file,*)    rp1, hp1      
    write(chk_file,*)    rp1, hp1
    factp1 = hp1/hp0
endif

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    nlinep, facthsml       
write(chk_file,*)    nlinep, facthsml  

write(chk_file,*) ' --- we are reading from cylider/hemisphere point generation -----'
write(chk_file,*) ' give factor for smoothing length. Total is 2*facthsml'
write(chk_file,*)  facthsml

! cgra05  = const(1)/2.0
                                        ! nlinep  = (-3 + sqrt (9. + 12.*real(npoin-1)))/6
npoin_s = 1 + 3*nlinep*(nlinep+1)
mss     = pi*rp0*rp0*hp0/npoin_s
if (nlinep.ge.1) then  ! changed MP 28.12.2018 use 1 point only
    drr = (rp0/nlinep)  ! rr=rp0 and hsml = diameter   
else
    drr = rp0
endif

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

ipoin        = 1
x_s (1,1)    = xp0
x_s (2,1)    = yp0
vx_s(:,1)    = 0.0
itype_s (1)  = 1
rho_s   (1)  = hp0*factp1
mass_s  (1)  = mss*factp1
hsml_s  (1)  = 2.*drr
u_s     (1)  = 0.0       

DO ilinep = 1,nlinep
   delth  = 2.*pi/(6.*ilinep)
   rpp    = (ilinep)*drr
   if ( (rpp/rp0).gt.1.00 ) rpp=0.9999999*rp0 ! avoids -ve root below
   factp1 = 1.0
   if (icdens.eq.2.and.rpp.le.rp1) factp1 = hp1/hp0
    
   hfact  = sqrt(1.0-(rpp/rp0)*(rpp/rp0))
   if (hfact.lt.0.1) hfact=0.1
   DO ith = 1, 6*ilinep
      ipoin = ipoin + 1
      if (ipoin.gt.npoin_s) exit !       ------  take care of roundoffs
      theta = ith*delth
      rx    = rpp*cos(theta)
      ry    = rpp*sin(theta)
      x_s   (1,ipoin)   = xp0 + rx
      x_s   (2,ipoin)   = yp0 + ry
      vx_s  (:,ipoin)   = 0.0
      itype_s (ipoin)   = 1
      if (icdens.eq.0) then
         rho_s  (ipoin)  = hp0   
         mass_s (ipoin)  = mss
      elseif (icdens.eq.1) then
         rho_s (ipoin)   = hp0*hfact  
         mass_s(ipoin)   = mss*hfact
      elseif (icdens.eq.2) then
         rho_s (ipoin)   = hp0*factp1  
         mass_s(ipoin)   = mss*factp1
      endif
     hsml_s  (ipoin) = 2.*drr*facthsml
     u_s     (ipoin) = 0.0
   ENDDO
ENDDO   

1005     format(a60)  


END SUBROUTINE  Read_Cylinder_SW


! ----------------------------------------------------------------------   

        subroutine Read_Multilinear_X_SW

! ---------------------------------------------------------------------- 

implicit none

integer(ink) nptsh , iptsh, ipoin_s, i
real   (irk) x0ptsh, x1ptsh, xxptsh, f0ptsh, f1ptsh, hhptsh, deltptsh
real   (irk) cgra05
real   (irk) facthsml
real   (irk), allocatable:: xptsh(:), hptsh(:)
       
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_s, facthsml ,nptsh     
write(chk_file,*)    npoin_s, facthsml ,nptsh 

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) ) 
allocate ( xptsh(nptsh),        hptsh(nptsh) ) ! allocate p_s(npoin_s) deleted

read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*) (xptsh(i),i=1,nptsh)
write(chk_file,*) (xptsh(i),i=1,nptsh)
read (dat_file,*) (hptsh(i),i=1,nptsh)
write(chk_file,*) (hptsh(i),i=1,nptsh)

cgra05       = 9.8/2.0       !------ we do it here because properties have not been read!!

if (npoin_s.gt.1) then
   deltptsh     = (xptsh(nptsh)-xptsh(1))/(npoin_s - 1)!!
   do ipoin_s     = 1, npoin_s
      xxptsh    = xptsh(1) + (ipoin_s-1)*deltptsh
      do iptsh  = 1,nptsh-1
         x0ptsh = xptsh(iptsh)
         x1ptsh = xptsh(iptsh+1)
         if (xxptsh.ge.x0ptsh.and.xxptsh.le.x1ptsh) then
             f0ptsh = (xxptsh-x1ptsh)/(x0ptsh-x1ptsh)
             f1ptsh = 1.-f0ptsh
             hhptsh = f0ptsh*hptsh(iptsh) + f1ptsh*hptsh(iptsh+1)
!             if (hhptsh.le.1.e-3) hhptsh = 1.e-3   *** bug killed 7th march 16
         endif
      enddo
      x_s  (1,ipoin_s) = xxptsh
      vx_s (1,ipoin_s) = 0.0
      itype_s(ipoin_s) = 1
      rho_s  (ipoin_s) = hhptsh
      mass_s (ipoin_s) = hhptsh*deltptsh
      if ((ipoin_s.eq.1).or.(ipoin_s.eq.npoin_s)) mass_s(ipoin_s)=mass_s(ipoin_s)/2.0
      hsml_s (ipoin_s) = 2.*deltptsh*facthsml
      u_s    (ipoin_s) = 0.0
   enddo
elseif (npoin_s.eq.1) then
      x_s  (1,1) = xptsh(1)
      vx_s (1,1) = 0.0
      itype_s(1) = 1
      rho_s  (1) = hptsh(1)
      mass_s (1) = hptsh(1) 
      hsml_s (1) = 2.*facthsml
      u_s    (1) = 0.0
endif      

1005     format(a60) 

deallocate(xptsh)
deallocate(hptsh)

END subroutine Read_Multilinear_X_SW


! ----------------------------------------------------------------------   

        subroutine Read_Multilinear_X_SWOLD

! ---------------------------------------------------------------------- 

implicit none

integer(ink) nptsh , iptsh, ipoin_s, i
real   (irk) x0ptsh, x1ptsh, xxptsh, f0ptsh, f1ptsh, hhptsh, deltptsh
real   (irk) cgra05
real   (irk) facthsml
real   (irk), allocatable:: xptsh(:), hptsh(:)
       
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_s, facthsml ,nptsh     
write(chk_file,*)    npoin_s, facthsml ,nptsh 

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) ) 

                                                        ! allocate p_s(npoin_s) deleted
allocate ( xptsh(nptsh),        hptsh(nptsh) ) 

read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*) (xptsh(i),i=1,nptsh)
write(chk_file,*) (xptsh(i),i=1,nptsh)
read (dat_file,*) (hptsh(i),i=1,nptsh)
write(chk_file,*) (hptsh(i),i=1,nptsh)

deltptsh     = (xptsh(nptsh)-xptsh(1))/(npoin_s - 1)
cgra05       = 9.8/2.0       !                  ------ we do it here because properties have not been read!!!!
do ipoin_s     = 1, npoin_s
   xxptsh    = xptsh(1) + (ipoin_s-1)*deltptsh
   do iptsh  = 1,nptsh-1
      x0ptsh = xptsh(iptsh)
      x1ptsh = xptsh(iptsh+1)
      if (xxptsh.ge.x0ptsh.and.xxptsh.le.x1ptsh) then
          f0ptsh = (xxptsh-x1ptsh)/(x0ptsh-x1ptsh)
          f1ptsh = 1.-f0ptsh
          hhptsh = f0ptsh*hptsh(iptsh) + f1ptsh*hptsh(iptsh+1)
          if (hhptsh.le.1.e-3) hhptsh = 1.e-3
      endif
   enddo
   x_s  (1,ipoin_s) = xxptsh
   vx_s (1,ipoin_s) = 0.0
   itype_s(ipoin_s) = 1
   rho_s  (ipoin_s) = hhptsh
   mass_s (ipoin_s) = hhptsh*deltptsh
   if ((ipoin_s.eq.1).or.(ipoin_s.eq.npoin_s)) mass_s(ipoin_s)=mass_s(ipoin_s)/2.0
   hsml_s (ipoin_s) = 2.*deltptsh*facthsml
   u_s    (ipoin_s) = 0.0
enddo

1005     format(a60) 

deallocate(xptsh)
deallocate(hptsh)

END subroutine Read_Multilinear_X_SWOLD


! ----------------------------------------------------------------------   

        subroutine Read_Multilinear_X_SW_trams

! ---------------------------------------------------------------------- 

implicit none


integer(ink) nptsh , ipoin_s, i, idiv
integer(ink) ipoin, nseg, iseg
integer(ink), allocatable:: ndivsh(:)

real   (irk) x0ptsh, x1ptsh, xxptsh, f0ptsh, f1ptsh, hhptsh, deltptsh
real   (irk) facthsml
real   (irk), allocatable:: xptsh(:), hptsh(:)   
       
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_s, facthsml ,nptsh     
write(chk_file,*)    npoin_s, facthsml ,nptsh

nseg = nptsh -1

allocate ( ndivsh (nseg) )  
allocate ( xptsh    (nptsh), hptsh(nptsh) )

read (dat_file,1005) text
write(chk_file,1005) text 
read (dat_file,*) (xptsh(i) ,i=1,nptsh)
write(chk_file,*) (xptsh(i) ,i=1,nptsh)
read (dat_file,*) (hptsh(i) ,i=1,nptsh)
write(chk_file,*) (hptsh(i) ,i=1,nptsh)
read (dat_file,*) (ndivsh(i),i=1,nseg )
write(chk_file,*) (ndivsh(i),i=1,nseg )

npoin_s = 1                     !    -- obtain total nr of points
do iseg = 1, nseg
   npoin_s = npoin_s + ndivsh(iseg)
enddo

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )    

ipoin_s = 0
DO iseg = 1, nseg
   deltptsh  = (xptsh(iseg+1)-xptsh(iseg))/ndivsh(iseg)
   x0ptsh    = xptsh(iseg)
   x1ptsh    = xptsh(iseg+1)
   DO idiv   = 1, ndivsh(iseg)
      xxptsh = x0ptsh + deltptsh*(idiv-1)
      f0ptsh = (xxptsh-x1ptsh)/(x0ptsh-x1ptsh)
      f1ptsh = 1.-f0ptsh
      hhptsh = f0ptsh*hptsh(iseg) + f1ptsh*hptsh(iseg+1)
      ipoin_s = ipoin_s +1   
      x_s  (1,ipoin_s) = xxptsh
      vx_s (1,ipoin_s) = 0.0
      itype_s(ipoin_s) = 1
      rho_s  (ipoin_s) = hhptsh
      mass_s (ipoin_s) = hhptsh*deltptsh
      hsml_s (ipoin_s) = 2.*deltptsh*facthsml
   ENDDO   
ENDDO 

x_s  (1,npoin_s) = xptsh(nptsh)
vx_s (1,npoin_s) = 0.0
itype_s(npoin_s) = 1
rho_s  (npoin_s) = hptsh(nptsh)
where (rho_s.lt.h_inf_SW) rho_s = h_inf_SW
mass_s (npoin_s) = hhptsh*deltptsh/2.0
mass_s (1)       = mass_s(1)/2.0
hsml_s(npoin_s)  = hsml_s(npoin_s-1)
u_s = 0.0


1005     format(a60) 

deallocate(xptsh)
deallocate(hptsh)
deallocate(ndivsh)

END subroutine Read_Multilinear_X_SW_trams

! ----------------------------------------------------------------------   

      subroutine Read_1Dbk_SW
      
! ----------------------------------------------------------------------   

!     This subroutine is used to generate initial data for the 
!     1 d break dam problem
!      
!     Reads:                Xl,Xc,Xr   (X0 is at the dam)
!                           hl and Hr  heights of water left and right of X0
!                           npoin      number of particles to be generated
! 
!      

implicit none  

integer(ink)  nl, nr, i , iomega  
real   (irk)  Xl,   Xc,   Xr,  hl,  hr, mss, pleft, pright
real   (irk)  xll, xrr, htot,  fl,  fr, dxl,   dxr, dxlr

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    npoin_s   
write(chk_file,*)    npoin_s 

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s(npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    Xl, Xc, Xr, hl, hr, iomega   
write(chk_file,*)    Xl, Xc, Xr, hl, hr, iomega 

!      ------       Initialize

xll   = (Xc-Xl)
xrr   = (Xr-Xc)
htot  = xll*hl + xrr*hr
fl    = xll*hl/htot
fr    = 1.-fl
nl    = npoin_s * fl
nr    = npoin_s - nl
dxl   = xll/(nl-1)
dxr   = xrr/nr

!      ------       Generate real particles 

If (iomega.eq.0) then

   x_s     = 0.
   vx_s    = 0.

   mss  = hl*xll/nl
   do i = 1,nl
      mass_s (i) = mss
      hsml_s (i) = dxl*2.4   ! ****2
      itype_s(i) = 1
      x_s  (1,i) = xl + dxl*(i-1)
   enddo 

   mss  = hr*xrr/nr
   do i = 1,nr
      mass_s (nl+i) = mss
      hsml_s (nl+i) = dxr*2.5   ! ***** changed 7 June 2013 from 2.0
      itype_s(nl+i) = 1
      x_s  (1,nl+i) = xc + dxr*(i)
   enddo 

!   pleft = (const(1)/2)*hl*hl 
!   pright= (const(1)/2)*hr*hr 
   do i = 1,npoin_s
      if (x_s(1,i).le.(xc + 1.e-6)) then
          u_s  (i)  = 0.0
          rho_s(i)  = hl
 !         p_s  (i)  = pleft
      else  
          u_s  (i)  = 0.0
          rho_s(i)  = hr
 !         p_s  (i)  = pright
      endif        
   enddo   
   
elseif (iomega.eq.1) then

   if (hr.le.1.e-6) then
       dxlr = (xc-xl)/(npoin_s - 1)
   else
       dxlr = (xr-xl)/(npoin_s - 1)
   endif
   do i = 1,npoin_s
      itype_s  (i) = 1
      hsml_s   (i) = 2.*dxlr
      x_s    (1,i) = xl + (i-1)*dxlr
      if (x_s(1,i).le.xc) then
          rho_s(i) = hl
      else
          rho_s(i) = hr
      endif
      mass_s(i)  = dxlr*rho_s(i)
      u_s   (i)  = 0.0
   enddo

else
   write(*,*) ' iomega is incorrect only 0 or 1 ', iomega
   pause
endif 

vx_s = 0.0

1005     format(a60) 



END subroutine Read_1Dbk_SW

! ----------------------------------------------------------------------   

  subroutine Read_rect_SW (xshiftg,yshiftg) 

! ---------------------------------------------------------------------- 

implicit none  
 
integer(ink) npoinx, npoiny
integer(ink) ipoinx, ipoiny, ipoin
   
real   (irk) x0, y0, x1, y1,  b, d0,  facthsml
real   (irk) xs, ys    
real   (irk) xc, yc,  a, vnxi, vneta, vvxi, vveta, mss 
real   (irk) xi, eta, area_Poin, cgra05, xshiftg,yshiftg  

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    npoinx, npoiny       
write(chk_file,*)    npoinx, npoiny 

read (dat_file,1005) text         
write(chk_file,1005) text        
read (dat_file,*)    x0,y0,x1,y1,b,d0,facthsml        
write(chk_file,*)    x0,y0,x1,y1,b,d0,facthsml

x0    = x0 - xshiftg
x1    = x1 - xshiftg
y0    = y0 - yshiftg
y1    = y1 - yshiftg

xc    = (x0+x1)/2.
yc    = (y0+y1)/2.
a     = ( (x1-xc)**2.+(y1-yc)**2 )**0.5
vnxi  = (x1-xc)/a
vneta = (y1-yc)/a
vvxi  = - vneta
vveta =   vnxi

xs = x0 + b*vneta
ys = y0 - b*vnxi

cgra05  = const(1)/2.0
                                         
npoin_s   = npoinx * npoiny
area_Poin = 4*a*b/npoin_s
mss       = area_Poin*d0 

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

                                                        ! allocate p_s(npoin_s) deleted

ipoin      = 0
do ipoiny = 1, npoiny
do ipoinx = 1, npoinx

   ipoin   = ipoin  + 1
   
   x_s   (1,ipoin) = xs + 2.*a*vnxi *(ipoinx-1)/(npoinx-1) + 2.*b*vvxi *(ipoiny-1)/(npoiny-1)
   x_s   (2,ipoin) = ys + 2.*a*vneta*(ipoinx-1)/(npoinx-1) + 2.*b*vveta*(ipoiny-1)/(npoiny-1)
   rho_s   (ipoin) = d0    
   vx_s  (:,ipoin) = 0.0
   itype_s (ipoin) = 1
   mass_s  (ipoin)  = mss
   hsml_s  (ipoin) = (area_Poin**0.5) * facthsml
   u_s     (ipoin) = 0.0     
   
enddo
enddo

1005     format(a60) 

END subroutine Read_rect_SW    

!------------------------------------------

subroutine Read_3Dcloud (xshift, yshift)

!------------------------------------------

!      ------  Gets Topographic data for SW 

implicit none
 
integer(ink)  In_3Dcloud, OUT_2Dcloud, Chk_3Dcloud
integer(ink)  np3D , np2D 

real   (irk)  xshift, yshift  
real   (irk), allocatable:: mass3D(:), x3D(:,:), v3D(:,:)     
real   (irk), allocatable:: h2D   (:), x2D(:,:), v2D(:,:), Z2D(:), x2D_tmp(:,:)
real   (irk), allocatable::          aux2D(:,:)

real   (irk)  xmin2D,  xmax2D,  delt2Dx, delt2D, facthsml2D 
real   (irk)  ymin2D,  ymax2D,  delt2Dy
integer(ink)  npoi2Dx, npoi2Dy 

integer(ink)  ip3D, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i
real   (irk)  r3D, dist2D, xx, yy, Zp, yg, xg, cgra05, dArea
real   (irk)  distx2, disty2, dist

character(16) problem_name            ! name of PFC3D file
character(60) text                    ! general purpose char string

In_3Dcloud    = 86
Chk_3Dcloud   = 88
OUT_2Dcloud   = 90

! -----------------------  changed ---------------

1005 format (a60)
1006 format (a16)

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 
len1 = len_trim(problem_name)

open (  In_3Dcloud ,  file = problem_name(1:len1)//'.cloud.in3D'      )
open (  Chk_3Dcloud,  file = problem_name(1:len1)//'.cloud.chk3D'      )
open (  OUT_2Dcloud,  file = problem_name(1:len1)//'.cloud.out2D'     )


!      ------  Read number of points in cloud and deltx 

read (In_3Dcloud ,1) text
write(Chk_3Dcloud,1) text
read (In_3Dcloud ,*) np3D, delt2D, facthsml2D
write(Chk_3Dcloud,*) np3D, delt2D, facthsml2D

allocate ( mass3D(np3D),x3D(3,np3D) , v3D(3,np3D))

if(allocated(const)) then    !!  MPJuly 2019 from Chuan
    cgra05  = const(1)/2.0
else
    cgra05  =9.81/2.0
endif

read (In_3Dcloud ,1) text
write(Chk_3Dcloud,1) text
1 format(a60)

xmin2D    =  1.e9
ymin2D    =  1.e9
xmax2D    = -1.e9
ymax2D    = -1.e9

DO ip3D =1, np3D
   read (In_3Dcloud ,*) r3d, mass3D(ip3D), (x3D(idimn,ip3D), idimn=1,3), (v3D(idimn,ip3D), idimn=1,3)  
   write(Chk_3Dcloud,*) r3d, mass3D(ip3D), (x3D(idimn,ip3D), idimn=1,3), (v3D(idimn,ip3D), idimn=1,3)
   
   x3D(1,ip3D) = x3D(1,ip3D) - xshift                   !      ------  normalize coords 
   x3D(2,ip3D) = x3D(2,ip3D) - yshift  
   if (x3D(1,ip3D).le.xmin2D) xmin2D = x3D(1,ip3D)
   if (x3D(1,ip3D).ge.xmax2D) xmax2D = x3D(1,ip3D)
   if (x3D(2,ip3D).le.ymin2D) ymin2D = x3D(2,ip3D)
   if (x3D(2,ip3D).ge.ymax2D) ymax2D = x3D(2,ip3D)
ENDDO

delt2Dx = delt2D
delt2Dy = delt2D
npoi2Dx = (xmax2D-xmin2D)/delt2Dx + 1
npoi2Dy = (ymax2D-ymin2D)/delt2Dy + 1 
delt2Dx = (xmax2D-xmin2D)/(npoi2Dx-1)
delt2Dy = (ymax2D-ymin2D)/(npoi2Dy-1)
 
np2D    = npoi2Dx * npoi2Dy

allocate ( aux2D(6,np2D), h2D   (np2D), x2D(2,np2D), v2D(2,np2D), Z2D(np2D) )
allocate ( x2D_tmp (2,np2D) )
aux2D   = 0.0;  h2D = 0.0; x2D = 0.0; v2D = 0.0; x2D_tmp = 0.0

dist2D = (delt2Dx*delt2Dx + delt2Dy*delt2Dy)**0.5

!  Get coordinates of 2D structured mesh npoi2dx x npoi2dy (divs)

ip2D  = 0
DO ipoi2Dy = 1, npoi2Dy 
    yy     = ymin2D + (ipoi2Dy-1)*delt2Dy
    DO ipoi2Dx = 1, npoi2Dx 
       ip2D    = ip2D + 1
       xx      = xmin2D + (ipoi2Dx-1)*delt2Dx
       x2D(1,ip2D) = xx
       x2D(2,ip2D) = yy
       call Get_Z_Topo (ip2D,xx, yy, Zp)
       Z2D(  ip2D) = Zp
    ENDDO
ENDDO


DO ip3D = 1, np3D
   xx = x3D(1,ip3D)
   yy = x3D(2,ip3D)
   call Get_Z_Topo (ip3D,xx, yy, Zp)
   ipoi2Dx    = (xx-xmin2D)/delt2Dx + 1
   ipoi2Dy    = (yy-ymin2D)/delt2Dy + 1
   inx0       = max(1,ipoi2Dx-1)
   iny0       = max(1,ipoi2Dy-1)
   inx1       = min(npoi2Dx,ipoi2Dx + 2)
   iny1       = min(npoi2Dy,ipoi2Dy + 2)
   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoi2Dx + inx
      distx2  = ( xx - x2D(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - x2D(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = max ( aux2D(1,ipoi2D),  x3D(3,ip3D)-Zp )
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +    v3D(1,ip3D)*mass3D(ip3D)/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) +    v3D(2,ip3D)*mass3D(ip3D)/dist
      aux2D (4,ipoi2D) = aux2D (4,ipoi2D) + mass3D  (ip3D)/dist
      aux2D (6,ipoi2D) = aux2D (6,ipoi2D) + 1      
   enddo
   enddo
ENDDO

!      ------  Get values at nodes

ip   = 0
ip2D = 0
DO ipoi2Dy = 1, npoi2Dy
DO ipoi2Dx = 1, npoi2Dx
   ip2D    = ip2D + 1
   if ( aux2D(1,ip2D).ge.0.1) then
        ip          = ip + 1
        v2D    (1,ip) = aux2D (2,ip2D)/aux2D (4,ip2D)
        v2D    (2,ip) = aux2D (3,ip2D)/aux2D (4,ip2D)
        h2D    (  ip) = aux2D (1,ip2D)
        x2D_tmp(:,ip) = x2D(:,ip2D)
                write(OUT_2Dcloud,*) x2D(1,ip2D), x2D(2,ip2D), v2D(1,ip2D), v2D(2,ip2D), h2D(ip2D)        
   endif   
enddo
enddo

!      ------  Transfer variables to soil arrays after allocating them. 
!              Only 2d pts of grid with h>0.1

npoin_s = ip
allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

dArea   = delt2Dx*delt2Dy
do ip   = 1, npoin_s
    x_s (:,ip) = x2D_tmp(:,ip)
    vx_s(:,ip) = v2D(:,ip)
    rho_s (ip) = h2D  (ip)
    mass_s(ip) = dArea*rho_s(ip)
enddo

u_s     = 0.0
itype_s = 1
hsml_s  = delt2D*facthsml2D

                

!      Check output

 write(OUT_2Dcloud,*) ' ======================================='
 write(OUT_2Dcloud,*) '         ip    rho_s  x_s   vx_s        '  
 write(OUT_2Dcloud,*) ' ======================================='
 do i = 1,npoin_s
    write(OUT_2Dcloud,*)  rho_s(i), x_s(1,i),x_s(2,i),vx_s(1,i), vx_s(2,i)
 enddo
 
close (  In_3Dcloud  )
close (  Chk_3Dcloud )
close (  OUT_2Dcloud )
 
deallocate ( mass3D, x3D, v3D )
deallocate ( h2D, x2D, v2D, Z2D, aux2D, x2D_tmp )


 
END subroutine Read_3Dcloud




!------------------------------------------

subroutine Read_Source (xshift, yshift)

!------------------------------------------

!      ------  Gets Topographic data for SW 

implicit none
 
integer(ink)  In_Source, Chk_source
integer(ink)  ips, idimn, len1

real   (irk)  xshift, yshift 
real   (irk)  xmin2D,  xmax2D, ymin2D,  ymax2D
real   (irk)  delt2D, facthsml2D 
real   (irk)  cgra05, dArea
real   (irk)  total_mass

character(16) problem_name            ! name of PFC3D file
character(60) text                    ! general purpose char string

In_Source    = 86
Chk_Source   = 88

print *,'Input the Source file ?'
read  *, problem_name
len1 = len_trim(problem_name)

open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'  )

!      ------  Read  

read (In_Source ,1) text
write(Chk_Source,1) text
read (In_Source ,*) npoin_s, delt2D, facthsml2D
write(Chk_Source,*) npoin_s, delt2D, facthsml2D

allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

                                                        ! allocate p_s(npoin_s) deleted

cgra05  = const(1)/2.0

read (In_Source ,1) text
write(Chk_Source,1) text
1 format(a60)

xmin2D    =  1.e9
ymin2D    =  1.e9
xmax2D    = -1.e9
ymax2D    = -1.e9
dArea   = delt2D * delt2D 

total_mass = 0.0
DO ips =1, npoin_s
   read (In_Source ,*) (x_s(idimn,ips), idimn=1,ndimn), rho_s(ips)  
   write(Chk_Source,*) (x_s(idimn,ips), idimn=1,ndimn), rho_s(ips)

   if (rho_s(ips).lt.h_inf_SW)  rho_s(ips)=h_inf_SW ! added 4 dec 2011
     
   x_s(1,ips) = x_s(1,ips) - xshift                     !      ------  normalize coords 
   x_s(2,ips) = x_s(2,ips) - yshift  
   if (x_s(1,ips).le.xmin2D) xmin2D = x_s(1,ips)
   if (x_s(1,ips).ge.xmax2D) xmax2D = x_s(1,ips)
   if (x_s(2,ips).le.ymin2D) ymin2D = x_s(2,ips)
   if (x_s(2,ips).ge.ymax2D) ymax2D = x_s(2,ips)
   vx_s(:,ips) = 0.0
   mass_s(ips) = dArea  * rho_s(ips)
   total_mass  = total_mass + mass_s(ips)
ENDDO

write(Chk_Source,*) '  TOTAL INITIAL VOLUME = ',total_mass

u_s     = 0.0
itype_s = 1
hsml_s  = delt2D*facthsml2D

close (  In_Source  )
close (  Chk_Source )

 
END subroutine Read_Source



!------------------------------------------

subroutine Read_Source_plus  (xshift, yshift, icunk)

!------------------------------------------

!      ------  Gets Topographic data for SW 

implicit none
 
integer(ink), intent(IN), OPTIONAL :: icunk     ! we use a -ve value to read dx and dy
logical   ic_dxdy                               ! when we provide icunk this var is YES

integer(ink)  In_Source, Chk_Source
integer(ink)  np_Src , np2D 

real   (irk)  xshift, yshift  
real   (irk), allocatable::  x_Src(:,:), h_Src(:)     
real   (irk), allocatable::    x2D(:,:), h2D  (:), x2D_tmp(:,:) 
real   (irk), allocatable::  aux2D(:,:)

real   (irk)  xmin2D,  xmax2D,  delt2Dx, delt2D, facthsml2D 
real   (irk)  ymin2D,  ymax2D,  delt2Dy
integer(ink)  npoi2Dx, npoi2Dy 

integer(ink)  ip_Src, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i
real   (irk)  dist2D, xx, yy, Zp, yg, xg, cgra05, dArea
real   (irk)  distx2, disty2, dist
real   (irk)  total_mass

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

ic_dxdy = PRESENT (icunk)

In_Source    = 86
Chk_Source   = 88

!print *,'Input the source file ?'  ** changed
!read  *, problem_name

! -----------------------  changed ---------------

1005 format (a60)
1006 format (a16)
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 

! -----------------------  changed ---------------


len1 = len_trim(problem_name)

open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'     )

!      ------  Read number of points in cloud and deltx 

read (In_Source ,1) text
write(Chk_Source,1) text
if (ic_dxdy) then
   read (In_Source ,*) np_Src, delt2Dx, delt2Dy, facthsml2D
   write(Chk_Source,*) np_Src, delt2Dx, delt2Dy, facthsml2D
   delt2D = max (delt2Dx, delt2Dy)
else
   read (In_Source ,*) np_Src, delt2D, facthsml2D
   write(Chk_Source,*) np_Src, delt2D, facthsml2D
   delt2Dx = delt2D
   delt2Dy = delt2D
endif 

allocate ( x_Src(ndimn,np_Src) , h_Src(np_Src))

! cgra05  = const(1)/2.0

read (In_Source ,1) text
write(Chk_Source,1) text
1 format(a60)

xmin2D    =  1.e9
ymin2D    =  1.e9
xmax2D    = -1.e9
ymax2D    = -1.e9

DO ip_Src =1, np_Src
   read (In_Source ,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src)  
   write(Chk_Source,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src) 
   x_Src(1,ip_Src) = x_Src(1,ip_Src) - xshift                   !      ------  normalize coords 
   x_Src(2,ip_Src) = x_Src(2,ip_Src) - yshift  
   if (x_Src(1,ip_Src).le.xmin2D) xmin2D = x_Src(1,ip_Src)
   if (x_Src(1,ip_Src).ge.xmax2D) xmax2D = x_Src(1,ip_Src)
   if (x_Src(2,ip_Src).le.ymin2D) ymin2D = x_Src(2,ip_Src)
   if (x_Src(2,ip_Src).ge.ymax2D) ymax2D = x_Src(2,ip_Src)
ENDDO

npoi2Dx = (xmax2D-xmin2D)/delt2Dx + 1
npoi2Dy = (ymax2D-ymin2D)/delt2Dy + 1 
delt2Dx = (xmax2D-xmin2D)/(npoi2Dx-1)
delt2Dy = (ymax2D-ymin2D)/(npoi2Dy-1)
 
np2D    = npoi2Dx * npoi2Dy

allocate ( x2D (ndimn,np2D), h2D(np2D),aux2D(3,np2D)  )
allocate ( x2D_tmp (2,np2D) )
aux2D   = 0.0;  h2D = 0.0; x2D = 0.0 ; x2D_tmp = 0.0

dist2D = (delt2Dx*delt2Dx + delt2Dy*delt2Dy)**0.5

!  Get coordinates of 2D structured mesh npoi2dx x npoi2dy (divs)

ip2D  = 0
DO ipoi2Dy = 1, npoi2Dy 
    yy     = ymin2D + (ipoi2Dy-1)*delt2Dy
    DO ipoi2Dx = 1, npoi2Dx 
       ip2D    = ip2D + 1
       xx      = xmin2D + (ipoi2Dx-1)*delt2Dx
       x2D(1,ip2D) = xx
       x2D(2,ip2D) = yy
    ENDDO
ENDDO

DO ip_Src = 1, np_Src
   xx = x_Src(1,ip_Src)
   yy = x_Src(2,ip_Src)
   ipoi2Dx    = (xx-xmin2D)/delt2Dx + 1
   ipoi2Dy    = (yy-ymin2D)/delt2Dy + 1
   inx0       = max(1,ipoi2Dx-1)
   iny0       = max(1,ipoi2Dy-1)
   inx1       = min(npoi2Dx,ipoi2Dx + 2)
   iny1       = min(npoi2Dy,ipoi2Dy + 2)
   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoi2Dx + inx
      distx2  = ( xx - x2D(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - x2D(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = aux2D (1,ipoi2D) + 1 
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +        1.0/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) + h_Src(ip_Src)/dist      
   enddo
   enddo
ENDDO

!      ------  Get values at nodes

ip   = 0
ip2D = 0
DO ipoi2Dy = 1, npoi2Dy
DO ipoi2Dx = 1, npoi2Dx
   ip2D    = ip2D + 1
   if ( aux2D(2,ip2D).ge.1.e-6) then
        ip        = ip + 1
        h2D  (ip) = aux2D (3,ip2D)/aux2D (2,ip2D)  
        x2D_tmp(:,ip) = x2D(:,ip2D)     
   endif   
enddo
enddo

!      ------  Transfer variables to soil arrays after allocating them. 
!              Only 2d pts of grid with h>0.1

npoin_s = ip
allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

                                                        ! allocate p_s(npoin_s) deleted


total_mass = 0.0
dArea   = delt2Dx*delt2Dy
do ip   = 1, npoin_s
    x_s (:,ip) = x2D_tmp(:,ip)
    vx_s(:,ip) = 0.0
    rho_s (ip) = h2D  (ip)
    if (rho_s(ip).lt.h_inf_SW) rho_s(ip)=h_inf_SW ! added 4 dec 2011
    mass_s(ip) = dArea*rho_s(ip)
    total_mass  = total_mass + mass_s(ip)
enddo

write(Chk_Source,*) '  TOTAL INITIAL VOLUME = ',total_mass


u_s     = 0.0
itype_s = 1
hsml_s  = delt2D*facthsml2D

close (  In_Source  )
close (  Chk_Source )
 
deallocate ( x_Src, h_Src )
deallocate ( h2D, x2D, aux2D,x2D_tmp)

END subroutine Read_Source_plus




!------------------------------------------

subroutine Read_Source_plus_New (icunk_s, xshift, yshift)

!------------------------------------------

!      ------  Gets Topographic data for SW 

implicit none
 
integer(ink)  In_Source, Chk_Source, icunk_s
integer(ink)  np_Src , np_Src0, np2D    ! 
integer(ink)  ndivsx, ndivsy            ! for each grip point we will explore ndivs divs of
                                        ! the resulting mesh
real   (irk)  xshift, yshift  
real   (irk), allocatable::  x_Src  (:,:), h_Src (:)    
real   (irk), allocatable::  x_Src0 (:,:), h_Src0(:)      
real   (irk), allocatable::  x2D    (:,:), h2D   (:), x2D_tmp(:,:) 
real   (irk), allocatable::  aux2D  (:,:)
real   (irk), allocatable::  fct_area (:)

real   (irk)  xmin2D,  xmax2D,  delt2Dx, delt2D, facthsml2D 
real   (irk)  ymin2D,  ymax2D,  delt2Dy
integer(ink)  npoi2Dx, npoi2Dy  !  

integer(ink)  ip_Src, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i
integer(ink)  n1, n2, n3, n4, ipoig, np_s
 
real   (irk)  dist2D, xx, yy, Zp, yg, xg, cgra05, dArea
real   (irk)  distx2, disty2, dist
real   (irk)  dx_s
real   (irk)  total_mass

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

In_Source    = 86
Chk_Source   = 88

! new vars:     dx_s............representative spacing of points (read)
!               ndivsx, ndivsy..for each grip point we will explore ndivs divs of
!                               the resulting mesh 

!               ndivSx = min (2,(dx_S/delt2Dx))
!               ndivSy = min (2,(dx_S/delt2Dy))
!
!               search indexes in the output grid are now:
!                 
!               inx0 = max(1,      ipoi2Dx - ndivSx + 1 )
!               iny0 = max(1,      ipoi2Dy - ndivSy     )
!               inx1 = min(npoi2Dx,ipoi2Dx + ndivSx + 1 )
!               iny1 = min(npoi2Dy,ipoi2Dy + ndivSy     )

1005 format (a60)
1006 format (a60)    !   BUG when using platform MP Nov  2023
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 

len1 = len_trim(problem_name_Global) !   BUG in platform excel file output MP 7 Nov 2023  ************
open (  In_Source ,  file = problem_name_Global(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name_Global(1:len1)//'.pts.chk'     )    

!len1 = len_trim(problem_name)
!open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
!open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'     )

!      ------  Read number of points in cloud and deltx 

read (In_Source ,1) text
write(Chk_Source,1) text

read (In_Source ,*) np_Src, dx_s, delt2Dx, delt2Dy, facthsml2D
write(Chk_Source,*) np_Src, dx_s, delt2Dx, delt2Dy, facthsml2D
delt2D = max (delt2Dx, delt2Dy)

1 format(a60)

IF (np_Src.gt.0 ) then
    
    allocate ( x_Src(ndimn,np_Src) , h_Src(np_Src))

    read (In_Source ,1) text
    write(Chk_Source,1) text

    xmin2D    =  1.e9
    ymin2D    =  1.e9
    xmax2D    = -1.e9
    ymax2D    = -1.e9

    DO ip_Src =1, np_Src
       read (In_Source ,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src)  
       write(Chk_Source,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src) 
       x_Src(1,ip_Src) = x_Src(1,ip_Src) - xshift                   !      ------  normalize coords 
       x_Src(2,ip_Src) = x_Src(2,ip_Src) - yshift  
       if (x_Src(1,ip_Src).le.xmin2D) xmin2D = x_Src(1,ip_Src)
       if (x_Src(1,ip_Src).ge.xmax2D) xmax2D = x_Src(1,ip_Src)
       if (x_Src(2,ip_Src).le.ymin2D) ymin2D = x_Src(2,ip_Src)
       if (x_Src(2,ip_Src).ge.ymax2D) ymax2D = x_Src(2,ip_Src)
    ENDDO

ELSE    !   we have given np_Src = -999, points are stored in coorg, topol(15:)
    
    np_Src0 = npoig
    np_s    = 0

    xmin2D    =  1.e9
    ymin2D    =  1.e9
    xmax2D    = -1.e9
    ymax2D    = -1.e9

    allocate ( x_Src0 (ndimn,np_Src0) , h_Src0(np_Src0))
    x_Src0 = 0.0 ; h_Src0 = 0.0
    Do ipoig = 1,npoig
       if (topol(15,ipoig).gt.1.e-3) then
           np_s = np_s + 1
           x_Src0 (1,np_s) = coorg (1,ipoig) 
           if (ndimn.eq.2) then
               x_Src0 (2,np_s) = coorg (2,ipoig)  
           endif
           h_Src0 (np_s) = topol (15,ipoig)
           if (x_Src0(1,np_s).le.xmin2D) xmin2D = x_Src0(1,np_s)
           if (x_Src0(1,np_s).ge.xmax2D) xmax2D = x_Src0(1,np_s)
           if (x_Src0(2,np_s).le.ymin2D) ymin2D = x_Src0(2,np_s)
           if (x_Src0(2,np_s).ge.ymax2D) ymax2D = x_Src0(2,np_s)
       endif
    enddo

    np_Src = np_s
    allocate ( x_Src(ndimn,np_Src) , h_Src(np_Src))
    x_Src (:,1:np_Src) = x_Src0 (:,1:np_Src)
    h_Src (1:np_Src)   = h_Src0 (1:np_Src)       
    deallocate ( x_Src0 , h_Src0 )
ENDIF

npoi2Dx = (xmax2D-xmin2D)/delt2Dx + 1
npoi2Dy = (ymax2D-ymin2D)/delt2Dy + 1 
delt2Dx = (xmax2D-xmin2D)/(npoi2Dx-1)
delt2Dy = (ymax2D-ymin2D)/(npoi2Dy-1)
 
np2D    = npoi2Dx * npoi2Dy

allocate ( x2D (ndimn,np2D), h2D(np2D),aux2D(4,np2D)  )
allocate ( x2D_tmp (2,np2D) )
allocate ( fct_Area (np2D) )

aux2D   = 0.0;  h2D = 0.0; x2D = 0.0 ; x2D_tmp = 0.0

dist2D = (delt2Dx*delt2Dx + delt2Dy*delt2Dy)**0.5

!  Get coordinates of 2D structured mesh npoi2dx x npoi2dy (divs)

ip2D  = 0
DO ipoi2Dy = 1, npoi2Dy 
    yy     = ymin2D + (ipoi2Dy-1)*delt2Dy
    DO ipoi2Dx = 1, npoi2Dx 
       ip2D    = ip2D + 1
       xx      = xmin2D + (ipoi2Dx-1)*delt2Dx
       x2D(1,ip2D) = xx
       x2D(2,ip2D) = yy
    ENDDO
ENDDO

ndivSx = max (2.,(dx_S/delt2Dx))
ndivSy = max (2.,(dx_S/delt2Dy))

DO ip_Src = 1, np_Src
   xx = x_Src(1,ip_Src)
   yy = x_Src(2,ip_Src)
   ipoi2Dx    = (xx-xmin2D)/delt2Dx + 1
   ipoi2Dy    = (yy-ymin2D)/delt2Dy + 1
   inx0       = max(1,      ipoi2Dx - ndivSx + 1 )
   iny0       = max(1,      ipoi2Dy - ndivSy     )
   inx1       = min(npoi2Dx,ipoi2Dx + ndivSx + 1 )
   iny1       = min(npoi2Dy,ipoi2Dy + ndivSy     )
   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoi2Dx + inx
      distx2  = ( xx - x2D(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - x2D(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = aux2D (1,ipoi2D) + 1 
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +        1.0/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) + h_Src(ip_Src)/dist      
   enddo
   enddo
ENDDO

!      ------  Obtain aux2D(4,ipoi2D) code at nodes of output mesh
!                     4 ... interior node, connected to 4 rects. full dA
!                     3 ... boundary corner, convex, 3 rects..........dA/3
!                     2 ... boundary side, plane,    2 rects          dA/2
!                     1 ... corner node, 1 rect only..                dA/4


DO ipoi2Dy = 1, npoi2Dy-1
DO ipoi2Dx = 1, npoi2Dx-1
   ipoi2D  = (ipoi2Dy-1)*npoi2dx + ipoi2Dx
   n1      = ipoi2D
   n2      = n1 + 1 
   n3      = n1 + npoi2dx  
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

!      ------  Get values at nodes

ip   = 0
ip2D = 0
DO ipoi2Dy = 1, npoi2Dy
DO ipoi2Dx = 1, npoi2Dx
   ip2D    = ip2D + 1
   if ( aux2D(2,ip2D).ge.1.e-6) then
        ip        = ip + 1
        h2D  (ip) = aux2D (3,ip2D)/aux2D (2,ip2D)  
        x2D_tmp(:,ip) = x2D(:,ip2D)
        fct_Area (ip) = aux2D (4,ip2D) ! 4 for inner to 1 at corners
   endif   
enddo
enddo

!      ------  Transfer variables to soil arrays after allocating them. 
!              Only 2d pts of grid with h>0.1

npoin_s = ip
allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )
                                                        ! allocate p_s(npoin_s) deleted

total_mass = 0.0
dArea   = delt2Dx*delt2Dy/4   ! after we will multiply for 4,3,2 or 1
do ip   = 1, npoin_s
    x_s (:,ip) = x2D_tmp(:,ip)
    vx_s(:,ip) = 0.0
    rho_s (ip) = h2D  (ip)
    if (rho_s(ip).lt.h_inf_SW) rho_s(ip)=h_inf_SW ! added 4 dec 2011
    mass_s(ip) = dArea*rho_s(ip)*fct_Area(ip)
    total_mass  = total_mass + mass_s(ip)
enddo

write(Chk_Source,*) '  TOTAL INITIAL VOLUME = ',total_mass


u_s     = 0.0
itype_s = 1
hsml_s  = delt2D*facthsml2D

if (icunk_s.eq.8) then  !      ------  Read expansion and sg¡hift of soil pts 

read (In_Source ,1) text
write(Chk_Source,1) text

read (In_Source ,*) h_expans, mass_expans, (mass_shift (idimn), idimn= 1,ndimn )
write(Chk_Source,*) h_expans, mass_expans, (mass_shift (idimn), idimn= 1,ndimn ) 
    
endif

close (  In_Source  )
close (  Chk_Source )
 
deallocate ( x_Src, h_Src )
deallocate ( h2D, x2D, aux2D,x2D_tmp)

END subroutine Read_Source_plus_New


!------------------------------------------

subroutine Read_Source_plus_New_WATER (icunk_w, xshift, yshift)

!------------------------------------------

!      ------  Gets Topographic data for SW: water for DFs 

implicit none
 
integer(ink)  In_Source, Chk_Source, icunk_w
integer(ink)  np_Src , np2D
integer(ink)  ndivsx, ndivsy            ! for each grip point we will explore ndivs divs of
                                        ! the resulting mesh
real   (irk)  xshift, yshift  
real   (irk), allocatable::  x_Src(:,:), h_Src(:)     
real   (irk), allocatable::    x2D(:,:), h2D  (:), x2D_tmp(:,:) 
real   (irk), allocatable::  aux2D(:,:)
real   (irk), allocatable::  fct_area(:)

real   (irk)  xmin2D,  xmax2D,  delt2Dx, delt2D, facthsml2D 
real   (irk)  ymin2D,  ymax2D,  delt2Dy
integer(ink)  npoi2Dx, npoi2Dy  !  

integer(ink)  ip_Src, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i
integer(ink)  n1, n2, n3, n4
 
real   (irk)  dist2D, xx, yy, Zp, yg, xg, cgra05, dArea
real   (irk)  distx2, disty2, dist
real   (irk)  dx_s
real   (irk)  total_mass

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

In_Source    = 86
Chk_Source   = 88

! new vars:     dx_s............representative spacing of points (read)
!               ndivsx, ndivsy..for each grip point we will explore ndivs divs of
!                               the resulting mesh 

!               ndivSx = min (2,(dx_S/delt2Dx))
!               ndivSy = min (2,(dx_S/delt2Dy))
!
!               search indexes in the output grid are now:
!                 
!               inx0 = max(1,      ipoi2Dx - ndivSx + 1 )
!               iny0 = max(1,      ipoi2Dy - ndivSy     )
!               inx1 = min(npoi2Dx,ipoi2Dx + ndivSx + 1 )
!               iny1 = min(npoi2Dy,ipoi2Dy + ndivSy     )

1005 format (a60)
1006 format (a16)
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 

len1 = len_trim(problem_name)

open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'     )

!      ------  Read number of points in cloud and deltx 

read (In_Source ,1) text
write(Chk_Source,1) text

read (In_Source ,*) np_Src, dx_s, delt2Dx, delt2Dy, facthsml2D
write(Chk_Source,*) np_Src, dx_s, delt2Dx, delt2Dy, facthsml2D
delt2D = max (delt2Dx, delt2Dy)

allocate ( x_Src(ndimn,np_Src) , h_Src(np_Src))

read (In_Source ,1) text
write(Chk_Source,1) text
1 format(a60)

xmin2D    =  1.e9
ymin2D    =  1.e9
xmax2D    = -1.e9
ymax2D    = -1.e9

DO ip_Src =1, np_Src
   read (In_Source ,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src)  
   write(Chk_Source,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src) 
   x_Src(1,ip_Src) = x_Src(1,ip_Src) - xshift                   !      ------  normalize coords 
   x_Src(2,ip_Src) = x_Src(2,ip_Src) - yshift  
   if (x_Src(1,ip_Src).le.xmin2D) xmin2D = x_Src(1,ip_Src)
   if (x_Src(1,ip_Src).ge.xmax2D) xmax2D = x_Src(1,ip_Src)
   if (x_Src(2,ip_Src).le.ymin2D) ymin2D = x_Src(2,ip_Src)
   if (x_Src(2,ip_Src).ge.ymax2D) ymax2D = x_Src(2,ip_Src)
ENDDO

npoi2Dx = (xmax2D-xmin2D)/delt2Dx + 1
npoi2Dy = (ymax2D-ymin2D)/delt2Dy + 1 
delt2Dx = (xmax2D-xmin2D)/(npoi2Dx-1)
delt2Dy = (ymax2D-ymin2D)/(npoi2Dy-1)
 
np2D    = npoi2Dx * npoi2Dy

allocate ( x2D (ndimn,np2D), h2D(np2D),aux2D(4,np2D)  )
allocate ( x2D_tmp (2,np2D) )
allocate ( fct_Area (np2D) )

aux2D   = 0.0;  h2D = 0.0; x2D = 0.0 ; x2D_tmp = 0.0

dist2D = (delt2Dx*delt2Dx + delt2Dy*delt2Dy)**0.5

!  Get coordinates of 2D structured mesh npoi2dx x npoi2dy (divs)

ip2D  = 0
DO ipoi2Dy = 1, npoi2Dy 
    yy     = ymin2D + (ipoi2Dy-1)*delt2Dy
    DO ipoi2Dx = 1, npoi2Dx 
       ip2D    = ip2D + 1
       xx      = xmin2D + (ipoi2Dx-1)*delt2Dx
       x2D(1,ip2D) = xx
       x2D(2,ip2D) = yy
    ENDDO
ENDDO

ndivSx = max (2.,(dx_S/delt2Dx))
ndivSy = max (2.,(dx_S/delt2Dy))

DO ip_Src = 1, np_Src
   xx = x_Src(1,ip_Src)
   yy = x_Src(2,ip_Src)
   ipoi2Dx    = (xx-xmin2D)/delt2Dx + 1
   ipoi2Dy    = (yy-ymin2D)/delt2Dy + 1
   inx0       = max(1,      ipoi2Dx - ndivSx + 1 )
   iny0       = max(1,      ipoi2Dy - ndivSy     )
   inx1       = min(npoi2Dx,ipoi2Dx + ndivSx + 1 )
   iny1       = min(npoi2Dy,ipoi2Dy + ndivSy     )
   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoi2Dx + inx
      distx2  = ( xx - x2D(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - x2D(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = aux2D (1,ipoi2D) + 1 
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +        1.0/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) + h_Src(ip_Src)/dist      
   enddo
   enddo
ENDDO

!      ------  Obtain aux2D(4,ipoi2D) code at nodes of output mesh
!                     4 ... interior node, connected to 4 rects. full dA
!                     3 ... boundary corner, convex, 3 rects..........dA/3
!                     2 ... boundary side, plane,    2 rects          dA/2
!                     1 ... corner node, 1 rect only..                dA/4


DO ipoi2Dy = 1, npoi2Dy-1
DO ipoi2Dx = 1, npoi2Dx-1
   ipoi2D  = (ipoi2Dy-1)*npoi2dx + ipoi2Dx
   n1      = ipoi2D
   n2      = n1 + 1 
   n3      = n1 + npoi2dx  
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

!      ------  Get values at nodes

ip   = 0
ip2D = 0
DO ipoi2Dy = 1, npoi2Dy
DO ipoi2Dx = 1, npoi2Dx
   ip2D    = ip2D + 1
   if ( aux2D(2,ip2D).ge.1.e-6) then
        ip        = ip + 1
        h2D  (ip) = aux2D (3,ip2D)/aux2D (2,ip2D)  
        x2D_tmp(:,ip) = x2D(:,ip2D)
        fct_Area (ip) = aux2D (4,ip2D) ! 4 for inner to 1 at corners
   endif   
enddo
enddo

!      ------  Transfer variables to soil arrays after allocating them. 
!              Only 2d pts of grid with h>0.1

npoin_w = ip
allocate ( x_w(ndimn,npoin_w),  rho_w(npoin_w), mass_w (npoin_w), hsml_w (npoin_w) )
allocate (vx_w(ndimn,npoin_w),  u_w  (npoin_w), itype_w(npoin_w) )
                                                        ! allocate p_s(npoin_s) deleted

total_mass = 0.0
dArea   = delt2Dx*delt2Dy/4   ! after we will multiply for 4,3,2 or 1
do ip   = 1, npoin_w
    x_w (:,ip) = x2D_tmp(:,ip)
    vx_w(:,ip) = 0.0
    rho_w (ip) = h2D  (ip)
    if (rho_w(ip).lt.h_inf_SW) rho_w(ip)=h_inf_SW ! added 4 dec 2011
    mass_w(ip) = dArea*rho_w(ip)*fct_Area(ip)
    total_mass  = total_mass + mass_w(ip)
enddo

write(Chk_Source,*) '  TOTAL INITIAL VOLUME = ',total_mass


u_w     = 0.0
itype_w = 2
hsml_w  = delt2D*facthsml2D

close (  In_Source  )
close (  Chk_Source )
 
deallocate ( x_Src, h_Src )
deallocate ( h2D, x2D, aux2D,x2D_tmp)

END subroutine Read_Source_plus_New_WATER


!-----------------------------------------------------

subroutine Read_Source_plus66  (xshift, yshift, icunk)

!-----------------------------------------------------

!      ------  Gets Topographic data for SW 

implicit none
 
integer(ink), intent(IN), OPTIONAL :: icunk     ! we use a -ve value to read dx and dy
logical   ic_dxdy                               ! when we provide icunk this var is YES

integer(ink)  In_Source, Chk_Source
integer(ink)  np_Src , np2D 

real   (irk)  xshift, yshift  
real   (irk), allocatable::  x_Src(:,:), h_Src(:)     
real   (irk), allocatable::    x2D(:,:), h2D  (:), x2D_tmp(:,:) 
real   (irk), allocatable::  aux2D(:,:)

real   (irk)  xmin2D,  xmax2D,  delt2Dx, delt2D, facthsml2D 
real   (irk)  ymin2D,  ymax2D,  delt2Dy
integer(ink)  npoi2Dx, npoi2Dy 

integer(ink)  ip_Src, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i
real   (irk)  dist2D, xx, yy, Zp, yg, xg, cgra05, dArea
real   (irk)  distx2, disty2, dist
real   (irk)  total_mass
real   (irk)  htmp

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

ic_dxdy = PRESENT (icunk)

In_Source    = 86
Chk_Source   = 88

1005 format (a60)
1006 format (a16)
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 

len1 = len_trim(problem_name)

open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'     )

!      ------  Read number of points in cloud and deltx 

read (In_Source ,1) text
write(Chk_Source,1) text
if (ic_dxdy) then
   read (In_Source ,*) np_Src, delt2Dx, delt2Dy, facthsml2D
   write(Chk_Source,*) np_Src, delt2Dx, delt2Dy, facthsml2D
   delt2D = max (delt2Dx, delt2Dy)
else
   read (In_Source ,*) np_Src, delt2D, facthsml2D
   write(Chk_Source,*) np_Src, delt2D, facthsml2D
   delt2Dx = delt2D
   delt2Dy = delt2D
endif 

allocate ( x_Src(ndimn,np_Src) , h_Src(np_Src))

! cgra05  = const(1)/2.0

read (In_Source ,1) text
write(Chk_Source,1) text
1 format(a60)

xmin2D    =  1.e9
ymin2D    =  1.e9
xmax2D    = -1.e9
ymax2D    = -1.e9

DO ip_Src =1, np_Src
   read (In_Source ,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src)  
   write(Chk_Source,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src) 
   x_Src(1,ip_Src) = x_Src(1,ip_Src) - xshift                   !      ------  normalize coords 
   x_Src(2,ip_Src) = x_Src(2,ip_Src) - yshift  
   if (x_Src(1,ip_Src).le.xmin2D) xmin2D = x_Src(1,ip_Src)
   if (x_Src(1,ip_Src).ge.xmax2D) xmax2D = x_Src(1,ip_Src)
   if (x_Src(2,ip_Src).le.ymin2D) ymin2D = x_Src(2,ip_Src)
   if (x_Src(2,ip_Src).ge.ymax2D) ymax2D = x_Src(2,ip_Src)
ENDDO

npoi2Dx = (xmax2D-xmin2D)/delt2Dx + 1
npoi2Dy = (ymax2D-ymin2D)/delt2Dy + 1 
delt2Dx = (xmax2D-xmin2D)/(npoi2Dx-1)
delt2Dy = (ymax2D-ymin2D)/(npoi2Dy-1)
 
np2D    = npoi2Dx * npoi2Dy

allocate ( x2D (ndimn,np2D), h2D(np2D),aux2D(3,np2D)  )
allocate ( x2D_tmp (2,np2D) )
aux2D   = 0.0;  h2D = 0.0; x2D = 0.0 ; x2D_tmp = 0.0

dist2D = (delt2Dx*delt2Dx + delt2Dy*delt2Dy)**0.5

!  Get coordinates of 2D structured mesh npoi2dx x npoi2dy (divs)

ip2D  = 0
DO ipoi2Dy = 1, npoi2Dy 
    yy     = ymin2D + (ipoi2Dy-1)*delt2Dy
    DO ipoi2Dx = 1, npoi2Dx 
       ip2D    = ip2D + 1
       xx      = xmin2D + (ipoi2Dx-1)*delt2Dx
       x2D(1,ip2D) = xx
       x2D(2,ip2D) = yy
    ENDDO
ENDDO

DO ip_Src = 1, np_Src
   xx = x_Src(1,ip_Src)
   yy = x_Src(2,ip_Src)
   ipoi2Dx    = (xx-xmin2D)/delt2Dx + 1
   ipoi2Dy    = (yy-ymin2D)/delt2Dy + 1
   inx0       = max(1,ipoi2Dx-1)
   iny0       = max(1,ipoi2Dy-1)
   inx1       = min(npoi2Dx,ipoi2Dx + 2)
   iny1       = min(npoi2Dy,ipoi2Dy + 2)
   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoi2Dx + inx
      distx2  = ( xx - x2D(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - x2D(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = aux2D (1,ipoi2D) + 1 
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +        1.0/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) + h_Src(ip_Src)/dist      
   enddo
   enddo
ENDDO

!      ------  Get values at nodes

ip   = 0
ip2D = 0
DO ipoi2Dy = 1, npoi2Dy
DO ipoi2Dx = 1, npoi2Dx
   ip2D    = ip2D + 1
   if ( aux2D(2,ip2D).ge.1.e-6) then
        htmp = aux2D (3,ip2D)/aux2D (2,ip2D)
        if (htmp.gt.h_inf_SW) then
           ip            = ip + 1
           h2D (ip)      = htmp  
           x2D_tmp(:,ip) = x2D(:,ip2D)
        endif   
   endif   
enddo
enddo

!      ------  Transfer variables to soil arrays after allocating them. 
!               

npoin_s = ip
allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )

total_mass = 0.0
dArea   = delt2Dx*delt2Dy
do ip   = 1, npoin_s
    x_s (:,ip) = x2D_tmp(:,ip)
    vx_s(:,ip) = 0.0
    rho_s (ip) = h2D  (ip) 
    mass_s(ip) = dArea*rho_s(ip)
    total_mass  = total_mass + mass_s(ip)
enddo

write(Chk_Source,*) '  TOTAL INITIAL VOLUME = ',total_mass

u_s     = 0.0
itype_s = 1
hsml_s  = delt2D*facthsml2D

close (  In_Source  )
close (  Chk_Source )
 
deallocate ( x_Src, h_Src )
deallocate ( h2D, x2D, aux2D,x2D_tmp)

END subroutine Read_Source_plus66


!------------------------------------------

subroutine Read_Source_plus_New66 (xshift, yshift)

!------------------------------------------

!      ------  Gets Topographic data for SW 

implicit none
 
integer(ink)  In_Source, Chk_Source
integer(ink)  np_Src , np2D
integer(ink)  ndivsx, ndivsy            ! for each grip point we will explore ndivs divs of
                                        ! the resulting mesh
real   (irk)  xshift, yshift  
real   (irk), allocatable::  x_Src(:,:), h_Src(:)     
real   (irk), allocatable::    x2D(:,:), h2D  (:), x2D_tmp(:,:) 
real   (irk), allocatable::  aux2D(:,:)
real   (irk), allocatable::  fct_area(:)

real   (irk)  xmin2D,  xmax2D,  delt2Dx, delt2D, facthsml2D 
real   (irk)  ymin2D,  ymax2D,  delt2Dy
integer(ink)  npoi2Dx, npoi2Dy  !  

integer(ink)  ip_Src, idimn, ip2D, ipoi2Dx, ipoi2Dy, ipoi2D   
integer(ink)  len1, inx0, iny0, inx1, iny1, inx, iny, ip, i
integer(ink)  n1, n2, n3, n4
 
real   (irk)  dist2D, xx, yy, Zp, yg, xg, cgra05, dArea
real   (irk)  distx2, disty2, dist
real   (irk)  dx_s
real   (irk)  total_mass
real   (irk)  htmp

character(16) problem_name            ! name of source file
character(60) text                    ! general purpose char string

In_Source    = 86
Chk_Source   = 88

! new vars:     dx_s............representative spacing of points (read)
!               ndivsx, ndivsy..for each grip point we will explore ndivs divs of
!                               the resulting mesh 

!               ndivSx = min (2,(dx_S/delt2Dx))
!               ndivSy = min (2,(dx_S/delt2Dy))
!
!               search indexes in the output grid are now:
!                 
!               inx0 = max(1,      ipoi2Dx - ndivSx + 1 )
!               iny0 = max(1,      ipoi2Dy - ndivSy     )
!               inx1 = min(npoi2Dx,ipoi2Dx + ndivSx + 1 )
!               iny1 = min(npoi2Dy,ipoi2Dy + ndivSy     )

1005 format (a60)
1006 format (a16)
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,1006) problem_name
write(chk_file,1005) problem_name 

len1 = len_trim(problem_name)

open (  In_Source ,  file = problem_name(1:len1)//'.pts'      )
open (  Chk_Source,  file = problem_name(1:len1)//'.pts.chk'     )

!      ------  Read number of points in cloud and deltx 

read (In_Source ,1) text
write(Chk_Source,1) text

read (In_Source ,*) np_Src, dx_s, delt2Dx, delt2Dy, facthsml2D
write(Chk_Source,*) np_Src, dx_s, delt2Dx, delt2Dy, facthsml2D
delt2D = max (delt2Dx, delt2Dy)

allocate ( x_Src(ndimn,np_Src) , h_Src(np_Src))

read (In_Source ,1) text
write(Chk_Source,1) text
1 format(a60)

xmin2D    =  1.e9
ymin2D    =  1.e9
xmax2D    = -1.e9
ymax2D    = -1.e9

DO ip_Src =1, np_Src
   read (In_Source ,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src)  
   write(Chk_Source,*) (x_Src(idimn,ip_Src), idimn=1,ndimn), h_Src(ip_Src) 
   x_Src(1,ip_Src) = x_Src(1,ip_Src) - xshift                   !      ------  normalize coords 
   x_Src(2,ip_Src) = x_Src(2,ip_Src) - yshift  
   if (x_Src(1,ip_Src).le.xmin2D) xmin2D = x_Src(1,ip_Src)
   if (x_Src(1,ip_Src).ge.xmax2D) xmax2D = x_Src(1,ip_Src)
   if (x_Src(2,ip_Src).le.ymin2D) ymin2D = x_Src(2,ip_Src)
   if (x_Src(2,ip_Src).ge.ymax2D) ymax2D = x_Src(2,ip_Src)
ENDDO

npoi2Dx = (xmax2D-xmin2D)/delt2Dx + 1
npoi2Dy = (ymax2D-ymin2D)/delt2Dy + 1 
delt2Dx = (xmax2D-xmin2D)/(npoi2Dx-1)
delt2Dy = (ymax2D-ymin2D)/(npoi2Dy-1)
 
np2D    = npoi2Dx * npoi2Dy

allocate ( x2D (ndimn,np2D), h2D(np2D),aux2D(4,np2D)  )
allocate ( x2D_tmp (2,np2D) )
allocate ( fct_Area (np2D) )

aux2D   = 0.0;  h2D = 0.0; x2D = 0.0 ; x2D_tmp = 0.0

dist2D = (delt2Dx*delt2Dx + delt2Dy*delt2Dy)**0.5

!  Get coordinates of 2D structured mesh npoi2dx x npoi2dy (divs)

ip2D  = 0
DO ipoi2Dy = 1, npoi2Dy 
    yy     = ymin2D + (ipoi2Dy-1)*delt2Dy
    DO ipoi2Dx = 1, npoi2Dx 
       ip2D    = ip2D + 1
       xx      = xmin2D + (ipoi2Dx-1)*delt2Dx
       x2D(1,ip2D) = xx
       x2D(2,ip2D) = yy
    ENDDO
ENDDO

ndivSx = max (2.,(dx_S/delt2Dx))
ndivSy = max (2.,(dx_S/delt2Dy))

DO ip_Src = 1, np_Src
   xx = x_Src(1,ip_Src)
   yy = x_Src(2,ip_Src)
   ipoi2Dx    = (xx-xmin2D)/delt2Dx + 1
   ipoi2Dy    = (yy-ymin2D)/delt2Dy + 1
   inx0       = max(1,      ipoi2Dx - ndivSx + 1 )
   iny0       = max(1,      ipoi2Dy - ndivSy     )
   inx1       = min(npoi2Dx,ipoi2Dx + ndivSx + 1 )
   iny1       = min(npoi2Dy,ipoi2Dy + ndivSy     )
   DO inx     = inx0,inx1
   DO iny     = iny0,iny1
      ipoi2D  = (iny-1)*npoi2Dx + inx
      distx2  = ( xx - x2D(1,ipoi2D) ) 
      distx2  = distx2*distx2
      disty2  = ( yy - x2D(2,ipoi2D) )
      disty2  = disty2*disty2
      dist    = max( 1.e-6, (distx2+disty2)**0.5 )
      dist    = dist/dist2D
      aux2D (1,ipoi2D) = aux2D (1,ipoi2D) + 1 
      aux2D (2,ipoi2D) = aux2D (2,ipoi2D) +        1.0/dist
      aux2D (3,ipoi2D) = aux2D (3,ipoi2D) + h_Src(ip_Src)/dist      
   enddo
   enddo
ENDDO

!      ------  Obtain aux2D(4,ipoi2D) code at nodes of output mesh
!                     4 ... interior node, connected to 4 rects. full dA
!                     3 ... boundary corner, convex, 3 rects..........dA/3
!                     2 ... boundary side, plane,    2 rects          dA/2
!                     1 ... corner node, 1 rect only..                dA/4


DO ipoi2Dy = 1, npoi2Dy-1
DO ipoi2Dx = 1, npoi2Dx-1
   ipoi2D  = (ipoi2Dy-1)*npoi2dx + ipoi2Dx
   n1      = ipoi2D
   n2      = n1 + 1 
   n3      = n1 + npoi2dx  
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

!      ------  Get values at nodes

ip   = 0
ip2D = 0
DO ipoi2Dy = 1, npoi2Dy
DO ipoi2Dx = 1, npoi2Dx
   ip2D    = ip2D + 1
   if ( aux2D(2,ip2D).ge.1.e-6) then
        htmp = aux2D (3,ip2D)/aux2D (2,ip2D)
        if (htmp.gt.h_inf_SW) then
           ip            = ip + 1
           h2D (ip)      = htmp  
           x2D_tmp(:,ip) = x2D(:,ip2D)
        endif   
   endif   
enddo
enddo

!      ------  Transfer variables to soil arrays after allocating them. 
!              Only 2d pts of grid with h>0.1

npoin_s = ip
allocate ( x_s(ndimn,npoin_s),  rho_s(npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
allocate (vx_s(ndimn,npoin_s),  u_s  (npoin_s), itype_s(npoin_s) )
                                                        ! allocate p_s(npoin_s) deleted

total_mass = 0.0
dArea   = delt2Dx*delt2Dy/4   ! after we will multiply for 4,3,2 or 1
do ip   = 1, npoin_s
    x_s (:,ip) = x2D_tmp(:,ip)
    vx_s(:,ip) = 0.0
    rho_s (ip) = h2D  (ip) 
    mass_s(ip) = dArea*rho_s(ip)*fct_Area(ip)
    total_mass  = total_mass + mass_s(ip)
enddo

write(Chk_Source,*) '  TOTAL INITIAL VOLUME = ',total_mass


u_s     = 0.0
itype_s = 1
hsml_s  = delt2D*facthsml2D

close (  In_Source  )
close (  Chk_Source )
 
deallocate ( x_Src, h_Src )
deallocate ( h2D, x2D, aux2D,x2D_tmp)

END subroutine Read_Source_plus_New66


! ----------------------------------------------------------------------   

      subroutine Read_reservoir_SW 

! ----------------------------------------------------------------------   

!     This routine gets reservoir water nodes info from topo routine... 
!     Reservoir is defined by a dam (we give x0,y0 and xf,yf and Z of water surface 
!     nx and ny_dam are components of vector normal to dam, pointing inwards.
!     to check that a point is in the reservoir, its Z topo should be below water surface Z, Z_swl, 
!     and the product (x-x0)*v >0

implicit none  

integer(ink) ngx, ngy, igx, igy 
integer(ink) ipoin_s

real   (irk) Z_swl, dgx, dgy, facthsml, lgx, lgy, xx, yy, RR 
real   (irk) Zp, dArea
real   (irk) x0_dam, y0_dam, xf_dam, yf_dam, nx_dam, ny_dam             !   initial and final point of dam, and normal vector
real   (irk) xmin_rsv, xmax_rsv, ymin_rsv, ymax_rsv     

real   (irk), allocatable:: xstmp(:,:), rhostmp(:)  ! tmp vars          !   the reservoir is within rectangle defined by xmin xmax, ymin and ymax_dam

read (dat_file,1005) text 
write(chk_file,1005) text
read (dat_file,*)    xmin_rsv, xmax_rsv, ymin_rsv, ymax_rsv, dgx, dgy  ! dgx is deltax for spacing nodes in X 
write(chk_file,*)    xmin_rsv, xmax_rsv, ymin_rsv, ymax_rsv, dgx, dgy  
read (dat_file,1005) text 
write(chk_file,1005) text
read (dat_file,*)    Z_swl, x0_dam, y0_dam, xf_dam, yf_dam  
write(chk_file,*)    Z_swl, x0_dam, y0_dam, xf_dam, yf_dam 

1005     format(a60) 

if (ndimn.eq.1) then
    write(*,*) ' routine read_reservoir only for 2D problems '
    write(*,*) ' hit any key '
    read (*,*) text
    STOP
endif

nx_dam = -(yf_dam-y0_dam)       ! this is the normal pointing inwards reservoir
ny_dam =   xf_dam-x0_dam

lgx   = (xmax_rsv-xmin_rsv-dgx)
lgy   = (ymax_rsv-ymin_rsv-dgy)

ngx   = (lgx/dgx) + 1
ngy   = (lgy/dgy) + 1

dgx   = lgx/(ngx-1)
dgy   = lgy/(ngy-1)

darea = dgx*dgy

npoin_s = ngx*ngy
allocate (xstmp(ndimn,npoin_s), rhostmp(npoin_s) ) 
                                                        ! allocate p_s(npoin_s) deleted
ipoin_s = 0
DO igy  = 1,ngy
DO igx  = 1,ngx

   xx     = xmin_rsv + (igx-1)*dgx + dgx/2.     !      ------  Avoid points running out of the domain
   yy     = ymin_rsv + (igy-1)*dgy + dgy/2.
   
   RR = (xx-x0_dam)*nx_dam + (yy-y0_dam)*ny_dam
   
   facthsml = 1.1
   call Get_Z_Topo (ipoin_s,xx,yy,Zp)
   
!          ipoin keeps the number of points of this grid element which belongs to computational domain
   
   if (Zp.lt.Z_swl.and.RR.ge.0.0.and.Zp.gt.-1.e6) then
      ipoin_s             = ipoin_s + 1
      xstmp   (1,ipoin_s) = xx
      xstmp   (2,ipoin_s) = yy
      rhostmp   (ipoin_s) = Z_swl - Zp
   endif
   
ENDDO
ENDDO   

npoin_s = ipoin_s

allocate (x_s(ndimn,npoin_s), vx_s(ndimn,npoin_s), rho_s (npoin_s), mass_s(npoin_s) )
allocate ( hsml_s  (npoin_s), itype_s   (npoin_s), u_s   (npoin_s) ) 

DO ipoin_s = 1, npoin_s

      x_s  (1,ipoin_s) = xstmp   (1,ipoin_s)
      x_s  (2,ipoin_s) = xstmp   (2,ipoin_s)
      rho_s  (ipoin_s) = rhostmp   (ipoin_s)
      mass_s (ipoin_s) = rho_s (ipoin_s)*darea
      hsml_s (ipoin_s) = facthsml * amax1(dgx,dgy)
      itype_s(ipoin_s) = 1
      vx_s (:,ipoin_s) = 0.0
      u_s    (ipoin_s) = 0.0
      
enddo

deallocate (xstmp, rhostmp)

 
END SUBROUTINE  Read_reservoir_SW


! ----------------------------------------------------------------------   

      subroutine Get_vps_points 
      
! ----------------------------------------------------------------------   

implicit none

integer(ink) icunk_vps  

!      ------   

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    icunk_vps   
write(chk_file,*)    icunk_vps 

1005     format(a60)

if (icunk_vps.eq.1)  then 
    call Read_vps_points_SW
elseif (icunk_vps.eq.2) then
    call Read_Segments_vps
else
    write(*,*) ' incorrect icunk_s option ', icunk_vps
endif
 
END SUBROUTINE  Get_vps_points 



! ----------------------------------------------------------------------   

      subroutine Read_vps_points_SW

! ---------------------------------------------------------------------- 

!     This subroutine is used to read virtual particles 

implicit none  

integer(ink)  ngroups_vp, igvp, Np_vp, ip_vp, ivp0, idimn, ivirt 
real   (irk)  Area_vp, Height_vp, hsml_vp, mass_vp,  twopi 

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    ngroups_vp, npoin_vps
write(chk_file,*)    ngroups_vp, npoin_vps

allocate ( x_vps(ndimn,npoin_s),  rho_vps(npoin_s), mass_vps(npoin_s), hsml_vps (npoin_s) )
allocate (vx_vps(ndimn,npoin_s),  p_vps  (npoin_s), u_vps   (npoin_s), itype_vps(npoin_s) )

twopi = 4.0*atan(1.0)
ivp0  = 0

DO igvp = 1,ngroups_vp
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    Np_vp, Area_vp, Height_vp, hsml_vp 
   write(chk_file,*)    Np_vp, Area_vp, Height_vp, hsml_vp 
   mass_vp = Area_vp*Height_vp/Np_vp
   read (dat_file,1005) text
   write(chk_file,1005) text
   if ((ivp0 + Np_vp).gt.npoin_vps ) Then
      write(*,*) ' Read_vps_points_SW: check number of vps. total exceeds npoin_vps '
      PAUSE
      STOP
   endif
   DO ip_vp  = 1, Np_vp
      ivirt  = ivp0 + ip_vp
      read (dat_file,*) (x_vps (idimn,ivirt),idimn=1,ndimn) 
      write(chk_file,*) (x_vps (idimn,ivirt),idimn=1,ndimn) 
      rho_vps    (ivirt) = Height_vp
      mass_vps   (ivirt) = mass_vp
!      p_vps      (ivirt) = (const(1)/2.)*Height_vp*Height_vp
      hsml_vps   (ivirt) = hsml_vp
      u_vps      (ivirt) = 0.0
      itype_vps  (ivirt) = -3
   enddo
   ivp0 = ivp0 + Np_vp
Enddo

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    vp_r0, vp_D, vp_n1, vp_n2
write(chk_file,*)    vp_r0, vp_D, vp_n1, vp_n2
   
1005     format(a60)  

END subroutine Read_vps_points_SW



!----------------------------------------------------------------------------------

Subroutine Read_Segments_vps

!----------------------------------------------------------------------------------

implicit none  

integer(ink) Nseg, iseg, ndiv_arc,  idiv_arc,  idimn, iv_ps
integer(ink) ipoin, idest , nvps_max, nvps_temp

real   (irk) ds_arc, d_arc, d_arc2, h0_arc, hf_arc, hh, cgra05  

real   (irk), allocatable:: x_temp(:,:), rho_temp(:), hsml_temp(:) 
real   (irk), allocatable:: x0_arc(:),   xf_arc(:),    vt_arc(:) 

nvps_max = npoig

allocate ( x_temp(ndimn,nvps_max), rho_temp(nvps_max), hsml_temp(nvps_max) )
allocate ( x0_arc(ndimn), xf_arc(ndimn), vt_arc(ndimn) )

x_temp = 0; rho_temp=0.; hsml_temp=0.; x0_arc =0.; xf_arc =0.; vt_arc =0.

nvps_temp = 0

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    Nseg    
write(chk_file,*)    Nseg

   
DO iseg = 1, Nseg               !      ------  Begin loop around Rectilinear segments (a-clockwise)

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)    !  --  Initial and final points
   write(chk_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)    
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    ds_arc                                                  !  --  Sep between points, area associated
   write(chk_file,*)    ds_arc  
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    h0_arc, hf_arc                                                  !  --   heights at initial and end points  
   write(chk_file,*)    h0_arc, hf_arc 
   
   vt_arc    = (xf_arc - x0_arc)
   d_arc2    = ( vt_arc(1)*vt_arc(1) + vt_arc(2)*vt_arc(2) )**0.5
   ndiv_arc  = d_arc2/ds_arc
   ds_arc    = d_arc2/ndiv_arc
   vt_arc    = vt_arc / d_arc2 
   
   DO idiv_arc = 1, ndiv_arc + 1
      idest    = nvps_temp + idiv_arc
      if (nvps_temp + idiv_arc.GT.nvps_max) then
         write(*,*) ' in routine read  vps max nr of pts exceeded...modify code to 2*npoig '
         write(*,*) nvps_max, nvps_temp+idiv_arc
         stop
      endif            
      x_temp (:,idest) = x0_arc(:) + (ds_arc/2.)*vt_arc(:) + (idiv_arc-1)*ds_arc*vt_arc(:)
      rho_temp (idest) = h0_arc + (idiv_arc -1)*(hf_arc - h0_arc)/ndiv_arc
      hsml_temp(idest) = ds_arc
   Enddo
      
   nvps_temp  = nvps_temp + ndiv_arc +  1
      
ENDDO  

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)    vp_r0, vp_D, vp_n1, vp_n2   
write(chk_file,*)    vp_r0, vp_D, vp_n1, vp_n2
   
npoin_vps = nvps_temp

allocate ( x_vps(ndimn,npoin_vps),  rho_vps(npoin_vps), mass_vps(npoin_vps), hsml_vps (npoin_vps) )
allocate (vx_vps(ndimn,npoin_vps),  p_vps  (npoin_vps), u_vps   (npoin_vps), itype_vps(npoin_vps) )

cgra05  = 9.8/2.0   

DO iv_ps = 1, npoin_vps
   x_vps  (:,iv_ps) = x_temp(:,iv_ps)
   vx_vps (:,iv_ps) = 0.0
   rho_vps  (iv_ps) = rho_temp(iv_ps)
   p_vps    (iv_ps) = cgra05*rho_vps(iv_ps)*rho_vps(iv_ps)
   mass_vps (iv_ps) = hsml_temp(iv_ps)*hsml_temp(iv_ps)*rho_vps(iv_ps)
   hsml_vps (iv_ps) = 2.* hsml_temp(iv_ps)
   u_vps    (iv_ps) = 0.0
   itype_vps(iv_ps) = -3
ENDDO

deallocate( x_temp, rho_temp, hsml_temp,x0_arc, xf_arc, vt_arc ) 

1005 format (a60)

END Subroutine Read_Segments_vps
  

!----------------------------------------------------------------------------------

     Subroutine Setup_Global_Arrays_SW   

!----------------------------------------------------------------------------------

!      ------  Allocate global arrays, pass local to global, deallocate, initialize

    implicit none

    integer(ink) itotv, ipoin, idimn,iiabs, iivn0
    integer(ink) it1, it2, it3
     
    if (.NOT.allocated (x00) ) then
      allocate (     x00(ndimn,ntotal),    x0 (ndimn,ntotal),     x (ndimn,ntotal),   vx0 (ndimn,ntotal),   vx (ndimn,ntotal) )
      allocate (     Z00(ntotal), rho00(ntotal) )
      allocate (     dx (ndimn,ntotal),   dvx (ndimn,ntotal),    av (ndimn,ntotal) )
      allocate ( indvxdt(ndimn,ntotal),exdvxdt(ndimn,ntotal),ardvxdt(ndimn,ntotal) )
      allocate ( mass   (ntotal), rho(ntotal), p(ntotal), u(ntotal), hsml(ntotal), c(ntotal) )
      allocate ( mass0  (ntotal) )
      allocate ( u0     (ntotal), rho0(ntotal),  du(ntotal),drho(ntotal) )
      allocate ( itype(ntotal), countiac(ntotal) )
      allocate (If_Out_Domain(ntotal))
      allocate (Xmin_Domain(ndimn))
      allocate (Xmax_Domain(ndimn))
    endif 

!    allocate ( mass   (ntotal), rho(ntotal), p(ntotal), u(ntotal), hsml(ntotal), c(ntotal), s(ntotal), e(ntotal) ) 
!    allocate ( ahdudt (ntotal),avdudt(ntotal),eta(ntotal) ) 
!    allocate ( u0     (ntotal), rho0(ntotal),  du(ntotal),drho(ntotal), ds  (ntotal), t(ntotal), tdsdt(ntotal)   ) 
    
    
    if (ic_ws_Interact.eq.2) then  ! for DFs type
        if (.NOT.allocated (h_DF) ) then
           allocate (h_DF(ntotal), n_DF(ntotal), h_sw (ntotal), v_sw (ndimn,ntotal) )
        endif
        h_DF  = 0.0
        n_DF  = 0.0
        h_sw  = 0.0
        v_sw  = 0.0
        
        if (icpwp.GE.1) then                !    Saeed 31st May  MP may 2020 cases
           if (.NOT.allocated (auxS1)) then
              allocate ( auxS1(ntotal)    , auxS10(ntotal), dauxS1(ntotal))
              allocate ( Pwp_avgd (ntotal))
           endif
           auxS1    = 0.0; auxS10 = 0.0; dauxS1 = 0.0
           Pwp_avgd = 0.0
        endif  

    endif
    
    mass0 = 0.0
    
    if_out_domain = 0
    
    npoin0 = 0
        
!      ------  From local soil and water arrays to global array, then deallocate
    
    itotv = 0
    if (npoin_s.gt.0) then                              !      ------  soil comes first
       do ipoin = 1, npoin_s
          itotv = itotv + 1
          do idimn = 1, ndimn
             x  (idimn,itotv) = x_s (idimn,ipoin)
             vx (idimn,itotv) = vx_s(idimn,ipoin)
          enddo
          itype (itotv) = itype_s(ipoin)
          rho   (itotv) = rho_s  (ipoin)
          mass  (itotv) = mass_s (ipoin)
          mass0 (itotv) = mass_s (ipoin)
          hsml  (itotv) = hsml_s (ipoin)
          u     (itotv) = u_s    (ipoin)
          p     (itotv) = 0.0
!          p     (itotv) = p_s    (ipoin)       !  *** take care. p_s  has not been defined
       enddo
       deallocate (x_s, vx_s, itype_s, rho_s, mass_s, hsml_s, u_s )
       if (ic_histogram.eq.1) then                      !      ------  if histogram, all particles are out at t=0
          allocate (if_out_domain_LAST(ntotal))
          if_out_domain     (1:np_bakG_Total)         = 0
          if_out_domain     (np_bakG_Total+1:npoin_s) = 1
          if_out_domain_LAST(1:npoin_s) = if_out_domain (1:npoin_s)         
          allocate ( list_avail_nodes(npoin) )
          list_avail_nodes = 0
          do ipoin = 0, npoin-1
             list_avail_nodes(npoin-ipoin) = ipoin + 1
          enddo
          n_avail_nodes = npoin
          do ipoin = 0, np_bakG_Total-1
             list_avail_nodes(npoin-ipoin) = 0    ! because these nodes are in the model always since the beginning             
          enddo
          n_avail_nodes = npoin - np_bakG_Total
          npoin0 = np_bakG_Total
       endif
    endif
    
    if (npoin_w.gt.0) then                              !      ------  then water
        do ipoin    = 1, npoin_w
          itotv    = itotv + 1
          do idimn = 1, ndimn
             x  (idimn,itotv) = x_w (idimn,ipoin)
             vx (idimn,itotv) = vx_w(idimn,ipoin)
          enddo
          itype (itotv) = itype_w(ipoin)
          rho   (itotv) = rho_w  (ipoin)
          mass  (itotv) = mass_w (ipoin)
          hsml  (itotv) = hsml_w (ipoin)
          u     (itotv) = u_w    (ipoin)
          p     (itotv) = 0.0
!          p     (itotv) = p_w    (ipoin)
       enddo  
       deallocate (x_w, vx_w, itype_w, rho_w, mass_w, hsml_w, u_w )   ! , p_w
    endif
    
    if (npoin_vps.gt.0) then                            !      ------  then vps
        do ipoin    = 1, npoin_vps
           itotv    = itotv + 1
           do idimn = 1, ndimn
                 x  (idimn,itotv) = x_vps (idimn,ipoin)
                 vx (idimn,itotv) = vx_vps(idimn,ipoin)
           enddo
           itype (itotv) = itype_vps(ipoin)
           rho   (itotv) = rho_vps  (ipoin)
           mass  (itotv) = mass_vps (ipoin)
           hsml  (itotv) = hsml_vps (ipoin)
           u     (itotv) = u_vps    (ipoin)
           p     (itotv) = p_vps    (ipoin)
        enddo
        deallocate (x_vps, vx_vps, itype_vps, rho_vps, mass_vps, hsml_vps, u_vps, p_vps )
    endif
    
    if (ic_abs.eq.1) then
        DO iiabs = 1, nabs
           x  (:,npoin+nvirt+iiabs) = x_abs (:,iiabs)
           hsml (npoin+nvirt+iiabs) = hsml_abs(iiabs)
           itype(npoin+nvirt+iiabs) = -12
           rho  (npoin+nvirt+iiabs) = 999.
           mass (npoin+nvirt+iiabs) = 999.
           p    (npoin+nvirt+iiabs) = 999.
           u    (npoin+nvirt+iiabs) = 999.
           labs             (iiabs) = npoin + nvirt + iiabs             !   we do not use it...
        enddo
    endif   
    
    if (ic_vn0.GE.1) then    !  mp 20th dec 2018 ic_vn0=2 for seepage trough wall
        
        ic2_vn0 = ic_vn0
        if(.NOT.allocated(if_correct) ) allocate ( if_correct(npoin) )
        if_correct = 0
        it1 = nvn0
        DO iivn0 = 1, nvn0
           x  (:,npoin+nvirt+nabs+iivn0) = x_vn0 (:,iivn0)
           hsml (npoin+nvirt+nabs+iivn0) = 2.* hsml_vn0(iivn0)     !  Chuan
           itype(npoin+nvirt+nabs+iivn0) = -14
           rho  (npoin+nvirt+nabs+iivn0) = 999.
           mass (npoin+nvirt+nabs+iivn0) = 999.
           p    (npoin+nvirt+nabs+iivn0) = 999.
           u    (npoin+nvirt+nabs+iivn0) = 999.
           lvn0                  (iivn0) = npoin + nvirt + nabs +iivn0  !   we do not use it...
        enddo
    endif       
    
    
!      ------  Initialize (changed, because we have to initialize x,vx,itype,rho,mass,hsml,u and p)

!      tmp initial vx set to 10 m/s
!       vx = 10.0 
    
        x0       = x  ;  x00      =  x0;  vx0     = vx ;  Z00      = 0.0 
        u0       = u  ;  rho0     = rho;  rho00   = rho;   
        
        dx       = 0.0;  dvx      = 0.0;  du      = 0.0;  drho     = 0.0                
        indvxdt  = 0.0;  exdvxdt  = 0.0;  ardvxdt = 0.0;  av       = 0.0
        c        = 0.0


        mass0 = mass  !  CH MP  8 April 2019

!        c        = 0.0;  s        = 0.0;  e       = 0.0
!        ds       = 0.0;  t        = 0.0;  tdsdt   = 0.0
!        ahdudt   = 0.0;  avdudt   = 0.0;  eta     = 0.0

END SUBROUTINE Setup_Global_Arrays_SW 


! -----------------------------------------------------------------
 
       SUBROUTINE Pint_Update_SW
 
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

END subroutine Pint_Update_SW



! -----------------------------------------------------------------
 
       SUBROUTINE Adaptive_dt_sph_SW 
 
! -----------------------------------------------------------------

!   Subroutine to calculate the critical time step using only c = sqrt(g.h) + v and hsml

implicit none

integer(ink) ipoin, idimn, itotv, jtotv , i, j, ielem
real   (irk) cc, vv, dtmin, dtmin2, v_diff,hsml_av
real   (irk) dt_aux, dxg_mod
real   (irk) cgra, Le 

dtmin    = 1.e10; dtmin2 = 1.e10
cgra     = const(1)

if (ic_adapt_dt_sph.eq.1 ) then                         ! dt = min(hsml/c)
   do ipoin = np_bakG_Total + 1, npoin                  ! see that we skip bakg points from gates
      if (if_out_domain(ipoin).eq.1) cycle
      if(rho(ipoin).lt.h_inf_SW) rho(ipoin) = h_inf_SW  !      old T08 limitation
      c(ipoin) = sqrt (cgra*rho(ipoin))
      itotv =  ipoin
      cc    = MAX ( 1.e-6, c(ipoin) )
      cc    = hsml(itotv)/(3*cc)
      dtmin = MIN (cc,dtmin)
   enddo
   cgra = const(1)
   ipoin    = ipoin     !  *** debug only
else if (ic_adapt_dt_sph.eq.2) then             !      ------  dt = min(hsml/ (c+v)) 
   do ipoin = np_bakG_Total + 1, npoin 
      if (if_out_domain(ipoin).eq.1) cycle
      if(rho(ipoin).lt.h_inf_SW) rho(ipoin) = h_inf_SW  !      old T08 limitation
      c(ipoin) = sqrt (cgra*rho(ipoin)) 
      if (if_out_domain(ipoin).eq.1) cycle
      itotv =  ipoin
      vv    = 0.0
      do idimn = 1,ndimn
         vv    = vv + vx(idimn,itotv)**2.
      enddo         
      vv    = sqrt (vv)
      cc    = MAX ( 1.e-6, c(itotv) ) + vv
      cc    = hsml(itotv)/(3*cc)
      dtmin = MIN (cc,dtmin)   
   enddo
else if (ic_adapt_dt_sph.eq.4) then             !      ------  dt = min(dxg/ (v)) 
   dxg_mod = (deltxg*deltxg + deltyg*deltyg)**0.5
   do ipoin = np_bakG_Total + 1, npoin 
      if (if_out_domain(ipoin).eq.1) cycle   
      itotv =  ipoin
      vv    = 0.0
      do idimn = 1,ndimn
         vv    = vv + vx(idimn,itotv)**2.
      enddo         
      vv    = sqrt (vv)
      cc    = dxg_mod/vv
      dtmin = MIN (cc,dtmin)   
   enddo
else if (ic_adapt_dt_sph.eq.3) then             !      ------  dt = min(hsml/ (c+v))  plus hsml/diffV 
   do ipoin = np_bakG_Total + 1, npoin 
      if (if_out_domain(ipoin).eq.1) cycle
      if(rho(ipoin).lt.h_inf_SW) rho(ipoin) = h_inf_SW  !      old T08 limitation
      c(ipoin) = sqrt (cgra*rho(ipoin))
      if (if_out_domain(ipoin).eq.1) cycle
      itotv =  ipoin
      vv    = 0.0
      do idimn = 1,ndimn
         vv  =  vv + vx(idimn,itotv)**2.
      enddo      
      vv    = sqrt (vv)
      cc    = MAX ( 1.e-6, c(itotv) ) + vv
      cc    = hsml(itotv)/(3*cc)
      dtmin = MIN (cc,dtmin)   
   enddo
   
   current=> last          !    Initializes linked list to last element 
   DO WHILE (associated(current))
      itotv = current%pair_i 
      jtotv = current%pair_j 
      v_diff = 0.0
      do idimn = 1,ndimn
         v_diff = v_diff + ( vx(idimn,itotv)-vx(idimn,jtotv) )**2.
      enddo
      v_diff = sqrt(v_diff)
      if (v_diff.le.1.e-7) v_diff=1.e-7
      hsml_av = ( hsml(itotv) + hsml(jtotv) )/2.
      dtmin2 = amin1( dtmin2, hsml_av/v_diff)
      dtmin2 = dtmin2/2.
      current=>current%next                
   ENDDO
   
   dtmin = MIN (dtmin,dtmin2)

else if (ic_adapt_dt_sph.eq.7) then             !      --FS  dt = Le/v
!   if (.not.allocated(intmat_FS)) then          ! 1st tstep intmat_FS not computed yet
      do ipoin = 1, npoin-1
         Le    = abs ( x(1,ipoin+1) - x(1, ipoin))
         vv    = abs (vx(1,ipoin+1) + vx(1,ipoin))/2.
         if (vv.gt.1.e-6 ) then                  ! if not, keep dtmin& do nothing
             cc    = 0.5*Le/vv
             dtmin = min (cc, dtmin)
         endif
      enddo
!   else                                         !  we v¡can use it now
!      do ielem = 1, nelem_FS
!         Le = geome_FS  (3, ielem)
!         i  = intmat_FS (1,ielem)
!         j  = intmat_FS (2,ielem)
!         vv = abs (vx(1,i)+vx(1,j))/2.
!         if (vv.gt.1.e-6 ) then                  ! if not, keep dtmin& do nothing
!             cc    = 0.5*Le/vv
!             dtmin = min (cc, dtmin)
!         endif
!      enddo
!   endif       

elseif (ic_adapt_dt_sph.eq.10 ) then                    ! dt = min(hsml/c, dxg/v)!       
   dxg_mod = (deltxg*deltxg + deltyg*deltyg)**0.5       !  MP June 2023 for Picacho jump
   cgra = const(1)
   do ipoin = np_bakG_Total + 1, npoin                  ! see that we skip bakg points from gates
      if (if_out_domain(ipoin).eq.1) cycle
      if(rho(ipoin).lt.h_inf_SW) rho(ipoin) = h_inf_SW  !      old T08 limitation
      c(ipoin) = sqrt (cgra*rho(ipoin))
      itotv =  ipoin
      cc    = MAX ( 1.e-6, c(ipoin) )
      cc    = hsml(itotv)/(3*cc)
      dtmin = MIN (cc,dtmin)
      vv    = 0.0
      do idimn = 1,ndimn
         vv  =  vv + vx(idimn,itotv)**2.
      enddo      
      vv    = sqrt (vv)
      cc    = dxg_mod/vv
      dtmin = MIN (cc,dtmin) 
   enddo
   ipoin    = ipoin     !  *** debug only
endif   

if (dtmin.lt.dt_sph) then
   dt_sph = dtmin
else
   if (ic_adapt_dt_sph.ne.0) dt_sph = amin1 (1.1*dt_sph, dtmin)
endif

END SUBROUTINE Adaptive_dt_sph_SW


! -----------------------------------------------------------------
 
       SUBROUTINE IntForces_SW ( p_only )
 
! -----------------------------------------------------------------

!   Subroutine to calculate the internal forces on the right hand side 
!              of the SW equations, i.e. the pressure gradient 
!   New argument is if_p_only, if set to 1 compute p(:) only    ! ** MP for FS
!
!   Sept 2021 CHGD MP Unsaturated materials See old routine just at the end _OLD  
implicit none

integer(ink), intent (IN), optional :: p_only                ! ** MP for FS
logical  if_p_only                                           ! ** MP for FS

integer(ink) ipoin, i, j, idimn
real   (irk) cgra, dens, cgra05, he, rhoij, h
real   (irk) xx, yy, zp, Zpx, Zpy
real   (irk) hxx, hxy, hyy,  kin_visc
real   (irk) denss, densw, facts, factw, fact_pwpavgd, factss, factww 
real   (irk) alphaS, alph, fact_alph, fact_alph2, alphI, alphJ, alphI_inv, alphJ_inv
real   (irk) mI, mJ

real   (irk) vcc
real   (irk), allocatable:: dvx(:),   divvK0(:)
real   (irk), allocatable:: gradn (:,:)                 !  grad n
real   (irk), allocatable:: tauxx(:), tauxy(:), tauyy(:)
real   (irk), allocatable:: wxy(:)  ! used for gradient normalization tests ******

integer (ink) ic_cosTheta         ! C16..1  use g projected along normal to terrain 
real    (irk) deng

real    (irk) fact_s2w     ! auxiliar for solid behaviour

!      ------  Initialization 

if_p_only = PRESENT (p_only)           ! ** MP for FS

if (.not.allocated (dvx) ) then 
    allocate (  dvx(ndimn) )
    dvx = 0.0
endif

if (.not.allocated (gradn).AND. ic_ws_interact.eq.2) then
    allocate (gradn (ndimn, npoin))
    gradn = 0.0
endif

if (.not.allocated (tauxy) ) then !      ------  tauxy will be tau*h/dens
    allocate ( tauxx(ntotal))
    allocate ( tauxy(ntotal))
    allocate ( tauyy(ntotal))
    allocate (   wxy(ntotal))
endif

tauxx   = 0.0 ; tauxy  = 0.0 ; tauyy  = 0.0   !  we will not worry about interactions
wxy     = 0.0                                 ! for the test *****

indvxdt = 0.0                                 !     other than 0 or -21 
dvx     = 0.0                                

cgra        = const(1)            ! gravz already defined for sloped plane. Default is g
cgra05      = gravz*0.5           ! for pressure, cgra05 uses gravz
dens        = const  (2)            
denss       = const (17)             
densw       = const (18)
ic_cosTheta = const(16) + 0.01    !  1 ... use g projected along normal to terrain
                                  
if ( ic_K0.eq.1) then             
   allocate (divvK0(ntotal))
   divvK0 = 0.0
   current=> last
   Do while (associated(current))
      i = current%pair_i 
      j = current%pair_j 
      if ( iabs(itype(i)).eq.1.AND.iabs(itype(j)).eq.1 ) THEN
         do idimn = 1,ndimn
            dvx(idimn) = vx(idimn,i) - vx(idimn,j) 
         enddo        
         vcc = dvx(1)*current%dwdx(1)        
         do idimn = 2,ndimn
            vcc = vcc + dvx(idimn)*current%dwdx(idimn)
         enddo     
         divvK0(i) = divvK0(i) + mass(j)*vcc
         divvK0(j) = divvK0(j) + mass(i)*vcc    
      endif
      current=>current%next
   enddo
 
endif   

!      ------  Get pressures, i.e.  g h*h/2  and multiply by K0 if required

!  to use h H alg with WIR

if      (ic_ws_Interact.eq.1) then     !  we need to compute hs  
         call Get_hs_at_w_SW           !  and its gradient at water nodes
elseif (ic_ws_Interact.eq.2) then      !   DFs get hs & vs at w nodes
         call Get_s_at_w_DF_SW         !           hw & vw at s noes
endif                                  !           h_DF = hs+hw & n_DF

where(rho.lt.h_inf_SW) rho=h_inf_SW    ! new filter 04 december 2011   

do i = 1,npoin

   if ( if_Out_Domain(i).eq.1) CYCLE
                            
   if ( ic_cosTheta.eq.1) then          !  we will use  gz if ic_cosTheta=1
        xx    = x(1,i)
        yy    = (yming+ymaxg)/2
        if (ndimn.eq.2) yy = x(2,i)
        call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)   
        deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
        gravz =  cgra/deng
        cgra05 = 0.5*gravz                                      
   endif

   if (itype(i).eq.2) then              ! 26th May 2014 use densw if read.
       if (const(18).gt.1.e-6) then
          dens = const(18)
       else
          dens = 1000.
       endif
   endif    
   
   IF (ic_ws_Interact.eq.0) THEN        ! No interaction. Use 0.5g*h^2 or 
                                        ! 0.5g*(h^2-Z^2)  SWL=0 !!!
      if (ic_SWalg.eq.0)  then
      
         p(i) = cgra05*rho(i)*rho(i)    ! classical landslide formulation
         
      elseif (ic_SWalg.eq.1)  then      ! classical landslide formulation
      
         xx = x(1,i)                    ! get xx and yy -> Zp from topo module
         yy = (yming+ymaxg)/2   
         if (ndimn.eq.2) yy = x(2,i)  
         call Get_Z_Topo (i, xx, yy, Zp )
         p(i) = cgra05*(rho(i)*rho(i)-Zp*Zp)
         
      endif  
   
   ELSEIF (ic_ws_Interact.eq.1) THEN    ! Interaction case Use 0.5g*h^2 
                                        ! or 0.5g*(h^2-(Z+hs)^2)  SWL=0 !!!   
      if (ic_SWalg.eq.0)  then
      
         p(i) = cgra05*rho(i)*rho(i)    ! classical landslide formulation 
         
      elseif (ic_SWalg.eq.1)  then      ! soil as in landslides, 
                                        ! water as water bodies. Zp = hs+Z
         if (itype(i).eq.1) then
         
            p(i) = cgra05*rho(i)*rho(i)
            
         elseif(itype(i).eq.2) then

            Zp   = hs_plus_Z_w (i)
            p(i) = cgra05*(rho(i)*rho(i)-Zp*Zp)
           
         endif
         
      endif  
   
   ELSEIF (ic_ws_Interact.eq.2) THEN           ! Sept 2021 CHGD MP Major MArch 2014 DFs
         
     alphaS = const(14)                        ! CH MP 4 March 2019                         
     alph   = alphaS                           ! alph*h saturated layer
     if (ic_hrelSat.EQ.1) then
         alph = hrelSat_DF(i)                
     endif    
     fact_alph =  1.0 
     fact_alph2 = 1.0  
     if (i.gt.npoin_s) then                    !  just for water points 
         fact_alph =  alph 
         fact_alph2 = alph*alph                 
     endif      
                                               !  solid:  Ps = grad(0.5*gz*(hs*h)
     p(i) = cgra05*h_DF(i)*rho(i)*fact_alph    !  water:  Pw = grad(0.5*gz*hw*hsat) hsat=h*alph   
                                                  
                                               !  pwp contribution...hw is n*h-df*alph
     if (icpwp.GE.1) then                      !  s: -(pwp/denss)*n*h_DF*alph  because of -(1/hs)grad *                
         if (i.le.npoin_s) then                !  w: +(pwp/densw)*n*h_DF*alph  
             p(i) = p(i) - (pwp_avgd(i)/denss)*n_DF(i)*h_DF(i)*alph
         elseif (i.gt.npoin_s.and.i.le.npoin) then
             p(i) = p(i) + (pwp_avgd(i)/densw)*n_DF(i)*h_DF(i)*alph           
         endif    
     endif    
      
   ENDIF
   
   if (ic_K0.eq.1) then                  ! Here we multiply sigmazz or pressure by K0
      if ( divvK0(i).ge.0.AND.i.le.npoin_s ) then
         p(i) = p(i)*K0_active
      else
         p(i) = p(i)*K0_pasive
      endif
   endif
enddo                                    ! ------  end loop ipoin
      
if (if_p_only) then                 ! ** MP for FS 
    return
endif

if (ic_end_solid.eq.1) then        !  solid behaviour via decreasing p 
    if (law_end_solid.eq.0) then   !  MP 30th Jan 2024
       fact_s2w = 0.0  
       if (time_sph.GT.t1_end_solid) then
           fact_s2w = 1.0
       endif
    elseif ( law_end_solid.eq.1 ) then
       fact_s2w = 1.0
       if  (time_sph.lt.t1_end_solid) then
           fact_s2w = (time_sph/t1_end_solid)**fct_exp_end_solid 
       endif
    elseif( law_end_solid.eq.2) then
       fact_s2w = 0.0
       if  (time_sph.gt.t1_end_solid) then
          fact_s2w = 1.0 - exp( -fct_exp_end_solid *  (time_sph-t1_end_solid) )
       endif 
   endif       
    p (1:npoin_s)= p(1:npoin_s)*fact_s2w
endif

!      ------  Calculate SPH sum for pressure force -p,a/rho

current=> last          !    Initializes linked list to last element 

DO WHILE (associated(current))

   i = current%pair_i 
   j = current%pair_j 

   IF ( current%Pint_type.eq.0 ) THEN

      he = 0.e0 
      mI = mass(i)                  !  --- chgd MP Sept 2021                                 
      mJ = mass(j)                  !   compute (-1/hs)*grap Ps  same for w

      if(pa_sph.eq.1) then  
         rhoij = 1.e0/(rho(i)*rho(j))                     !   ----  For SPH algorithm 1 
         do idimn = 1,ndimn
            h = -(p(i) + p(j))*current%dwdx(idimn)        !   ----  Pressure part
            h = h*rhoij
            indvxdt(idimn,i) = indvxdt(idimn,i) + mJ*h
            indvxdt(idimn,j) = indvxdt(idimn,j) - mI*h
         enddo  
      
      else if (pa_sph.eq.2) then                           !      ----- For SPH algorithm 2  
   
         do idimn = 1,ndimn 
            h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*current%dwdx(idimn) 
            indvxdt(idimn,i) = indvxdt(idimn,i) + mJ*h
            indvxdt(idimn,j) = indvxdt(idimn,j) - mI*h
         enddo
      
      else if (pa_sph.eq.3) then                           !      ----- For SPH algorithm 3  
                                                           !            ALPHA correction factor JB  
         do idimn = 1,ndimn 
            h = -(p(i)/(alphaJB(i)*rho(i)**2) + p(j)/(alphaJB(j)*rho(j)**2) )*current%dwdx(idimn) 
            he = he + (vx(idimn,j) - vx(idimn,i))*h
            indvxdt(idimn,i) = indvxdt(idimn,i) + mJ*h
            indvxdt(idimn,j) = indvxdt(idimn,j) - mI*h
         enddo
         
      else if (pa_sph.eq.4) then                           !      ----- For SPH  4  corr grad, 
                                                           !           alpha for convenience 1D ONLY  
         do idimn = 1,ndimn 
            indvxdt(idimn,i) = indvxdt(idimn,i) - ((p(j) - p(i))*current%dwdx(idimn)*mJ/rho(j)) 
            indvxdt(idimn,j) = indvxdt(idimn,j) - ((p(j) - p(i))*current%dwdx(idimn)*mI/rho(i))
         enddo
      
      endif 
   ENDIF
   current=>current%next                
ENDDO

IF (ic_ws_interact.eq.2.AND.ic_gradn.eq.1 ) THEN       !      ------ Obtain (1/h)*grad n ...correction Shiomi
                                                       !  Change MP 30th June 2019 gradn only in general case
                                                       !    Not in shiomi approx

    current=> last          !    Initializes linked list to last element 
    
    DO WHILE (associated(current))

       i = current%pair_i 
       j = current%pair_j 

       IF ( current%Pint_type.eq.0 ) THEN
     
          if(pa_sph.eq.1) then                  !   ----  For SPH algorithm 1  
             rhoij = 1.e0/(rho(i)*rho(j))                      
             do idimn = 1,ndimn
                h = (n_DF(i) + n_DF(j))*current%dwdx(idimn)        
                h = h*rhoij
                gradn(idimn,i) = gradn(idimn,i) + mass(j)*h
                gradn(idimn,j) = gradn(idimn,j) - mass(i)*h
             enddo  
                 
          else if (pa_sph.eq.2) then            !      ----- For SPH algorithm 2     
             do idimn = 1,ndimn 
                h = (n_DF(i)/rho(i)**2 + n_DF(j)/rho(j)**2)*current%dwdx(idimn)  
                gradn(idimn,i) = gradn(idimn,i) + mass(j) * h
                gradn(idimn,j) = gradn(idimn,j) - mass(i) * h
             enddo
          endif
   
       ENDIF
       
       current=>current%next 
       
    ENDDO
    !
    !       (-1/hs)*(grad n) facts (-1/hw)*(grad n) factw
    !       facts =  0.5*(densw/denss)*gz*(h_df*alph)**2 + (h_df*alph)* pwp/denss
    !       factw = - 0.5*gz*(h_df*alph)**2              - (h_df*alph)* pwp/densw

    factss = 0.0; factww = 0.0
         
    do i = 1, npoin_s            ! ----------------  SOIL contribution
       factss  =  cgra05*densw/denss 
       alph   = const(14)                       
       if (ic_hrelSat.EQ.1) then
           alph = hrelSat_DF(i)                
       endif 
       factss = factss * (h_DF(i)*alph)**2. 
       fact_pwpavgd = 0.0
       if (icpwp.gt.0) then
            fact_pwpavgd = (h_DF(i)*alph)*(pwp_avgd(i)/denss)
       endif
       factss = factss + fact_pwpavgd
       if (h_sw(i).gt.h_inf_SW) then         !   MP changed 13 may 19, to avoid gradn in w w/o s and viceversa
           gradn(:,i) =  gradn(:,i)*( - factss - fact_pwpavgd  ) 
       else
           gradn(:,i) = 0.0
       endif            
    enddo
   
                                 ! ----------------  WATER contribution        
     do i = npoin_s +1, npoin  
        factww  = cgra05                         
        if (ic_hrelSat.EQ.1) then
            alph = hrelSat_DF(i)                
        endif 
        factww = factww * (h_DF(i)*alph)**2. 
        fact_pwpavgd = 0.0
        if (icpwp.gt.0) then 
            fact_pwpavgd = (h_DF(i)*alph)*(pwp_avgd(i)/densw)
        endif
        factww = factww + fact_pwpavgd
        if (h_sw(i).gt.h_inf_SW) then         !   MP changed 13 may 19, to avoid gradn in w w/o s and viceversa
            gradn(:,i) =  gradn(:,i)*( factww  + fact_pwpavgd  )
        else
           gradn(:,i) = 0.0
        endif    
     enddo

    indvxdt = indvxdt + gradn
    
ENDIF    

!IF (const(12).ge.1.e-10.AND.ndimn.eq.2) then !      ------  Only if viscosity is present
IF (1.eq.0) then
   
! ***********
  write(*,*) ' stopped. This is for body visco. Deactivate const(12) makeit zero '
  write (*,*) ' hit any key '
  read  (*,*) i 
  pause
  ! rho = 1.0
  ! mass = 1.0
  vx (1,:) = 0.0
  do i = 1,npoin
     vx(2,i) =  x(1,i) 
  enddo   
  
! ************+

  kin_visc = const(12) / dens    
    
   current=> last            !    Initializes linked list to last element 
   DO WHILE (associated(current))
   
      i = current%pair_i 
      j = current%pair_j  

      IF ( current%Pint_type.eq.0.OR.current%Pint_type.eq.-1 ) THEN

         if (current%Pint_type.eq.0) then
            do idimn = 1,ndimn
               dvx(idimn) = vx(idimn,j) - vx(idimn,i)
            enddo
         elseif (current%Pint_type.eq.-21) then 
            do idimn = 1,ndimn
               dvx(idimn) =  - vx(idimn,i)   !  part. j belongs to wall vj=0
            enddo     
         endif
                
         hxx = 2.e0*dvx(1)*current%dwdx(1) - dvx(2)*current%dwdx(2) 
         hxy =      dvx(1)*current%dwdx(2) + dvx(2)*current%dwdx(1)
         hyy = 2.e0*dvx(2)*current%dwdx(2) - dvx(1)*current%dwdx(1)  
         hxx = 2.e0/3.e0*hxx
         hyy = 2.e0/3.e0*hyy           
         tauxx(i) = tauxx(i) + mass(j)*hxx/rho(j)
         tauxx(j) = tauxx(j) + mass(i)*hxx/rho(i)   
         tauxy(i) = tauxy(i) + mass(j)*hxy/rho(j)
         tauxy(j) = tauxy(j) + mass(i)*hxy/rho(i)            
         tauyy(i) = tauyy(i) + mass(j)*hyy/rho(j)
         tauyy(j) = tauyy(j) + mass(i)*hyy/rho(i) 
                   
         current=>current%next
      ENDIF
     
   ENDDO
   
   current=> last            !    Initializes linked list to last element 
   DO WHILE (associated(current))
   
      i = current%pair_i 
      j = current%pair_j  

      IF ( current%Pint_type.eq.0.OR.current%Pint_type.eq.-1 ) THEN

         if (current%Pint_type.eq.0) then
            do idimn = 1,ndimn
               dvx(idimn) = x(idimn,j) - x(idimn,i)
            enddo      
         endif 
         
         hxy =      dvx(1)*current%dwdx(1)  
         wxy(i) = wxy(i) + mass(j)*hxy/rho(j)
         wxy(j) = wxy(j) + mass(i)*hxy/rho(i) 
                   
         current=>current%next
      ENDIF
     
   ENDDO
   
   do i = 1, npoin
      tauxy(i) = tauxy(i)/wxy(i)   !  ******
      tauxx(i) = kin_visc * tauxx(i)*rho(i)  ! dens*g*h/dens *h
      tauxy(i) = kin_visc * tauxy(i)*rho(i)  
      tauyy(i) = kin_visc * tauyy(i)*rho(i)
   enddo
 
   
! ******************

   do i = 1, npoin
      write(res_file2,10001) i, x(1,i), x(2,i), tauxy(i), tauxx(i), tauyy(i)
   enddo
   i = i
10001 format (i8,2x, 5(g11.4,2x))
! ******************   

!      ------  Calculate SPH sum for viscous force (eta Tab),b/ grav 

   
! ******************

   tauxy = 0.0
   tauxx = 0.0
   tauyy = 0.0
   do i = 1, npoin
      tauxx(i) =   x(1,i) 
   enddo
   indvxdt = 0.0

! ******************   


   current=> last          !    Initializes linked list to last element 
   
   DO WHILE (associated(current))

      i = current%pair_i 
      j = current%pair_j 
   
      IF ( current%Pint_type.eq.0.OR.current%Pint_type.eq.-1 ) THEN

         if (pa_sph.eq.1) then             !   ----  For SPH algorithm 1     
            rhoij = 1.e0/(rho(i)*rho(j))  
            do idimn = 1,ndimn
               h = 0.0
               if (idimn.eq.1) then        
                  h = h + ( tauxx(i) +  tauxx(j) ) * current%dwdx(1) 
                  h = h + ( tauxy(i) +  tauxy(j) ) * current%dwdx(2)
               elseif (idimn.eq.2) then                     
                  h = h + ( tauxy(i) +  tauxy(j) ) * current%dwdx(1) +  &
                          ( tauyy(i) +  tauyy(j) ) * current%dwdx(2)
               endif  
               h = h*rhoij
               indvxdt(idimn,i) = indvxdt(idimn,i) + mass(j)*h
               indvxdt(idimn,j) = indvxdt(idimn,j) - mass(i)*h
            enddo   
   
         else if (pa_sph.eq.2) then             !   ----  For SPH algorithm 2
            do idimn = 1,ndimn 
               h = 0.0 
               if (idimn.eq.1) then                            
                   h = h + ( tauxx(i)/rho(i)**2 +  tauxx(j)/rho(j)**2 )*current%dwdx(1)
                   h = h + ( tauxy(i)/rho(i)**2 +  tauxy(j)/rho(j)**2 )*current%dwdx(2)
               elseif (idimn.eq.2) then                     
                   h = h + ( tauxy(i)/rho(i)**2                    &
                         +   tauxy(j)/rho(j)**2 )*current%dwdx(1)  &
                         + ( tauyy(i)/rho(i)**2                    &  
                         +   tauyy(j)/rho(j)**2 )*current%dwdx(2)  
               endif 
               indvxdt(idimn,i) = indvxdt(idimn,i) + mass(j)*h
               indvxdt(idimn,j) = indvxdt(idimn,j) - mass(i)*h
            enddo
         endif

         current=>current%next

      ENDIF    
                     
   ENDDO
! ******************

   write(res_file2,*) ' -------------------------------------------'
   do i = 1, npoin
      write(res_file2,10001) i, x(1,i), x(2,i) , indvxdt (1,i), indvxdt (2,i) 
   enddo
   stop
! ******************  

   
ENDIF                                 !      -----  end of viscous terms computations

if (pa_sph.eq.4) then                 !      ----- For SPH  4  corr grad, alpha for convenience 1D ONLY   
    do i = 1,npoin
       indvxdt(:,i) = indvxdt(:,i)/alphaJB(i)
    enddo
endif 

if (ic_K0.eq.1) then
   deallocate (divvK0)
endif   

if (ic_histogram.eq.1.and.np_bakg_Total.gt.0) then
   indvxdt(:,1:np_bakg_Total) = 0.0
endif   
 
END SUBROUTINE   IntForces_SW 



! -----------------------------------------------------------------
 
       SUBROUTINE IntForces_SW_OLD ( p_only )
 
! -----------------------------------------------------------------

!   Subroutine to calculate the internal forces on the right hand side 
!              of the SW equations, i.e. the pressure gradient 
!   New argument is if_p_only, if set to 1 compute p(:) only    ! ** MP for FS
 
implicit none

integer(ink), intent (IN), optional :: p_only                ! ** MP for FS
logical  if_p_only                                           ! ** MP for FS

integer(ink) ipoin, i, j, idimn
real   (irk) cgra, dens, cgra05, he, rhoij, h
real   (irk) xx, yy, zp, Zpx, Zpy
real   (irk) hxx, hxy, hyy,  kin_visc
real   (irk) denss, densw, facts, factw, fact_pwpavgd 
real   (irk) alphaS, alph, fact_alph, fact_alph2, alphI, alphJ, alphI_inv, alphJ_inv
real   (irk) mI, mJ

real   (irk) vcc
real   (irk), allocatable:: dvx(:),   divvK0(:)
real   (irk), allocatable:: gradn (:,:)                 !  grad n
real   (irk), allocatable:: tauxx(:), tauxy(:), tauyy(:)
real   (irk), allocatable:: wxy(:)  ! used for gradient normalization tests ******

integer (ink) ic_cosTheta         ! C16..1  use g projected along normal to terrain 
real    (irk) deng

real    (irk) fact_s2w     ! auxiliar for solid behaviour

!      ------  Initialization 

if_p_only = PRESENT (p_only)           ! ** MP for FS

if (.not.allocated (dvx) ) then 
    allocate (  dvx(ndimn) )
    dvx = 0.0
endif

if (.not.allocated (gradn).AND. ic_ws_interact.eq.2) then
    allocate (gradn (ndimn, npoin))
    gradn = 0.0
endif

if (.not.allocated (tauxy) ) then !      ------  tauxy will be tau*h/dens
    allocate ( tauxx(ntotal))
    allocate ( tauxy(ntotal))
    allocate ( tauyy(ntotal))
    allocate (   wxy(ntotal))
endif

tauxx   = 0.0 ; tauxy  = 0.0 ; tauyy  = 0.0   !  we will not worry about interactions
wxy     = 0.0                                 ! for the test *****

indvxdt = 0.0                                 !     other than 0 or -21 
dvx     = 0.0                                

cgra        = const(1)            ! gravz already defined for sloped plane. Default is g
cgra05      = gravz*0.5           ! for pressure, cgra05 uses gravz
dens        = const  (2)            
denss       = const (17)             
densw       = const (18)
ic_cosTheta = const(16) + 0.01    !  1 ... use g projected along normal to terrain
                                  
if ( ic_K0.eq.1) then             
   allocate (divvK0(ntotal))
   divvK0 = 0.0
   current=> last
   Do while (associated(current))
      i = current%pair_i 
      j = current%pair_j 
      if ( iabs(itype(i)).eq.1.AND.iabs(itype(j)).eq.1 ) THEN
         do idimn = 1,ndimn
            dvx(idimn) = vx(idimn,i) - vx(idimn,j) 
         enddo        
         vcc = dvx(1)*current%dwdx(1)        
         do idimn = 2,ndimn
            vcc = vcc + dvx(idimn)*current%dwdx(idimn)
         enddo     
         divvK0(i) = divvK0(i) + mass(j)*vcc
         divvK0(j) = divvK0(j) + mass(i)*vcc    
      endif
      current=>current%next
   enddo
 
endif   

!      ------  Get pressures, i.e.  g h*h/2  and multiply by K0 if required

!  to use h H alg with WIR

if      (ic_ws_Interact.eq.1) then     !  we need to compute hs  
         call Get_hs_at_w_SW           !  and its gradient at water nodes
elseif (ic_ws_Interact.eq.2) then      !   DFs get hs & vs at w nodes
         call Get_s_at_w_DF_SW         !           hw & vw at s noes
endif                                  !           h_DF = hs+hw & n_DF

where(rho.lt.h_inf_SW) rho=h_inf_SW    ! new filter 04 december 2011   

do i = 1,npoin

   if ( if_Out_Domain(i).eq.1) CYCLE
                            
   if ( ic_cosTheta.eq.1) then          !  we will use  gz if ic_cosTheta=1
        xx    = x(1,i)
        yy    = (yming+ymaxg)/2
        if (ndimn.eq.2) yy = x(2,i)
        call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)   
        deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
        gravz =  cgra/deng
        cgra05 = 0.5*gravz                                      
   endif

   if (itype(i).eq.2) then              ! 26th May 2014 use densw if read.
       if (const(18).gt.1.e-6) then
          dens = const(18)
       else
          dens = 1000.
       endif
   endif    
   
   IF (ic_ws_Interact.eq.0) THEN        ! No interaction. Use 0.5g*h^2 or 
                                        ! 0.5g*(h^2-Z^2)  SWL=0 !!!
      if (ic_SWalg.eq.0)  then
      
         p(i) = cgra05*rho(i)*rho(i)    ! classical landslide formulation
         
      elseif (ic_SWalg.eq.1)  then      ! classical landslide formulation
      
         xx = x(1,i)                    ! get xx and yy -> Zp from topo module
         yy = (yming+ymaxg)/2   
         if (ndimn.eq.2) yy = x(2,i)  
         call Get_Z_Topo (i, xx, yy, Zp )
         p(i) = cgra05*(rho(i)*rho(i)-Zp*Zp)
         
      endif  
   
   ELSEIF (ic_ws_Interact.eq.1) THEN    ! Interaction case Use 0.5g*h^2 
                                        ! or 0.5g*(h^2-(Z+hs)^2)  SWL=0 !!!   
      if (ic_SWalg.eq.0)  then
      
         p(i) = cgra05*rho(i)*rho(i)    ! classical landslide formulation 
         
      elseif (ic_SWalg.eq.1)  then      ! soil as in landslides, 
                                        ! water as water bodies. Zp = hs+Z
         if (itype(i).eq.1) then
         
            p(i) = cgra05*rho(i)*rho(i)
            
         elseif(itype(i).eq.2) then

            Zp   = hs_plus_Z_w (i)
            p(i) = cgra05*(rho(i)*rho(i)-Zp*Zp)
           
         endif
         
      endif  
   
   ELSEIF (ic_ws_Interact.eq.2) THEN           ! MArch 2014 DFs
         
     alphaS = const(14)                        ! CH MP 4 March 2019                         
     alph   = alphaS                           ! alph*h saturated layer
     if (ic_hrelSat.EQ.1) then
         alph = hrelSat_DF(i)                
     endif    
     fact_alph =  1.0 
     fact_alph2 = 1.0  
     if (i.gt.npoin_s) then                    ! correct only water pressure terms
         fact_alph =  alph 
         fact_alph2 = alph*alph                 
     endif                                     !  hs*dvs/dt = grad(0.5*((1-n).h^2.gz)
     p(i) = cgra05*h_DF(i)*rho(i)*fact_alph2   !  same for both cases (4th May 14)
                                               !  hw*dvw/dt = grad(0.5*((  n  h^2.gz)                                         
     if (icpwp.GE.1) then                      !  add the Dpwp contribution from dN_DF
         if (i.le.npoin_s) then                !  GE 1 MP 30th May 2017
             p(i) = p(i) - pwp_avgd(i)*n_DF(i)*h_DF(i)/denss
         elseif (i.gt.npoin_s.and.i.le.npoin) then
             p(i) = p(i) + alph * pwp_avgd(i)*n_DF(i)*h_DF(i)/densw
         endif    
     endif    
      
   ENDIF
   
   if (ic_K0.eq.1) then                  ! Here we multiply sigmazz or pressure by K0
      if ( divvK0(i).ge.0.AND.i.le.npoin_s ) then
         p(i) = p(i)*K0_active
      else
         p(i) = p(i)*K0_pasive
      endif
   endif
enddo                                    ! ------  end loop ipoin
      
if (if_p_only) then                 ! ** MP for FS
    return
endif

if (ic_end_solid.eq.1) then        !  solid behaviour via decreasing p
    if (law_end_solid.eq.0) then
        fact_s2w = 1.0
        if (time_sph.lt.t1_end_solid) fact_s2w = 0.0
    elseif (law_end_solid.eq.1) then
        fact_s2w = (time_sph/t1_end_solid)**fct_exp_end_solid
        if (fact_s2w.gt.1.0) fact_s2w= 1.0
    endif    
    p (1:npoin_s)= p(1:npoin_s)*fact_s2w
endif

!      ------  Calculate SPH sum for pressure force -p,a/rho

current=> last          !    Initializes linked list to last element 

DO WHILE (associated(current))

   i = current%pair_i 
   j = current%pair_j 

   IF ( current%Pint_type.eq.0 ) THEN
      
      alphI = 1.0                        !  Default for soil particles
      alphJ = 1.0
      alphI_inv = 1.0
      alphJ_inv = 1.0
      if (i.gt.npoin_s) then             ! I is water, so J is water too
         if (ic_hrelSat.EQ.1) then
             alphI = hrelSat_DF(i) 
             if (alphI.lt.1.e-6) then
                 alphI_inv = 0.0
             else
                 alphI_inv = 1.0/alphI
             endif
             alphJ = hrelSat_DF(j) 
             if (alphJ.lt.1.e-6) then
                 alphJ_inv = 0.0
             else
                 alphJ_inv = 1.0/alphJ
             endif
         endif
      else
         alphI = alphI 
      endif

      he = 0.e0
      if (ic_hrelSat.eq.1) then
         mI = mass(i) * alphI_inv
         mJ = mass(j) * alphJ_inv
      else
         mI = mass(i)  
         mJ = mass(j)  
      endif   

      if(pa_sph.eq.1) then  
         rhoij = 1.e0/(rho(i)*rho(j))                     !   ----  For SPH algorithm 1 
         do idimn = 1,ndimn
            h = -(p(i) + p(j))*current%dwdx(idimn)        !   ----  Pressure part
            h = h*rhoij
            indvxdt(idimn,i) = indvxdt(idimn,i) + mJ*h
            indvxdt(idimn,j) = indvxdt(idimn,j) - mI*h
         enddo  
      
      else if (pa_sph.eq.2) then                           !      ----- For SPH algorithm 2  
   
         do idimn = 1,ndimn 
            h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*current%dwdx(idimn) 
            indvxdt(idimn,i) = indvxdt(idimn,i) + mJ*h
            indvxdt(idimn,j) = indvxdt(idimn,j) - mI*h
         enddo
      
      else if (pa_sph.eq.3) then                           !      ----- For SPH algorithm 3  
                                                           !            ALPHA correction factor JB  
         do idimn = 1,ndimn 
            h = -(p(i)/(alphaJB(i)*rho(i)**2) + p(j)/(alphaJB(j)*rho(j)**2) )*current%dwdx(idimn) 
            he = he + (vx(idimn,j) - vx(idimn,i))*h
            indvxdt(idimn,i) = indvxdt(idimn,i) + mJ*h
            indvxdt(idimn,j) = indvxdt(idimn,j) - mI*h
         enddo
         
      else if (pa_sph.eq.4) then                           !      ----- For SPH  4  corr grad, 
                                                           !           alpha for convenience 1D ONLY  
         do idimn = 1,ndimn 
            indvxdt(idimn,i) = indvxdt(idimn,i) - ((p(j) - p(i))*current%dwdx(idimn)*mJ/rho(j)) 
            indvxdt(idimn,j) = indvxdt(idimn,j) - ((p(j) - p(i))*current%dwdx(idimn)*mI/rho(i))
         enddo
      
      endif 
   ENDIF
   current=>current%next                
ENDDO

IF (ic_ws_interact.eq.2.AND.ic_gradn.eq.1 ) THEN       !      ------ Obtain (1/h)*grad n ...correction Shiomi
                                                       !  Change MP 30th June 2019 gradn only in general case
                                                       !    Not in shiomi approx

    current=> last          !    Initializes linked list to last element 
    
    DO WHILE (associated(current))

       i = current%pair_i 
       j = current%pair_j 

       IF ( current%Pint_type.eq.0 ) THEN
     
          if(pa_sph.eq.1) then                  !   ----  For SPH algorithm 1  
             rhoij = 1.e0/(rho(i)*rho(j))                      
             do idimn = 1,ndimn
                h = (n_DF(i) + n_DF(j))*current%dwdx(idimn)        
                h = h*rhoij
                gradn(idimn,i) = gradn(idimn,i) + mass(j)*h
                gradn(idimn,j) = gradn(idimn,j) - mass(i)*h
             enddo  
                 
          else if (pa_sph.eq.2) then            !      ----- For SPH algorithm 2     
             do idimn = 1,ndimn 
                h = (n_DF(i)/rho(i)**2 + n_DF(j)/rho(j)**2)*current%dwdx(idimn)  
                gradn(idimn,i) = gradn(idimn,i) + mass(j) * h
                gradn(idimn,j) = gradn(idimn,j) - mass(i) * h
             enddo
          endif
   
       ENDIF
       
       current=>current%next 
       
    ENDDO
    
    facts  =  cgra05*densw/denss     
    do i = 1, npoin_s
       fact_pwpavgd = 0.0
       if (icpwp.gt.0) then
            fact_pwpavgd = h_DF(i)*pwp_avgd(i)/denss
       endif
       if (h_sw(i).gt.h_inf_SW) then         !   MP changed 13 may 19, to avoid gradn in w w/o s and viceversa
           gradn(:,i) =  gradn(:,i)*( - facts * h_DF(i)*h_DF(i) - fact_pwpavgd  ) 
       else
           gradn(:,i) = 0.0
       endif            
    enddo

     alph   = const(14) 
     facts  = cgra05      
     do i = npoin_s +1, npoin                         
        if (ic_hrelSat.EQ.1) then
            alph = hrelSat_DF(i)                
        endif 
        fact_pwpavgd = 0.0
        if (icpwp.gt.0) then 
            fact_pwpavgd = h_DF(i)*pwp_avgd(i)/densw
        endif
       if (h_sw(i).gt.h_inf_SW) then         !   MP changed 13 may 19, to avoid gradn in w w/o s and viceversa
           gradn(:,i) =  gradn(:,i)*( facts * h_DF(i)* h_DF(i)*alph + fact_pwpavgd  ) 
       else
           gradn(:,i) = 0.0
       endif    
     enddo
!    do i = npoin_s +1, npoin
!       facts = cgra05
!       fact_pwpavgd = 0.0
!       if (icpwp.gt.0) then 
!           fact_pwpavgd = h_DF(i)*pwp_avgd(i)/densw
!       endif
!       gradn(:,i) =  gradn(:,i)*( facts * h_DF(i)* h_DF(i) + fact_pwpavgd  )  ! chk sign

!    GRADN = 0.0   !  ***********************************************

    indvxdt(:, 1:npoin) = indvxdt(:, 1:npoin) + gradn                          ! Saeidmt 22MAy2021  Not considered
    
ENDIF    

!indvxdt(:, 1:npoin) = 0.0   !  WORK19 STOPS SOIL & water
!indvxdt(:, 1:npoin_s) = 0.0   !  WORK19 STOPS SOIL only

!      ------   Deal with body viscosity iff  const(12) is not zero
!               Calculate SPH sum for shear tensor Tauxy  

!IF (const(12).ge.1.e-10.AND.ndimn.eq.2) then !      ------  Only if viscosity is present
IF (1.eq.0) then
   
! ***********
  write(*,*) ' stopped. This is for body visco. Deactivate const(12) makeit zero '
  write (*,*) ' hit any key '
  read  (*,*) i 
  pause
  ! rho = 1.0
  ! mass = 1.0
  vx (1,:) = 0.0
  do i = 1,npoin
     vx(2,i) =  x(1,i) 
  enddo   
  
! ************+

  kin_visc = const(12) / dens    
    
   current=> last            !    Initializes linked list to last element 
   DO WHILE (associated(current))
   
      i = current%pair_i 
      j = current%pair_j  

      IF ( current%Pint_type.eq.0.OR.current%Pint_type.eq.-1 ) THEN

         if (current%Pint_type.eq.0) then
            do idimn = 1,ndimn
               dvx(idimn) = vx(idimn,j) - vx(idimn,i)
            enddo
         elseif (current%Pint_type.eq.-21) then 
            do idimn = 1,ndimn
               dvx(idimn) =  - vx(idimn,i)   !  part. j belongs to wall vj=0
            enddo     
         endif
                
         hxx = 2.e0*dvx(1)*current%dwdx(1) - dvx(2)*current%dwdx(2) 
         hxy =      dvx(1)*current%dwdx(2) + dvx(2)*current%dwdx(1)
         hyy = 2.e0*dvx(2)*current%dwdx(2) - dvx(1)*current%dwdx(1)  
         hxx = 2.e0/3.e0*hxx
         hyy = 2.e0/3.e0*hyy           
         tauxx(i) = tauxx(i) + mass(j)*hxx/rho(j)
         tauxx(j) = tauxx(j) + mass(i)*hxx/rho(i)   
         tauxy(i) = tauxy(i) + mass(j)*hxy/rho(j)
         tauxy(j) = tauxy(j) + mass(i)*hxy/rho(i)            
         tauyy(i) = tauyy(i) + mass(j)*hyy/rho(j)
         tauyy(j) = tauyy(j) + mass(i)*hyy/rho(i) 
                   
         current=>current%next
      ENDIF
     
   ENDDO
   
   current=> last            !    Initializes linked list to last element 
   DO WHILE (associated(current))
   
      i = current%pair_i 
      j = current%pair_j  

      IF ( current%Pint_type.eq.0.OR.current%Pint_type.eq.-1 ) THEN

         if (current%Pint_type.eq.0) then
            do idimn = 1,ndimn
               dvx(idimn) = x(idimn,j) - x(idimn,i)
            enddo      
         endif 
         
         hxy =      dvx(1)*current%dwdx(1)  
         wxy(i) = wxy(i) + mass(j)*hxy/rho(j)
         wxy(j) = wxy(j) + mass(i)*hxy/rho(i) 
                   
         current=>current%next
      ENDIF
     
   ENDDO
   
   do i = 1, npoin
      tauxy(i) = tauxy(i)/wxy(i)   !  ******
      tauxx(i) = kin_visc * tauxx(i)*rho(i)  ! dens*g*h/dens *h
      tauxy(i) = kin_visc * tauxy(i)*rho(i)  
      tauyy(i) = kin_visc * tauyy(i)*rho(i)
   enddo
 
   
! ******************

   do i = 1, npoin
      write(res_file2,10001) i, x(1,i), x(2,i), tauxy(i), tauxx(i), tauyy(i)
   enddo
   i = i
10001 format (i8,2x, 5(g11.4,2x))
! ******************   

!      ------  Calculate SPH sum for viscous force (eta Tab),b/ grav 

   
! ******************

   tauxy = 0.0
   tauxx = 0.0
   tauyy = 0.0
   do i = 1, npoin
      tauxx(i) =   x(1,i) 
   enddo
   indvxdt = 0.0

! ******************   


   current=> last          !    Initializes linked list to last element 
   
   DO WHILE (associated(current))

      i = current%pair_i 
      j = current%pair_j 
   
      IF ( current%Pint_type.eq.0.OR.current%Pint_type.eq.-1 ) THEN

         if (pa_sph.eq.1) then             !   ----  For SPH algorithm 1     
            rhoij = 1.e0/(rho(i)*rho(j))  
            do idimn = 1,ndimn
               h = 0.0
               if (idimn.eq.1) then        
                  h = h + ( tauxx(i) +  tauxx(j) ) * current%dwdx(1) 
                  h = h + ( tauxy(i) +  tauxy(j) ) * current%dwdx(2)
               elseif (idimn.eq.2) then                     
                  h = h + ( tauxy(i) +  tauxy(j) ) * current%dwdx(1) +  &
                          ( tauyy(i) +  tauyy(j) ) * current%dwdx(2)
               endif  
               h = h*rhoij
               indvxdt(idimn,i) = indvxdt(idimn,i) + mass(j)*h
               indvxdt(idimn,j) = indvxdt(idimn,j) - mass(i)*h
            enddo   
   
         else if (pa_sph.eq.2) then             !   ----  For SPH algorithm 2
            do idimn = 1,ndimn 
               h = 0.0 
               if (idimn.eq.1) then                            
                   h = h + ( tauxx(i)/rho(i)**2 +  tauxx(j)/rho(j)**2 )*current%dwdx(1)
                   h = h + ( tauxy(i)/rho(i)**2 +  tauxy(j)/rho(j)**2 )*current%dwdx(2)
               elseif (idimn.eq.2) then                     
                   h = h + ( tauxy(i)/rho(i)**2                    &
                         +   tauxy(j)/rho(j)**2 )*current%dwdx(1)  &
                         + ( tauyy(i)/rho(i)**2                    &  
                         +   tauyy(j)/rho(j)**2 )*current%dwdx(2)  
               endif 
               indvxdt(idimn,i) = indvxdt(idimn,i) + mass(j)*h
               indvxdt(idimn,j) = indvxdt(idimn,j) - mass(i)*h
            enddo
         endif

         current=>current%next

      ENDIF    
                     
   ENDDO
! ******************

   write(res_file2,*) ' -------------------------------------------'
   do i = 1, npoin
      write(res_file2,10001) i, x(1,i), x(2,i) , indvxdt (1,i), indvxdt (2,i) 
   enddo
   stop
! ******************  

   
ENDIF                                 !      -----  end of viscous terms computations

if (pa_sph.eq.4) then                 !      ----- For SPH  4  corr grad, alpha for convenience 1D ONLY   
    do i = 1,npoin
       indvxdt(:,i) = indvxdt(:,i)/alphaJB(i)
    enddo
endif 

if (ic_K0.eq.1) then
   deallocate (divvK0)
endif   

if (ic_histogram.eq.1.and.np_bakg_Total.gt.0) then
   indvxdt(:,1:np_bakg_Total) = 0.0
endif   
 
END SUBROUTINE   IntForces_SW_OLD 


!----------------------------------------------------------------------

      subroutine Trigger_FOSM_SW  

!----------------------------------------------------------------------
 
implicit none

integer (ink) ipoin, iSlope, ipoig, iproblem, icase, indx, idest
integer (ink) na,  nh,  ia , ih 

real    (irk) xx, yy, Zp, Zpx, Zpy, ZpGrad
real    (irk) Slope_max, Slope_min
Real    (irk) cohesion  
Real    (irk) tanfi, gtanfi , factpw
real    (irk) cosTh, sinTh, xden, xnum1, xnum21, xnum22, xnum23

real    (irk) h, hinv, pwp, pwprel , frictx 

real    (irk) deng, xtmp
real    (irk) cgra, cgra05              ! g and g/2  

integer (ink) ic_cosTheta               ! 
                                        ! ------------- DENSITIES ------------- 
real    (irk) dens, denss, densw, densd !   adde dens dry. This is the input !!!  
real    (irk) porosity, Sr              !   Porosity and degree of saturation
real    (irk) dens1, dens2              !   mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    !   average and effective average dens of 2 layers
real    (irk) alphaS                    !   ratio hsat/h, from hrelpw. If below 0.1 -> 0

integer (ink) ic_fct_dens               !   0  classic  1 dens correction 2 DFs 
real    (irk) fact_dens   ! not used    !   if 1 and icpwp=0 computed as  fact_dens = (dens-densw)/dens
                                        !           icpwp=  computed using hs/h, porosity...
                                        !   if 2  for DF   fact_dens = (denss-densw)/denss                                 
real    (irk) pi
real    (irk) PGA, pAi, phi
real    (irk) dZpgrad , FoS, Fos_old
real    (irk) FoS_mean, FoS_sigma, sigmai, dXi, FoS1, FoS2, dFdX, FoS_CoeffVar
real    (irk) PoF, aa, LgN_sigma, LgN_mean

pi        = 4.*atan(1.0)
iproblem  = nproblems_index             ! index of problem

!  ......TASk  1....... INITIALIZATION  .............................

ipoin      = 1
xx         = x (1,1) 
rho (1)    = const (75)  
h          = rho (1)
PGA        = const (76)
pwprel     = const (77)  !  new 28th april mp 2022
cgra       = const  (1)  !   we will use gravz, x &y if ic_cosTheta=1
                         !   cgra used for drag forces, is g and not gz
cohesion   = const(6)
tanfi      = const(9)

!      ------ TASk  2  densities ---------------------------------

densd         = const  (2)       !  important: this is the data user will provide
denss         = const (17)
densw         = const (18)
Sr            = const (19)
alphas        = const (14)
    
if ( (denss.LT.1.e-3).OR.(densw.lt.1.e-3) ) then         !  Filters zero solid or water densities  
   write (*,*) ' denss or densw wrong...', denss, densw
   write (*,*) ' strike a num key to acknowledge code will stop '
   read  (*,*) ipoin
   STOP
endif
    
porosity     = (denss-densd)/ denss                    !     new April 2022
dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
dens_bar_eff = dens_bar - densw*alphaS                 !     modified MP 10.01.2019 
fact_dens    = dens_bar_eff/dens_bar                   !     multiplies tanFi  NOT USED

!      ------ TASk  3  generate list of slopes Slope (N_Slopes)

N_Slopes      = const (71)
Slope_max_cut = const (72)
Slope_min_cut = const (73)
n_stats       = const (74)

if (.NOT.allocated(Slope)) then    ! allocate at first problem
   allocate (Slope (N_slopes) )                         ; Slope = 0.0
   allocate (FoS_Stats  (N_Slopes, N_stats*N_Gzones))   ; FoS_Stats = 0.0
   allocate (FoS_DEM    (npoig,    N_stats*N_Gzones))   ; FoS_DEM   = 0.0
endif

if (nproblems_index.eq.1) then  ! avoid unnecessary recalculation

   Slope_max = -999.
   Slope_min = 1.e3
   DO ipoig = 1, npoig
      xx = coorg(1,ipoig)
      yy = (yming+ymaxg)/2    
      if (ndimn.eq.2) then
          yy = coorg(2,ipoig)
      endif
      call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
      ZpGrad = (Zpx*Zpx + Zpy*Zpy)
      If (Zpgrad.gt.1.e-6) ZpGrad = ZpGrad**0.5
      if (Zpgrad.GT.Slope_Max_cut) Zpgrad = Slope_Max_cut
      if (Zpgrad.LT.Slope_Min_cut) Zpgrad = Slope_Min_cut
      if (Zpgrad.GT.Slope_max) slope_max = ZpGrad
      if (Zpgrad.LT.Slope_min) slope_min = ZpGrad
   ENDDO
   dZpgrad = (slope_max-slope_min)/(N_slopes)
   do islope = 1, N_slopes
      Slope(islope) = slope_min + (dZpgrad/2.0) + (islope-1)*dZpgrad
   enddo 
ENDIF 

!   TASK 4 ............  compute FoS for iproblem, islope ...............
!                        store it in FoS_Slope Ç(iproblem, islope)

DO islope = 1, N_Slopes

!      ------  Obtain gravities ------------------------------
        
   Zpgrad = Slope (islope)
   Zpx    = Zpgrad           ! This is a local x along steepes descent
   Zpy    = 0.0                                                            
       
   deng   = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
   gravz  =  cgra/deng
   gravx  = abs ( -cgra*Zpx/deng )                  !  is -gravz*Zpx indeed. Left for pedagogical..
   cgra05 = 0.5*gravz 
   cosTh  = gravz/cgra                              ! bug corrected MP 25 April 2022
   sinTh  = gravx/cgra

!      ------  Obtain friction ------------------------------ 
   
   h     = rho(1)
   if (ic_ws_interact.eq.2)  h = h_DF(ipoin)            !   CH MP 22 Feb 2019
      
   factpw = 0.0
   if (icpwp.GE.1) then           !  pwp is always the relative value at the bottom  
      pwprel = Const (77)
   else
      pwprel = 0.0
   endif 
   factpw = 1. - pwprel           !   new  we use relative values
  
   xden   = dens_bar * h * (gravx  + PGA*cosTh)
   xnum1  = cohesion
   xnum21 = dens_bar_eff * gravz * h 
   xnum22 = factpw - (PGA * dens_bar *sinTh)/( gravz * dens_bar_eff) !  MP 11 May 2022     
   xnum23 = tanFi

   Fos_old  = (xnum1 + xnum21*xnum22*xnum23)/xden

   if (ma_xtp.eq.0) then    !  correct loop index
       na = 1
   else
       na = ma_xtp
   endif
   if (mh_xtp.eq.0) then
       nh = 1
   else
       nh = mh_xtp
   endif

   DO ia = 1, na                       !  ---  Loop  Accelerations
      
      if (allocated (Acc_xtp)) then 
         PGA = Acc_xtp(ia)
      endif
      
      DO ih = 1, nh                    !  ---  Loop  h rel sat
         
         if (allocated (hrs_xtp)) then 
            alphas = hrs_xtp(ia)
         endif

         FoS_mean = Get_FOS (densd,  densw,    denss,   Sr, alphas, &    !      -------- FUNCTION Get_FOS
                              h, pwprel, cohesion, tanFi, PGA) 
         FoS_sigma= 0.0                                                  ! compute sum dfdx*sigma
         DO icase  = 1, ncases_mat
            indx   = CCases_mat (icase,1)
            idest  = (iproblem-1)*2 +1 
            Sigmai = FOSMms (idest+2, icase)
            dXi    = Sigmai / 10.
            if (indx.EQ.9) then                !  ---------------------tanFi
               Fos1 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi-dXi, PGA)
               Fos2 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi+dXi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.6) then                !  ------------------cohesion
               Fos1 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion-dXi, tanFi, PGA)
               Fos2 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion+dXi, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.2) then                !  ------------------densd
               Fos1 = Get_FOS (densd-dXi,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd+dXi,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.17) then                !  ------------------denss
               Fos1 = Get_FOS (densd,densw, denss-dXi, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd,densw, denss+dXi, Sr, alphas, &     
                            h, pwprel, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.18) then                !  ------------------densw
               Fos1 = Get_FOS (densd,densw-dXi, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd,densw+dXi, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.19) then                !  -------------------Sr
               Fos1 = Get_FOS (densd,densw, denss, Sr-dXi, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd,densw, denss, Sr+dXi, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.14) then                ! -------------------- alphas
               Fos1 = Get_FOS (densd,densw, denss, Sr, alphas-dXi, &     
                               h, pwprel, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd,densw, denss, Sr, alphas+dXi, &     
                               h, pwprel, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.75) then                ! -------------------- h
               Fos1 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h-dXi, pwprel, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h+dXi, pwprel, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.76) then                ! -------------------- PGA
               Fos1 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA-dXi)
               Fos2 = Get_FOS (densd,densw, denss, Sr, alphas, &     
                               h, pwprel, cohesion, tanFi, PGA+dXi)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            elseif (indx.EQ.77) then                ! --------------------pwp
               Fos1 = Get_FOS (densd,densw      ,    denss,    Sr, alphas, &     
                               h    , pwprel-dXi, cohesion, tanFi, PGA)
               Fos2 = Get_FOS (densd,      densw,    denss, Sr, alphas, &     
                               h    , pwprel+dXi, cohesion, tanFi, PGA)
               dFdX = (FoS2 - FoS1)/(2.0*dXi)
               FoS_sigma = FoS_sigma + (dFdX*sigmai)**2.
            endif
         ENDDO
         FoS_sigma = FoS_sigma**0.5
         if (FoS_mean.GT.1.e-6) then 
             FoS_CoeffVar = FoS_sigma / FoS_mean   
         else
             FoS_CoeffVar = 999.
         endif 
         idest     = (iproblem-1)*N_stats
                                                    ! Assume LogNormal MP 2 May 2022
         aa = 81./130.    !  alternatives 2/pi or 5/8
         LgN_sigma =  ( log (1.+(FoS_sigma/FoS_mean)**2.) )**0.5
         LgN_mean  =  log (Fos_mean) - 0.5*LgN_sigma*LgN_sigma
         xtmp      =  - LgN_mean / LgN_sigma
         if (xtmp.lt.0) then
              PoF =  1. - 0.5 * (  1. + sqrt (1.- exp (-aa* xtmp * xtmp ))  )
         else
              PoF =       0.5 * (  1. + sqrt (1.- exp (-aa* xtmp * xtmp ))  )
         endif

         FoS_stats (islope, idest+1 ) = PoF
         FoS_stats (islope, idest+2 ) = FoS_mean
         FoS_stats (islope, idest+3 ) = FoS_sigma
         FoS_stats (islope, idest+4 ) = FoS_CoeffVar
         FoS_stats (islope, idest+5 ) = (FoS_mean-1.)/FoS_sigma
         if (allocated(PAi_xtp)) then
            pAi = PAi_xtp(ia)
         else
            pAi = 1.0
         endif
         if (allocated(Phi_xtp)) then
            phi = Phi_xtp(ih)
         else
            phi = 1.0
         endif
         FoS_stats (islope, idest+6 ) = FoS_stats (islope, idest+6 ) +  PoF * PAi *Phi 
      ENDDO                      !  PGA
   ENDDO                         !  --- a lpha                       
ENDDO                            !  ----  islopes

contains
! --------------------------------------------------------------------------------------------

Function Get_FOS (densd0,densw0, denss0, Sr0, alphas0, &         !      -------- FUNCTION Get_FOS
                  h0, pwp0, cohesion0, tanFi0, PGA0) 
! --------------------------------------------------------------------------------------------
   implicit none

   real   (irk) Get_FoS
   real   (irk) densd0,densw0, denss0, Sr0, alphas0, h0, pwp0, cohesion0, tanFi0, PGA0
   real   (irk) porosity0, dens20, dens10, dens_bar0, dens_bar_eff0
   real   (irk) xnum10, xden0, xnum210, xnum220, xnum230
   real   (irk) factpw0, hinv0 
       
   porosity0     = (denss0-densd0)/ denss0                       !     new April 2022
   dens20        = (1.0-porosity0)*denss0 + porosity0*densw0*Sr0  !  -- dens of upper layer (can be unsat)
   dens10        = (1.0-porosity0)*denss0 + porosity0*densw0     !     lower layer (sat)   Sr=alphaS
   dens_bar0     = (1.0-alphaS0)*dens20 + alphaS0*dens10          !     average dens
   dens_bar_eff0 = dens_bar0 - densw0*alphaS0                     !     modified MP 10.01.2019 
      
   factpw0 = 0.0
   if (icpwp.eq.1.or.icpwp.eq.2.or.icpwp.eq.3) then     !  pwp is always the relative value at the bottom  
      factpw0 = 1. - pwp0                            !   pwp  *hinv  now we use relative values
   endif 
   
   xden0   = dens_bar0 * h0 * (gravx  + PGA0*cosTh)
   xnum10  = cohesion0
   xnum210 = dens_bar_eff0 * gravz * h0 
   xnum220 = factpw0 - (PGA0 * dens_bar0 *sinTh)/( gravz * dens_bar_eff0)  ! mp 11 may 2022     
   xnum230 = tanFi0

   Get_Fos = (xnum10 + xnum210*xnum220*xnum230)/xden0
   
end function Get_FoS


end subroutine Trigger_FOSM_SW 


!----------------------------------------------------------------------

      subroutine Trigger_FoS_SW  

!----------------------------------------------------------------------
 
implicit none

integer (ink) ipoin, iSlope, ipoig, iproblem

real    (irk) xx, yy, Zp, Zpx, Zpy, ZpGrad
real    (irk) Slope_max, Slope_min
Real    (irk) cohesion  
Real    (irk) tanfi, gtanfi , factpw
real    (irk) cosTh, sinTh, xden, xnum1, xnum21, xnum22, xnum23

real    (irk) h, hinv, pwp, pwprel, frictx 

real    (irk) deng
real    (irk) cgra, cgra05              ! g and g/2  

integer (ink) ic_cosTheta               ! 
                                        ! ------------- DENSITIES ------------- 
real    (irk) dens, denss, densw, densd !   adde dens dry. This is the input !!!  
real    (irk) porosity, Sr              !   Porosity and degree of saturation
real    (irk) dens1, dens2              !   mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    !   average and effective average dens of 2 layers
real    (irk) alphaS                    !   ratio hsat/h, from hrelpw. If below 0.1 -> 0

integer (ink) ic_fct_dens               !   0  classic  1 dens correction 2 DFs 
real    (irk) fact_dens   ! not used    !   if 1 and icpwp=0 computed as  fact_dens = (dens-densw)/dens
                                        !           icpwp=  computed using hs/h, porosity...
                                        !   if 2  for DF   fact_dens = (denss-densw)/denss                                 
real    (irk) pi
real    (irk) PGA
real    (irk) dZpgrad , FoS

pi        = 4.*atan(1.0)
iproblem  = nproblems_index             ! index of problem

!  ......TASk  1....... INITIALIZATION  .............................

N_Slopes      = const (71)
Slope_max_cut = const (72)
Slope_min_cut = const (73)
n_stats       = const (74)

if (.NOT.allocated(Slope)) then    ! allocate at first problem
   allocate (Slope (N_slopes) )                         ; Slope = 0.0
   allocate (FoS_Slope  (nproblems, N_Slopes))          ; FoS_Slope = 0.0
   allocate (FoS_Stats  (N_Slopes, N_stats*N_Gzones))   ; FoS_Stats = 0.0
   allocate (FoS_DEM    (npoig,    N_stats*N_Gzones))   ; FoS_DEM   = 0.0
endif

ipoin      = 1
rho (1)    = const (75)  
h          = rho (1)
xx         = x (1,1) 
PGA        = const (76)
pwprel     = const (77)  !  new 28th april mp 2022
cgra       = const(1)    !   we will use gravz, x &y if ic_cosTheta=1
                         !   cgra used for drag forces, is g and not gz
cohesion   = const(6)
tanfi      = const(9)

ic_cosTheta = const(16) + 0.01        ! 1 ... use g projected along normal to terrain

!      ------ TASk  2  densities ---------------------------------

ic_fct_dens   = const (20) + 0.1
dens          = const  (2)
densd         = const  (2)       !  important: this is the data user will provide
denss         = const (17)
densw         = const (18)
Sr            = const (19)
alphas        = const (14)
ic_hrelSat    = 0         !   const (34) + 0.01
    
if ( (denss.LT.1.e-3).OR.(densw.lt.1.e-3) ) then         !  Filters zero solid or water densities
   write (*,*) ' denss or densw wrong...', denss, densw
   write (*,*) ' strike a num key to acknowledge code will stop '
   read  (*,*) ipoin
   STOP
endif
    
! porosity     = (denss-dens)/(denss-densw)
porosity     = (denss-densd)/ denss                    !     new April 2022
dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
dens_bar_eff = dens_bar - densw*alphaS                 !     modified MP 10.01.2019 
fact_dens    = dens_bar_eff/dens_bar                   !     multiplies tanFi  NOT USED


!      ------ TASk  3  generate list of slopes Slope (N_Slopes)

if (nproblems_index.eq.1) then  ! avoid unnecessary recalculation

   Slope_max = -999.
   Slope_min = 1.e3
   DO ipoig = 1, npoig
      xx = coorg(1,ipoig)
      yy = (yming+ymaxg)/2    
      if (ndimn.eq.2) then
          yy = coorg(2,ipoig)
      endif
      call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
      ZpGrad = (Zpx*Zpx + Zpy*Zpy)
      If (Zpgrad.gt.1.e-6) ZpGrad = ZpGrad**0.5
      if (Zpgrad.GT.Slope_Max_cut) Zpgrad = Slope_Max_cut
      if (Zpgrad.LT.Slope_Min_cut) Zpgrad = Slope_Min_cut
      if (Zpgrad.GT.Slope_max) slope_max = ZpGrad
      if (Zpgrad.LT.Slope_min) slope_min = ZpGrad
   ENDDO
   dZpgrad = (slope_max-slope_min)/(N_slopes)
   do islope = 1, N_slopes
      Slope(islope) = slope_min + (dZpgrad/2.0) + (islope-1)*dZpgrad
   enddo 
ENDIF 

!   TASK 4 ............  compute FoS for iproblem, islope ...............
!                        store it in FoS_Slope Ç(iproblem, islope)


DO islope = 1, N_Slopes

!      ------  Obtain gravities ------------------------------
        
   Zpgrad = Slope (islope)
   Zpx    = Zpgrad           ! This is a local x along steepes descent
   Zpy    = 0.0                                                            
       
   deng   = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
   gravz  =  cgra/deng
   gravx  = abs ( -cgra*Zpx/deng )                  !  is -gravz*Zpx indeed. Left for pedagogical..
   cgra05 = 0.5*gravz 
   cosTh  = gravz/cgra                              ! bug corrected MP 25 April 2022
   sinTh  = gravx/cgra

!      ------  Obtain friction ------------------------------ 
   
   h     = rho(1)
   if (ic_ws_interact.eq.2)  h = h_DF(ipoin)     !   CH MP 22 Feb 2019
      
   factpw = 0.0                   ! changed MP 28th April 2022  We use C77 = peprel
   if (icpwp.GE.1) then           !  pwp is always the relative value at the bottom  
      pwprel = Const (77)
   else
      pwprel = 0.0
   endif 
   factpw = 1. - pwprel           !   new  we use relative values
   
   xden   = dens_bar * h * (gravx  + PGA*cosTh)
   xnum1  = cohesion
   xnum21 = dens_bar_eff * gravz * h 
   xnum22 = factpw - (PGA * dens_bar *sinTh)/( gravz * dens_bar_eff)   ! factpw    
   xnum23 = tanFi

   Fos    = (xnum1 + xnum21*xnum22*xnum23)/xden

   FoS_Slope (iproblem, islope) = FoS             ! Store in FoS_Slope matrix

ENDDO   ! islopes

END subroutine  Trigger_FoS_SW  



! -----------------------------------------------------------------
 
       SUBROUTINE ExtForces_SW (ic_RKiter)        
 
! -----------------------------------------------------------------

!   SLOW version. Uses Get_basal_VB_Ext_SW Jan 2023
!
!   Subroutine to calculate the external forces: Slope and Bottom friction
!
!  changed for topo dependent properties 21st march 2014

!   Added also virtual parts... 
!   02 July Voellmy uses xi->cmanning tau = frict + rho*g*V*V/xi
!           To make it work use a frictional material plus manning not zero   
!   27 July: nfrict=3 for Bham, properties depend on mass fraction of soil
!            tauy = tauy0*(1-exp(-Cs*s)) where
!                   Cs is a coefficient stored in position Manning const(3)
!                   s = (mass-mass0)/mass= (mass soil)/(mass soil + mass water)
!
!   2016 Dec: Bug in Bingham, it was not 
!             Added HB 23, and 24+25 simple frict+coh+viscous laws (dv/dz)^m
!             New routines for slope and basal friction. Moved from main to
!                          contained routines


!  ......................... DECLARATIONS .............................

implicit none

integer(ink) ic_RKiter     ! argument. If 0, RK loop, 1, after it

integer(ink) ipoigx, ipoigy, ipoig
integer(ink) ipoin,idimn, i,j
integer(ink) inx0, inx1, inx
integer(ink) iny0, iny1, iny
real   (irk) distx2, disty2, dist, auxx, tmp_mass, fact_s2w

real   (irk), SAVE, allocatable:: auxv(:,:)     ! MP CH 30th May 2019
real   (irk), allocatable:: tmp_extf(:)         ! MP 1 may 18, to use ic_end_solid option

! integer(ink) nfric  not local 7th dec 2017

real   (irk) xx, yy, Zp, Zpx, Zpy, ZpGrad
real   (irk) top4, top5, top6, top7, top8, top9

Real   (irk) cmanning, xi_voellmy, xi_Voellmy_expf  !  MP July 2023
real   (irk) tauy0, constK, dens_ws  
real   (irk) visco, tanfi8, hfrict0, c73, tanfi0, Bfact, Cv, comp  
Real   (irk) tanfi, gtanfi, gtanfi8
real   (irk) c8_Muf, C6_n1, C7_n2, C67_N, fact9               !  ---  Perzyna
real   (irk) Cm
real   (irk) nsoilB, csoilB, constKB, tauy0B, densB, fsoilB   !  ---  Constants for Bingham + erosion
real   (irk) K_Bagnold, Mu_CF, constMu_CF                     !          new cf model  
real   (irk) m_vexp, tb_abs

real   (irk) h, vxx, vyy, hu, hv, hinv, vmod,hvmod,  hvmodinv, pwp, frictx, fricty, frict_abs
real   (irk) cf1, cf2a, a_hat, a2, b2, c2, xh, Tau0,cf3a, cf3b,   cf4b
real   (irk) cf5, factpw, cv2byR, cosa, sina, cos2, sin2, cbs, xnum, xeden, rinv, factaux
real   (irk) cf8, hrelpw, hrelpw2, hfrict, xden
real   (irk) erosion_base, erosion_rate 
real   (irk) sw_erosion_rate
real   (irk) tmp_avail_eros

real   (irk) frictxx, frictyy
real   (irk) tau_ref            !  we use it in new Bingham corrected routine  MP 19 feb 2008
real   (irk) deng
real   (irk) ZZmin, ZZmax, ZZdiscr, Tandelta
integer(ink) itopo_zones, iexit

real   (irk) velx, vely, velinv 
real   (irk) tempx, tempy, temp, vcomp 
real   (irk) sinth, costh
integer(ink) iitime
                                       !  ------  This is for virtual particles, just in case
real   (irk) rr, f, rr0, dd, p1, p2
real   (irk), allocatable:: dx(:)                               !      

real    (irk) C4                        ! old for  erosion constant (either topo or constant)
real    (irk) hhsZ                      ! stands for h + hs + Z

real    (irk) cgra, cgra05              ! g and g/2  

integer (ink) ic_cosTheta               ! C16..1  use g projected along normal to terrain 

                                        ! ------------- DENSITIES ------------- 
                                        !
real    (irk) dens, denss, densw        !   ic_fct_dens = C20  0 no correction  1.. correct 
real    (irk) porosity, Sr              !   Porosity and degree of saturation
real    (irk) dens1, dens2              !   mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    !   average and effective average dens of 2 layers
real    (irk) fact_pwp_max,fact_pwp_min !   f  actors limiting (1-dpwp/h) in pwp computa
real    (irk) fpwp_max, fpwp_min        !   fact. limiting pwp:  h*fact
real    (irk) alphaS                    !   ratio hsat/h, from hrelpw. If below 0.1 -> 0

integer (ink) ic_fct_dens               !   0  classic  1 dens correction 2 DFs 
real    (irk) fact_dens                 !   if 1 and icpwp=0 computed as  fact_dens = (dens-densw)/dens
                                        !           icpwp=  computed using hs/h, porosity...
                                        !   if 2  for DF   fact_dens = (denss-densw)/denss 
real    (irk) betaw                     !   auxiliar for d Pwp                                   
 
real    (irk) V_T, expm                 !  constants for sw law Stored in C27 and C28 
integer (ink) darcinian_law             !     1 for Darcy, 2 for new  C26

real    (irk) alpha_slip, slip_fact     ! change MP 6 april 2018. Vbottom = alfa*Vtop  = 2a/(1+a) *V averaged
                               
real    (irk) pi
integer (ink) icm
real    (irk) hn, hns, hnw, dhns, dhnw, n0b 

pi        = 4.*atan(1.0)
c73       = 7./3.

if(.NOT.allocated(dx) ) allocate (dx(ndimn))

!  ......................... INITIALIZATION  .............................

cgra       = const(1)    !   we will use gravz, x &y if ic_cosTheta=1
                         !   cgra used for drag forces, is g and not gz
dens       = const(2)
dens_ws    = 1000./dens  !  Default used for interaction type 1 so far

!   nfric      = const(5) + 0.01 7th dec 2017 MP is a global var, defined in input_sw

cmanning   = const(3)   !   C3 controls several rheological laws

xi_voellmy = const(3)   !   default : Manning    June 06 2016

if (ic_end_solid.eq.1) then        !  solid behaviour via decreasing p 
    if (law_end_solid.eq.0) then   !  MP 30th Jan 2024
       fact_s2w = 0.0  
       if (time_sph.GT.t1_end_solid) then
           fact_s2w = 1.0
       endif
    elseif ( law_end_solid.eq.1 ) then
       fact_s2w = 1.0
       if  (time_sph.lt.t1_end_solid) then
           fact_s2w = (time_sph/t1_end_solid)**fct_exp_end_solid 
       endif
    elseif( law_end_solid.eq.2) then
       fact_s2w = 0.0
       if  (time_sph.gt.t1_end_solid) then
          fact_s2w = 1.0 - exp( -fct_exp_end_solid *  (time_sph-t1_end_solid) )
       endif 
   endif       
    if (const (41).gt.1.e-6) then
         xi_Voellmy =  const(3) + (const(41)-const(3)) * fact_s2w
    endif
endif

 if (time_sph.GT.t1_end_solid) then   !  MP andrei rosas case 25 Jan 2024
    if (const (41).gt.1.e-6) then
       xi_voellmy = const(41)
    endif
 endif

if (xi_Voellmy.LT.0) then
    constMu_CF = -(25.*const(3)/4.)
    Mu_CF      = - const(3)
endif

if (ic_ws_Interact.eq.2) then         ! *** DFs cannot use Voellmy 'cause is for water
   xi_Voellmy = const(12)              !     MP 24 may 18
   if (xi_Voellmy.LT.0) then
       constMu_CF = -(25.*const(12)/4.)
       Mu_CF      = - const(12)
   endif
endif
K_Bagnold  = const(3)   ! if negative, is the C-F model Mu_CF 
                        !  Mu_CF      = const(3)   Nov 2016
tauy0      = const(6)
constK     = const(7)
visco      = const(8)
tanfi8     = const(9)
hfrict0    = const(10)
tanfi0     = const(12)
Bfact      = const(13)
comp       = const(15)
C6_n1      = const(6)  
C7_n2      = const(7)  
C67_N      = C6_n1 + C7_n2
C8_Muf     = const(8)

ic_cosTheta = const(16) + 0.01        ! 1 ... use g projected along normal to terrain

!      ------  density factor and pwp 
!              Obtains (i) all densities and fact_dens (ii) pwp factors (iii) DF vars

call DFs_FSs_dens_and_pwp_factors     !  MP March 2023
!
!     The results are :   
!          (i) accesible from everywhere using module SW Vars
!              fact_dens_DF(:), 
!              dens_DF (:), densw_DF(:), dens1_DF(:), dens2_DF(:), dens_bar_DF(:), dens_bar_eff_DF(:)
!              fact_pwp_max_DF (:), fact_pwp_min_DF (:), fpwp_max_DF(:), fpwp_min_DF(:)
!              fact_dens_DF(:),
!         (ii) accesible from extforces only (vars declared here, sub contained in)
!              dens, densw, dens1, dens2, dens_bar, dens_bar_eff
!              fact_pwp_max, fact_pwp_min, fpwp_max, fpwp_min

Cv = 4*Bfact/(pi*pi)                  !   *** Pwp previous
hrelpw2  = const(14)*const(14) 

IF (ic_RKiter.eq.0)  then
    exdvxdt = 0.0    !!!  ** changed 02 July 2014
    du      = 0.0    !!!             15 July
ENDIF    

If ((const(21).gt.1.e-9).or.(ic_basal_erosion.eq.1)) then    ! 1st condition body erosion
    if (.not.allocated (auxv)) then                          ! 2nd topo dependent erosion
        allocate (auxv(2,npoig))                             ! auxv updates topo
!        if (topol(15,1).le.1.e-6) then                      ! it will be done only once
!            topol(15,:) = 999.0   
!        endif 
!        if ( ALL(topol(15,:)< 1.e-6) ) then                !  MP CH 27 May 2019 
!            topol(15,:) = 999.0                            ! deactivated 21st June Moved to input
!        endif            
    endif
    auxv = 0.0 
endif

!      ------   interaction water soil (when called from RK loop, 0)

if (ic_ws_Interact.eq.1.and.ic_RKiter.eq.0) then     
   current=>last
   do while (associated(current))
      i = current%pair_i
      j = current%pair_j
      if (current%Pint_type.eq.1) then
         do idimn = 1,ndimn                     ! will never pass trough here CHK
                                                ! note: cgra-> gravz         
            if    (ic_SwAlg.eq.0) then          ! we do not use the stabilized formulation for water
                                                ! soil chk costheta
                exdvxdt(idimn,i)  = exdvxdt(idimn,i)  - dens_ws*cgra*mass(j)*current%dwdx(idimn)   
                exdvxdt(idimn,j)  = exdvxdt(idimn,j)  +         cgra*mass(i)*current%dwdx(idimn)  ! w 
                
            elseif (ic_SwAlg.eq.1) then         ! stab for water, included in ext forces before here
            
                exdvxdt(idimn,i)  = exdvxdt(idimn,i)  - dens_ws*cgra*mass(j)*current%dwdx(idimn)  !   soil
                
            endif
         enddo
      endif
   current=>current%next
   enddo

   if (ic_ws_erosion.eq.1) then   ! interfase erosion w-s compute hs at w, hw at s, vws, h_df
      call Get_s_at_w_wir_SW
   endif
endif

Call Check_Out_Domain_SW     ! (a2) comment: not necessary, was done in RK4 after updating nodes)

DO ipoin = 1,npoin  ! =================================  LOOP in POINTS =========
    
!      ------   Obtain nfric and
!                if_frictional  (= 1 for frictional laws: 5 to 10, 26 and 27
!                if_CV2byR      ( =1  for 7,8,10, 26 and 27)
    
    nfric = const(5) + 0.01
    if_frictional = 0
    if_V2byR      = 0
    if ( ((nfric.ge.5).AND.(nfric.le.10)).OR.(nfric.eq.26).OR.(nfric.eq.27) ) then
        if_frictional = 1
        if ( (nfric.eq.7).OR.(nfric.eq.8).OR.(nfric.eq.10) &
             .OR.(nfric.eq.26).OR. (nfric.eq.27)) then
             if_V2byR = 1
        endif
   
   elseif (nfric.eq.202.OR.nfric.eq.203) then   !  MP 2023 Nov new Bagnold model
        if_frictional = 1; if_V2byR = 1
        xi_Voellmy      = const (102)
        xi_Voellmy_expf = const (103)
   endif  

   if (ic_hrelSat.EQ.1) then                     ! use factors from DFs_dens_and_pwp_factors
        fact_pwp_max = fact_pwp_max_DF (ipoin)   ! needed in basal friction subroutine
        fact_pwp_min = fact_pwp_min_DF (ipoin)   ! vars defined here, dub included in ext forces
        fpwp_max     = fpwp_max_DF(ipoin)
        fpwp_min     = fpwp_min_DF(ipoin)       
        fact_dens    = fact_dens_DF(ipoin)
   endif

   if ( if_Out_Domain(ipoin).eq.1) CYCLE  
   
   xx  = x(1,ipoin)            !  fixing bug in Hungr: before the rk if 
   yy  = (yming+ymaxg)/2     
   if (ndimn.eq.2) yy = x(2,ipoin)
   
   vxx = vx(1,ipoin)
   vyy = 0.0
   if (ndimn.eq.2) vyy = vx(2,ipoin) 

   h     = rho(ipoin)
   if (ic_ws_interact.eq.2)  h = h_DF(ipoin)     !   CH MP 22 Feb 2019
   hu    = h*vx(1,ipoin)
   hv    = 0.0
   if (ndimn.eq.2) hv = h*vx(2,ipoin)

   if (h.lt.comp) then   !  ** see here comp which limits 1/h
       h        = comp
       hu       = 0.0
       hv       = 0.0
       hinv     = 0.0
       vmod     = 0.0
       hvmod    = 0.0
       hvmodinv = 0.0
   else
       hinv     = 1/h
       hvmod    = sqrt(hu*hu+hv*hv)
       vmod     = hvmod*hinv
       hvmod    = sqrt(hu*hu+hv*hv)
       if (hvmod.ge.1.e-5) then
           hvmodinv = 1./hvmod
       else
           hvmodinv = 0.0
       endif 
   endif
                                                       
   IF (ic_RKiter.eq.0) THEN    !  *** only for extforces **                                                                                               
                                   
       if ( ic_ws_Interact.eq.1.AND.itype(ipoin).eq.2) THEN     ! ----  WIR Interact., type 2 is water  
           dens_ws = densw/dens     
           dens    = densw            
           nfric   = 1               
           if_frictional = 0                                    ! Avoid water being frictional
           if_V2byR      = 0                                        
       elseif  (ic_ws_Interact.eq.2.AND.itype(ipoin).eq.2) THEN ! ---- DFs  
           dens_ws  = densw/dens 
           dens     = densw
           nfric    = 1               
           if_frictional = 0                                    ! Avoid water being fictional
           if_V2byR      = 0                           
       else                                   
           dens   = const(2)                                    !  otherwise: no interaction
           nfric  = const(5) + 0.01                             !    use dens and nfric  
       endif                             
   
       Call Get_Grav_Ext_SW		       !  Obtain gravity, slope and forces 
                 
!      ------    Friction OR VB basal in slow

       If (ic_Slow.eq.1) THEN ! changed MP May 2023 then(SPH_t_Integ_Alg.eq.10) then
          if  (ic_hrelSat.eq.1.AND.itype(ipoin).eq.1) THEN  ! visco only forsolid
             Call Get_basal_VB_Ext_SW      ! in slow landslides obtains VB using ex+int 
          endif
       else
          Call basal_friction_Ext_SW
       endif
                 
!      ------    interfase Friction

     if ( ic_ws_Interact.eq.1.AND.ic_ws_friction.EQ.1) THEN   ! Interaction
        call get_ws_interfase_friction                        ! exdvxdt computed inside subroutine
        if (itype(ipoin).eq.1) then
            exdvxdt(1,ipoin) = frictx/(dens*rho(ipoin))
            if (ndimn.eq.2) then
               exdvxdt(2,ipoin) = fricty/(dens*rho(ipoin))
            endif
        elseif (itype(ipoin).eq.2) then
            exdvxdt(1,ipoin) = frictx/(densw*rho(ipoin))
            if (ndimn.eq.2) then
               exdvxdt(2,ipoin) = fricty/(densw*rho(ipoin))
            endif
        endif
     endif 

!      ------  PWP contribution    
 
     if (icpwp.eq.1) then                                 !      ------ PWP
!        du(ipoin) = -Bfact*hinv*hinv*u(ipoin)/hrelpw2 
         betaw     = Bfact*hinv*hinv/hrelpw2 
         du(ipoin) = - (u(ipoin)/dt_sph)*(1-exp(-betaw*dt_sph))
     endif    

     If (SPH_t_Integ_Alg.NE.10) then
        exdvxdt(1,ipoin) = exdvxdt(1,ipoin) - frictx/h     
        if (ndimn.eq.2) exdvxdt(2,ipoin) = exdvxdt(2,ipoin) - fricty/h  
     endif
     iitime = time/dt_sph

1000 format  (5(2x,g12.5)) 
 
!      ------  Erosion Z dependent March 2014 April 2018

     If (Ic_erosion.eq.1) then  !!!   MP 21st June 2019 ((const(21).gt.1.e-9).or.(ic_basal_erosion.eq.1)) then

        erosion_rate = 0.0
        call Get_Erosion_Rate_SW
        if (ic_Stef_eros.eq.1) then 
            erosion_rate_V(ipoin) = erosion_rate   ! MP march 2021 
        endif

        if (nconst.gt.45) then  ! change MP 6 april 2018. Vbottom = alfa*Vtop  = 2a/(1+a) *V averaged
           alpha_slip = 0.0
           slip_fact  = 1.0
           if (const(47).gt.1.e-6) then
              alpha_slip = const(47)
              slip_fact  = 2*alpha_slip/(1.+alpha_slip)
           endif
        endif

        if (ic_ws_interact.eq.0) then
            do idimn = 1,ndimn                  
               exdvxdt(idimn,ipoin) = exdvxdt(idimn,ipoin) - erosion_rate * vx(idimn,ipoin)*hinv*slip_fact 
            enddo
        elseif (ic_ws_interact.eq.2) then    !  DFs  
            do idimn = 1,ndimn              ! change MP 6 april 2018. Vbottom = alfa*Vtop  = 2a/(1+a) *V averaged
               exdvxdt(idimn,ipoin) = exdvxdt(idimn,ipoin) - erosion_rate * vx(idimn,ipoin)*hinv*slip_fact 
            enddo
        elseif (ic_ws_interact.eq.1) then    !  WIRs  only basal erosion for soil
            if (itype(ipoin).eq.1) then
               do idimn = 1,ndimn               !  
                  exdvxdt(idimn,ipoin) = exdvxdt(idimn,ipoin) - erosion_rate * vx(idimn,ipoin)*hinv*slip_fact
               enddo
            endif
        endif
              
     endif   ! basal erosiom 
!   
!    Interface erosion   -----------------------------
!
     if (ic_ws_interact.eq.1.AND.ic_ws_erosion.eq.1) then
        sw_erosion_rate = 0.0
        call  Get_sw_interface_erosion_SW
        if (itype(ipoin).eq.1) then
            do idimn = 1,ndimn              
               exdvxdt(idimn,ipoin) = exdvxdt(idimn,ipoin) &
                     - sw_erosion_rate * ( vx(1,ipoin)-v_sw(1,ipoin) )*hinv
            enddo
        elseif (itype(ipoin).eq.2) then
            do idimn = 1,ndimn              
               exdvxdt(idimn,ipoin) = exdvxdt(idimn,ipoin) &
                      + sw_erosion_rate * ( vx(1,ipoin)-v_sw(1,ipoin) )*hinv
            enddo
        endif
     endif
          
   ELSEIF (ic_RKiter.eq.1) then   
          
      if (ic_crush.eq.1) then
          Call basal_friction_Ext_SW     
          call Get_r_crush_SW          !  January 18  crushing grains
      endif       
      
      If  (Ic_erosion.eq.1) then       !!!   MP 21st June 2019 ((const(21).gt.1.e-9).or.(ic_basal_erosion.eq.1)) then

          Call Get_Grav_Ext_SW            !  we need slopes and frictx
          Call basal_friction_Ext_SW      !  and frictx.. in erosion subr  

          erosion_rate = 0.0
          call Get_Erosion_Rate_SW        ! changed January 2016

          erosion_base = dt_sph*erosion_rate 
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
             auxv(1,ipoig) = auxv(1,ipoig) + erosion_base/dist  
             auxv(2,ipoig) = auxv(2,ipoig) + 1.0/dist 
          enddo
          enddo

          icm = 1
          if (ic_ws_interact.eq.1.AND.itype(ipoin).eq.2) icm = 0    !  do not do for water in WIRs
          if (icm.eq.1) then
             mass(ipoin)  = mass(ipoin) + dt_sph * erosion_rate * mass(ipoin)/rho(ipoin)
             if (npoin_s.eq.1.OR.SPH_t_integ_Alg.eq.110.OR.SPH_t_integ_Alg.eq.100)  then   ! MP April 2021 running PLK01 1 block erosion
                 rho (ipoin)  = rho(ipoin)  + dt_sph * erosion_rate                        ! and August 2022 
             endif 
          endif

      endif 
      
      if  (ic_ws_interact.eq.1.AND.ic_ws_erosion.eq.1) then !   ***
         sw_erosion_rate = 0.0
         call Get_sw_interface_erosion_SW
         if (itype(ipoin).eq.2) then      ! ---water
            mass(ipoin) = mass(ipoin) - dt_sph * sw_erosion_rate*mass(ipoin)/rho(ipoin)
            rho (ipoin) = rho(ipoin)  + dt_sph * sw_erosion_rate
         elseif (itype(ipoin).eq.1) then  ! --- soil
            hn  = rho(ipoin)
            hns = rho(ipoin)*(1-n_wir(ipoin))
            hnw = rho(ipoin)* n_wir(ipoin)
            n0b = const(46)            ! this is base porosity
            dhns= dt_sph * (1-n0b)*erosion_rate  
            dhnw= dt_sph * (n0b *erosion_rate + sw_erosion_rate)
            hns = hns + dhns
            hnw = hnw + dhnw
            n_wir(ipoin)     = hnw/(hnw+hns)
            dens_wir (ipoin) = denss*(1-n_wir(ipoin)) + densw*n_wir(ipoin)
            mass(ipoin) = mass(ipoin) + (dhns + dhnw)*mass(ipoin)/rho(ipoin)
            rho (ipoin) = rho(ipoin)  + dt_sph * (hns + hnw)
         endif
      endif

   ENDIF   !   ***** ic_RKiter cases   
  
 
ENDDO   ! ==== here ends the loop in sph nodes ====================================    

if (ic_RKiter.eq.1) then

!                   ------  If erosion,  decrease Z in topol  ** UPDATED MARCH 2014
!                           when ic_RKiter=1

   If (Ic_erosion.eq.1) then       !!!   MP 21st June 2019 ((const(21).gt.1.e-9).or.(ic_basal_erosion.eq.1)) then
      Do ipoig =1,npoig
         if (auxv(2,ipoig).gt.1.e-6) then 
            xnum            = auxv(1,ipoig)/auxv(2,ipoig) 
            if ( xnum .gt. topol(15,ipoig))  xnum = topol(15,ipoig)
            if (topol(15,ipoig).gt.0.) then
                topol(14,ipoig) = topol(14,ipoig) +  xnum 
                topol (1,ipoig) = topol (1,ipoig) -  xnum
                topol(15,ipoig) = topol(15,ipoig) -  xnum  ! changed January 2016 
            endif
         endif
      enddo
   endif 

!                   ------  If Darcinian interaction...
   
   darcinian_law = const(26) + 0.001
   
   if (ic_ws_Interact.eq.2.and.darcinian_law.gt.0 ) then   ! MP 8th march 2023 !!! 
      if (ic_slow.eq.1.AND.ic_DF.EQ.1) then 
         call Get_s_at_w_DF_SW
         do ipoin = npoin_s + 1, npoin_s + npoin_w
            vx(:,ipoin) = v_sw(:,ipoin)
         enddo 
      else
         call Get_Drag_SW
      endif
   endif    
   
!   do ipoin = npoin_s+1, npoin                                  !  WORK19  WE DECREASE hw only
!      fact_s2w    = rho(ipoin)
!      rho (ipoin) = rho(ipoin) - dt_sph*(2.*9.8*rho(ipoin))**0.5
!      if(rho(ipoin).LT.const(15)) rho(ipoin) = const(15)
!      mass(ipoin) = mass(ipoin)*rho(ipoin)/fact_s2w
!   enddo 

elseif (ic_RKiter.eq.0) then  ! ---------------  
 !  darcinian_law = const(26) + 0.001   
 !  if (ic_ws_Interact.eq.2.and.darcinian_law.gt.0.and.SPH_t_integ_Alg.eq.10  ) then   ! *** MP March 2023 slow 
 !      call Get_s_at_w_DF_SW 
 !      call Get_Drag_SW
 !  endif   
   
   if (ic_ws_Interact.eq.2.and.darcinian_law.gt.0 ) then   ! MP 8th march 2023 !!! 
      if (ic_slow.eq.1.AND.ic_DF.EQ.1) then 
         call Get_s_at_w_DF_SW
         do ipoin = npoin_s + 1, npoin_s + npoin_w
            vx(:,ipoin) = v_sw(:,ipoin)
         enddo 
      endif
   endif    

   if (ic_end_solid.eq.1) then            !     MP 18 may 2019  new routine in TOOLs module
       call solid_slab
   endif

   call get_virtual_parts_EXT_SW  

ENDIF

! exdvxdt(:, 1:npoin_s) = 0.0 !  WORK19 STOPS SOIL

CONTAINS 

! -------------------------------------------------------------------

       Subroutine Get_basal_VB_Ext_SW
       
!-------------------------------------------------------------------

!      ----   Slow 

implicit none

real    (irk), allocatable:: Taubx(:), nb(:), Vbx(:)
real    (irk) Taub_mod, Vb_mod, xnum, xden, frictmp, Tandelta 
real    (irk) tanFi, cohesion, n_exp, B_vb, Ru, sB, B_temp, h_shear

integer(ink) idimn 

if (.NOT.allocated (Taubx)) then
   allocate (Taubx(ndimn), nb(ndimn), Vbx (ndimn))
endif

if (icpwp.eq.1) then         !      ------  PWP computation
    pwp = u(ipoin)
elseif (icpwp.eq.2) then
    pwp = Uaux(1,ipoin)      !      ------  point 1 is at the bottom
elseif (icpwp.eq.3) then
    pwp = Uaux(1,ipoin)      !      ------  The only point 1 is at the bottom
else
    pwp = 0.0
endif 

if ( ts_Rain % type_Val.EQ.2  ) then   ! added JMS 14th March 2024
   pwp =  ts_Rain %  pwp_presc_basal
endif

hfrict   = max(h,hfrict0)   ! ** hfrict  limits friction 1/h

nfric = const(5) + 0.01

Taubx (:) = (exdvxdt(:,ipoin) + indvxdt(:,ipoin))*h*dens_bar  ! -----modify for DFs 2 phases
Taub_mod  = Taubx (1) * Taubx (1)
if (ndimn.eq.2) Taub_mod  = Taub_mod + Taubx (2) * Taubx (2)
Taub_mod  = Taub_mod**0.5

if (Taub_mod.GT.1.e-5) then

   nb (1) = Taubx (1)/Taub_mod
   if (ndimn.eq.2)  nb (2) = Taubx (2)/Taub_mod

   B_vb     = const (8)
   cohesion = const (6)
   TanFi    = const (9)
                          
   if (ic_basal_friction.eq.1) then               !!  added to modify friction in the slow model too   
      call Get_Topo_Basal_Friction_SW (Tandelta)  !  JwSM 17th  march 2024
        if (Tandelta.ge.0) then
            tanfi = Tandelta
        endif
    endif  

  !   if (x(1,ipoin).lt.37.) TanFi = 2*TanFi   !!!!!   ONLY FOR TEGUS G2  MP 17th march 2024
   n_exp    = const (7)
   h_shear  = const(12)
   Ru       = (1. - hinv*min(pwp,h*0.9))    !   15th march  jw SM
   ! frictmp  = tanFi*Ru*fact_dens
   frictmp  = tanFi*Ru*dens_bar_eff_DF(ipoin)/dens_bar_DF(ipoin) !*** temporal
   ! sB     = cohesion + frictmp * dens_bar * gravz * h
   sB       = cohesion + dens_bar_eff_DF(ipoin) * gravz * h * tanFi * Ru 
   if (nfric.eq.101) then
      Vb_mod   = B_vb * (Taub_mod / sB)** n_exp
   elseif (nfric.eq.102) then
      if (taub_mod.gt.sB) then
         Vb_mod = B_vb * ( (Taub_mod -sB) /sB)** n_exp 
      else
         Vb_mod = 0.0
      endif
   elseif (nfric.eq.103) then
      B_temp = (b_vb/h_shear) * min(h_shear,h) 
      Vb_mod   = B_vb * (Taub_mod / sB)** n_exp
   elseif (nfric.eq.104) then
      if (taub_mod.gt.sB) then
         B_temp = (b_vb/h_shear) * min(h_shear,h)
         Vb_mod = B_temp * ( (Taub_mod -sB) /sB)** n_exp 
      else
         Vb_mod = 0.0
      endif
   endif

   Vx (1,ipoin) = Vb_mod*nb(1)
   if (ndimn.eq.2) Vx (2,ipoin) = Vb_mod*nb(2)

    

ELSE
   Vx (:,ipoin) = 0.0
ENDIF

end subroutine Get_basal_VB_Ext_SW

!-------------------------------------------------------------------

       Subroutine DFs_FSs_dens_and_pwp_factors     !  MP March 2023
       
!-------------------------------------------------------------------

!      ----   

implicit none

integer (ink) ipoin

! real   (irk) dens, denss, densw, Sr, alphas
! real   (irk) porosity, dens1, dens2, dens_bar, dens_bar_eff
! real   (irk) fact_dens
! real   (irk) fact_pwp_max, fact_pwp_min, fpwp_max, fpwp_min

!      ------  Default values     

ic_fct_dens   = const (20) + 0.1
dens          = const  (2)
denss         = const (17)
densw         = const (18)
Sr            = const (19)
alphas        = const (14)
ic_hrelSat    = const (34) + 0.01
fact_pwp_max  = 1.0
fact_pwp_min  = 0.0
fpwp_max      = 1.0
fpwp_min      = 0.0

if ( (ic_hrelSat.EQ.1).AND.(.NOT.allocated (fact_dens_DF)) ) then
      allocate  ( fact_dens_DF(npoin) )   
      allocate  ( dens1_DF(npoin), dens2_DF   (npoin), densw_DF       (npoin) ) 
      allocate  ( dens_DF (npoin), dens_bar_DF(npoin), dens_bar_eff_DF(npoin) )
      allocate  ( fact_pwp_max_DF (npoin), fact_pwp_min_DF (npoin) )
      allocate  ( fpwp_max_DF     (npoin), fpwp_min_DF     (npoin) )
      fact_dens_DF = 0.0; dens1_DF=0.0; dens2_DF=0.0; densw_DF=0.0
      dens_DF = 0.0; dens_bar_DF=0.0; dens_bar_eff_DF=0.0
      fact_pwp_max_DF= 0.0; fact_pwp_min_DF=0.0; fpwp_max_DF=0.0;fpwp_min_DF=0.0
endif

IF (ic_fct_dens.EQ.0) THEN                           !  1 phase, w or w/o pwp, correction in tanfi

    if ( densw.lt.1.e-3) densw = 1000.0
    dens_ws = densw/dens
    if ( alphas.LT.0.01) then                         !  means alphas ignored, so it is assumed to be 1
        alphas  = 1.0
    endif
    fact_dens    = 1.0                                        ! Default
    dens_bar     = dens                                       !  needed for slow  Jan 2023
    dens_bar_eff = dens - densw

ELSEIF (ic_fct_dens.EQ.1) then                             !  1 phase, w or w/o pwp, correction OK
      
    IF (ic_hrelSat.eq.0) then       !  all dens  constants to compute factors						  
      
	   if (alphaS.le.0.01) then                                  !  neglect saturated part
          fact_dens = 1.0                                        !  (dens-densw)/dens   
       endif
       porosity     = (denss-dens)/(denss-densw)              !
       dens2        = (1-porosity)*denss                      !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !       
       fact_dens    = dens_bar_eff/dens_bar                   !    will multipliy tanFi
       if (icpwp.gt.0) then                                   !     If only affects pwp inf      
           fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
           fact_pwp_min = 0.
           fpwp_max     = 1.0
           fpwp_min     = -alphaS*densw/dens_bar_eff
       endif 

    ELSEIF (ic_hrelSat.eq.1) then 
	   
       DO ipoin = 1, npoin
           alphaS       = hrelSat_DF (ipoin)                      !  MP chgd march 2023
           porosity     = n_DF (ipoin)            
           dens2        = (1-porosity)*denss                      !  dens of upper layer (can be unsat)
           dens1        = (1-porosity)*denss + porosity*densw     !  ower layer (sat)   Sr=alphaS
           dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
           dens_bar_eff = dens_bar - densw*alphaS                 !       
           fact_dens    = dens_bar_eff/dens_bar                   !     multiplies tanFi
           dens_bar_DF     (ipoin) = dens_bar
           dens_bar_eff_DF (ipoin) = dens_bar_eff
           fact_dens_DF    (ipoin) = fact_dens
           if (icpwp.gt.0) then                                   !     If only affects pwp inf      
               fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
               fact_pwp_min = 0.
               fpwp_max     = 1.0
               fpwp_min     = -alphaS*densw/dens_bar_eff
               fact_pwp_max_DF (ipoin) = fact_pwp_max   
               fact_pwp_min_DF (ipoin) = fact_pwp_min
               fpwp_max_DF     (ipoin) = fpwp_max
               fpwp_min_DF     (ipoin) = fpwp_min
           endif
       ENDDO 
  
	ENDIF
    
ELSEIF (ic_fct_dens.EQ.2) then                              !  2 phases DF, Hsat variable and pwp optionals

    IF (ic_hrelSat.eq.0) then                                !  all densities constants to compute factors
       porosity     =  (denss-dens)/(denss-densw)            !  mp march 2023
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !       
       fact_dens    = dens_bar_eff/(denss*(1-porosity))       !     multiplies tanFi
       if (icpwp.gt.0) then                                   !     If only affects pwp inf      
           fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
           fact_pwp_min = 0.
           fpwp_max     = 1.0
           fpwp_min     = -alphaS*densw/dens_bar_eff
       endif 

    ELSEIF (ic_hrelSat.eq.1) then                   !  all dens and factors depend on ipoin   

       DO ipoin = 1, npoin
           alphas       = hrelSat_DF (ipoin)                      !  MP chgd Aug 2021
           porosity     = n_DF (ipoin)           !
           dens2        = (1-porosity)*denss + porosity*densw*Sr  !  dens of upper layer (can be unsat)
           dens1        = (1-porosity)*denss + porosity*densw     !  ower layer (sat)   Sr=alphaS
           dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
           dens_bar_eff = dens_bar - densw*alphaS                 !       
           fact_dens    = dens_bar_eff/(denss*(1-porosity))       !     multiplies tanFi
           dens1_DF        (ipoin) = dens1
           dens2_DF        (ipoin) = dens2
           dens_bar_DF     (ipoin) = dens_bar
           dens_bar_eff_DF (ipoin) = dens_bar_eff
           fact_dens_DF    (ipoin) = fact_dens
           if (icpwp.gt.0) then                                   !     If only affects pwp inf      
               fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
               fact_pwp_min = 0.
               fpwp_max     = 1.0
               fpwp_min     = -alphaS*densw/dens_bar_eff
               fact_pwp_max_DF (ipoin) = fact_pwp_max   
               fact_pwp_min_DF (ipoin) = fact_pwp_min
               fpwp_max_DF     (ipoin) = fpwp_max
               fpwp_min_DF     (ipoin) = fpwp_min
           endif
       ENDDO 

    ENDIF

ELSE 
    write (*,*) ' ic_fct_dens given is wrong...', ic_fct_dens
    write (*,*) ' strike a num key to acknowledge '
    read  (*,*) ipoin
    STOP
ENDIF

END Subroutine DFs_FSs_dens_and_pwp_factors     !  MP March 2023

!-------------------------------------------------------------------

       Subroutine get_ws_interfase_friction

!-------------------------------------------------------------------

!
!     Located inside ExtForces_SW 
!     computes interfase friction w and s in submarine landslides
!     frictx units are rho*g*h. Divide by dens*h after output
          
implicit none

integer(ink) nfric_ws  
real   (irk) Mu_sw  ! this is the visco constant C56 (either topo or constant)
real   (irk) Lmix   ! either mixing length C57 or hs, hw
real   (irk) MubyL, tmp

nfric_ws = const(55) + 0.001
if (nfric_ws.eq.1) then
    Mu_sw    = const(56)
    Lmix     = const(57)
    if (Lmix.le.1.e-6) Lmix = h
    MubyL    = Mu_sw/Lmix
    frictx   =  - (vx(1,ipoin)-v_sw(1,ipoin))*MubyL   ! for soils vs-vw, for w vw-vs
    if (ndimn.eq.2) then
       fricty   =  - (vx(2,ipoin)-v_sw(2,ipoin))*MubyL
    endif
else
    write(*,*) ' wrong interfase fricti law  Hit any key  '
    read (*,*) tmp
endif


end Subroutine get_ws_interfase_friction

!-------------------------------------------------------------------

       Subroutine Get_sw_interface_Erosion_SW   
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes erosion_rate dh/dt for w and s in submarine landslides
!        Gets erosion as a rate of change of mass dm=(m/rho)*drho
!         drho = Es*h*abs(vx s- vx w)   where Es is [L]-1  around 1.e-3, 5.e-3...
!         Erosion rate: dh/dt = Es*h*abs(vx s- vx w)
!     using: 
!        xx,yy coordinates of SPH node
!        Zp,Zpx, Zpy topo Z and slope
!        h & vmod   depth and mod. of velocity
!     In the main region, we use const
!        const(51)=1 & const(52) is Hungr constant
!             (17) denss (soil grains) 
!             (18) densw (pore fluid)
!             (21) type of erosion law 1 Hungr
!             (22) Hungr constant
!       

implicit none

integer(ink) sw_eros_law  
real   (irk) vmod_diff  
real   (irk) Keros  ! this is the erosion constant C52 (either topo or constant)
real   (irk) Lmix   ! either mixing length C53 or hs, hw

real   (irk) tmp, tmp1, auxT
real   (irk) Cstar, denss, densw, vol_conc  
real   (irk) cosa, sina, gradZv, theta, thetae, xnum, xden, tanthetae
real   (irk) C_basal, tanfi_basal, beta_basal
real   (irk) Taub, Tauf
real   (irk) n0_eros, fact_n0         ! DFs  mp 9 feb 18

!      ------  Obtain erosion in main region (only Hungr & TB so far)

sw_erosion_rate = 0.0
sw_eros_law = const (50) + 0.1
if (sw_eros_law.eq.1) then
    Keros    = const (51)
    Lmix     = const (52)
    if (Lmix.le.1.e-6) Lmix = h
    vmod_diff   = (vx(1,ipoin)-v_sw(1,ipoin))**2
    if (ndimn.eq.2) then
       vmod_diff   = vmod_diff + (vx(2,ipoin)-v_sw(2,ipoin))**2
    endif
    vmod_diff = vmod_diff**0.5
    sw_erosion_rate = Keros * Lmix * vmod_diff
else
   write(*,*) ' entered get eros sw w/o proper eros law', sw_eros_law
   PAUSE
endif

End Subroutine Get_sw_interface_Erosion_SW  


!-------------------------------------------------------------------

       Subroutine Get_r_crush_SW  
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes new grsin diameter
!         
!     uses fron ext forces: xx, yy, vxx, vyy, 
!          frictx, fricty, vmod, ipoin, dens,  h  

implicit none

integer(ink) law_crush 
real   (irk) r0_crush, rf_crush, B_crush , Taub_crush

law_crush = const(36)

r0_crush = const(37)
rf_crush = const(38)
B_crush  = const(39)

Taub_crush =  frictx*frictx
if (ndimn.eq.2) Taub_crush =  Taub_crush +  fricty*fricty
Taub_crush = dens*h*Taub_crush **0.5

w_crush(ipoin) = w_crush(ipoin) + dt_sph*vmod*Taub_crush
r_crush(ipoin) = r0_crush * B_crush / ( w_crush(ipoin) + B_crush)

END SUBROUTINE Get_r_crush_SW 


!-------------------------------------------------------------------

       Subroutine Get_Grav_Ext_SW
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes grav
!     using: 
!         
!         
!         

implicit none 

!      ------  Obtain gravities  

!     write(*,*) ' ipoin, xx', ipoin,xx
        
     if (if_frictional) then               !  nfric.eq.7.or.nfric.eq.8.or.nfric.eq.10.or.nfric.eq.25
        call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy, nfric, top4, top5, top6, top7, top8, top9)
     else 
        call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
     endif 
                                                             
     if ( ic_cosTheta.eq.1) then                          !  if we use the theta formulation   
        call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)     !  we will define them again              
        deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
        gravz =  cgra/deng
        gravx = -cgra*Zpx/deng                            !  is -gravz*Zpx indeed. Left for pedagogical..
        gravy = -cgra*Zpy/deng
        cgra05 = 0.5*gravz                                ! *** think of joining with nfrict 7 
     elseif ( ic_cosTheta.eq.0) then 
        call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
     
  ! ******************************  DEBUG MP **June 2022
  !      if (x(1,1).le.10.00) then
  !         Zpx = -1
  !      else
  !        Zpx = 0.0
  !     endif
    ! *********************************                  
        deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
        gravz =  cgra 
        gravx = -cgra*Zpx                               !  see the diff with the intrinsic coord system  
        gravy = -cgra*Zpy 
        cgra05 = 0.5*gravz
     elseif ( ic_cosTheta.eq.2) then                    !  do nothing, we have gx gy gz 
        deng  = deng
     endif
   
     ZpGrad = (Zpx*Zpx + Zpy*Zpy)
     If (Zpgrad.gt.1.e-6) ZpGrad = ZpGrad**0.5
   
     IF (ic_ws_Interact.eq.0)  THEN    
                                          ! No interaction. Use g.Zp,x or (h+Z)gradZ  SWL=0 !!
                                          ! Only done within RK loop   
        if (ic_SWalg.eq.0)  then          ! classical landslide formulation
   
          exdvxdt(1,ipoin)                 = exdvxdt(1,ipoin) + gravx !    **  changed  - gravz*Zpx
          if (ndimn.eq.2) exdvxdt(2,ipoin) = exdvxdt(2,ipoin) + gravy !    **  changed  - gravz*Zpy 
         
        elseif (ic_SWalg.eq.1)  then      !  water bodies stable formulation
                   
          exdvxdt(1,ipoin)     = exdvxdt(1,ipoin)   &
                               + ( (rho(ipoin) + Zp) * gravx )/rho(ipoin)
          if (ndimn.eq.2) then   
              exdvxdt(2,ipoin) = exdvxdt(2,ipoin)   &                      ! Chuan 7 May spotted bug
                               + ( (rho(ipoin) + Zp) * gravy )/rho(ipoin)
          endif 
          
        endif 
   
     ELSEIF (ic_ws_Interact.eq.1) THEN    ! Interaction case 
                                                             ! Only don within RK loop 
   
       if     (ic_SWalg.eq.0)  then      !  classical landslide formulation for both. 
                                         !  We deal with interaction later 

              exdvxdt(1,ipoin)                 = exdvxdt(1,ipoin) + gravx !    **  changed  - gravz*Zpx - gravz*Zpx
              if (ndimn.eq.2) exdvxdt(2,ipoin) = exdvxdt(2,ipoin) + gravy !    **  changed  - gravz*Zpy - gravz*Zpy 
      
       elseif (ic_SWalg.eq.1)  then      !  

         if (itype(ipoin).eq.2) then    ! Water bodies interacting with soil *** check gravx use
        
             hhsZ = rho(ipoin) + hs_plus_Z_w(ipoin)        
             exdvxdt(1,ipoin)     = exdvxdt(1,ipoin)   &                   ! leave it ***
                                  - gravz* ( hhsZ * gradhs_plus_Z_w(1,ipoin) )/rho(ipoin)
             if (ndimn.eq.2) then   
                 exdvxdt(2,ipoin) = exdvxdt(2,ipoin)   &                    ! BUG corrected MP after Chuan 4.07.2019
                                  - gravz* ( hhsZ * gradhs_plus_Z_w(2,ipoin) )/rho(ipoin)
             endif                               

         elseif (itype(ipoin).eq.1) then   ! Soil interacting with water bodies
   
                exdvxdt(1,ipoin)                 = exdvxdt(1,ipoin)  + gravx !**  changed - gravz*Zpx - gravz*Zpx
                if (ndimn.eq.2) exdvxdt(2,ipoin) = exdvxdt(2,ipoin)  + gravy !**  changed - gravz*Zpy - gravz*Zpy
             
         endif 
         
       endif
   
     ELSEIF (ic_ws_Interact.eq.2) THEN    ! DF case
   
      exdvxdt(1,ipoin)                 = exdvxdt(1,ipoin)  + gravx !    **  changed  - gravz*Zpx - gravz*Zpx
      if (ndimn.eq.2) exdvxdt(2,ipoin) = exdvxdt(2,ipoin)  + gravy !    **  changed  - gravz*Zpy - gravz*Zpy
      
     ENDIF

END SUBROUTINE   Get_Grav_Ext_SW



!-------------------------------------------------------------------

       Subroutine basal_friction_Ext_SW
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes basal friction
!     using: 
!         
!     Moved to the end filtering of velocities, max taub    
!         

implicit none 

real   (irk) tanfi_mu, b_muI, r_muI, n_muI, Froude, ff   !    used for mu(I) rheo MP 25 nov 17
real   (irk) mu2, Inert, Inert0                          !  I=inertia nr.. mu2 and I0 are params 2nd mu(I) law

real   (irk) vdot, deccelT, deccelN, accelT, frictx0, fricty0, fricT, TotFrict, hhh, ctf
real   (irk) alpha_slip, fa                              ! MP april 11th 2018

real   (irk) Tauy, Mu, Tauy0, Tauy_expf, Mu0, Mu_expf, Cvv   !  MP 21st July 2023  model 201
real   (irk) xi_v, xi_expf 
real   (irk) lambda_B , Cvv_max                     !  MP 19 nov 2023 Bagnold       
real   (irk), allocatable:: deccel(:), accel (:), vstar(:), vcurr (:)  ! used for filter tau
real   (irk), allocatable:: Rot (:,:), v_t   (:), v_n  (:)

            ! h     = rho(ipoin)    obtained in extforces before calling this subroutine
            ! hu    = h*vx(1,ipoin)
            ! hv    = 0.0
            ! if (ndimn.eq.2) hv = h*vx(2,ipoin)

            ! if (h.lt.comp) then   !  ** see here comp which limits 1/h
            !   h        = comp
            !   hu       = 0.0
            !   hv       = 0.0
            !   hinv     = 0.0
            !   vmod     = 0.0
            !   hvmod    = 0.0
            !   hvmodinv = 0.0
            !else
            !   hinv     = 1/h
            !   hvmod    = sqrt(hu*hu+hv*hv)
            !   vmod     = hvmod*hinv
            !   hvmod    = sqrt(hu*hu+hv*hv)
            !   if (hvmod.ge.1.e-5) then
            !      hvmodinv = 1./hvmod
            !   else
            !      hvmodinv = 0.0
            !   endif  
            ! endif

if (R_Drum.gt.1.e-3) then    !!  *** DRUM   
        sinth = xx/R_drum
        costh = sqrt (1-sinth*sinth)
        hu    = - w_drum_radps * R_drum * costh*rho(ipoin) + hu
endif   

if (icpwp.eq.1) then         !      ------  PWP computation
    pwp = u(ipoin)
elseif (icpwp.eq.2) then
    pwp = Uaux(1,ipoin)      !      ------  point 1 is at the bottom
elseif (icpwp.eq.3) then
    pwp = Uaux(1,ipoin)      !      ------  The only point 1 is at the bottom
else
    pwp = 0.0
endif 

hfrict   = max(h,hfrict0)   ! ** hfrict  limits friction 1/h

frictx = 0.0
fricty = 0.0 

if(nfric.eq.1) then                     !      ------  Newtonian fluid, turbulent
   cf1      = cgra*cmanning*cmanning
   frictx   = cf1*hu*hvmod/hfrict**c73 
   fricty   = cf1*hv*hvmod/hfrict**c73 
   TotFrict = (frictx*frictx + fricty*fricty)**0.5  
elseif (nfric.eq.2) then                !      ------ Newtonian fluid, laminar regime  
   cf2a     = constK/(8.*dens)
   frictx   = cf2a*hu/hfrict**2.
   fricty   = cf2a*hv/hfrict**2. 
   frictx0  = cf2a  
   TotFrict = (frictx*frictx + fricty*fricty)**0.5
    !elseif (nfric.eq.3) then           !      ------  Coussot law: deactivated 
    !   cf3a = 1.93*(constK**0.9)*(tauy0**0.1)
    !   cf3b   = (tauy0 + cf3a*(vmod/hfrict)**0.3)/dens
    !   frictx = cf3b*hu*hvmodinv
    !   fricty = cf3b*hv*hvmodinv
                       
elseif (nfric.eq.3) then                !      ------  Bingham erosion
  
   ! elseif (nfric.eq.3) then ! ------ Bingham erosion
    frictx  = 0.0
    fricty  = 0.0
    frictxx = 0.0
    frictyy = 0.0
    if (mass(ipoin).gt.1.e-6) then ! nsoilB varies between 0 (water) and ->
        nsoilB = (mass(ipoin)-mass0(ipoin))/mass(ipoin)
    else
        write(*,*) ' mass of point ', ipoin, ' is zero...code stopped'
        pause
    endif
    csoilB  = const(3)
    fsoilB  = 1. - exp(- csoilB * nsoilB ) !CHUAN
    constKB = constK * fsoilB ! viscosity
    tauy0B  = tauy0 * fsoilB ! cohesion
    densB   = 1000. + (dens - 1000.)*nsoilB
    if (nsoilB.gt.0.01) then ! *** cannot cycle here
        a_hat = (6.* constKB *hvmod)/(tauy0B *hfrict**2.) ! changed
        if (a_hat .gt. 50) then                !   viscous, xi-> 0
            cf4b   = 3*constKB/densB           !   dims are  accel*h 
            frictx = cf4b * hu/hfrict**2
            fricty = cf4b * hv/hfrict**2
        elseif (a_hat .lt.0.01) then           !   cohesive, xi-> 1
            cf4b   = tauy0B/densB              !   dims are  accel*h 
            frictx = cf4b * hu * hvmodinv
            fricty = cf4b * hv * hvmodinv
        else                                   ! intermediate
            a2     = 48.
            b2     = -32.*a_hat-114.
            c2     = 65.
            xh     = (-b2-(b2*b2-4.*a2*c2)**0.5)/(2.*a2)
            Tau0   = tauy0B/xh
            cf4b   = Tau0/densB                ! dims are  accel*h
            frictx = cf4b*hu*hvmodinv
            fricty = cf4b*hv*hvmodinv
       endif
    endif
    cmanning = const(8)
    cf1     = cgra*cmanning*cmanning
    frictxx = cf1*hu*hvmod/hfrict**c73
    frictyy = cf1*hv*hvmod/hfrict**c73
    if (abs(frictxx).gt.abs(frictx)) frictx = frictxx
    if (abs(frictyy).gt.abs(fricty)) fricty = frictyy 

    TotFrict = (frictx*frictx + fricty*fricty)**0.5

elseif (nfric.eq.201) then   !      ------  Bingham with tau and mu dep on (1-n_df(ipoin))
                             !              MP 21st July 2023
    frictx  = 0.0          
    fricty  = 0.0
    frictxx = 0.0
    frictyy = 0.0

    Tauy0      = const(102)
    Tauy_expf  = const(103)
    Mu0        = const(104)
    Mu_expf    = const(105)
    Cvv        = 1.-n_DF(ipoin)

    Tauy       = Tauy0 * exp (Tauy_expf * Cvv)
    Mu         = Mu0   * exp (Mu_expf   * Cvv)
        
    a_hat = (6.* Mu *hvmod)/(tauy  *hfrict**2.) ! changed
    if (a_hat .gt. 50) then                     !   viscous, xi-> 0
        cf4b   = 3*Mu/dens                      !   dims are  accel*h 
        frictx = cf4b * hu/hfrict**2
        fricty = cf4b * hv/hfrict**2
    elseif (a_hat .lt.0.01) then                !   cohesive, xi-> 1
        cf4b   = Tauy /dens                     !   dims are  accel*h 
        frictx = cf4b * hu * hvmodinv
        fricty = cf4b * hv * hvmodinv
     else                                       ! intermediate
        a2     = 48.
        b2     = -32.*a_hat-114.
        c2     = 65.
        xh     = (-b2-(b2*b2-4.*a2*c2)**0.5)/(2.*a2)
        Tau0   = tauy/xh
        cf4b   = Tau0/dens                ! dims are  accel*h
        frictx = cf4b*hu*hvmodinv
        fricty = cf4b*hv*hvmodinv 
    endif

    cmanning = const(8)
    cf1     = cgra*cmanning*cmanning
    frictxx = cf1*hu*hvmod/hfrict**c73
    frictyy = cf1*hv*hvmod/hfrict**c73
    if (abs(frictxx).gt.abs(frictx)) frictx = frictxx
    if (abs(frictyy).gt.abs(fricty)) fricty = frictyy 

    TotFrict = (frictx*frictx + fricty*fricty)**0.5
  
   !      ------  Bingham and HB   exact solution

elseif (nfric.eq.4.or.nfric.eq.23.or.nfric.eq.24.or.nfric.eq.25) then
       
    Tauy0      = const (6)
    constMu_CF = const (7)
    constK     = const (7)
    m_vexp     = const (8)
    Cm         = m_vexp/ ((m_vexp+1)*(2*m_vexp+1))
    a_hat      = (1./Cm)* ( (constK/tauy0)**(1./m_vexp) )*(hvmod/hfrict**2)

    if (tauy0.lt.1.e-32)then
      write(*,*)' Error entering Yield: must be greater than 1.e-32'
      stop
    endif
    
    if ( nfric.eq.4) then        !  Bingham 
         
        a_hat  = (6.*constK*hvmod)/(tauy0*hfrict**2.)
        if (a_hat .gt. 50) then                !   viscous, xi-> 0
            cf4b   = 3*constK/dens             !   dims are  accel*h 
            frictx = cf4b * hu/hfrict**2
            fricty = cf4b * hv/hfrict**2
        elseif (a_hat .lt.0.01) then           !   cohesive, xi-> 1
            cf4b   = tauy0/dens                !   dims are  accel*h 
            frictx = cf4b * hu * hvmodinv
            fricty = cf4b * hv * hvmodinv
        else                                   ! intermediate
            a2     = 48.
            b2     = -32.*a_hat-114.
            c2     = 65.
            xh     = (-b2-(b2*b2-4.*a2*c2)**0.5)/(2.*a2)
            Tau0   = tauy0/xh
            cf4b   = Tau0/dens                ! dims are  accel*h
            frictx = cf4b*hu*hvmodinv
            fricty = cf4b*hv*hvmodinv
       endif

    endif
  
    if (nfric.eq.23) then    !  HB models, including m=1 *** verify a_hat cuts ****  MP & Chuan 
       
       if (m_vexp.eq.1) then ! Bingham
          a2  = 48.
          b2  = -32.*a_hat-114.
          c2  = 65.
       elseif ( m_vexp.eq.2) then  
          a2  = 10.
          b2  = - a_hat- 19.2188
          c2  = 9.1875
       elseif ( abs(m_vexp-0.333333).le.1.e-3 ) then 
          a_hat = (4.5*constK*hvmod)/(tauy0*hfrict**2.)
          a2  = 48.
          b2  = - 32.* a_hat- 114.
          c2  = 65. 
       elseif ( m_vexp.eq.0.5 ) then 
          a_hat = (4.5*constK*hvmod)/(tauy0*hfrict**2.)
          a2  = 48.
          b2  = - 32.* a_hat- 114.
          c2  = 65. 
       else
         write(*,*)' Error entering m exp only 1,2,1/2 and 1/3'
         stop
       endif

       xh   = (-b2-(b2*b2-4.*a2*c2)**0.5)/(2.*a2)
       Tau0 = tauy0/xh
          
    endif
  
    if (nfric.eq.24.or.nfric.eq.25) then    !  approx HB models, including m=1 
       Tau0 = (constMu_CF*(vmod/hfrict)**m_vexp)
       if (nfric.eq.24) then 
           Cm = 1.   
       elseif (nfric.eq.25) then
           Cm = ((2*m_vexp+1)/m_vexp)**m_vexp
       endif    
       Tau0 = Tauy0 + Tau0*Cm
    endif
    
    cf4b   = Tau0/dens
    frictx = cf4b*hu*hvmodinv
    fricty = cf4b*hv*hvmodinv 
    TotFrict = (frictx*frictx + fricty*fricty)**0.5  
    
!      ------  Frictional fluid, modified friction law, V2/R  -----------------------
!              We have included here Pwp effects 

elseif ( if_frictional.eq.1)   then !  setup in extforces, 17th dec 17

    cv2byR = 0.0    
    
!      ------  Frictional fluid with curvature effect    
!              if (nfric.eq.7.or.nfric.eq.8.or.nfric.eq.10) then

    if (if_V2byR.eq.1) then       !  setup in extforces, 17th dec 17
        vxx   = hu*hinv                            
        vyy   = hv*hinv
        vmod  = (vxx*vxx + vyy*vyy)**0.5
        cv2byR=0.0
        if(vmod.ge.1e-5) then
           cosa    = vxx/vmod
           sina    = vyy/vmod
           cos2    = cosa*cosa
           sin2    = sina*sina
           cbs     = sina*cosa
           xnum    = top7*cos2 + top8*cbs + top9*sin2
           xden    = top4*cos2 + top5*cbs + top6*sin2
           rinv    = xnum/xden
           cv2byR  = vmod*vmod*rinv
        endif
    endif 

    tanfi = tanfi8                       !    Pore water pressures. Apparent friction angle

    if (nfric.eq.8.or.  &                   !  setup in extforces, 17th dec 17
        nfric.eq.26.or.nfric.eq.27) then    !    used for mu(I) rheo MP 25 nov 17
        xi_voellmy = 0.0               !    because C3 is also voellemy
        if (nfric.eq.8.or.nfric.eq.26) then
           b_muI      = const(3)
           r_muI      = const(6)
           n_muI      = const(7)
           Froude     = vmod/(gravz*h)**0.5
           ff         = (1 + b_muI*(Froude*r_muI/h)**n_muI)
           tanfi      = tanfi * ff
        elseif (nfric.eq.27) then
           Mu2        =  const(3)
           r_muI      =  const(6)
           Inert0     =  const(7)
           Froude     =  vmod/(gravz*h)**0.5
           Inert      =  Froude*r_muI/h
           tanfi      =  tanfi + (Mu2-tanfi)/(1.0 + (Inert0/Inert))  
        endif           
    endif

    hhh = rho(ipoin)                     !    MP June 2017 pwp normalized by hs+hw
    if (ic_ws_interact.eq.2) hhh = h_DF(ipoin)
    if (pwp.gt.fpwp_max*hhh) pwp = fpwp_max*hhh   !    filters out pwp
    if (pwp.lt.fpwp_min*hhh) pwp = fpwp_min*hhh
    u(ipoin) = pwp

    if (icpwp.eq.1.or.icpwp.eq.2.or.icpwp.eq.3) then     !  pwp is always the relative value at the bottom 
    
       factpw = 1. - pwp*hinv                            !   pwp  *hinv  now we use relative values
       if (factpw.lt.fact_pwp_min)   factpw = fact_pwp_min
       if (factpw.gt.fact_pwp_max )  factpw = fact_pwp_max
       
!       tanfi = tanfi8*factpw
       tanfi = tanfi *factpw       !  tanfi modified by mu rheo, or tanfi8
       
    elseif (icpwp.eq.0) then                             ! obsolete   
!       if (Bfact.ge.1.e-6) then
!           tanfi  = tanfi8 + (tanfi0-tanfi8)*exp(-time/Bfact)   !  MP Dec 2023 OBSOLETE
!       endif
    endif


!   ------ In the case there is basal friction associate to terrain 
!           check that ipoin belongs to any zone
!           and if so, obain tan delta

    if (ic_basal_friction.eq.1) then   
       call Get_Topo_Basal_Friction_SW (Tandelta)
!       tanfi = min(tanfi,Tandelta)                  ! MP 1st June 18, Roberto
        if (Tandelta.ge.0) then
            tanfi = Tandelta
        endif
    endif   
      
!    gtanfi = (cgra+cv2byR)*tanfi                ! changes: (i)   use gravz
    gtanfi = ( gravz + cv2byR)*tanfi*fact_dens  !          (ii) (denss-densw)/densw DFs
    if (gtanfi.lt.0.0) gtanfi=0.0               !          (iii) (dens-densw)/dens fully sat 
          
    cf5      = gtanfi                           ! error corrected, it was *h              
    frictx   = cf5*hu*hvmodinv* h
    fricty   = 0.0 
    frictx0  = cf5  
    if (ndimn.eq.2) then
       fricty = cf5*hv*hvmodinv* h
    endif                                       !       Dims are [g.h]
    
!      ------  Voellmy and new CF model: tau=tau frict + (25*Mu_CF/4)*(V2/h2)  Mu=50 approx

    if (Xi_Voellmy.gt.0.0) then                          ! Voellmy = frictional + turbulent
       
        vxx    = vx(1,ipoin)
        vyy    = 0.0
        if (ndimn.eq.2) vyy = vx(2,ipoin)
        vmod   = vxx*vxx + vyy*vyy
        if (vmod.gt.0.0) then
            vmod = sqrt(vmod)
        else
            vmod = 0.0
        endif
        if (nfric.eq.202) then               !  MP 23 July 2023 
            Cvv        = 1.-n_DF(ipoin)
            xi_v       = const (102)
            xi_expf    = const (103)
            xi_voellmy = xi_v * exp (-Cvv*xi_expf)
            ns_dep_props(1,ipoin) = xi_voellmy
        elseif (nfric.eq.203) then          ! MP 19 Nov 2023 Bagnold w/o h**2 term
                                            ! MP 19 Nov 2023 Bagnold w/o h**2 term
            Cvv        = 1.-n_DF(ipoin)
            Mu0        = const (102)
            Cvv_max    = const (103) 
            lambda_B   = ((Cvv_max/Cvv)**(1./3.)) - 1.00  !  MP and SM Nov 2023
            lambda_B   = max (lambda_B,0.01)      
            xi_voellmy = (dens*cgra*lambda_B**2.) /Mu0   ! Momentum eqn divides by dens; fric * cgra
             
            ns_dep_props(1,ipoin) = xi_voellmy
        endif

        frictx = frictx + cgra  *vmod*vxx/xi_voellmy 
        fricty = fricty + cgra  *vmod*vyy/xi_voellmy
          
    elseif (Xi_Voellmy.lt.0.) then                          !  This is new CF model
       
        vxx    = vx(1,ipoin)
        vyy    = 0.0
        if (ndimn.eq.2) vyy = vx(2,ipoin)
            vmod   = vxx*vxx + vyy*vyy
        if (vmod.gt.0.0) then
            vmod = sqrt(vmod)
        else
            vmod = 0.0
        endif
        frictx = frictx + constMu_CF*vmod*vxx/(hfrict*hfrict)
        fricty = fricty + constMu_CF*vmod*vyy/(hfrict*hfrict)
        
    endif              
            
    if (nfric.eq.9.or.nfric.eq.10) then 
   
       frict_abs  = (frictx*frictx + fricty*fricty)**0.5   
               
       !  correction MP 10-april 2018
       !       fact9      = 1. +  ( (2.*C8_Muf*vmod/hfrict)**(1./C67_N) )*(frict_abs**(-C6_n1/C67_N))

       fact9      = ( (2.*C8_Muf*vmod/hfrict)**(1./C67_N) )*(frict_abs**(-C6_n1/C67_N))
       
       if ((nconst.gt.45) .and. (const(47).GT.1.e-6)) then    !  ch MP 11 April 2018
          alpha_slip = const(47)
          fa         = (1-alpha_slip)/(1+alpha_slip)
          fact9      = fact9*fa
       endif
           
       frictx     = frictx * (1 + fact9 )
       if (ndimn.eq.2) then
          fricty = fricty * (1 + fact9 )
       endif 

    endif 
    
    TotFrict = (frictx*frictx + fricty*fricty)**0.5
            
elseif (nfric.eq.11) then 

 !      ------  Bagnold model : Taub    = (25/4)*muB*sinFiB*(v/h)**2
 !                              frictx  = KBagnold*(v/h)**2  KB=(1/rho)*(25/4)*muB*sinFiB 
 !                              exdvxdt = frictx/h  = Taub /(rho*h)    
 
   vxx    = hu*hinv
   vyy    = hv*hinv
   frictx = K_Bagnold*vxx*vmod/(hfrict*hfrict)
   fricty = K_Bagnold*vyy*vmod/(hfrict*hfrict)
   TotFrict = (frictx*frictx + fricty*fricty)**0.5

elseif(nfric.eq.22) then                            ! simplified Bham + non-newtonian
               
   Tauy0      = const (6)
   constMu_CF = const (7)
   m_vexp     = const (8)
   tb_abs     = Tauy0 + (constMu_CF*(vmod/hfrict)**m_vexp)
   tb_abs     = tb_abs/dens
   if (vxx.ge.1.e-4) then
      frictx     = tb_abs*vxx/vmod
   endif
   if (ndimn.eq.2.and.vyy.ge.1.e-4) then
      frictx     = frictx + tb_abs*vyy/vmod
   endif                                      
   TotFrict = (frictx*frictx + fricty*fricty)**0.5

elseif(nfric.gt.27)then
    write(*,*) ' Error: nfric not implemented'
    stop
endif                                                !      ------  End friction laws
 
!      ------  Limits to Taub, chk v reversal
!
!              vstar is an approximation of v at n+1
!                    computed as dt* (accel + deccel)
!                    accel  is the acceleratio of int +ext forces
!                    deccel is the basal friction deccel. Sign is minus frictx
!                                                         value = frictx/h
!                    v_t and v_n unit vectors along v and its normal
!                    Rot = ( 0, -1; 1,0)  v_n = Rot¨v_t
!
!              when velocity changes sign, it may happen that friction causes this reversal
!         
! if (nfric.ne.4) then
if(.not.allocated(accel)) then
   allocate ( accel(ndimn), deccel(ndimn), vstar(ndimn), vcurr(ndimn) )
   allocate ( rot(ndimn,ndimn), v_t(ndimn), v_n (ndimn)  )
endif

deccel    = 0.0
deccel(1) = - frictx*hinv
if (ndimn.eq.2) deccel(2) = - fricty*hinv

accel(:) = exdvxdt(:,ipoin) + indvxdt(:,ipoin)

IF (vmod.gt.1.e-5) THEN
   rot = 0.0
   if (ndimn.eq.2) then
      rot(1,2) = -1.
      rot(2,1) = 1.
   endif   
   vcurr(:) = vx(:,ipoin)
   vstar    = vcurr + dt_sph*(accel + deccel)
   vdot     = dot_product (vstar, vcurr) 
   if (vdot.lt.0 ) then
       v_t     = vcurr/vmod
       v_n     = matmul(Rot,v_t)
       accelT  = dot_product (accel , v_t)      !  Projection of  accel on t direction
       deccelT = dot_product (deccel, v_t)      !  Projection of deccel on t direction
       if ((accelT*deccelT).lt.0) then
           ! deccelT = -accelT                    ! it will provide 0 total accel +deccel a,ong t
           ! deccelN = dot_product (deccel, v_n)  ! not modified
           ! deccel  = deccelT*v_t + deccelN*v_n
           deccel  = -vcurr/dt_sph - accel        ! chnaged MP 11 May19 to avoid accel=0->v=ct
           frictx  = -deccel (1) *h
           if(ndimn.eq.2) fricty = -deccel(2)*h
       endif
   else                                         !  weird case, vx opposite sign than bot accel and deccel
       deccel = -vcurr/dt_sph                   !  we assume deccel will make v zero, then apply acel
!       frictx  = -deccel (1) *h                 !   MP  14 th january 2019
!       if(ndimn.eq.2) fricty = -deccel(2)*h   
   endif 
ELSE   !  --- vmod smaller than 1.e-5

   accelT = (dot_product (accel, accel))**0.5

   if (accelT.gt.1.e-5) then       ! MP 28th Feb 2017
       ctf = TotFrict
       if (if_frictional.eq.1) ctf = cf5
       frictx = ctf*accel(1)/accelT
       if (ndimn.eq.2) fricty = ctf*accel(2)/accelT 
       fricT = (frictx*frictx + fricty*fricty)**0.5
       if (accelT.le.fricT) then  
          frictx = accel(1)                  !  exdvxdt(1,ipoin)*accel(1)/accelT
          if (ndimn.eq.2)  fricty = accel(2) !  exdvxdt(2,ipoin)*accel(2)/accelT 
       endif
   else
      frictx = 0.0; fricty = 0.0  
   endif 

   frictx = frictx*h
   if (ndimn.eq.2) fricty = fricty*h    
           
ENDIF 

!   endif
!  ----- end of filter ------


END SUBROUTINE   basal_friction_Ext_SW



!-------------------------------------------------------------------

       Subroutine Get_Drag_SW   
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes new velocities vs and vw
!        
!        denss* (1-n)* h * vs/dt = ...+ alphas *h * R   we have ntegraed along saturated part only
!        densw*    n * h * vw/dt = ...- alphas *h * R 
!                           R    = Cd * (Vw - Vs)  Cd = C1
!        solving the ODEs dvs/dt = Cs*(Vw-Vs)  
!                         dvw/dt = Cw*(Vs-Vw)   
!        where
!              C1 is Cd, depends on law  
!                     2 Anderson.... n(n-1)(rhos-rhof)g/(VT*n**m)
!                     1 Darcy        n**2 g rhow/Kw
!                     3 karmanK      kwd  = kwd0 * nnn**3/(1- nnn*nnn  )  
!              Cw     = C1/(denss*(1-n_DF(ipoin))    
!              Cs     = C1/(densw*(  n_DF(ipoin))
!              a_drag = Cs/Cw = alpha_sat * (denss*(1-n_DF(ipoin)) /C1/(densw*(  n_DF(ipoin)) 
!              C1 depends on drag law
!
!     Several cases ....(general, fluid limit, solid limit) ................
!     
!     Case 1 general  ic_ss = 1
!          f1 = (Vs0 + a_drag*Vw0)  f2 = (Vs0-Vw0) expp = exp(-(1+a_drag*dt_sph))        
!          Vs = ( f1 + a_drag*t2*expp)/(1+a_drag)
!          Vw = ( f1 +        t2*expp)/(1+a_drag)
!        subcase 1.1 assymptotic, kw very small. Identif const(27)<0.0  ic_ss = 1
!          Vs = f1/(1+a_drag)   Vw = Vs
!        subcase 1.2 assymptotic kw very large  Identif const(27)>>1.e5  ic_ss = 2
!          Vs = Vs0   Vw = Vw0    
!     Case 2 general  fluid limit n-> 1    
!          auxCd depends on law
!            law 2 Anderson ((rhos-rhow)/rhos) * (g/VT)*(1/n**(m-1))
!            law 1 Darcy    (n**2/(1-n))*(rhow/rhos)*(g/kwd)  (1-n) use n0 C29
!            law 3 Karman   ((1+n)/n)*(rhow/rhos)*(g/kwd0)           
!          Vs = Vw0 + (Vs0-Vw0)*exp(-auxcd*dt)   
!          Vw = Vw0
!        subcase 1.1 assymptotic, kw very small. Identif const(27)<0.0  ic_ss = 1
!          Vs = Vw0   Vw = Vw0
!        subcase 1.2 assymptotic kw very large  Identif const(27)>>1.e5  ic_ss = 2
!          Vs = Vs0   Vw = Vw0  
!     Case 3 general  solid limit n-> 1  
!          Vs = Vs0   Vw = Vs0  
!
! .................................................................................
!     Note: vars inherited from ExternalForces_SW:
!           darcinian_law, const, denss, densw, dt_sph, grav
!
!        read ndf0 = const(29) but, if we have forgotten take it as 0.01  
!        if (ndf0.le.1.e-2) ndf0 = 1.e-2
!

implicit none

integer(ink) ipoin, idimn, ic_ss  
real   (irk) V_T, expm, C1, Cw, Cs, a_drag, Cn1, ndf0, C10
real   (irk) Vs0, Vw0, Vs, Vw, f1, f2, expp, g , nnn
real   (irk) kwd, kwd0, auxCd, alphas, alph

g    = abs(const (1))
ic_ss = 0                    ! index for steady state assympt solution

if (darcinian_law.eq.1) then
    kwd  = const(27) 
    ndf0 = const(29)
    if (ndf0.le.1.e-2) ndf0 = 1.e-2  ! changed MP 5 feb 18
    if (kwd.lt.1.e-10) then
       ic_ss = 1
    elseif (kwd.gt.1.e6) then
       ic_ss = 2
    endif
elseif (darcinian_law.eq.2) then
    V_T  = const(27)
    expm = const(28)  
    ndf0 = const(29)
    if (ndf0.le.1.e-2) ndf0 = 1.e-2
    if (V_T.lt.1.e-10) then
       ic_ss = 1
    elseif (V_T.gt.1.e6) then
       ic_ss = 2
    endif 
    if (V_T.lt.1.e-10) ic_ss = 1
elseif (darcinian_law.eq.3) then
    kwd0 = const(27)
    ndf0 = const(29)
    if (ndf0.le.1.e-2) ndf0 = 1.e-2 
    if (kwd0.lt.1.e-10) then
       ic_ss = 1
    elseif (kwd0.gt.1.e6) then
       ic_ss = 2
    endif
endif

do ipoin = 1, npoin  ! --------------------------------------IPOIN LOOP   

   nnn  = n_DF(ipoin)
         
   alphaS = const(14)                        ! CH MP Sept 2021                         
   alph   = alphaS                           ! alph*h saturated layer
   if (ic_hrelSat.EQ.1) then
       alphas = hrelSat_DF(ipoin)                
   endif    

   do idimn = 1, ndimn  ! ................................. IDIMN LOOP
         
      if  (itype(ipoin).eq.1) then
           Vs0    = vx  (idimn,ipoin)
           Vw0    = v_sw(idimn,ipoin)
      else if (itype(ipoin).eq.2) then
           Vw0    = vx    (idimn,ipoin) 
           Vs0    = v_sw  (idimn,ipoin)
      endif 
     
      if (alph.lt.1.e-3 ) then   !  almost i¡unsat soil, no interaction
         Vs  = Vs0               ! moved here MP March 2023, vs0 vw0 defined now
         Vw  = Vw0
         CYCLE
      endif  

      if ( (n_DF(ipoin).GT.ndf0).AND.(n_DF(ipoin).LT.(1.0-ndf0)) ) THEN ! general porosity case 

          if (ic_ss.eq.0) then                         !    general interaction case                    
              if (darcinian_law.eq.2) then             !      ......  Anderson Law
                  C1   = nnn*(1-nnn)*(denss-densw)*const(1)
                  C1   = C1 / (V_T * nnn**expm )
              elseif (darcinian_law.eq.1) then         !     .......  Darcy
                  C1   = nnn*nnn * densw*g /kwd
              elseif (darcinian_law.eq.3) then         !     ........  Darcy + Karman K
                  kwd  = kwd0 * nnn**3/(1- nnn*nnn  )
                  C1   = nnn*nnn * densw*g /kwd     
              endif 

              Cs   = C1/ (  denss*( 1-n_DF(ipoin) ) )
              Cw   = C1/ (  densw*(   n_DF(ipoin) ) ) 
              a_drag = alphas * (  densw*( n_DF(ipoin) ))/( denss*( 1-n_DF(ipoin) ))    ! Unsat think MARCH 2023 *****    
              f1     = (Vs0 + a_drag*Vw0)
              f2     = (Vs0-Vw0)  
              expp   = exp(-(1 + a_drag)* Cw *dt_sph)                          ! general case
              Vs  = ( f1 + a_drag*f2*expp)/(1.0 + a_drag)
              Vw  = ( f1 -        f2*expp)/(1.0 + a_drag)   !Chuan correction 1
           
          elseif (ic_ss.eq.1) then     ! very low permeability   direct formula steady state
              a_drag = alphas *(densw*( n_DF(ipoin) ))/( denss*( 1-n_DF(ipoin) ))       
              Vs     = (1/(1.+a_drag))*(vs0 + a_drag*vw0)       ! BUG march 2023
              Vw     = vs
!              if (ic_slow.eq.1) then    !   for slow cases MP March 2023
!                 Vs = Vs0
!                 Vw = Vs
!              endif
           elseif (ic_ss.eq.2) then
              Vs = vs0
              Vw = vw0  
           endif

       elseif (n_DF(ipoin).GE.(1.0-ndf0)) then                              !  -----------fluid limit  n->1 

          if (ic_ss.eq.0) then
             if (darcinian_law.eq.2) then                       !        Anderson Law 
                 auxCd  = ((denss-densw)/denss )*(g/V_T)*(1./nnn*(expm-1))
             elseif (darcinian_law.eq.1) then                   !        Darcy
                 auxCd  = (nnn*nnn/(1-ndf0)) * (densw/denss) * (g/kwd)
             elseif (darcinian_law.eq.3) then                   !        Darcy KK
                 auxCd  = ((1+nnn)/nnn)*(densw/denss)*(g/kwd0) 
             endif 
             Vs  = Vw0 + (Vs0-Vw0)*exp(-auxCd*dt_sph)   
             Vw  = Vw0

          elseif (ic_ss.eq.1) then   ! direct formula steady state
             vw = vw0
             vs = vw0                ! mp5 feb 18
          elseif (ic_ss.eq.2) then
             vw = vw0
             vs = vs0
          endif

       elseif     (n_DF(ipoin).LE.ndf0 ) then                               !  ----------------  solid limit
            Vs     = Vs0  
            Vw     = Vs0
       endif
         
       if  (itype(ipoin).eq.1) then
            vx  (idimn,ipoin)    = Vs 
            v_sw(idimn,ipoin)    = Vw      
       else if (itype(ipoin).eq.2) then
            vx    (idimn,ipoin)  = Vw   
            v_sw  (idimn,ipoin)  = Vs 
       endif 
         
    enddo               ! idimn ---------------------------- 
   
enddo                   ! ipoin            ----------------------------


END SUBROUTINE Get_Drag_SW  


!-------------------------------------------------------------------

       Subroutine Get_Drag_SW_OLD   
       
!-------------------------------------------------------------------
!     ****** before UNSAT sept 2021
!     Located inside ExtForces_SW  
!     computes new velocities vs and vw
!         
!        solving the ODEs dvs/dt = Cs*(Vw-Vs)  
!                         dvw/dt = Cw*(Vs-Vw)   
!        where
!              C1 is Cd, depends on law  
!                     2 Anderson.... n(n-1)(rhos-rhof)g/(VT*n**m)
!                     1 Darcy        n**2 g rhow/Kw
!                     3 karmanK      kwd  = kwd0 * nnn**3/(1- nnn*nnn  )  
!              Cw     = (Cd*alphasat)/(denss*(1-n_DF(ipoin))      
!              Cs     = Cd/(densw*(  n_DF(ipoin))     *** see the alphasat here Sept 2021 MP 
!                      we do not include the alph factor, simplified LHS and RHS
!                      but when alphasat -> 0 vw = vw0 and vs = vs0
!              a_drag = Cs/Cw    
!              Cd depends on drag law
!
!     Several cases ....(general, fluid limit, solid limit) ................
!     
!     Case 1 general  ic_ss = 1
!          f1 = (Vs0 + a_drag*Vw0)  f2 = (Vs0-Vw0) expp = exp(-(1+a_drag*dt_sph))        
!          Vs = ( f1 + a_drag*t2*expp)/(1+a_drag)
!          Vw = ( f1 +        t2*expp)/(1+a_drag)
!        subcase 1.1 assymptotic, kw very small. Identif const(27)<0.0  ic_ss = 1
!          Vs = f1/(1+a_drag)   Vw = Vs
!        subcase 1.2 assymptotic kw very large  Identif const(27)>>1.e5  ic_ss = 2
!          Vs = Vs0   Vw = Vw0    
!     Case 2 general  fluid limit n-> 1    
!          auxCd depends on law
!            law 2 Anderson ((rhos-rhow)/rhos) * (g/VT)*(1/n**(m-1))
!            law 1 Darcy    (n**2/(1-n))*(rhow/rhos)*(g/kwd)  (1-n) use n0 C29
!            law 3 Karman   ((1+n)/n)*(rhow/rhos)*(g/kwd0)           
!          Vs = Vw0 + (Vs0-Vw0)*exp(-auxcd*dt)   
!          Vw = Vw0
!        subcase 1.1 assymptotic, kw very small. Identif const(27)<0.0  ic_ss = 1
!          Vs = Vw0   Vw = Vw0
!        subcase 1.2 assymptotic kw very large  Identif const(27)>>1.e5  ic_ss = 2
!          Vs = Vs0   Vw = Vw0  
!     Case 3 general  solid limit n-> 1  
!          Vs = Vs0   Vw = Vs0  
!
! .................................................................................
!     Note: vars inherited from ExternalForces_SW:
!           darcinian_law, const, denss, densw, dt_sph, grav
!
!        read ndf0 = const(29) but, if we have forgotten take it as 0.01  
!        if (ndf0.le.1.e-2) ndf0 = 1.e-2
!

implicit none

integer(ink) ipoin, idimn, ic_ss  
real   (irk) V_T, expm, C1, Cw, Cs, a_drag, Cn1, ndf0, C10
real   (irk) Vs0, Vw0, Vs, Vw, f1, f2, expp, g , nnn
real   (irk) kwd, kwd0, auxCd

g    = abs(const (1))
ic_ss = 0                    ! index for steady state assympt solution

if (darcinian_law.eq.1) then
    kwd  = const(27) 
    ndf0 = const(29)
    if (ndf0.le.1.e-2) ndf0 = 1.e-2  ! changed MP 5 feb 18
    if (kwd.lt.1.e-10) then
       ic_ss = 1
    elseif (kwd.gt.1.e6) then
       ic_ss = 2
    endif
elseif (darcinian_law.eq.2) then
    V_T  = const(27)
    expm = const(28)  
    ndf0 = const(29)
    if (ndf0.le.1.e-2) ndf0 = 1.e-2
    if (V_T.lt.1.e-10) then
       ic_ss = 1
    elseif (V_T.gt.1.e6) then
       ic_ss = 2
    endif 
    if (V_T.lt.1.e-10) ic_ss = 1
elseif (darcinian_law.eq.3) then
    kwd0 = const(27)
    ndf0 = const(29)
    if (ndf0.le.1.e-2) ndf0 = 1.e-2 
    if (kwd0.lt.1.e-10) then
       ic_ss = 1
    elseif (kwd0.gt.1.e6) then
       ic_ss = 2
    endif
endif

do ipoin = 1, npoin  ! --------------------------------------IPOIN LOOP   

   nnn  = n_DF(ipoin)

   do idimn = 1, ndimn  ! ................................. IDIMN LOOP
      
      if  (itype(ipoin).eq.1) then
           Vs0    = vx  (idimn,ipoin)
           Vw0    = v_sw(idimn,ipoin)
      else if (itype(ipoin).eq.2) then
           Vw0    = vx    (idimn,ipoin) 
           Vs0    = v_sw  (idimn,ipoin)
      endif 

      if ( (n_DF(ipoin).GT.ndf0).AND.(n_DF(ipoin).LT.(1.0-ndf0)) ) THEN ! general porosity case 

          if (ic_ss.eq.0) then                         !    general interaction case                    
              if (darcinian_law.eq.2) then             !      ......  Anderson Law
                  C1   = nnn*(1-nnn)*(denss-densw)*const(1)
                  C1   = C1 / (V_T * nnn**expm )
              elseif (darcinian_law.eq.1) then         !     .......  Darcy
                  C1   = nnn*nnn * densw*g /kwd
              elseif (darcinian_law.eq.3) then         !     ........  Darcy + Karman K
                  kwd  = kwd0 * nnn**3/(1- nnn*nnn  )
                  C1   = nnn*nnn * densw*g /kwd     
              endif 

              Cs   = C1/ (  denss*( 1-n_DF(ipoin) ) )
              Cw   = C1/ (  densw*(   n_DF(ipoin) ) ) 
              a_drag = (  densw*( n_DF(ipoin) ))/( denss*( 1-n_DF(ipoin) ))       
              f1     = (Vs0 + a_drag*Vw0)
              f2     = (Vs0-Vw0)  
              expp   = exp(-(1 + a_drag)* Cw *dt_sph)                                             ! general case
              Vs  = ( f1 + a_drag*f2*expp)/(1.0 + a_drag)
              Vw  = ( f1 -        f2*expp)/(1.0 + a_drag)   !Chuan correction 1
           
          elseif (ic_ss.eq.1) then                              ! direct formula steady state
              a_drag = (  densw*( n_DF(ipoin) ))/( denss*( 1-n_DF(ipoin) )) 
              Vs     = (a_drag/(1.+a_drag))*(vs0 + a_drag*vw0)
              Vw     = vs
           elseif (ic_ss.eq.2) then
              Vs = vs0
              Vw = vw0  
           endif

       elseif (n_DF(ipoin).GE.(1.0-ndf0)) then                              !  -----------fluid limit  n->1 

          if (ic_ss.eq.0) then
             if (darcinian_law.eq.2) then                       !        Anderson Law 
                 auxCd  = ((denss-densw)/denss )*(g/V_T)*(1./nnn*(expm-1))
             elseif (darcinian_law.eq.1) then                   !        Darcy
                 auxCd  = (nnn*nnn/(1-ndf0)) * (densw/denss) * (g/kwd)
             elseif (darcinian_law.eq.3) then                   !        Darcy KK
                 auxCd  = ((1+nnn)/nnn)*(densw/denss)*(g/kwd0) 
             endif 
             Vs  = Vw0 + (Vs0-Vw0)*exp(-auxCd*dt_sph)   
             Vw  = Vw0

          elseif (ic_ss.eq.1) then   ! direct formula steady state
             vw = vw0
             vs = vw0                ! mp5 feb 18
          elseif (ic_ss.eq.2) then
             vw = vw0
             vs = vs0
          endif

       elseif     (n_DF(ipoin).LE.ndf0 ) then                               !  ----------------  solid limit
            Vs     = Vs0  
            Vw     = Vs0
       endif
         
       if  (itype(ipoin).eq.1) then
            vx  (idimn,ipoin)    = Vs 
            v_sw(idimn,ipoin)    = Vw      
       else if (itype(ipoin).eq.2) then
            vx    (idimn,ipoin)  = Vw   
            v_sw  (idimn,ipoin)  = Vs 
       endif 
         
    enddo               ! idimn ---------------------------- 
   
enddo                   ! ipoin            ----------------------------


END SUBROUTINE Get_Drag_SW_OLD  

!-------------------------------------------------------------------

       Subroutine Get_Erosion_Rate_SW   
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes erosion_rate dh/dt
!        Gets erosion as a rate of change of mass dm=(m/rho)*drho
!         drho = Es*h*vx   where Es is [L]-1  around 1.e-3, 5.e-3...
!         Erosion rate: dh/dt = Es*h*vx
!     using: 
!        xx,yy coordinates of SPH node
!        Zp,Zpx, Zpy topo Z and slope
!        h & vmod   depth and mod. of velocity
!     In the main region, we use const
!        const(21)=1 & const(22) is Hungr constant
!        (17) denss (soil grains) 
!        (18) densw (pore fluid)
!        (21) type of erosion law 1 Hungr
!        (22) Hungr constant
!     DFs identified by ic_ws_interact = 2
!        use const(25) = basal porosity n0_eros
!            basal_props(10, index_topo_zone) n0_eros

implicit none

integer(ink) eros_law, index_topo_zone, eros_geom
integer(ink) ipoi2Dx, ipoi2Dy, ipoi2D, n1, n2, n3, n4

real   (irk) Keros  ! this is the erosion constant (either topo or constant)
real   (irk) tmp, tmp1, auxT
real   (irk) Cstar, denss, densw, vol_conc  
real   (irk) cosa, sina, gradZv, theta, thetae, xnum, xden, tanthetae
real   (irk) C_basal, tanfi_basal, beta_basal
real   (irk) Taub, Tauf
real   (irk) n0_eros, fact_n0         ! DFs  mp 9 feb 18

!      ------  Obtain erosion in main region (only Hungr & TB so far)

erosion_rate = 0.0
eros_law = const (21) + 0.1
if (eros_law.eq.1) then
    Keros    = const (22)
    erosion_rate = Keros * h * vmod
elseif (eros_law.eq.2) then
    Keros    = const (22)
    Cstar    = const (23)
    denss    = const (17)
    densw    = const (18)
    if (denss.le.1.e-3.or.densw.le.1.e-3) then
       write (*,*) ' C17 denss or C18 densw are zero'
       PAUSE
    endif
    vol_conc = (dens-densw)/(denss-densw)
    erosion_rate = 0.0
    if (vmod.gt.1.e-5) then
       cosa      = vxx/vmod
       sina      = vyy/vmod
       gradZv    = Zpx*cosa + Zpy*sina
       theta     = atan2 (-gradZv,1. )
       xnum      = (denss-densw)*Vol_conc
       xden      = xnum + densw
       xnum      = xnum*tanfi
       tanthetae = xnum/xden
       thetae = atan(tanthetae) ! atan2 (xden, xnum)
       if (theta.gt.thetae) then
          erosion_rate = Keros*Cstar*vmod*tan(theta-thetae)
       endif
    endif
elseif (eros_law.eq.3) then   ! Saeid, April 2020 CHGD
    Keros       = const (22)
    C_basal     = const (23)
    tanfi_basal = const (24)  
    beta_basal  = const (25)  ! pwp factor
    denss       = const (17)
    densw       = const (18)
    if (denss.le.1.e-3.or.densw.le.1.e-3) then
       write (*,*) ' C17 denss or C18 densw are zero'
       PAUSE
    endif 
    erosion_rate = 0.0
    if (vmod.gt.1.e-5) then
       Taub = dens_bar*(frictx*frictx + fricty*fricty )**0.5
       Tauf = C_basal + dens_bar_eff*cgra*h*(1.-beta_basal)*tanfi_basal
       if (Taub.gt.Tauf.and.vmod.gt.1.e-3) then
          erosion_rate =  Keros*(Taub-Tauf)/dens_bar_eff
       endif
    endif
! else
!   write(*,*) ' entered get eros sw w/o proper eros law', eros_law
!   PAUSE
endif

!      -------  DFs modify erosion_rate using porosity of basal surface DFs only  mp 9 feb 18

if (ic_ws_interact.eq.2) then
   n0_eros = 0
   if (nconst.gt.45) then   ! MP 3rd April 2018
     n0_eros = const(46)  
   endif
   fact_n0 = 1
   if (ipoin.le.npoin_s) then
       if(n0_eros.lt.1.e-6)  then
          fact_n0 = ( 1- n_DF(ipoin) )
       else
          fact_n0 = 1-n0_eros
       endif
    elseif ((ipoin.gt.npoin_s).AND.(ipoin.le.(npoin_s+npoin_w))) then
       if(n0_eros.lt.1.e-6)  then
          fact_n0 = n_DF(ipoin)
       else
          fact_n0 = n0_eros
       endif
    endif
    erosion_rate = fact_n0 * erosion_rate
endif

!      ------  Obtain erosion in topo regions

IF (ic_basal_erosion.eq.1) THEN
!     erosion_rate = 0.0       MP CHANGED 23th June. Not necessary, done at beginning of subr
   call Which_Topo_Zone(xx, yy, 2, index_topo_zone)
   if (index_topo_zone.gt.0) then
      eros_law  = basal_law (index_topo_zone)
      eros_geom = basal_props (1,index_topo_zone) + 0.01
      if (eros_law.eq.0) then
          erosion_rate  = 0
          if (allocated(Activ_Pwp_Stef_Eros_Pt )) then
              Activ_Pwp_Stef_Eros_Pt (ipoin)= 0    !     ch MP 4th May 20121 to skip in not eros basal zones
          endif
      elseif (eros_law.eq.1) then
          Keros         = basal_props (7,index_topo_zone) 
          erosion_rate  = Keros * h * vmod
      elseif (eros_law.eq.2) then
          Keros    = basal_props (7,index_topo_zone) 
          Cstar    = basal_props (8,index_topo_zone) 
          denss    = const (17)
          densw    = const (18)
          if (denss.le.1.e-3.or.densw.le.1.e-3) then
             write (*,*) ' C17 denss or C18 densw are zero'
             PAUSE
          endif
          vol_conc = (dens-densw)/(denss-densw)
          erosion_rate = 0.0
          if (vmod.gt.1.e-5) then
             cosa = vxx/vmod
             sina = vyy/vmod
             gradZv = Zpx*cosa + Zpy*sina
             theta  = atan2 ( -gradZv,1.)
             xnum   = (denss-densw)*Vol_conc
             xden   = xnum + densw
             xnum   = xnum*tanfi
             thetae = atan2 (xden, xnum)
             if (theta.gt.thetae) then
                erosion_rate = Keros*Cstar*vmod*tan(theta-thetae)
             endif      
          else
             write(*,*) ' from get erosion rate SW, C21 should be 1'
             PAUSE
          endif
      elseif (eros_law.eq.3) then
         Keros       = basal_props (7, index_topo_zone)
         C_basal     = basal_props (8, index_topo_zone)
         tanfi_basal = basal_props (9, index_topo_zone)  
         beta_basal  = basal_props (10,index_topo_zone)  ! pwp factor
         denss       = const (17)
         densw       = const (18)
         if (denss.le.1.e-3.or.densw.le.1.e-3) then
            write (*,*) ' C17 denss or C18 densw are zero'
            PAUSE
         endif 
         erosion_rate = 0.0
         if (vmod.gt.1.e-5) then
            Taub = dens*h*(frictx*frictx + fricty*fricty )**0.5
            Tauf = C_basal + dens*cgra*h*(1.-beta_basal)*tanfi_basal
            if (Taub.gt.Tauf.and.vmod.gt.1.e-3) then
               erosion_rate = Keros*(Taub-Tauf)/dens
            endif
         endif
      endif

        
!      -------  modify erosion_rate using porosity of basal surface DFs only  mp 9 feb 18

      if (ic_ws_interact.eq.2) then
          n0_eros = basal_props(10, index_topo_zone) 
          fact_n0 = 1
          if (ipoin.le.npoin_s) then
              if(n0_eros.lt.1.e-6)  then
                 fact_n0 = ( 1- n_DF(ipoin) )
              else
                 fact_n0 = 1-n0_eros
              endif
           elseif ((ipoin.gt.npoin_s).AND.(ipoin.le.(npoin_s+npoin_w))) then
              if(n0_eros.lt.1.e-6)  then
                 fact_n0 = n_DF(ipoin)
              else
                 fact_n0 = n0_eros
              endif
           endif
           erosion_rate = erosion_rate * fact_n0  
      endif
             
      if ((eros_geom.eq.101).OR.(eros_geom.eq.102)) then       !   2019 MP CH 29 may         
         ipoi2Dx = ( xx - xming)/deltxg +1
         ipoi2Dy = ( yy - yming)/deltyg +1
         ipoi2D  = (ipoi2Dy-1)*npoigx + ipoi2Dx
         n1      = ipoi2D
         n2      = n1 + 1 
         n3      = n1 + npoigx  
         n4      = n3 + 1
         auxt    = topol(15,n1) + topol(15,n2)+topol(15,n3)+topol(15,n4)
         if (auxt.lt.1.e-6) then
             erosion_rate = 0.0
         endif  
      endif    
            
   endif                !   END LOOP index_topo_zone :GT.0

ENDIF                   !   END LOOP ic_basal_erosion
         

END SUBROUTINE   Get_Erosion_Rate_SW


!-------------------------------------------------------------------

       Subroutine Get_Topo_Basal_Friction_SW   (Tandelta)
       
!-------------------------------------------------------------------
!
!     Located inside ExtForces_SW 
!     computes tan delta 
!     using: 
!        xx,yy coordinates of SPH node
!        Zp,Zpx, Zpy topo Z and slope
!        h & vmod   depth and mod. of velocity

implicit none

integer(ink) frict_law, index_topo_zone

real   (irk) Tandelta  ! this is the basal frict constant  

!      ------  Obtain basal frict in topo regions

Tandelta = -999.
if (ic_basal_friction.eq.1) THEN
   call Which_Topo_Zone(xx, yy, 1, index_topo_zone)
   if (index_topo_zone.gt.0) then
      frict_law = basal_law    (index_topo_zone)
      Tandelta   = basal_props (7,index_topo_zone) 
   endif   
endif   
         

END SUBROUTINE   Get_Topo_Basal_Friction_SW



!-------------------------------------------------------------------

       Subroutine  get_virtual_parts_EXT_SW 
       
!-------------------------------------------------------------------
!
!      

implicit none

!      ------  Boundary particle force and penalty anti-penetration force.

   if (nvirt.gt.0) then    

       rr0 = vp_r0             !   1.25e-5
       dd  = vp_D              !   1.e-2
       p1  = vp_n1             !   12
       p2  = vp_n2             !   4

       current=> last          !    Initializes linked list to last element
      
       DO WHILE (associated(current))
          i = current%pair_i 
          j = current%pair_j 
          if (current%Pint_type.eq.-1) then
             rr = 0.      
             do idimn = 1,ndimn
                dx(idimn) =  x(idimn,i) -  x(idimn,j)
                rr        = rr + dx(idimn)*dx(idimn)
             enddo  
             rr = sqrt(rr)
             if(rr.lt.rr0) then
                f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
                do idimn = 1, ndimn
                   exdvxdt(idimn, i) = exdvxdt(idimn, i) + dd*dx(idimn)*f
                   exdvxdt(idimn, j) = exdvxdt(idimn, j) - dd*dx(idimn)*f
                enddo
             endif
          endif
          current=>current%next
       ENDDO   

   endif   


   if (ic_histogram.eq.1.and.np_bakg_Total.gt.0) then
      exdvxdt(:,1:np_bakg_Total) = 0.0
   endif

deallocate (dx)

END SUBROUTINE  get_virtual_parts_EXT_SW  


   
END SUBROUTINE   ExtForces_SW 


!-------------------------------------------------------------------

       Subroutine Get_Rain_SW
       
!-------------------------------------------------------------------

!      ------  MP Feb 2023 new version for Fs march 2023

implicit none

character(60) text, description 
1005 format (a60)

integer(ink) np_ts  , if_pwp, if_drain, if_rain, type_Val
real   (irk) unit_ts, h_basal,  L_basal ,  tanTh_basal,  kw_basal ,  poros_basal
real   (irk), allocatable::  time_ts (:), val_ts(:)

integer (ink) ic_fct_dens, ic_hrelSat, ipoin,  ipts
real    (irk) ddt, dt1, dt2, dtt, t_tmp, tfct, Intens, I_pwp, I_rain
real    (irk) dens , denss, densw , Sr , alphas, hrelpw, alph 
real    (irk) porosity, dens1, dens2, dens_bar, dens_bar_eff 
real    (irk) hw,   massw,  hDFw, hsatw 
real    (irk) dhw, dmassw, dhsat, dalpha  
real    (irk) hsat, dhSat_Rain, dhSat_Drain
real    (irk) pwp_presc_basal                ! added JMS 14y¡th March 2024

type_Val    = ts_Rain % type_Val
np_ts       = ts_Rain % np_ts
unit_ts     = ts_Rain % unit_ts
tscale_series = ts_Rain % unit_ts   ! 13th march Jw MP SM
if_pwp      = ts_Rain % if_pwp
if_drain    = ts_Rain % if_drain
h_basal     = ts_Rain % h_basal
L_basal     = ts_Rain % L_basal 
tanTh_basal = ts_Rain % tanTh_basal
poros_basal = ts_Rain % poros_basal
kw_basal    = ts_Rain % kw_basal     
pwp_presc_basal =  ts_Rain %  pwp_presc_basal   ! JMS March 2024
I_rain = 0.0; I_pwp = 0.0

if_rain = 0
if (type_Val.eq.1 ) if_rain = 1
if (type_val.eq.2 ) then
   !  if (if_pwp.eq.0) then
   !        write (*,*) ' if_pwp is wrong... ', if_pwp
   !        write (*,*) ' inconsistent with type_val 2 '
   !         write (*,*) ' strike a num key to acknowledge '
   !        write (*,*) ipoin
   !        STOP
   !     endif
endif

allocate (time_ts(np_ts), val_ts(np_ts))
 
DO ipts = 1, np_ts 
   time_ts (ipts) = ts_Rain % time_ts(ipts)
   val_ts  (ipts) = ts_Rain % val_ts  (ipts) 
ENDDO

!      ------  TASK 1  obtain rain intensity etc at time_sph

Intens = 0.0        ! default for rain or pwp intensity
tfct   = 100. 
ipts   = 0
DO while (tfct.GT.0.AND.ipts.lt.np_ts) 
   ipts  = ipts +1
   !t_tmp = time_sph * unit_ts
   t_tmp = time_sph * tscale_code / tscale_series
   dt1   = - t_tmp         + time_ts(ipts) 
   dt2   = time_ts(ipts+1) - t_tmp  
   dtt   =abs(dt1) +  abs(dt2)  
   if (dt1*dt2.LE.0.0) then
      Intens = (abs(dt1)*val_ts(ipts+1) + abs(dt2)*val_ts(ipts))/dtt 
      if (type_Val.eq.1 ) then       ! only for rain Pwp directly in m  MP chgd may 2023
         Intens = Intens/1000        !  rain intensity in m.unit time-13
         Intens = Intens * tscale_code / 3600        ! 
      endif
      tfct  =  dt1*dt2
   endif 
ENDDO

If (type_Val.EQ.1 ) then   !   -------   RAIN infiltration and draining

   !      ------  TASK 2  for all water pics update hw, hsat, alphas, mass
   !                      includes dewatering too (TASK5)

   !      ------  Default values     

   I_Rain = Intens

   ic_fct_dens   = const (20) + 0.1
   dens          = const  (2)
   denss         = const (17)
   densw         = const (18)
   Sr            = const (19)
   alphas        = const (14)
   ic_hrelSat    = const (34) + 0.01

   dhSat_Rain = 0.0 ;   dhSat_Drain= 0.0

   IF (ic_DF.EQ.1) then                !  2 phases DF, Hsat variable and pwp optionals
                                       !  all dens and factors depend on ipoin   
       DO ipoin = npoin_s +1, npoin_s + npoin_w
          alphaS       = hrelSat_DF (ipoin)
          porosity     = n_DF (ipoin)  
          massw        = mass (ipoin)
          hw           = rho  (ipoin)   
          hDFw         = h_DF (ipoin)
          dhSat_Rain   = I_Rain*dt_sph/porosity   !  Intensity in m/s per unit area
          if (if_drain.eq.1) then                 ! drainage trough upper layer (colluvium, etc)
             hsat          =  h_DF(ipoin) * alphaS
             dhSat_Drain   = -Kw_basal*TanTh_basal*hsat*porosity*dt_sph/L_basal
          elseif (if_drain.eq.2) then             ! drainage trough lower saturated layer 
             hsat         = h_basal
             dhSat_Drain  =  -Kw_basal*TanTh_basal*hsat*poros_basaL*dt_sph/L_basal
          endif
          dhsat        = dhSat_Drain + dhSat_Rain
          dalpha       = dhsat / hDFw
          dhw          = dalpha * porosity * hDFw
          dmassw       = massw * dhw /hw  
           
          rho(ipoin)         = rho(ipoin) +  dhw 
          mass (ipoin)       = mass(ipoin) + dmassw
          hrelSat_DF (ipoin) = hrelSat_DF (ipoin) + dalpha               
       ENDDO 
                    !      ------  TASK 3  for all soil pics 
                    !              update hw, hsat, alphas, mass
       call Get_s_at_w_DF_SW  !  note 2023 do not use Sr for hw and hs...only for densities 
   
   ELSEIF (ic_FS.EQ.1) then 

      IF (ic_hrelSat.eq.1) then           !  1 phase FS, Hsat variable and pwp optionals
                                          !  all dens and factors depend on ipoin   
         DO ipoin =  1, npoin
            alphaS       = hrelSat_DF (ipoin)                      !  
            porosity     = n_DF (ipoin) 
            dens2        = (1-porosity)*denss                      !  dens of upper layer (can be unsat)
            dens1        = (1-porosity)*denss + porosity*densw     !  ower layer (sat)   Sr=alphaS
            dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
            dens_bar_eff = dens_bar - densw*alphaS  
            dhSat_Rain   = I_Rain*dt_sph/porosity                  !  Intensity in m/s per unit area
            if (if_drain.eq.1) then                                ! drainage trough upper layer (colluvium, etc)
               hsat       =  rho(ipoin) * alphaS                   ! h_DF(ipoin) * alphaS  SM JW Mp 11th march
               dhSat_Drain= -Kw_basal*TanTh_basal*hsat*porosity*dt_sph/L_basal
            elseif (if_drain.eq.2) then                            ! drainage trough lower saturated layer 
               hsat         = h_basal
               dhSat_Drain  =  -Kw_basal*TanTh_basal*hsat*poros_basaL*dt_sph/L_basal
            endif
            dhsat        = dhSat_Drain + dhSat_Rain
            dalpha       = dhsat / rho(ipoin)
            hrelSat_DF (ipoin) = hrelSat_DF (ipoin) + dalpha               
          ENDDO  
      ENDIF

!      ------  TASK 4  at soil nodes, increase pwp 
!                      because of response of basal layer caused by d sigma zz total
!      ------  Default values     

   IF (If_Pwp.NE.0) THEN

      ic_fct_dens   = const (20) + 0.1
      dens          = const  (2)
      denss         = const (17)
      densw         = const (18)
      Sr            = const (19)
      alphas        = const (14)
      hrelpw        = alphas
      alph          = alphas
      ic_hrelSat    = const (34) + 0.01

      IF (ic_DF.EQ.1) then             !  2 phases DF, Hsat variable and pwp optionals
         DO ipoin = 1, npoin_s
           IF (ic_hrelSat.eq.0) then        !  all densities constants to compute factors
               porosity     = (denss-dens)/(denss-densw)              !
           ELSEIF (ic_hrelSat.eq.1) then 
               porosity     = n_DF (ipoin)
           ENDIF
           dens2        = (1-porosity)*denss                      !  -- dens of upper layer (can be unsat)
           dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
           alphas       = hrelSat_DF(ipoin)                       ! 
           dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
           dens_bar_eff = dens_bar - densw*alphaS                 !     
       
           Uaux (1,ipoin) = Uaux (1,ipoin) + ((dens1-dens2)/ dens_bar_eff) * dhsat          
         ENDDO         
      ELSEIF (ic_FS.eq.1) then              !  1 phase FS, Hsat variable and pwp optionals
         DO ipoin = 1, npoin
            alphaS       = hrelSat_DF (ipoin)                      !  
            porosity     = n_DF       (ipoin)  
            dens2        = (1-porosity)*denss                      !  -- dens of upper layer (can be unsat)
            dens1        = (1-porosity)*denss + porosity*densw     ! 
            dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
            dens_bar_eff = dens_bar - densw*alphaS                 !     
       
           Uaux (1,ipoin) = Uaux (1,ipoin) + ((dens1-dens2)/ dens_bar_eff) * dhsat          
         ENDDO  
      ENDIF
   ENDIF

ENDIF   !      -------------------- endif RAIN and DRAIN

ELSEIF (type_Val.EQ.2) Then
 
   ts_Rain %  pwp_presc_basal = Intens   !  changed JMS 14th March 2024
 
ENDIF


END Subroutine Get_Rain_SW


!-------------------------------------------------------------------

       Subroutine Get_Pwp_SW ( drho, dn )   
       
!-------------------------------------------------------------------

!    version MP March 2023 Fs with ic_fctdens = 1 and sat zone

implicit none

real   (irk), intent(IN), optional ::  dn   (npoin)
real   (irk), intent(IN)           ::  drho (npoin)

real,allocatable,SAVE:: pwp(:), pwp0(:)
real,allocatable     :: wi (:)

logical if_dn

integer(ink)  ipoin, index_topo_zone, i, j
integer(ink)  ntpwp, it, ip 
integer(ink)  LBc1, LBc2, npwp
integer(ink)  iexit, itopo_zones
integer(ink)  drag_law                ! 1..darcy  2 Anderson  3 karmanKozeny

real   (irk)  dtpwp, dxpwp, dxpwp2, local_t 
real   (irk)  hh,   dhh,  zz,  zzhh,    xi,    hpwp,  dhhstep
real   (irk)  Bfact, Cv,  hinv, comp
real   (irk)  VBc1, VBc2
real   (irk)  beta_inflow   !  when VBc1 is -989 then inflow from the bottom, P is b*liq MP Sept 2021
real   (irk)  xx,yy,Zp,Zpx,Zpy, ZZmin, ZZmax, ZZdiscrd
real   (irk)  dpwp_dn                                   ! aux it is pwp generated by dn
real   (irk)  nn, dnn, nn1, factpwp, dpwp_dn_step

real    (irk) cgra, cgra05              ! grav and 0.5 grav 

                                        ! ------------- DENSITIES ------------- 
                                        !
real    (irk) dens, denss, densw        ! ic_fct_dens = C20  0 no correction  1.. correct 
real    (irk) porosity, Sr              ! Porosity and degree of saturation
real    (irk) dens1, dens2              ! mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    ! average and effective average dens of 2 layers
real    (irk) fact_pwp_max,fact_pwp_min ! factors limiting (1-dpwp/h) in pwp computa
real    (irk) fpwp_max, fpwp_min        ! fact. limiting pwp:  h*fact
real    (irk) alphaS, alphaS2           ! ratio hsat/h, from alphaS. If below 0.1 -> 0

real    (irk) Kv, fact_Kv               !  volum stiffness skeleton for DFs pwp nad factor Xit
real    (irk) kw_bar                    !  permeability  coeff LT-1
real    (irk) ndf0, kwd, kwd0           ! n.GT.ndf0,  kwd darcy LT-1, kwdo initial kw Karman
real    (irk) V_T, exp_m, Cv0           ! V_T C27 Anderson, exp_m C28 Anderson, Cvo before* n's

integer (ink) ic_fct_dens               !  0  classic  1 dens correction 2 DFs 
real    (irk) fact_dens                 !  if 1 and icpwp=0 computed as  fact_dens = (dens-densw)/dens
                                        !           icpwp=  computed using hs/h, porosity...
                                        !  if 2  for DF   fact_dens = (denss-densw)/denss 
real    (irk) betaw                     ! aux for dPwp  
real    (irk) sumdhh, fact_a            ! aux vars    

real    (irk) btw, Cj, f1, f2, f3       ! stefan aux vars. Note btw is erosion_rate_V(ipoin) March 2021
real    (irk) dtpwpE                    ! dt for the advection problem  
real    (irk) VB_adv                    ! tmp pwp rel to liquef at eroded soil   

integer (ink) Flag_q_loc, Flag_q_glob   ! aux flags for flux MP CHGD Aug 2021   
integer (ink) iVBc1                     !  integer value of VBc1    

real   (irk)   fact_s2w 

!      -------  Jan 2924 solid slab, after breaking pwp=0
!               law 0, pwp is made zero and then return
!               otherwise, compute factor  to use it at the end

if (ic_end_solid.eq.1) then        !  solid behaviour via decreasing p 
    if (law_end_solid.eq.0) then   !  MP 30th Jan 2024
       fact_s2w = 0.0  
       if (time_sph.GT.t1_end_solid) then
           fact_s2w = 1.0
           if (allocated (Uaux)) Uaux = 0.0
           if (allocated (pwp))  pwp = 0.0
           RETURN
       endif
    elseif ( law_end_solid.eq.1 ) then
       fact_s2w = 1.0
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

!      ------  initialize and allocate

if(.NOT.allocated(wi) ) allocate (wi(npoin))
wi = 0.0

!      ------  Constants

cgra        = const (1)
cgra05      = 0.5*cgra

if (ic_FS.eq.1) then                 ! MP for FS slow March 2023
   call mini_FSs_initial 
else
   call mini_DFs_dens_and_pwp_factors 
endif 

Bfact    = const(13)
Cv       = 4*Bfact/(pi*pi)
alphaS   = const(14)
alphaS2  = const(14)*const(14) 
comp     = const(15)

npwp     = nUaux
if (.not.allocated (pwp)) then 
   allocate (pwp (npwp) )
   allocate (pwp0(npwp) )
   pwp  = 0.0
   pwp0 = 0.0
endif

drag_law = const(26)
if_dn       = present (dn)

if (drag_law.GT.0) then          ! if (if_dn) then MP march 2023 we could use laws for FS 
   Kv       = const(31)          ! Volum stiffness skeleton, used for pwp dfs from dn
   fact_Kv  = const(32)
   Kv       = Kv * fact_Kv 
   drag_law = const(26)          ! corrected Saeed 29 may 2018
   if (drag_law.EQ.1) then       !  Darcy
      Kwd   = const(27)
      Cv    = (Kwd*Kv)/(densw*cgra)
      Bfact = Cv*pi*pi/4.
   elseif (drag_law.EQ.2) then        !  Anderson
      V_T   = const(27)
      exp_m = const(28)
      ndf0  = const(29)
      Cv0   = Kv*V_T / ( (denss-densw)*cgra)
   elseif (drag_law.EQ.3) then        ! Karman Kozeny
      Kwd   = const(27)
      Cv0   = (Kwd*Kv)/(densw*cgra)
   endif    
endif 

if (ic_basal_flux.EQ.1) then
   basal_flux = 0.0
endif

IF (npwp.gt.1) THEN             !  icpwp=2, FD mesh   or npwp > 1  ------------------------------- 

   Flag_q_glob = 0              !  CHGD MP Aug 2021 Signals flux in rack at any point

   DO ipoin = 1, npoin_s

      Flag_q_loc = 0              !  CHGD MP Aug 2021 Signals flux in rack at this point

      if ( if_Out_Domain(ipoin).eq.1) then
          Uaux(:,ipoin) = 0.0
      endif
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      
      xx = x(1,ipoin)
      yy = (yming+ymaxg)/2   
      if (ndimn.eq.2) yy = x(2,ipoin)   

      LBc1     = BC_pwp_type(1)
      VBc1     = BC_pwp_val (1)
      LBc2     = BC_pwp_type(2)
      VBc2     = BC_pwp_val (2) 
                                          !  ------  TOPO dependent BCs        
      if (ic_basal_PwpBCs.eq.1) then      !  simple version. Update later     
                                          !  using basal_type                                         
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
          call Which_Topo_Zone(xx, yy, 3, index_topo_zone)
          if (index_topo_zone.gt.0) then    ! 25 sept 2018 MP change
              LBc1  = basal_props(7,index_topo_zone)
              VBc1  = basal_props(8,index_topo_zone)
              if (ic_basal_flux.EQ.1) then           !  Saeidmt 28Oct2021
                  Flag_q_loc = 1; Flag_q_glob = 1    !  CHGD MP Sept 2021
              endif
          endif
          
      endif
   
      hh       = rho (ipoin)
      if (ic_FS.EQ.1) then               !  MP March 2023 for FS slow
         alphaS = hrelSat_DF(ipoin)
      endif
      hpwp     = hh*alphaS               ! either defined by default or here FS MArch 2023 MP
      dhh      = drho(ipoin)     
      dpwp_dn  = 0.0                     ! CH mp 4th July 17

      IF (if_dn) THEN                    ! DFs computations ------------------------------------

           if (ic_hrelSat.EQ.1) then                  !  in sub mini we have obtained ct densities
              alphas        = hrelSat_DF (ipoin)      !  here we deal with variable saturated layer 
              dens1         = dens1_DF   (ipoin)
              dens2         = dens2_DF   (ipoin)
              dens_bar      = dens_bar_DF (ipoin)
              dens_bar_eff  = dens_bar_eff_DF (ipoin)
              fact_pwp_max  = fact_pwp_max_DF (ipoin)
              fact_pwp_min  = fact_pwp_min_DF (ipoin)
              fpwp_max      = fpwp_max_DF (ipoin)
              fpwp_min      = fpwp_min_DF (ipoin)
          endif

           dpwp_dn      = 0.0
           nn           = n_DF (ipoin)
           dnn          = dn (ipoin)
           nn1          = max ( (1-nn), 0.05)
           hh           = h_DF(ipoin)
           dhh          = drho(ipoin)/nn1
           hpwp         = hh*alphaS     !   CHGD MP aug 2021 hrelpw                               
           dpwp_dn      = - (Kv * dnn)/(nn1*cgra*dens_bar_eff)   ! if dn>0 then dpwpw<0
           
           if (allocated (if_correct)) then
               if (if_correct(ipoin).eq.1) dpwp_dn = 0.0  ! MP CHG for pwp near wall in 1D 3 March 2020
           endif                                          ! modified 23 April20, and MP Sept 2021

           if (drag_law.eq.1) then
                                    !  Cv   already obtained from C13
                                    !  Bfact = Cv*pi*pi/4 already computed by default
              kw_bar =  const(27)   !   Added MP sept 2021
           elseif (drag_law.eq.2) then
              Cv     = Cv0*(nn**(exp_m+1))/nn1   ! ***  wrong should be exp_m+1  MP Feb 2023
              Bfact  = Cv*pi*pi/4
              kw_bar = ((nn**exp_m)/nn1)*densw*V_T/(denss-densw)   !   Added MP sept 2021
           elseif (drag_law.eq.3) then 
              Cv    = Cv0 * (nn**3)/ ( (1+nn)*nn1)
              Bfact = Cv*pi*pi/4
              kw_bar = Const(27)*(nn**3)/ ( (1+nn)*nn1)            !   MP July 2023 BUG
           endif
        
      endif                             ! DFs computations

      if ((VBc1.LT.-900).AND.(npwp.GT.1)) then         ! condition changed of location MP CH 22 Feb19
         ! VBc1     = -(densw/dens_bar_eff)*alphas*hh  ! when TOTAl pwp = 0 at screen  MP changed 19th dec 2018
           VBc1     = -(densw/dens_bar_eff)            ! care, will be * by hh*alphas later CHGD MP Aug 2021
           fpwp_min = VBc1
           if (ic_basal_flux.EQ.1) then                !  MP CHGD Aug 2021 
               flag_q_loc =1
               flag_q_glob=1
           endif
      endif

      dxpwp    = hpwp/(npwp-1)
      dxpwp2   = dxpwp*dxpwp

      if (Cv.GT.1.e-20) then
         dtpwp    = dxpwp2/(2*Cv)                ! Consolidation explicit FTCS Fdiffs  MArch 2021
      else
         dtpwp    = dxpwp2/(2.e-20)
      endif
      
      if (ic_Stef_eros.eq.1) then
         btw      = erosion_rate_V (ipoin)
         if (btw.le.1.e-10) then
            dtpwpE   = dxpwp * 1.e10
         else
            dtpwpE   = dxpwp /btw              ! Advection LW or TG 1 iteration ML
         endif  
         dtpwp    = min (dtpwp, dtpwpE)
      endif
      dtpwp    = safetpwp*dtpwp

      ntpwp    = ceiling (dt_sph/dtpwp)
      dtpwp    = dt_sph / float(ntpwp)

      dhhstep      = dhh/ntpwp
      dpwp_dn_step = dpwp_dn/ntpwp   ! DFs case, divide by ntimestps local
      pwp(:)   = Uaux (:,ipoin)

!      ------  loop in dt's

      local_t  = 0.0

      DO it    = 1, ntpwp
         
         !      -------------  ADVECTION step
           
         if (ic_Activ_Pwp_Stef_eros.EQ.1) Then 
         if (allocated (Activ_Pwp_Stef_Eros_Pt)) then
         if ( Activ_Pwp_Stef_Eros_Pt (ipoin).EQ.1) then  
            VB_adv = BC_adv_Val (1)               !  default. Note, relative to liquef, * by hh
            pwp(1) = VB_adv * hh                  !  at the bottom we prescribe 
            pwp0   = pwp                          !     the pwp as the failing soil pwp
            do ip = 2, npwp-1
               zz   = dxpwp*(ip-1)
               zzhh = hh - zz                     ! dist to free surface
               xi   = zz/hh
               Cj =  btw * (1.0- xi)*dtpwp/dxpwp  !  ******  btw * dtpwp/dxpwp  ! at the bottom C  = 1
               f1 =  0.50 * Cj * (1.00 + Cj)
               f3 = -0.50 * Cj * (1.00 - Cj)
               f2 =  1.00 - Cj*Cj
               pwp(ip) = f1*pwp0(ip-1) + f2*pwp0(ip) + f3*pwp0(ip+1)   !  FD 1 step Lax Wendroff scheme
            enddo     
         endif
         endif
         endif

         !      -------------  CONSOLIDATION step       
      
         pwp0     = pwp
         local_t = local_t + dtpwp

         DO ip = 2, npwp - 1
            zz   = dxpwp*(ip-1)
            zzhh = hh - zz                     ! dist to free surface
            xi   = zz/hh
            
            sumdhh   = dhhstep * (1.-xi)
            if (ic_hrelSat.EQ.1) then
                if (abs(alphas-1.).LT.0.001) then
                     sumdhh = dhhstep * (1.-xi)
                else
                     sumdhh = dhhstep*((1.-alphas)*dens2/dens_bar_eff &                     ! Saeidmt Aug2021   **Added  dhhstep*()
                            + alphas *(dens1-densw)*(1.-xi/alphaS)/dens_bar_eff)
                endif  
            endif

            pwp (ip) = pwp0 (ip) + (Cv*dtpwp/dxpwp2)* (pwp0(ip-1)-2.*pwp0(ip)+pwp0(ip+1)) &
                       + sumdhh + dpwp_dn_step
            
            if     ( pwp(ip).gt.(zzhh*fpwp_max) ) then
                    pwp(ip) =  (zzhh*fpwp_max)                         
            elseif ( pwp(ip).lt.(zzhh*fpwp_min) ) then  !    
                    pwp(ip) =   (zzhh*fpwp_min)               
            endif
         enddo
   
        if (LBc1.eq.1) then             !  Boundary conditions 1..var  2..flux 1 bottom 2..top
           pwp(1) = VBc1 * hh * alphaS  !  CHGD  MP AUG 2021, it was OK only for val=0
        elseif (LBc1.eq.2) then
           pwp (1) = (4.*pwp (2)-pwp (3)-2.*VBc1*dxpwp ) /3.
        endif
        if (pwp(1).gt.(hh)) then   !   corrected
             pwp(1) = hh
        endif
   
        if (LBc2.eq.1) then
           pwp (npwp) = VBc2 * hh   !  changed by MP march 2021, it was OK only for val=0
        elseif (LBc2.eq.2) then
           pwp (npwp) = (4.*pwp (npwp-1)-pwp (npwp-2) + 2.*VBc2*dxpwp ) /3.
        endif
        if (pwp(npwp).gt.(0.0)) then
            pwp(npwp) = 0.0
        endif 
        pwp0 = pwp
     enddo
     
     Uaux(:,ipoin) = pwp(:)

     if (flag_q_loc.EQ.1) then    !  CHGD MP Aug 2021 Compute basal flux @ ipoin   
        basal_flux (ipoin) = -(kw_bar)  * ( -3.*pwp(1) + 4.*pwp(2) -pwp(3) )/ (2.*dxpwp)  ! Saeidmt 28Oct2021        kw_bar has a negative sign
     endif
   
   ENDDO    !  end loop ipoin
   
ELSEIF (npwp.eq.1) then         !  icpwp =3 or npwp =1 -----------------------------------

   Flag_q_glob = 0              !  CHGD MP Aug 2021 Signals flux in rack at any point

   DO ipoin = 1, npoin_s

      Flag_q_loc = 0              !  CHGD MP Aug 2021 Signals flux in rack at this point

      if ( if_Out_Domain(ipoin).eq.1) Uaux(1,ipoin) = 0.0
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      
      xx = x(1,ipoin)
      yy = (yming+ymaxg)/2   
      if (ndimn.eq.2) yy = x(2,ipoin)   

      LBc1     = BC_pwp_type(1)  
      VBc1     = BC_pwp_val (1)
      LBc2     = BC_pwp_type(2)
      VBc2     = BC_pwp_val (2) 
                                                 
      if (ic_basal_PwpBCs.eq.1.and.ic_basal_flux.EQ.1) then      !   ------  TOPO dependent BCs : grids ONLY     
                                                                 !           using basal_type                                         
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
          call Which_Topo_Zone(xx, yy, 3, index_topo_zone)
          if (index_topo_zone.gt.0) then    ! 25 sept 2018 MP change
              LBc1        = basal_props(7,index_topo_zone)
              VBc1        = basal_props(8,index_topo_zone)
              beta_inflow = basal_props(9,index_topo_zone)     !  MP sept 2021 basal pwp for inflow
              Flag_q_loc  = 1; Flag_q_glob = 1
          endif
          
      endif
   
      hh       = rho (ipoin)
      if (ic_FS.EQ.1) then               !  MP March 2023 for FS slow
         alphaS = hrelSat_DF(ipoin)
      endif
      hpwp     = hh*alphaS         !   CHGD MP Aug 2021 hrelpw  
      dhh      = drho(ipoin)       !   Classic case
      dpwp_dn  = 0.0
                        
      if (if_dn) then                     ! DFs computations

         if (ic_hrelSat.EQ.1) then                    !  in sub mini we have obtained ct densities
              alphas        = hrelSat_DF (ipoin)      !  here we deal with variable saturated layer 
              dens1         = dens1_DF   (ipoin)
              dens2         = dens2_DF   (ipoin)
              dens_bar      = dens_bar_DF (ipoin)
              dens_bar_eff  = dens_bar_eff_DF (ipoin)
              fact_pwp_max  = fact_pwp_max_DF (ipoin)
              fact_pwp_min  = fact_pwp_min_DF (ipoin)
              fpwp_max      = fpwp_max_DF (ipoin)
              fpwp_min      = fpwp_min_DF (ipoin)
         endif

         dpwp_dn      = 0.0
         nn           = n_DF (ipoin)
         dnn          = dn (ipoin)
         nn1          = max ( (1-nn), 0.05)
         hh           = h_DF(ipoin)
         dhh          = drho(ipoin)/nn1
         hpwp         = hh*alphaS         ! CHGD MP Aug 2021 
         dpwp_dn      = - (Kv * dnn)/(nn1*cgra*dens_bar_eff)       ! Saeidmt 28Oct2021   dpwp_dn should be negative      if dn>0 then dpwpw<0

         if (allocated (if_correct)) then
             if (if_correct(ipoin).eq.1) dpwp_dn = 0.0       ! MP CHG for pwp near wall in 1D 3 March 2020
         endif                                               ! modified 23 April
          
         if (ic_hrelSat.EQ.0) then                           ! Saeid included 31st March 2020
            dhh = dhh + dpwp_dn  !  DFs 
         else
            xi = 0.0                                         ! Saeid included 31st March 2020
            sumdhh = (1.-alphas)*dens2/dens_bar_eff &
                            + alphas *(dens1-densw)*(1.-xi/alphaS)/dens_bar_eff
            dhh = dhh*sumdhh + dpwp_dn
         endif

         if (drag_law.eq.1) then
            !  Cv   already obtained from C13
            !  Bfact = Cv*pi*pi/4 already computed by default
            kw_bar =  const(27)
         elseif (drag_law.eq.2) then
            Cv     = Cv0*(nn**exp_m)/nn1
            Bfact  = Cv*pi*pi/4
            kw_bar = ((nn**exp_m)/nn1)*densw*V_T/(denss-densw)   !  MP chgd Aug 2021
         elseif (drag_law.eq.3) then 
            Cv    = Cv0 * (nn**3)/ ( (1+nn)*nn1)
            Bfact = Cv*pi*pi/4
            kw_bar = Const(13)*(nn**3)/ ( (1+nn)*nn1)            ! chk
         endif
      endif 
         
      if (hh.lt.comp) then   !  ** see here comp which limits 1/h
         hh       = comp
         hinv     = 1.0/comp !  ** changed 15th July 2014
      else
         hinv     = 1/hh
      endif 

      iVBc1 = int(VBc1)
      if (iVBc1.EQ.-999) then                           ! condition changed of location MP CH 22 Feb19
         ! VBc1     = -(densw/dens_bar_eff)*alphas*hh  ! when TOTAl pwp = 0 at screen  MP changed 19th dec 2018
           VBc1     = -(densw/dens_bar_eff)            ! care, will be * by hh*alphas later CHGD MP Aug 2021
           fpwp_min = VBc1
           if (ic_basal_flux.EQ.1) then                !  MP CHGD Aug 2021 
               flag_q_loc =1
               flag_q_glob=1
               Uaux (1,ipoin)     = VBc1*hh*alphaS      !  we apply the BC as a 1st step
               basal_flux (ipoin) = -(kw_bar) 
           endif
      else if (iVBc1.EQ.-979) then     
           if (ic_basal_flux.EQ.1) then                 
               flag_q_loc =1
               flag_q_glob=1
               Uaux (1,ipoin)     = beta_inflow*hh*alphaS    !  MP Sept 2021 we apply the BC as a 1st step
               basal_flux (ipoin) =  beta_inflow*(kw_bar) 
           endif
      endif     
      
!      Uaux (1,ipoin) = Uaux (1,ipoin) - dt_sph*Bfact*hinv*hinv*Uaux (1,ipoin)/hrelpw2    + dhh
      
      if (ic_basal_flux.NE.1) then                        ! MP Aug 2021 CHGD 
         betaw = Bfact*hinv*hinv/alphas2                  ! BUG MP CHGD 20 Feb 2023  ***** 
         Uaux (1,ipoin) = Uaux(1,ipoin)*exp(-betaw*dt_sph) + (dhh/dt_sph)*(1.-exp(-betaw*dt_sph))/betaw 
      endif
      
      if (Uaux(1,ipoin).gt.hh) then
              Uaux(1,ipoin) = hh
      endif
      if (Uaux(1,ipoin).lt.hh*alphaS*fpwp_min) then
              Uaux(1,ipoin) = hh*alphaS*fpwp_min
      endif


   ENDDO    !   ends ipoin LOOP 
   
ENDIF       !   ends npwp type cases -----------------------------------------------------------------

if (ic_end_solid.eq.1) then        !  solid behaviour via decreasing p   MP 29 01 2024
   Uaux = Uaux * (1.0 - fact_s2w )        ! the factor has been computed at the beginning of subroutine   
endif

if (if_dn) then            !  INTEGRATE AND AVERAGE PWP  -----------------------

    if (npwp.eq.1) then    !  int. along depth Triangle
       do ipoin =1,npoin_s
          nn = n_DF(ipoin)
          nn1 = max ( (1-nn), 0.05)
          dens_bar_eff = (denss-densw)*nn1
          pwp_avgd (ipoin) = Uaux(1,ipoin)*cgra05*dens_bar_eff/ h_DF(ipoin)                  !  Saeidmt 28Oct2021   Add   / h_DF(ipoin)
       enddo
    elseif (npwp.gt.1) then
       do ipoin =1,npoin_s
          nn = n_DF(ipoin)
          nn1 = max ( (1-nn), 0.05)
          dens_bar_eff = (denss-densw)*nn1
          factpwp = cgra
          pwp_avgd (ipoin) = 0.5*(Uaux(1,ipoin)+ Uaux(npwp,ipoin))*dxpwp
          do ip = 1,npwp 
             pwp_avgd(ipoin) = pwp_avgd(ipoin) + Uaux(  ip,ipoin) *dxpwp 
          enddo
          pwp_avgd(ipoin) = pwp_avgd(ipoin) * gravz*dens_bar_eff/ h_DF(ipoin)  !  cgra-> gravz
       enddo
   endif
   
   pwp_avgd (npoin_s+1:npoin_s+npoin_w) = 0.0

   wi = 0.0                ! CHGD MP Aug 2021
   current=> last          !    Initializes linked list to last element 

   DO WHILE ( associated(current) )
         
      IF ( current%Pint_type.eq.1 ) THEN   ! check this is solid - water interact

         i = current%pair_i 
         j = current%pair_j 
         pwp_avgd (j) = pwp_avgd (j) + mass(i)*current%w*(pwp_avgd(i)/rho(i)) 
         wi       (j) = wi(j)        + mass(i)*current%w /rho(i)  

      endif      
      current => current%next            !
  enddo 
    
  do i = npoin_s +1, npoin              ! Saeed, when water is not interacting with soil
     if (wi(i).gt.0) then 
        pwp_avgd (i) = pwp_avgd (i) / wi(i)
     endif
  enddo
  
endif  

if (Flag_q_glob.EQ.1) then
   call Get_hw_from_basal_flux
endif

CONTAINS

!---------------------------------------------------------------

     subroutine Get_hw_from_basal_flux 

!---------------------------------------------------------------
 
!   Subroutine to calculate 
!
!           (a) h_sw (ntotal): hw at soil pts and hs at water points
!   We use:         
!            mass(   ntotv) : Particle masses  
!            current%Pint_type: 
!                 0 for particles of same kind
!                 1 pic I soil pic J water
!            wi(ntotal) used for normalization
!
! --------------------------------------------------------------

implicit none

integer(ink)  I, J , int_sum, iw, ic_input

real   (irk) tmpp 
real   (irk), allocatable:: wi(:)

  
if(.NOT.allocated(wi)) allocate ( wi(ntotal) )

!      ------  Task 1: transfer fuxes from solid in rack to fluid nodes

wi = 0.0
current=> last          !    Initializes linked list to last element 

DO WHILE ( associated(current) )      
   if ( current%Pint_type.eq.1 ) THEN   !  this is solid - water interact
      i = current%pair_i 
      j = current%pair_j 
      basal_flux (j) = basal_flux (j) + mass(i)*current%w*(basal_flux(i)/rho(i)) 
      wi         (j) = wi         (j) + mass(i)*current%w /rho(i)  
   endif      
   current => current%next            !
ENDDO 
    
do i = npoin_s +1, npoin              ! water points only
   if (wi(i).gt.0) then 
      basal_flux (i) = basal_flux (i) / wi(i)
   endif
enddo

!      ------  Task 2: Use flux to update hw at fluid nodes

Do i = npoin_s + 1, npoin_s + npoin_w
   nn      = n_DF (i)
   if ( rho(i).GE. 0.03) then            !  avoid draining almost dry mix
      tmpp    = basal_flux(i) * dt_sph * nn  
      mass(i) = mass(i) + tmpp * mass(i)/rho(i)
      rho (i) = rho (i) + tmpp
   else
      ipoin = ipoin
   endif
enddo

!      ------  Task 3: update data from soil and water

call Get_s_at_w_DF_SW (2)
  
end subroutine Get_hw_from_basal_flux

!-------------------------------------------------------------------

       Subroutine mini_FSs_initial 
       
!-------------------------------------------------------------------

!      ------  Flow slides with vriable degree of sat height MP March 2023 

implicit none

integer(ink) ipoin       !  note  dens, etc defined  in container subrouine
                         !  
!      ------  FSs values.

dens       = const   (2) 
denss      = const  (17)
densw      = const  (18)
Sr         = const  (19)
alphas     = const  (14)
ic_hrelSat = const  (34)

dens_bar     = dens
dens_bar_eff = dens_bar - alphas*densw
fact_dens   = dens_bar_eff/dens_bar

if (.NOT.allocated (hrelSat_DF) )  then
   allocate (hrelSat_DF  (npoin))
   allocate (fact_dens_DF(npoin)) 
   allocate (fact_pwp_max_DF (npoin) , fact_pwp_min_DF (npoin)) 
   allocate (fpwp_max_DF     (npoin) , fpwp_min_DF   (npoin))  
   fact_pwp_max_DF (:) = 0.0 
   fact_pwp_min_DF (:) = 0.0  
   fpwp_max_DF     (:) = 0.0 
   fpwp_min_DF     (:) = 0.0                                          
endif
       
DO ipoin = 1, npoin
    alphaS       = hrelSat_DF (ipoin)                      !   
    porosity     = (denss-dens)/(denss-alphaS*densw)       !  for Fs n does not change. Compute n0
    dens2        = (1-porosity)*denss                      !  dens of upper layer (can be unsat)
    dens1        = (1-porosity)*denss + porosity*densw     !  ower layer (sat)   Sr=alphaS
    dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
    dens_bar_eff = dens_bar - densw*alphaS                 !       
    fact_dens    = dens_bar_eff/dens_bar                   !     multiplies tanFi
    fact_dens_DF    (ipoin) = fact_dens
    if (icpwp.gt.0) then                                   !     If only affects pwp inf      
        fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
        fact_pwp_min = 0.
        fpwp_max     = 1.0
        fpwp_min     = -alphaS*densw/dens_bar_eff
        fact_pwp_max_DF (ipoin) = fact_pwp_max   
        fact_pwp_min_DF (ipoin) = fact_pwp_min
        fpwp_max_DF     (ipoin) = fpwp_max
        fpwp_min_DF     (ipoin) = fpwp_min
    endif
ENDDO 

END Subroutine mini_FSs_initial

!-------------------------------------------------------------------

       Subroutine mini_DFs_dens_and_pwp_factors
       
!-------------------------------------------------------------------

!      ----   

implicit none

integer (ink) ipoin

real    (irk) dens_ws

! real   (irk) dens, denss, densw, Sr, hrelpw, alphas
! real   (irk) porosity, dens1, dens2, dens_bar, dens_bar_eff
! real   (irk) fact_dens
! real   (irk) fact_pwp_max, fact_pwp_min, fpwp_max, fpwp_min

!      ------  Default values     

ic_fct_dens   = const (20) + 0.1
dens          = const  (2)
denss         = const (17)
densw         = const (18)
Sr            = const (19)
alphas        = const (14)
ic_hrelSat    = const (34) + 0.01
fact_pwp_max  = 1.0
fact_pwp_min  = 0.0
fpwp_max      = 1.0
fpwp_min      = 0.0

if ( (ic_hrelSat.EQ.1).AND.(.NOT.allocated (fact_dens_DF)) ) then
      allocate  ( fact_dens_DF(npoin) )   
      allocate  ( dens1_DF(npoin), dens2_DF   (npoin), densw_DF       (npoin) ) 
      allocate  ( dens_DF (npoin), dens_bar_DF(npoin), dens_bar_eff_DF(npoin) )
      allocate  ( fact_pwp_max_DF (npoin), fact_pwp_min_DF (npoin) )
      allocate  ( fpwp_max_DF     (npoin), fpwp_min_DF     (npoin) )
      fact_dens_DF = 0.0; dens1_DF=0.0; dens2_DF=0.0; densw_DF=0.0
      dens_DF = 0.0; dens_bar_DF=0.0; dens_bar_eff_DF=0.0
      fact_pwp_max_DF= 0.0; fact_pwp_min_DF=0.0; fpwp_max_DF=0.0;fpwp_min_DF=0.0
endif

IF (ic_fct_dens.EQ.0) THEN                           !  1 phase, w or w/o pwp, correction in tanfi

    if ( densw.lt.1.e-3) densw = 1000.0
    dens_ws = densw/dens
    if ( alphas.LT.0.01) then                         !  means hrelpw ignored, so it is assumed to be 1
        alphas  = 1.0
    endif
    fact_dens = 1.0                                   ! Default

ELSEIF (ic_fct_dens.EQ.1) then                                  !  1 phase, w or w/o pwp, correction OK
    
    if ( (denss.LT.1.e-3).OR.(densw.lt.1.e-3) ) then       !  Filters zero solid or water densities
       write (*,*) ' denss or densw wrong...', denss, densw
       write (*,*) ' strike a num key to acknowledge '
       read  (*,*) ipoin
       STOP
    endif
    
    if (alphaS.le.0.01) then                                  !  neglect saturated part
       fact_dens = 1.0                                        !  (dens-densw)/dens   
    else
       porosity     = (denss-dens)/(denss-densw)
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !     modified MP 10.01.2019 
       fact_dens    = dens_bar_eff/dens_bar                   !     multiplies tanFi
    endif

    if (icpwp.gt.0) then                                   !     If only affects pwp inf      
       fact_pwp_max = 1. + alphaS*densw/dens_bar_eff       !     ( 1-dpwp_bar/h)   
       fact_pwp_min = 0.
       fpwp_max     = 1.0
       fpwp_min     = -alphaS*densw/dens_bar_eff
    endif 
       
ELSEIF (ic_fct_dens.EQ.2) then                              !  2 phases DF, Hsat variable and pwp optionals

    IF (ic_hrelSat.eq.0) then                       !  all densities constants to compute factors
       porosity     = (denss-dens)/(denss-densw)              !
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !       
       fact_dens    = dens_bar_eff/(denss*(1-porosity))       !     multiplies tanFi
       if (icpwp.gt.0) then                                   !     If only affects pwp inf      
           fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
           fact_pwp_min = 0.
           fpwp_max     = 1.0
           fpwp_min     = -alphaS*densw/dens_bar_eff
       endif 

    ELSEIF (ic_hrelSat.eq.1) then                   !  all dens and factors depend on ipoin   

       DO ipoin = 1, npoin
           alphaS       = hrelSat_DF (ipoin)                      !  MP chgd bug aug 2021 alphas was wrong
           porosity     = n_DF (ipoin)                            !
           dens2        = (1-porosity)*denss + porosity*densw*Sr  !  dens of upper layer (can be unsat)
           dens1        = (1-porosity)*denss + porosity*densw     !  ower layer (sat)   Sr=alphaS
           dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
           dens_bar_eff = dens_bar - densw*alphaS                 !       
           fact_dens    = dens_bar_eff/(denss*(1-porosity))       !     multiplies tanFi
           dens1_DF        (ipoin) = dens1
           dens2_DF        (ipoin) = dens2
           dens_bar_DF     (ipoin) = dens_bar
           dens_bar_eff_DF (ipoin) = dens_bar_eff
           fact_dens_DF    (ipoin) = fact_dens
           if (icpwp.gt.0) then                                   !     If only affects pwp inf      
               fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
               fact_pwp_min = 0.
               fpwp_max     = 1.0
               fpwp_min     = -alphaS*densw/dens_bar_eff
               fact_pwp_max_DF (ipoin) = fact_pwp_max   
               fact_pwp_min_DF (ipoin) = fact_pwp_min
               fpwp_max_DF     (ipoin) = fpwp_max
               fpwp_min_DF     (ipoin) = fpwp_min
           endif
       ENDDO 

    ENDIF

ELSE 
    write (*,*) ' ic_fct_dens given is wrong...', ic_fct_dens
    write (*,*) ' strike a num key to acknowledge '
    read  (*,*) ipoin
    STOP
ENDIF

END Subroutine mini_DFs_dens_and_pwp_factors



END SUBROUTINE   Get_Pwp_SW   

!-------------------------------------------------------------------

       Subroutine Get_Pwp_SW_OLD ( drho, dn )   
       
!-------------------------------------------------------------------


implicit none

real   (irk), intent(IN), optional ::  dn   (npoin)
real   (irk), intent(IN)           ::  drho (npoin)

real,allocatable,SAVE:: pwp(:), pwp0(:)
real,allocatable     :: wi (:)

logical if_dn

integer(ink)  ipoin, index_topo_zone, i, j
integer(ink)  ntpwp, it, ip 
integer(ink)  LBc1, LBc2, npwp
integer(ink)  iexit, itopo_zones
integer(ink)  drag_law                ! 1..darcy  2 Anderson  3 karmanKozeny

real   (irk)  dtpwp, dxpwp, dxpwp2, local_t 
real   (irk)  hh,   dhh,  zz,  zzhh,    xi,    hpwp,  dhhstep
real   (irk)  Bfact, Cv, hrelpw, hrelpw2, hinv, comp
real   (irk)  VBc1, VBc2
real   (irk)  beta_inflow   !  when VBc1 is -989 then inflow from the bottom, P is b*liq MP Sept 2021
real   (irk)  xx,yy,Zp,Zpx,Zpy, ZZmin, ZZmax, ZZdiscrd
real   (irk)  dpwp_dn                                   ! aux it is pwp generated by dn
real   (irk)  nn, dnn, nn1, factpwp, dpwp_dn_step

real    (irk) cgra, cgra05              ! grav and 0.5 grav 

                                        ! ------------- DENSITIES ------------- 
                                        !
real    (irk) dens, denss, densw        ! ic_fct_dens = C20  0 no correction  1.. correct 
real    (irk) porosity, Sr              ! Porosity and degree of saturation
real    (irk) dens1, dens2              ! mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    ! average and effective average dens of 2 layers
real    (irk) fact_pwp_max,fact_pwp_min ! factors limiting (1-dpwp/h) in pwp computa
real    (irk) fpwp_max, fpwp_min        ! fact. limiting pwp:  h*fact
real    (irk) alphaS                    ! ratio hsat/h, from hrelpw. If below 0.1 -> 0

real    (irk) Kv, fact_Kv               !  volum stiffness skeleton for DFs pwp nad factor Xit
real    (irk) kw_bar                    !  permeability  coeff LT-1
real    (irk) ndf0, kwd, kwd0           ! n.GT.ndf0,  kwd darcy LT-1, kwdo initial kw Karman
real    (irk) V_T, exp_m, Cv0           ! V_T C27 Anderson, exp_m C28 Anderson, Cvo before* n's

integer (ink) ic_fct_dens               !  0  classic  1 dens correction 2 DFs 
real    (irk) fact_dens                 !  if 1 and icpwp=0 computed as  fact_dens = (dens-densw)/dens
                                        !           icpwp=  computed using hs/h, porosity...
                                        !  if 2  for DF   fact_dens = (denss-densw)/denss 
real    (irk) betaw                     ! aux for dPwp  
real    (irk) sumdhh, fact_a            ! aux vars    

real    (irk) btw, Cj, f1, f2, f3       ! stefan aux vars. Note btw is erosion_rate_V(ipoin) March 2021
real    (irk) dtpwpE                    ! dt for the advection problem  
real    (irk) VB_adv                    ! tmp pwp rel to liquef at eroded soil   

integer (ink) Flag_q_loc, Flag_q_glob   ! aux flags for flux MP CHGD Aug 2021   
integer (ink) iVBc1                     !  integer value of VBc1    

!      ------  initialize and allocate

if(.NOT.allocated(wi) ) allocate (wi(npoin))
wi = 0.0

!      ------  Constants

cgra        = const (1)
cgra05      = 0.5*cgra

if (ic_FS.eq.1) then                 ! MP for FS slow March 2023
   call mini_FSs_initial_OLD
else
   call mini_DFs_dens_and_pwp_factors_OLD 
endif 

Bfact    = const(13)
Cv       = 4*Bfact/(pi*pi)
hrelpw   = const(14)
hrelpw2  = const(14)*const(14) 
comp     = const(15)

npwp     = nUaux
if (.not.allocated (pwp)) then 
   allocate (pwp (npwp) )
   allocate (pwp0(npwp) )
   pwp  = 0.0
   pwp0 = 0.0
endif

if_dn       = present (dn)
if (if_dn) then
   Kv       = const(31)               ! Volum stiffness skeleton, used for pwp dfs from dn
   fact_Kv  = const(32)
   Kv       = Kv * fact_Kv 
   drag_law = const(26)               ! corrected Saeed 29 may 2018
   if (drag_law.EQ.1) then            !  Darcy
      Kwd   = const(27)
      Cv    = (Kwd*Kv)/(densw*cgra)
      Bfact = Cv*pi*pi/4.
   elseif (drag_law.EQ.2) then        !  Anderson
      V_T   = const(27)
      exp_m = const(28)
      ndf0  = const(29)
      Cv0   = Kv*V_T / ( (denss-densw)*cgra)
   elseif (drag_law.EQ.3) then        ! Karman Kozeny
      Kwd   = const(27)
      Cv0   = (Kwd*Kv)/(densw*cgra)
   endif    
endif 

if (ic_basal_flux.EQ.1) then
   basal_flux = 0.0
endif

IF (npwp.gt.1) THEN             !  icpwp=2, FD mesh   or npwp > 1  ------------------------------- 

   Flag_q_glob = 0              !  CHGD MP Aug 2021 Signals flux in rack at any point

   DO ipoin = 1, npoin_s

      Flag_q_loc = 0              !  CHGD MP Aug 2021 Signals flux in rack at this point

      if ( if_Out_Domain(ipoin).eq.1) then
          Uaux(:,ipoin) = 0.0
      endif
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      
      xx = x(1,ipoin)
      yy = (yming+ymaxg)/2   
      if (ndimn.eq.2) yy = x(2,ipoin)   

      LBc1     = BC_pwp_type(1)
      VBc1     = BC_pwp_val (1)
      LBc2     = BC_pwp_type(2)
      VBc2     = BC_pwp_val (2) 
                                          !  ------  TOPO dependent BCs        
      if (ic_basal_PwpBCs.eq.1) then      !  simple version. Update later     
                                          !  using basal_type                                         
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
          call Which_Topo_Zone(xx, yy, 3, index_topo_zone)
          if (index_topo_zone.gt.0) then    ! 25 sept 2018 MP change
              LBc1  = basal_props(7,index_topo_zone)
              VBc1  = basal_props(8,index_topo_zone)
              if (ic_basal_flux.EQ.1) then           !  Saeidmt 28Oct2021
                  Flag_q_loc = 1; Flag_q_glob = 1    !  CHGD MP Sept 2021
              endif
          endif
          
      endif
   
      hh       = rho (ipoin)
      if (ic_FS.EQ.1.AND.ic_hrelsat.eq.1) then  !  MP March 2023 for FS slow
         alphaS = hrelSat_DF(ipoin)
      endif
      hpwp     = hh*alphaS               ! old CHGD MP Aug 2021 hh*hrelpw
      dhh      = drho(ipoin)     
      dpwp_dn  = 0.0                     ! CH mp 4th July 17

      if (if_dn) then                    ! DFs computations

           if (ic_hrelSat.EQ.1) then                  !  in sub mini we have obtained ct densities
              alphas        = hrelSat_DF (ipoin)      !  here we deal with variable saturated layer 
              dens1         = dens1_DF   (ipoin)
              dens2         = dens2_DF   (ipoin)
              dens_bar      = dens_bar_DF (ipoin)
              dens_bar_eff  = dens_bar_eff_DF (ipoin)
              fact_pwp_max  = fact_pwp_max_DF (ipoin)
              fact_pwp_min  = fact_pwp_min_DF (ipoin)
              fpwp_max      = fpwp_max_DF (ipoin)
              fpwp_min      = fpwp_min_DF (ipoin)
          endif

           dpwp_dn      = 0.0
           nn           = n_DF (ipoin)
           dnn          = dn (ipoin)
           nn1          = max ( (1-nn), 0.05)
           hh           = h_DF(ipoin)
           dhh          = drho(ipoin)/nn1
           hpwp         = hh*alphaS     !   CHGD MP aug 2021 hrelpw                               
           dpwp_dn      = - (Kv * dnn)/(nn1*cgra*dens_bar_eff)   ! if dn>0 then dpwpw<0
           
           if (allocated (if_correct)) then
               if (if_correct(ipoin).eq.1) dpwp_dn = 0.0  ! MP CHG for pwp near wall in 1D 3 March 2020
           endif                                          ! modified 23 April20, and MP Sept 2021

           if (drag_law.eq.1) then
                                    !  Cv   already obtained from C13
                                    !  Bfact = Cv*pi*pi/4 already computed by default
              kw_bar =  const(27)   !   Added MP sept 2021
           elseif (drag_law.eq.2) then
              Cv     = Cv0*(nn**(exp_m+1))/nn1   ! ***  wrong should be exp_m+1  MP Feb 2023
              Bfact  = Cv*pi*pi/4
              kw_bar = ((nn**exp_m)/nn1)*densw*V_T/(denss-densw)   !   Added MP sept 2021
           elseif (drag_law.eq.3) then 
              Cv    = Cv0 * (nn**3)/ ( (1+nn)*nn1)
              Bfact = Cv*pi*pi/4
              kw_bar = Const(13)*(nn**3)/ ( (1+nn)*nn1)            !   Added MP sept 2021
           endif
        
      endif                             ! DFs computations

      if ((VBc1.LT.-900).AND.(npwp.GT.1)) then         ! condition changed of location MP CH 22 Feb19
         ! VBc1     = -(densw/dens_bar_eff)*alphas*hh  ! when TOTAl pwp = 0 at screen  MP changed 19th dec 2018
           VBc1     = -(densw/dens_bar_eff)            ! care, will be * by hh*alphas later CHGD MP Aug 2021
           fpwp_min = VBc1
           if (ic_basal_flux.EQ.1) then                !  MP CHGD Aug 2021 
               flag_q_loc =1
               flag_q_glob=1
           endif
      endif

      dxpwp    = hpwp/(npwp-1)
      dxpwp2   = dxpwp*dxpwp
         ! ***************************************
         !  Cv  = 1.0   ONLY for debuggu¡ing using PLK01
         ! ***************************************

      if (Cv.GT.1.e-20) then
         dtpwp    = dxpwp2/(2*Cv)                ! Consolidation explicit FTCS Fdiffs  MArch 2021
      else
         dtpwp    = dxpwp2/(2.e-20)
      endif
      
      if (ic_Stef_eros.eq.1) then
         btw      = erosion_rate_V (ipoin)
         ! ***************************************
         !  btw = 0.0  ONLY for debuggu¡ing using PLK01
         ! ***************************************
         if (btw.le.1.e-10) then
            dtpwpE   = dxpwp * 1.e10
         else
            dtpwpE   = dxpwp /btw              ! Advection LW or TG 1 iteration ML
         endif  
         dtpwp    = min (dtpwp, dtpwpE)
      endif
      dtpwp    = safetpwp*dtpwp

      ntpwp    = ceiling (dt_sph/dtpwp)
      dtpwp    = dt_sph / float(ntpwp)

      dhhstep      = dhh/ntpwp
      dpwp_dn_step = dpwp_dn/ntpwp   ! DFs case, divide by ntimestps local
      pwp(:)   = Uaux (:,ipoin)

!      ------  loop in dt's

      local_t  = 0.0

      DO it    = 1, ntpwp
         
         !      -------------  ADVECTION step
           
         if (ic_Activ_Pwp_Stef_eros.EQ.1) Then 
         if (allocated (Activ_Pwp_Stef_Eros_Pt)) then
         if ( Activ_Pwp_Stef_Eros_Pt (ipoin).EQ.1) then  
            VB_adv = BC_adv_Val (1)               !  default. Note, relative to liquef, * by hh
            pwp(1) = VB_adv * hh                  !  at the bottom we prescribe 
            pwp0   = pwp                          !     the pwp as the failing soil pwp
            do ip = 2, npwp-1
               zz   = dxpwp*(ip-1)
               zzhh = hh - zz                     ! dist to free surface
               xi   = zz/hh
               Cj =  btw * (1.0- xi)*dtpwp/dxpwp  !  ******  btw * dtpwp/dxpwp  ! at the bottom C  = 1
               f1 =  0.50 * Cj * (1.00 + Cj)
               f3 = -0.50 * Cj * (1.00 - Cj)
               f2 =  1.00 - Cj*Cj
               pwp(ip) = f1*pwp0(ip-1) + f2*pwp0(ip) + f3*pwp0(ip+1)   !  FD 1 step Lax Wendroff scheme
            enddo     
         endif
         endif
         endif

         !      -------------  CONSOLIDATION step       
      
         pwp0     = pwp
         local_t = local_t + dtpwp

         DO ip = 2, npwp - 1
            zz   = dxpwp*(ip-1)
            zzhh = hh - zz                     ! dist to free surface
            xi   = zz/hh
            
            sumdhh   = dhhstep * (1.-xi)
            if (ic_hrelSat.EQ.1) then
                if (abs(alphas-1.).LT.0.001) then
                     sumdhh = dhhstep * (1.-xi)
                else
                     sumdhh = dhhstep*((1.-alphas)*dens2/dens_bar_eff &                     ! Saeidmt Aug2021   **Added  dhhstep*()
                            + alphas *(dens1-densw)*(1.-xi/alphaS)/dens_bar_eff)
                endif  
            endif

            pwp (ip) = pwp0 (ip) + (Cv*dtpwp/dxpwp2)* (pwp0(ip-1)-2.*pwp0(ip)+pwp0(ip+1)) &
                       + sumdhh + dpwp_dn_step
            
            if     ( pwp(ip).gt.(zzhh*fpwp_max) ) then
                    pwp(ip) =  (zzhh*fpwp_max)                         
            elseif ( pwp(ip).lt.(zzhh*fpwp_min) ) then  !    
                    pwp(ip) =   (zzhh*fpwp_min)               
            endif
         enddo
   
        if (LBc1.eq.1) then             !  Boundary conditions 1..var  2..flux 1 bottom 2..top
           pwp(1) = VBc1 * hh * alphaS  !  CHGD  MP AUG 2021, it was OK only for val=0
        elseif (LBc1.eq.2) then
           pwp (1) = (4.*pwp (2)-pwp (3)-2.*VBc1*dxpwp ) /3.
        endif
        if (pwp(1).gt.(hh)) then   !   corrected
             pwp(1) = hh
        endif
   
        if (LBc2.eq.1) then
           pwp (npwp) = VBc2 * hh   !  changed by MP march 2021, it was OK only for val=0
        elseif (LBc2.eq.2) then
           pwp (npwp) = (4.*pwp (npwp-1)-pwp (npwp-2) + 2.*VBc2*dxpwp ) /3.
        endif
        if (pwp(npwp).gt.(0.0)) then
            pwp(npwp) = 0.0
        endif 
        pwp0 = pwp
     enddo
     
     Uaux(:,ipoin) = pwp(:)

     if (flag_q_loc.EQ.1) then    !  CHGD MP Aug 2021 Compute basal flux @ ipoin   
        basal_flux (ipoin) = -(kw_bar)  * ( -3.*pwp(1) + 4.*pwp(2) -pwp(3) )/ (2.*dxpwp)  ! Saeidmt 28Oct2021        kw_bar has a negative sign
     endif
   
   ENDDO    !  end loop ipoin
   
ELSEIF (npwp.eq.1) then         !  icpwp =3 or npwp =1 -----------------------------------

   Flag_q_glob = 0              !  CHGD MP Aug 2021 Signals flux in rack at any point

   DO ipoin = 1, npoin_s

      Flag_q_loc = 0              !  CHGD MP Aug 2021 Signals flux in rack at this point

      if ( if_Out_Domain(ipoin).eq.1) Uaux(1,ipoin) = 0.0
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      
      xx = x(1,ipoin)
      yy = (yming+ymaxg)/2   
      if (ndimn.eq.2) yy = x(2,ipoin)   

      LBc1     = BC_pwp_type(1)  
      VBc1     = BC_pwp_val (1)
      LBc2     = BC_pwp_type(2)
      VBc2     = BC_pwp_val (2) 
                                                 
      if (ic_basal_PwpBCs.eq.1.and.ic_basal_flux.EQ.1) then      !   ------  TOPO dependent BCs : grids ONLY     
                                                                 !           using basal_type                                         
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)
          call Which_Topo_Zone(xx, yy, 3, index_topo_zone)
          if (index_topo_zone.gt.0) then    ! 25 sept 2018 MP change
              LBc1        = basal_props(7,index_topo_zone)
              VBc1        = basal_props(8,index_topo_zone)
              beta_inflow = basal_props(9,index_topo_zone)     !  MP sept 2021 basal pwp for inflow
              Flag_q_loc  = 1; Flag_q_glob = 1
          endif
          
      endif
   
      hh       = rho (ipoin)
      hpwp     = hh*alphaS         !   CHGD MP Aug 2021 hrelpw  
      dhh      = drho(ipoin)       !   Classic case
      dpwp_dn  = 0.0
                        
      if (if_dn) then                     ! DFs computations

         if (ic_hrelSat.EQ.1) then                    !  in sub mini we have obtained ct densities
              alphas        = hrelSat_DF (ipoin)      !  here we deal with variable saturated layer 
              dens1         = dens1_DF   (ipoin)
              dens2         = dens2_DF   (ipoin)
              dens_bar      = dens_bar_DF (ipoin)
              dens_bar_eff  = dens_bar_eff_DF (ipoin)
              fact_pwp_max  = fact_pwp_max_DF (ipoin)
              fact_pwp_min  = fact_pwp_min_DF (ipoin)
              fpwp_max      = fpwp_max_DF (ipoin)
              fpwp_min      = fpwp_min_DF (ipoin)
         endif

         dpwp_dn      = 0.0
         nn           = n_DF (ipoin)
         dnn          = dn (ipoin)
         nn1          = max ( (1-nn), 0.05)
         hh           = h_DF(ipoin)
         dhh          = drho(ipoin)/nn1
         hpwp         = hh*alphaS         ! CHGD MP Aug 2021 
         dpwp_dn      = - (Kv * dnn)/(nn1*cgra*dens_bar_eff)       ! Saeidmt 28Oct2021   dpwp_dn should be negative      if dn>0 then dpwpw<0

         if (allocated (if_correct)) then
             if (if_correct(ipoin).eq.1) dpwp_dn = 0.0       ! MP CHG for pwp near wall in 1D 3 March 2020
         endif                                               ! modified 23 April
          
         if (ic_hrelSat.EQ.0) then                           ! Saeid included 31st March 2020
            dhh = dhh + dpwp_dn  !  DFs 
         else
            xi = 0.0                                         ! Saeid included 31st March 2020
            sumdhh = (1.-alphas)*dens2/dens_bar_eff &
                            + alphas *(dens1-densw)*(1.-xi/alphaS)/dens_bar_eff
            dhh = dhh*sumdhh + dpwp_dn
         endif

         if (drag_law.eq.1) then
            !  Cv   already obtained from C13
            !  Bfact = Cv*pi*pi/4 already computed by default
            kw_bar =  const(27)
         elseif (drag_law.eq.2) then
            Cv     = Cv0*(nn**exp_m)/nn1
            Bfact  = Cv*pi*pi/4
            kw_bar = ((nn**exp_m)/nn1)*densw*V_T/(denss-densw)   !  MP chgd Aug 2021
         elseif (drag_law.eq.3) then 
            Cv    = Cv0 * (nn**3)/ ( (1+nn)*nn1)
            Bfact = Cv*pi*pi/4
            kw_bar = Const(13)*(nn**3)/ ( (1+nn)*nn1)            ! chk
         endif
      endif 
         
      if (hh.lt.comp) then   !  ** see here comp which limits 1/h
         hh       = comp
         hinv     = 1.0/comp !  ** changed 15th July 2014
      else
         hinv     = 1/hh
      endif 

      iVBc1 = int(VBc1)
      if (iVBc1.EQ.-999) then                           ! condition changed of location MP CH 22 Feb19
         ! VBc1     = -(densw/dens_bar_eff)*alphas*hh  ! when TOTAl pwp = 0 at screen  MP changed 19th dec 2018
           VBc1     = -(densw/dens_bar_eff)            ! care, will be * by hh*alphas later CHGD MP Aug 2021
           fpwp_min = VBc1
           if (ic_basal_flux.EQ.1) then                !  MP CHGD Aug 2021 
               flag_q_loc =1
               flag_q_glob=1
               Uaux (1,ipoin)     = VBc1*hh*alphaS      !  we apply the BC as a 1st step
               basal_flux (ipoin) = -(kw_bar) 
           endif
      else if (iVBc1.EQ.-979) then     
           if (ic_basal_flux.EQ.1) then                 
               flag_q_loc =1
               flag_q_glob=1
               Uaux (1,ipoin)     = beta_inflow*hh*alphaS    !  MP Sept 2021 we apply the BC as a 1st step
               basal_flux (ipoin) =  beta_inflow*(kw_bar) 
           endif
      endif     
      
!      Uaux (1,ipoin) = Uaux (1,ipoin) - dt_sph*Bfact*hinv*hinv*Uaux (1,ipoin)/hrelpw2    + dhh
      
      if (ic_basal_flux.NE.1) then                        ! MP Aug 2021 CHGD
 !        betaw = Bfact*hinv*hinv/hrelpw2                 ! if NOT rack, we use this 
         betaw = Bfact*hinv*hinv/(alphas*alphas)          ! BUG MP CHGD 20 Feb 2023  ***** 
         Uaux (1,ipoin) = Uaux(1,ipoin)*exp(-betaw*dt_sph) + (dhh/dt_sph)*(1.-exp(-betaw*dt_sph))/betaw 
      endif
      
      if (Uaux(1,ipoin).gt.hh) then
              Uaux(1,ipoin) = hh
      endif
      if (Uaux(1,ipoin).lt.hh*alphaS*fpwp_min) then
              Uaux(1,ipoin) = hh*alphaS*fpwp_min
      endif


   ENDDO    !   ends ipoin LOOP 
   
ENDIF       !   ends npwp type cases -----------------------------------------------------------------

if (if_dn) then            !  INTEGRATE AND AVERAGE PWP  -----------------------

    if (npwp.eq.1) then    !  int. along depth Triangle
       do ipoin =1,npoin_s
          nn = n_DF(ipoin)
          nn1 = max ( (1-nn), 0.05)
          dens_bar_eff = (denss-densw)*nn1
          pwp_avgd (ipoin) = Uaux(1,ipoin)*cgra05*dens_bar_eff/ h_DF(ipoin)                  !  Saeidmt 28Oct2021   Add   / h_DF(ipoin)
       enddo
    elseif (npwp.gt.1) then
       do ipoin =1,npoin_s
          nn = n_DF(ipoin)
          nn1 = max ( (1-nn), 0.05)
          dens_bar_eff = (denss-densw)*nn1
          factpwp = cgra
          pwp_avgd (ipoin) = 0.5*(Uaux(1,ipoin)+ Uaux(npwp,ipoin))*dxpwp
          do ip = 1,npwp 
             pwp_avgd(ipoin) = pwp_avgd(ipoin) + Uaux(  ip,ipoin) *dxpwp 
          enddo
          pwp_avgd(ipoin) = pwp_avgd(ipoin) * gravz*dens_bar_eff/ h_DF(ipoin)  !  cgra-> gravz
       enddo
   endif
   
   pwp_avgd (npoin_s+1:npoin_s+npoin_w) = 0.0

   wi = 0.0                ! CHGD MP Aug 2021
   current=> last          !    Initializes linked list to last element 

   DO WHILE ( associated(current) )
         
      IF ( current%Pint_type.eq.1 ) THEN   ! check this is solid - water interact

         i = current%pair_i 
         j = current%pair_j 
         pwp_avgd (j) = pwp_avgd (j) + mass(i)*current%w*(pwp_avgd(i)/rho(i)) 
         wi       (j) = wi(j)        + mass(i)*current%w /rho(i)  

      endif      
      current => current%next            !
  enddo 
    
  do i = npoin_s +1, npoin              ! Saeed, when water is not interacting with soil
     if (wi(i).gt.0) then 
        pwp_avgd (i) = pwp_avgd (i) / wi(i)
     endif
  enddo
  
endif  

if (Flag_q_glob.EQ.1) then
   call Get_hw_from_basal_flux_OLD
endif

CONTAINS


!---------------------------------------------------------------

     subroutine Get_hw_from_basal_flux_OLD

!---------------------------------------------------------------
 
!   Subroutine to calculate 
!
!           (a) h_sw (ntotal): hw at soil pts and hs at water points
!   We use:         
!            mass(   ntotv) : Particle masses  
!            current%Pint_type: 
!                 0 for particles of same kind
!                 1 pic I soil pic J water
!            wi(ntotal) used for normalization
!
! --------------------------------------------------------------

implicit none

integer(ink)  I, J , int_sum, iw, ic_input

real   (irk) tmpp 
real   (irk), allocatable:: wi(:)

  
if(.NOT.allocated(wi)) allocate ( wi(ntotal) )

!      ------  Task 1: transfer fuxes from solid in rack to fluid nodes

wi = 0.0
current=> last          !    Initializes linked list to last element 

DO WHILE ( associated(current) )      
   if ( current%Pint_type.eq.1 ) THEN   !  this is solid - water interact
      i = current%pair_i 
      j = current%pair_j 
      basal_flux (j) = basal_flux (j) + mass(i)*current%w*(basal_flux(i)/rho(i)) 
      wi         (j) = wi         (j) + mass(i)*current%w /rho(i)  
   endif      
   current => current%next            !
ENDDO 
    
do i = npoin_s +1, npoin              ! water points only
   if (wi(i).gt.0) then 
      basal_flux (i) = basal_flux (i) / wi(i)
   endif
enddo

!      ------  Task 2: Use flux to update hw at fluid nodes

Do i = npoin_s + 1, npoin_s + npoin_w
   nn      = n_DF (i)
   if ( rho(i).GE. 0.03) then            !  avoid draining almost dry mix
      tmpp    = basal_flux(i) * dt_sph * nn  
      mass(i) = mass(i) + tmpp * mass(i)/rho(i)
      rho (i) = rho (i) + tmpp
   else
      ipoin = ipoin
   endif
enddo

!      ------  Task 3: update data from soil and water

call Get_s_at_w_DF_SW (2)
  
end subroutine Get_hw_from_basal_flux_OLD

!-------------------------------------------------------------------

       Subroutine mini_FSs_initial_OLD 
       
!-------------------------------------------------------------------

!      ------  Flow slides with vriable degree of sat height MP March 2023 

implicit none

integer(ink) ipoin 
real   (irk) Sr, alphas, dens, denss, densw, dens_bar , dens_bar_eff 
real   (irk) fact_dens       
!      ------  FSs values.

dens       = const   (2) 
denss      = const  (17)
densw      = const  (18)
Sr         = const  (19)
alphas     = const  (14)
ic_hrelSat = const  (34)

dens_bar     = dens
dens_bar_eff = dens_bar - alphas*densw
fact_dens   = dens_bar_eff/dens_bar

if (.NOT.allocated (hrelSat_DF) )  then
   allocate (hrelSat_DF  (npoin))
   allocate (fact_dens_DF(npoin)) 
   allocate (fact_pwp_max_DF (npoin) , fact_pwp_min_DF (npoin)) 
   allocate (fpwp_max_DF     (npoin) , fpwp_min_DF   (npoin))  
   fact_pwp_max_DF (:) = 0.0 
   fact_pwp_min_DF (:) = 0.0  
   fpwp_max_DF     (:) = 0.0 
   fpwp_min_DF     (:) = 0.0                                          
endif

hrelSat_DF  (:) = alphas
fact_dens_DF(:) = fact_dens

END Subroutine mini_FSs_initial_OLD

!-------------------------------------------------------------------

       Subroutine mini_DFs_dens_and_pwp_factors_OLD
       
!-------------------------------------------------------------------

!      ----   

implicit none

integer (ink) ipoin

real    (irk) alph, dens_ws

! real   (irk) dens, denss, densw, Sr, hrelpw, alphas
! real   (irk) porosity, dens1, dens2, dens_bar, dens_bar_eff
! real   (irk) fact_dens
! real   (irk) fact_pwp_max, fact_pwp_min, fpwp_max, fpwp_min

!      ------  Default values     

ic_fct_dens   = const (20) + 0.1
dens          = const  (2)
denss         = const (17)
densw         = const (18)
Sr            = const (19)
alphas        = const (14)
hrelpw        = alphas
alph          = alphas
ic_hrelSat    = const (34) + 0.01
fact_pwp_max  = 1.0
fact_pwp_min  = 0.0
fpwp_max      = 1.0
fpwp_min      = 0.0

if ( (ic_hrelSat.EQ.1).AND.(.NOT.allocated (fact_dens_DF)) ) then
      allocate  ( fact_dens_DF(npoin) )   
      allocate  ( dens1_DF(npoin), dens2_DF   (npoin), densw_DF       (npoin) ) 
      allocate  ( dens_DF (npoin), dens_bar_DF(npoin), dens_bar_eff_DF(npoin) )
      allocate  ( fact_pwp_max_DF (npoin), fact_pwp_min_DF (npoin) )
      allocate  ( fpwp_max_DF     (npoin), fpwp_min_DF     (npoin) )
      fact_dens_DF = 0.0; dens1_DF=0.0; dens2_DF=0.0; densw_DF=0.0
      dens_DF = 0.0; dens_bar_DF=0.0; dens_bar_eff_DF=0.0
      fact_pwp_max_DF= 0.0; fact_pwp_min_DF=0.0; fpwp_max_DF=0.0;fpwp_min_DF=0.0
endif

IF (ic_fct_dens.EQ.0) THEN                           !  1 phase, w or w/o pwp, correction in tanfi

    if ( densw.lt.1.e-3) densw = 1000.0
    dens_ws = densw/dens
    if ( alphas.LT.0.01) then                         !  means hrelpw ignored, so it is assumed to be 1
        alphas  = 1.0
        hrelpw  = 1.0
    endif
    fact_dens = 1.0                                   ! Default

ELSEIF (ic_fct_dens.EQ.1) then                                  !  1 phase, w or w/o pwp, correction OK
    
    if ( (denss.LT.1.e-3).OR.(densw.lt.1.e-3) ) then       !  Filters zero solid or water densities
       write (*,*) ' denss or densw wrong...', denss, densw
       write (*,*) ' strike a num key to acknowledge '
       read  (*,*) ipoin
       STOP
    endif
    
    if (alphaS.le.0.01) then                                  !  neglect saturated part
       fact_dens = 1.0                                        !  (dens-densw)/dens   
    else
       porosity     = (denss-dens)/(denss-densw)
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !     modified MP 10.01.2019 
       fact_dens    = dens_bar_eff/dens_bar                   !     multiplies tanFi
    endif

    if (icpwp.gt.0) then                                   !     If only affects pwp inf      
       fact_pwp_max = 1. + alphaS*densw/dens_bar_eff       !     ( 1-dpwp_bar/h)   
       fact_pwp_min = 0.
       fpwp_max     = 1.0
       fpwp_min     = -alphaS*densw/dens_bar_eff
    endif 
       
ELSEIF (ic_fct_dens.EQ.2) then                              !  2 phases DF, Hsat variable and pwp optionals

    IF (ic_hrelSat.eq.0) then                       !  all densities constants to compute factors
       porosity     = (denss-dens)/(denss-densw)              !
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !       
       fact_dens    = dens_bar_eff/(denss*(1-porosity))       !     multiplies tanFi
       if (icpwp.gt.0) then                                   !     If only affects pwp inf      
           fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
           fact_pwp_min = 0.
           fpwp_max     = 1.0
           fpwp_min     = -alphaS*densw/dens_bar_eff
       endif 

    ELSEIF (ic_hrelSat.eq.1) then                   !  all dens and factors depend on ipoin   

       DO ipoin = 1, npoin
           alph         = hrelSat_DF (ipoin); alphaS = alph       !  MP chgd bug aug 2021 alphas was wrong
           porosity     = n_DF (ipoin)           !
           dens2        = (1-porosity)*denss + porosity*densw*Sr  !  dens of upper layer (can be unsat)
           dens1        = (1-porosity)*denss + porosity*densw     !  ower layer (sat)   Sr=alphaS
           dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
           dens_bar_eff = dens_bar - densw*alphaS                 !       
           fact_dens    = dens_bar_eff/(denss*(1-porosity))       !     multiplies tanFi
           dens1_DF        (ipoin) = dens1
           dens2_DF        (ipoin) = dens2
           dens_bar_DF     (ipoin) = dens_bar
           dens_bar_eff_DF (ipoin) = dens_bar_eff
           fact_dens_DF    (ipoin) = fact_dens
           if (icpwp.gt.0) then                                   !     If only affects pwp inf      
               fact_pwp_max = 1. + alphaS*densw/dens_bar_eff      !     ( 1-dpwp_bar/h)   
               fact_pwp_min = 0.
               fpwp_max     = 1.0
               fpwp_min     = -alphaS*densw/dens_bar_eff
               fact_pwp_max_DF (ipoin) = fact_pwp_max   
               fact_pwp_min_DF (ipoin) = fact_pwp_min
               fpwp_max_DF     (ipoin) = fpwp_max
               fpwp_min_DF     (ipoin) = fpwp_min
           endif
       ENDDO 

    ENDIF

ELSE 
    write (*,*) ' ic_fct_dens given is wrong...', ic_fct_dens
    write (*,*) ' strike a num key to acknowledge '
    read  (*,*) ipoin
    STOP
ENDIF

END Subroutine mini_DFs_dens_and_pwp_factors_OLD


END SUBROUTINE   Get_Pwp_SW_OLD

!-------------------------------------------------------------------

       Subroutine OutputSave_SW  
       
!-------------------------------------------------------------------

implicit none

!      ------   

integer(ink)  i, id, len1
character(60) text                    ! general purpose char string 

99 format (i7,2x,i5,9(2x,g14.7))
98 format (i7,2x, 9(2x,g14.7))

len1 = len_trim(SPH_SW_problem_name)
open ( restart_file, file = SPH_SW_problem_name(1:len1)//'.restart.dat' )

write( restart_file, * )  ' ----- npoin_s      time at restart ----- ' , time  
write( restart_file, *)    npoin_s, time 
if (ndimn.eq.2) then 
   write( restart_file, * )  ' ipoin...itype....x ............   y ..........   vx............. vy ..........   h........   mass..........   hsml ..........   pwp   '
elseif (ndimn.eq.1) then 
   write( restart_file, * )  ' ipoin...itype....x ..........   vx..........   h........   mass..........   hsml or hrel ..........   pwp   '
endif
DO i = 1, npoin_s 
   if (ic_hrelSat.eq.1) then
       write(restart_file,99) i, itype(i), (x(id,i),id=1,ndimn), (vx(id,i),id=1,ndimn), &
                              rho (i), mass(i), hrelSat_DF(i), u (i) 
   else   
       write(restart_file,99) i, itype(i), (x(id,i),id=1,ndimn), (vx(id,i),id=1,ndimn), &
                              rho (i), mass(i), hsml(i), u (i) 
   endif 
ENDDO

IF (npoin_w.gt.0) THEN
    write( restart_file, * )  ' ----- npoin_w      time at restart ----- ' , time  
    write( restart_file, *)    npoin_w, time 
    if (ndimn.eq.2) then 
       write( restart_file, * )  ' ipoin...itype....x ............   y ..........   vx............. vy ..........   h........   mass..........   hsml ..........   pwp   '
    elseif (ndimn.eq.1) then 
       write( restart_file, * )  ' ipoin...itype....x ..........   vx..........   h........   mass..........   hsml ..........   pwp   '
    endif
    DO i = 1 + npoin_s, npoin_w + npoin_s  
       if (ic_hrelSat.eq.1) then
           write(restart_file,99) i, itype(i), (x(id,i),id=1,ndimn), (vx(id,i),id=1,ndimn), &
                                  rho (i), mass(i), hrelSat_DF(i), u (i) 
       else   
           write(restart_file,99) i, itype(i), (x(id,i),id=1,ndimn), (vx(id,i),id=1,ndimn), &
                              rho (i), mass(i), hsml(i), u (i) 
       endif     
    ENDDO
ENDIF 
close (restart_file)


IF (ic_erosion.GT.0) THEN   

   len1 = len_trim(SPH_SW_problem_name)
   open ( restart_top_file, file = SPH_SW_problem_name(1:len1)//'.restart.top' )
   write( restart_top_file, * )  ' ----- npoig      time at restart ----- ' , time  
   write( restart_top_file, *)    npoig, time 
   write( restart_top_file, * )  ' ipoig...Topol (1,:)....Topol (14,:) ..........Topol (16,:)  '
   DO i = 1, npoig 
      write(restart_top_file,98) i, Topol (1,i) , Topol (14,i) , Topol (15,i)  
   ENDDO
   close (restart_top_file)

ENDIF

End Subroutine OutputSave_SW


!-------------------------------------------------------------------

       Subroutine OutputMesh_SW  
       
!-------------------------------------------------------------------

    implicit none

!      ------  OUT for GID post.mesh  depends on IC_wir

if (ic_WIR.eq.0) then
    call OutputMesh_Classic
else if (ic_WIR.eq.1) then
    if (if_coarse_plt.eq.1) then
        call OutputMesh_coarse_SW
    else
       if (ic_TRIGGER.EQ.1) then
           call OutputMesh_WIR_TRIGGER
       else
           call OutputMesh_WIR
       endif
    endif
endif    

End Subroutine OutputMesh_SW



!-------------------------------------------------------------------

       Subroutine OutputMesh_Classic
       
!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  and init x0 --------

implicit none

integer  (ink) ipoin, ipoig, ipoigx, ipoigy, ieleg, n1, n2, n3, n4, i, idimn, inodg
integer  (ink) neleg
integer  (ink), allocatable:: intmag4(:,:)
integer  (ink) keleg, plteleg 

real     (irk) xx, yy, zp, xref, yref, xi, eta, sh1, sh2, sh3, sh4 

!      ------  generates background grid, plotting purposes only

neleg =  (npoigx-1)*(npoigy-1)
if(.NOT.allocated(intmag4) ) allocate ( intmag4(4,neleg) )
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

      
write (gid_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'

do ipoig = 1, npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)
   write(gid_msh,*) ipoig, xx, yy, zp
enddo

Call Check_Out_Domain_SW

do ipoin  = 1, npoin
   if ( if_Out_Domain(ipoin).eq.1) CYCLE             
   i      = ipoin + npoig
   xx     = x(1,ipoin)
   if (ndimn.eq.1) then
      yy  = (yming + ymaxg)/2.
   else
      yy  = x(2,ipoin)
   endif 
      
   call Get_Z_Topo ( ipoin,xx, yy, Zp )  
   Z00(ipoin) = Zp                         !      ------  Here we store the initial Z of particles. Used for plotting
   write(gid_msh,*) i, xx, yy, zp 
enddo
write (gid_msh,*) ' End coordinates'      
    
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg             ! **** chg
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gid_msh,*) ' End elements' 

write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '

write (gid_msh,*) ' Elements'
do i = 1, npoin                                         !       ------  ch 
   write(gid_msh,*) i, npoig+i, iabs(itype(i)) 
enddo
write (gid_msh,*) ' End elements'      
      
! write (gid_res,*) 'GiD Post Results File 1.0'   !  MP 2020 05 08 Goes after opening of file
      
close(gid_msh)
deallocate (intmag4)

End Subroutine OutputMesh_Classic


!-------------------------------------------------------------------

       Subroutine OutputMesh_WIR
       
!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  and init x0 --------
!              added Q3 mesh for 1D full pwp

implicit none

integer  (ink) ipoin, i,is_end, idimn, inodg
integer  (ink) ipoig, ipoigx, ipoigy, ieleg, neleg
integer  (ink), allocatable:: intmag4(:,:)
integer  (ink) keleg, plteleg 
integer  (ink) iaux, is_water

real     (irk) xx, yy, zz,  zp
real     (irk) hh, dhh          ! for 1d pwp meshes h and increm of h

integer  (ink), allocatable:: intmat_pwp(:,:)   ! 1dim, full pwp mesh
integer  (ink)  nel_pwp                         ! nr of elems in pwp mesh

!      ------  generates background grid, plotting purposes only

neleg =  (npoigx-1)*(npoigy-1)
if(.NOT.allocated(intmag4) ) allocate ( intmag4(4,neleg) )
intmag4 = 0
keleg   = 0   ! chg
do ipoigy = 1,npoigy-1
   do ipoigx  = 1, npoigx-1
      ieleg   = (ipoigy-1)*(npoigx-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoigx + ipoigx
      plteleg = poig_in_model(ipoig)         * poig_in_model(ipoig+1)       &
               *poig_in_model(ipoig+npoigx+1)* poig_in_model(ipoig+npoigx)   
      if (plteleg.eq.1) then!  
          keleg = keleg+1
          intmag4(1,keleg) = ipoig   
          intmag4(2,keleg) = ipoig  + 1
          intmag4(3,keleg) = ipoig  + npoigx + 1
          intmag4(4,keleg) = ipoig  + npoigx 
      endif
   enddo
enddo

write (gid_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)       !  + topol(14,ipoig)     !     Just for restart of topo TAke care Zp (ipoin) will be WRONG
   write(gid_msh,*) ipoig, xx, yy, zp           !     MP 06 06 2019 for topo.restart
enddo

do ipoig = 1, npoig                             !      ------  Q2 This is 2nd set of g nodes npoig+1:2*npoig
   i     = npoig + ipoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)                       !     + topol(14,ipoig)
   write(gid_msh,*) i, xx, yy, zp
enddo

Call Check_Out_Domain_SW

i = 2*npoig
do ipoin  = 1, npoin                            !      ------  3rd set: coordinates of SPH nodes
   i      = i + 1                               !   
   if ( if_Out_Domain(ipoin).eq.0) then            
      xx     = x(1,ipoin)
      if (ndimn.eq.1) then
         yy  = (yming + ymaxg)/2.
      else
         yy  = x(2,ipoin)
      endif 
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      Z00(ipoin) = Zp                              !      ------  Here we store the initial Z of particles. Used for plotting
      write(gid_msh,*) i, xx, yy, zp 
   elseif ( if_Out_Domain(ipoin).eq.1) then         ! Modified 15 April 2016 for restart + windows          
     ! xx     = xming                                ! outside pics go to corner :-)
     ! yy     = yming
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      Z00(ipoin) = Zp                              !      ------  Here we store the initial Z of particles. Used for plotting
      write(gid_msh,*) i, xx, yy, zp 
   endif
enddo

i        = 2*npoig + npoin
   
if (ndimn.eq.1.and.nUaux.gt.1.and.iabs(Gid_Mask_SW(12)).eq.1) THEN  !   ------  Q3 This is 4rd set of nodes
                                                                    !           ** MP iabs added Jan 20 2019
   nel_pwp  = (npoin_s-1)*(nUaux-1)                                 !           ** MP npoin-> npoin_s 20th Jan 19
   if(.NOT.allocated(intmat_pwp) ) allocate (intmat_pwp(4,nel_pwp) )   ! bakground mesh
   intmat_pwp = 0
   ieleg      = 0
   do ipoin = 1, npoin_s-1                                          !           ** MP npoin-> npoin_s 20th Jan 19
      do iaux  = 1, nUaux-1
         ieleg = ieleg + 1
         intmat_pwp (1,ieleg) = ieleg + (ipoin-1) 
         intmat_pwp (2,ieleg) = intmat_pwp (1,ieleg) + nUaux
         intmat_pwp (3,ieleg) = intmat_pwp (2,ieleg) + 1
         intmat_pwp (4,ieleg) = intmat_pwp (1,ieleg) + 1
      enddo
   enddo
         
   do ipoin = 1, npoin_s                             !             ** MP npoin-> npoin_s 20th Jan 19 
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      xx    = x  (1,ipoin)
      yy    = (yming + ymaxg)/2.
      call Get_Z_Topo ( ipoin, xx, yy, Zp)     
      if (allocated(h_df) ) then                     !             ** MP npoin-> npoin_s 20th Jan 19 
          hh = h_df(ipoin)
      else
           hh    = rho(ipoin)
      endif
      dhh   = hh/(nUaux-1)
      do iaux  = 1, nUaux 
         i     = i + 1
         zz    = Zp + (iaux-1)*dhh
         write(gid_msh,*) i, xx, yy, zz
      enddo
   enddo
   
endif   

write (gid_msh,*) ' End coordinates'
      
                                                !      ========= ELEMS  =========================== 
                                                
write (gid_msh,*) 'MESH SPH_Q1   dimension   3    ElemType Quadrilateral  Nnode 4 '  !------- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg               
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gid_msh,*) ' End elements' 

write (gid_msh,*) 'MESH SPH_Q2   dimension   3    ElemType Quadrilateral  Nnode 4 '  !----- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg)+npoig, inodg = 1,4), '    1 '
enddo
write (gid_msh,*) ' End elements'

is_water = 0
write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  SOIL ----- 
write (gid_msh,*) ' Elements'
do i = 1, npoin                                          
   is_end = i
   if (iabs(itype(i)).eq.2) is_water = 1
   if (iabs(itype(i)).eq.2) EXIT
   write(gid_msh,*) i, 2*npoig+i, iabs(itype(i))
enddo
write (gid_msh,*) ' End elements' 

if (is_water.eq.1) then
   write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  WATER  --- 
   write (gid_msh,*) ' Elements'
   do i = is_end, npoin                                             
      write(gid_msh,*) i-is_end+1, 2*npoig+i, iabs(itype(i))
   enddo
   write (gid_msh,*) ' End elements' 
endif

if (ndimn.eq.1.and.nUaux.gt.1.and.Gid_Mask_SW(12).eq.1) THEN  ! 
   i = 2*npoig + npoin
   write (gid_msh,*) 'MESH SPH_Q_pwp   dimension   3    ElemType Quadrilateral  Nnode 4 '  !----- 
   write (gid_msh,*) ' Elements'
   do ieleg = 1, nel_pwp
      write (gid_msh,*) ieleg, (intmat_pwp(inodg,ieleg)+ i, inodg = 1,4), '    1 '
   enddo
   write (gid_msh,*) ' End elements'
endif   
      
! write (gid_res,*) 'GiD Post Results File 1.0'   !  MP 2020 05 08 Goes after opening of file
      
close(gid_msh)
deallocate (intmag4)
if (allocated (intmat_pwp) ) then
    deallocate (intmat_pwp)
endif

End Subroutine OutputMesh_WIR

!-------------------------------------------------------------------

       Subroutine OutputMesh_WIR_TRIGGER
       
!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  and init x0 --------
!              added Q3 mesh for 1D full pwp

implicit none

integer  (ink) ipoin, i,is_end, idimn, inodg
integer  (ink) ipoig, ipoigx, ipoigy, ieleg, neleg
integer  (ink), allocatable:: intmag4(:,:)
integer  (ink) keleg, plteleg 
integer  (ink) iaux, is_water

real     (irk) xx, yy, zz,  zp
real     (irk) hh, dhh          ! for 1d pwp meshes h and increm of h

integer  (ink), allocatable:: intmat_pwp(:,:)   ! 1dim, full pwp mesh
integer  (ink)  nel_pwp                         ! nr of elems in pwp mesh

!      ------  generates background grid, plotting purposes only

neleg =  (npoigx-1)*(npoigy-1)
if(.NOT.allocated(intmag4) ) allocate ( intmag4(4,neleg) )
intmag4 = 0
keleg   = 0   ! chg
do ipoigy = 1,npoigy-1
   do ipoigx  = 1, npoigx-1
      ieleg   = (ipoigy-1)*(npoigx-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoigx + ipoigx
      plteleg = poig_in_model(ipoig)         * poig_in_model(ipoig+1)       &
               *poig_in_model(ipoig+npoigx+1)* poig_in_model(ipoig+npoigx)   
      if (plteleg.eq.1) then!  
          keleg = keleg+1
          intmag4(1,keleg) = ipoig   
          intmag4(2,keleg) = ipoig  + 1
          intmag4(3,keleg) = ipoig  + npoigx + 1
          intmag4(4,keleg) = ipoig  + npoigx 
      endif
   enddo
enddo

write (gid_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)       !  + topol(14,ipoig)     !     Just for restart of topo TAke care Zp (ipoin) will be WRONG
   write(gid_msh,*) ipoig, xx, yy, zp           !     MP 06 06 2019 for topo.restart
enddo

Call Check_Out_Domain_SW

write (gid_msh,*) ' End coordinates'
      
                                                !      ========= ELEMS  =========================== 
                                                
write (gid_msh,*) 'MESH SPH_Q1   dimension   3    ElemType Quadrilateral  Nnode 4 '  !------- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg               
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gid_msh,*) ' End elements' 
      
close(gid_msh)
deallocate (intmag4)
if (allocated (intmat_pwp) ) then
    deallocate (intmat_pwp)
endif

End Subroutine OutputMesh_WIR_TRIGGER

!-------------------------------------------------------------------

       Subroutine OutputMesh_WIR_TRIGGER_BAK
       
!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  and init x0 -------- PROPOSED CHANGES TO
!              added Q3 mesh for 1D full pwp               AVOID PLOTTING
!                                                          big or small slopes
implicit none

integer  (ink) ipoin, i,is_end, idimn, inodg
integer  (ink) ipoig, ipoigx, ipoigy, ieleg, neleg
integer  (ink), allocatable:: intmag4(:,:)
integer  (ink) keleg, plteleg 
integer  (ink) iaux, is_water

real     (irk) xx, yy, zz,  zp, Zx, Zy, Zgrad
real     (irk) hh, dhh          ! for 1d pwp meshes h and increm of h

integer  (ink), allocatable:: intmat_pwp(:,:)   ! 1dim, full pwp mesh
integer  (ink)  nel_pwp                         ! nr of elems in pwp mesh

write (gid_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)       !  + topol(14,ipoig)     !     Just for restart of topo TAke care Zp (ipoin) will be WRONG
   Zx = topol (2,ipoig)
   Zy = topol (3,ipoig)
   Zgrad = (Zx * Zx + Zy*Zy)
   if (Zgrad.GT.0.001) then
      Zgrad = Zgrad ** 0.5
   else
      Zgrad = 0.0
   endif
   Slope_Max_cut = const(72); Slope_Min_cut = const(73) 
   if ( Zgrad.GT.Slope_Max_cut.OR. Zgrad.LT.Slope_Min_cut) THEN  ! avoid plotting nodes with slopes
      poig_in_model(ipoig) = 0                                  !  outside range
   endif                                                        ! MP CHGD 02 April 2022
   write(gid_msh,*) ipoig, xx, yy, zp           !     MP 06 06 2019 for topo.restart
enddo

Call Check_Out_Domain_SW

write (gid_msh,*) ' End coordinates'

!      ------  generates background grid, plotting purposes only

neleg =  (npoigx-1)*(npoigy-1)
if(.NOT.allocated(intmag4) ) allocate ( intmag4(4,neleg) )
intmag4 = 0
keleg   = 0   ! chg
do ipoigy = 1,npoigy-1
   do ipoigx  = 1, npoigx-1
      ieleg   = (ipoigy-1)*(npoigx-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoigx + ipoigx
      plteleg = poig_in_model(ipoig)         * poig_in_model(ipoig+1)       &
               *poig_in_model(ipoig+npoigx+1)* poig_in_model(ipoig+npoigx)   
      if (plteleg.eq.1) then!  
          keleg = keleg+1
          intmag4(1,keleg) = ipoig   
          intmag4(2,keleg) = ipoig  + 1
          intmag4(3,keleg) = ipoig  + npoigx + 1
          intmag4(4,keleg) = ipoig  + npoigx 
      endif
   enddo
enddo      
                                                !      ========= ELEMS  =========================== 
                                                
write (gid_msh,*) 'MESH SPH_Q1   dimension   3    ElemType Quadrilateral  Nnode 4 '  !------- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg               
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gid_msh,*) ' End elements' 
      
close(gid_msh)
deallocate (intmag4)
if (allocated (intmat_pwp) ) then
    deallocate (intmat_pwp)
endif

End Subroutine OutputMesh_WIR_TRIGGER_BAK


!-------------------------------------------------------------------

       Subroutine OutputMesh_WIR_SAEID_bug
       
!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  and init x0 --------
!              added Q3 mesh for 1D full pwp

implicit none

integer  (ink) ipoin, i,is_end, idimn, inodg
integer  (ink) ipoig, ipoigx, ipoigy, ieleg, neleg
integer  (ink), allocatable:: intmag4(:,:)
integer  (ink) keleg, plteleg 
integer  (ink) iaux, is_water

real     (irk) xx, yy, zz,  zp
real     (irk) hh, dhh          ! for 1d pwp meshes h and increm of h

integer  (ink), allocatable:: intmat_pwp(:,:)   ! 1dim, full pwp mesh
integer  (ink)  nel_pwp                         ! nr of elems in pwp mesh

!      ------  generates background grid, plotting purposes only

neleg =  (npoigx-1)*(npoigy-1)
if(.NOT.allocated(intmag4) ) allocate ( intmag4(4,neleg) )
intmag4 = 0
keleg   = 0   ! chg
do ipoigy = 1,npoigy-1
   do ipoigx  = 1, npoigx-1
      ieleg   = (ipoigy-1)*(npoigx-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoigx + ipoigx
      plteleg = poig_in_model(ipoig)         * poig_in_model(ipoig+1)       &
               *poig_in_model(ipoig+npoigx+1)* poig_in_model(ipoig+npoigx)   
      if (plteleg.eq.1) then!  
          keleg = keleg+1
          intmag4(1,keleg) = ipoig   
          intmag4(2,keleg) = ipoig  + 1
          intmag4(3,keleg) = ipoig  + npoigx + 1
          intmag4(4,keleg) = ipoig  + npoigx 
      endif
   enddo
enddo

write (gid_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'                !      ------ COORDS  --------------------

do ipoig = 1, npoig                             !      ------  Q1 This is first set of g nodes 1:npoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)       !  + topol(14,ipoig)     !     Just for restart of topo TAke care Zp (ipoin) will be WRONG
   write(gid_msh,*) ipoig, xx, yy, zp           !     MP 06 06 2019 for topo.restart
enddo

do ipoig = 1, npoig                             !      ------  Q2 This is 2nd set of g nodes npoig+1:2*npoig
   i     = npoig + ipoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)                       !     + topol(14,ipoig)
   write(gid_msh,*) i, xx, yy, zp
enddo

do ipoig = 1, npoig                             !   Saeidmt 26May2021      ------  Q3 This is 3rd set of g nodes npoig+1:3*npoig
   i     = 2*npoig + ipoig
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)                       !     + topol(14,ipoig)
   write(gid_msh,*) i, xx, yy, zp
enddo

Call Check_Out_Domain_SW

i = 3*npoig                                     !   Saeidmt 26May2021    **old i = 2*npoig
do ipoin  = 1, npoin                            !      ------  3rd set: coordinates of SPH nodes
   i      = i + 1                               !   
   if ( if_Out_Domain(ipoin).eq.0) then            
      xx     = x(1,ipoin)
      if (ndimn.eq.1) then
         yy  = (yming + ymaxg)/2.
      else
         yy  = x(2,ipoin)
      endif 
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      Z00(ipoin) = Zp                              !      ------  Here we store the initial Z of particles. Used for plotting
      write(gid_msh,*) i, xx, yy, zp 
   elseif ( if_Out_Domain(ipoin).eq.1) then         ! Modified 15 April 2016 for restart + windows          
     ! xx     = xming                                ! outside pics go to corner :-)
     ! yy     = yming
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      Z00(ipoin) = Zp                              !      ------  Here we store the initial Z of particles. Used for plotting
      write(gid_msh,*) i, xx, yy, zp 
   endif
enddo

i        = 3*npoig + npoin                          !   Saeidmt 26May2021    **old    i        = 2*npoig + npoin
   
if (ndimn.eq.1.and.nUaux.gt.1.and.iabs(Gid_Mask_SW(12)).eq.1) THEN  !   ------  Q3 This is 4rd set of nodes
                                                                    !           ** MP iabs added Jan 20 2019
   nel_pwp  = (npoin_s-1)*(nUaux-1)                                 !           ** MP npoin-> npoin_s 20th Jan 19
   if(.NOT.allocated(intmat_pwp) ) allocate (intmat_pwp(4,nel_pwp) )   ! bakground mesh
   intmat_pwp = 0
   ieleg      = 0
   do ipoin = 1, npoin_s-1                                          !           ** MP npoin-> npoin_s 20th Jan 19
      do iaux  = 1, nUaux-1
         ieleg = ieleg + 1
         intmat_pwp (1,ieleg) = ieleg + (ipoin-1) 
         intmat_pwp (2,ieleg) = intmat_pwp (1,ieleg) + nUaux
         intmat_pwp (3,ieleg) = intmat_pwp (2,ieleg) + 1
         intmat_pwp (4,ieleg) = intmat_pwp (1,ieleg) + 1
      enddo
   enddo
         
   do ipoin = 1, npoin_s                             !             ** MP npoin-> npoin_s 20th Jan 19 
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      xx    = x  (1,ipoin)
      yy    = (yming + ymaxg)/2.
      call Get_Z_Topo ( ipoin, xx, yy, Zp)     
      if (allocated(h_df) ) then                     !             ** MP npoin-> npoin_s 20th Jan 19 
          hh = h_df(ipoin)
      else
           hh    = rho(ipoin)
      endif
      dhh   = hh/(nUaux-1)
      do iaux  = 1, nUaux 
         i     = i + 1
         zz    = Zp + (iaux-1)*dhh
         write(gid_msh,*) i, xx, yy, zz
      enddo
   enddo
   
endif   

write (gid_msh,*) ' End coordinates'
      
                                                !      ========= ELEMS  =========================== 
                                                
write (gid_msh,*) 'MESH SPH_Q1   dimension   3    ElemType Quadrilateral  Nnode 4 '  !------- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg               
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gid_msh,*) ' End elements' 

write (gid_msh,*) 'MESH SPH_Q2   dimension   3    ElemType Quadrilateral  Nnode 4 '  !----- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg)+npoig, inodg = 1,4), '    1 '
enddo
write (gid_msh,*) ' End elements'

write (gid_msh,*) 'MESH SPH_Q3   dimension   3    ElemType Quadrilateral  Nnode 4 '  !   Saeidmt 26May2021  **Add 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg)+2*npoig, inodg = 1,4), '    1 '
enddo
write (gid_msh,*) ' End elements'

is_water = 0
write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  SOIL ----- 
write (gid_msh,*) ' Elements'
do i = 1, npoin                                          
   is_end = i
   if (iabs(itype(i)).eq.2) is_water = 1
   if (iabs(itype(i)).eq.2) EXIT
   write(gid_msh,*) i, 3*npoig+i, iabs(itype(i))                           !   Saeidmt 26May2021  **old      write(gid_msh,*) i, 2*npoig+i, iabs(itype(i))     
enddo
write (gid_msh,*) ' End elements'  

if (is_water.eq.1) then
   write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  WATER  --- 
   write (gid_msh,*) ' Elements'
   do i = is_end, npoin                                             
      write(gid_msh,*) i-is_end+1, 3*npoig+i, iabs(itype(i))             !   Saeidmt 26May2021  **old      write(gid_msh,*) i-is_end+1, 2*npoig+i, iabs(itype(i))
   enddo
   write (gid_msh,*) ' End elements' 
endif

if (ndimn.eq.1.and.nUaux.gt.1.and.Gid_Mask_SW(12).eq.1) THEN  ! 
   i = 3*npoig + npoin                                                                                !   Saeidmt 26May2021  **old    i = 2*npoig + npoin  
   write (gid_msh,*) 'MESH SPH_Q_pwp   dimension   3    ElemType Quadrilateral  Nnode 4 '  !----- 
   write (gid_msh,*) ' Elements'
   do ieleg = 1, nel_pwp
      write (gid_msh,*) ieleg, (intmat_pwp(inodg,ieleg)+ i, inodg = 1,4), '    1 '
   enddo
   write (gid_msh,*) ' End elements'
endif   
      
! write (gid_res,*) 'GiD Post Results File 1.0'   !  MP 2020 05 08 Goes after opening of file
      
close(gid_msh)
deallocate (intmag4)
if (allocated (intmat_pwp) ) then
    deallocate (intmat_pwp)
endif

End Subroutine OutputMesh_WIR_SAEID_bug



!-------------------------------------------------------------------

       Subroutine OutputMesh_coarse_SW
       
!-------------------------------------------------------------------

!      ------  OUT for GID post.mesh  and init x0 --------

implicit none

integer  (ink) ipoin, i,is_end, idimn, inodg
integer  (ink) ipoig, ipoigx, ipoigy, ieleg, neleg
integer  (ink), allocatable:: intmag4(:,:)
integer  (ink) keleg, plteleg 

real     (irk) xx, yy, zp

!      ------  generates background grid, plotting purposes only

neleg =  (npoix_plt-1)*(npoiy_plt-1)
if(.NOT.allocated(intmag4) ) allocate ( intmag4(4,neleg) )
intmag4 = 0
keleg   = 0   ! chg
do ipoigy = 1,npoiy_plt-1
   do ipoigx  = 1, npoix_plt-1
      ieleg   = (ipoigy-1)*(npoix_plt-1) + ipoigx 
      ipoig   = (ipoigy-1)*npoix_plt + ipoigx
      plteleg = poi_in_model_plt(ipoig)*poi_in_model_plt(ipoig+1)*poi_in_model_plt(ipoig+npoix_plt+1)*poi_in_model_plt(ipoig+npoix_plt)  ! **** chg
      if (plteleg.eq.1) then! **** chg
          keleg = keleg+1
          intmag4(1,keleg) = ipoig   
          intmag4(2,keleg) = ipoig  + 1
          intmag4(3,keleg) = ipoig  + npoix_plt + 1
          intmag4(4,keleg) = ipoig  + npoix_plt 
      endif
   enddo
enddo

      
write (gid_msh,*) 'MESH SPH_quads   dimension   3    ElemType Quadrilateral  Nnode 4 '

write (gid_msh,*) ' coordinates'                !      ======== COORDS  ================

do ipoig = 1, npoi_plt                          !      ------  Q1 This is first set of g nodes 1:npoi_plt
   xx    = coor_plt(1,ipoig)
   yy    = coor_plt(2,ipoig)
   zp    = coor_plt(3,ipoig)
   write(gid_msh,*) ipoig, xx, yy, zp
enddo

do ipoig = 1, npoi_plt                          !      ------  Q2 This is 2nd set of g nodes npoig+1:2*npoig
   i     = npoi_plt + ipoig
   xx    = coor_plt(1,ipoig)
   yy    = coor_plt(2,ipoig)
   zp    = coor_plt(3,ipoig)
   write(gid_msh,*) i, xx, yy, zp
enddo

Call Check_Out_Domain_SW

do ipoin  = 1, npoin                            !      ------  3rd set: coordinates of SPH nodes
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   i      = ipoin + 2*npoi_plt
   xx     = x(1,ipoin)
   if (ndimn.eq.1) then
      yy  = (yming + ymaxg)/2.
   else
      yy  = x(2,ipoin)
   endif 
   call Get_Z_Topo ( ipoin, xx, yy, Zp)

   Z00(ipoin) = Zp                              !      ------  Here we store the initial Z of particles. Used for plotting
   write(gid_msh,*) i, xx, yy, zp 
enddo
write (gid_msh,*) ' End coordinates'      
                                                !      ======== ELEMS  ================ 
                                                
write (gid_msh,*) 'MESH SPH_Q1   dimension   3    ElemType Quadrilateral  Nnode 4 '  !----- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg               
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg), inodg=1,4), '    1 '
enddo
write (gid_msh,*) ' End elements' 


      
write (gid_msh,*) 'MESH SPH_Q2   dimension   3    ElemType Quadrilateral  Nnode 4 '  !------- 
write (gid_msh,*) ' Elements'
do ieleg = 1, keleg
   write (gid_msh,*) ieleg, (intmag4(inodg,ieleg)+npoi_plt, inodg = 1,4), '    1 '
enddo
write (gid_msh,*) ' End elements'



write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  SOIL ----- 
write (gid_msh,*) ' Elements'
do i = 1, npoin                                          
   is_end = i
   if (iabs(itype(i)).eq.2) EXIT
   write(gid_msh,*) i, 2*npoi_plt+i, iabs(itype(i))
enddo
write (gid_msh,*) ' End elements' 



write (gid_msh,*) 'MESH SPH_points dimension 3   ElemType Point  Nnode 1 '  ! ----  WATER  --------------------------------- 
write (gid_msh,*) ' Elements'
do i = is_end, npoin                                             
   write(gid_msh,*) i-is_end+1, 2*npoi_plt+i, iabs(itype(i))
enddo
write (gid_msh,*) ' End elements' 
     
! write (gid_res,*) 'GiD Post Results File 1.0'   !  MP 2020 05 08 Goes after opening of file
      
close(gid_msh)
deallocate (intmag4)

End Subroutine  OutputMesh_coarse_SW



!-------------------------------------------------------------------

       Subroutine OutputRes_SW  
       
!-------------------------------------------------------------------


    implicit none
    real    (irk) time_tmp
     
!      ------  OUT for GID post.res  depends on IC_wir


time_tmp = time                 !    MP 2020 05 08 for cases, print nr of problem instead of time
if (ic_cases_mat.GT.0.AND.nproblems.GT.1) then
   time     = nproblems_index
endif

if (ic_WIR.eq.0) then
    call OutputRes_Classic
else if (ic_WIR.eq.1) then
    if (if_coarse_plt.eq.1) then
        call OutputRes_coarse_SW
    else
        call OutputRes_WIR
    endif
endif    

if (ic_cases_mat.GT.0.AND.nproblems.GT.1) then
   time     = time_tmp
endif

End Subroutine OutputRes_SW


!-------------------------------------------------------------------

       Subroutine OutputRes_Classic  
       
!-------------------------------------------------------------------


!      ------  OUT for GID post.res  and init x0 --------

implicit none

integer  (ink)  i, ipoin, ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer  (ink)  inx0, inx1, inx, iny0, iny1, iny
real     (irk)  xx, yy, Zp 
real     (irk)  xdisp, ydisp, zdisp, vxx, vyy
real     (irk), allocatable::disp(:)
real     (irk), allocatable::auxv(:,:)
real     (irk)  distx2, disty2, dist, distg 

if(.NOT.allocated(disp) ) allocate ( disp(ndimn) ) 
if(.NOT.allocated(auxv) ) allocate ( auxv(5,npoig))

!      ------   First of all we get values at topo grid. Add for each cell value*(1/dist) and divide by sum (1/dist)

auxv = 0.0
distg = (deltxg*deltxg + deltyg*deltyg)**0.5

DO ipoin = 1,npoin                                      !      ------  ch
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
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
      auxv(1,ipoig) = auxv(1,ipoig) + rho(ipoin)/dist 
      auxv(5,ipoig) = auxv(5,ipoig) + 1.00/dist 
      auxv(4,ipoig) = auxv(4,ipoig) + (u(ipoin))/dist   ! pwp is now normalized /rho(ipoin)
      auxv(2,ipoig) = auxv(2,ipoig) + vx(1,ipoin)/dist 
      if(ndimn.eq.2) auxv(3,ipoig) = auxv(3,ipoig) + vx(2,ipoin)/dist
   enddo
   enddo
ENDDO

DO ipoig = 1, npoig
   if (auxv(5,ipoig).ge.1.e-6) then
       auxv(1:4,ipoig) = auxv(1:4,ipoig)/auxv(5,ipoig) 
   endif
enddo

!      ------   And now plot:
!               First, Position of avalanche


write(gid_res,*) 'Result "height" "Height" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoig                                    !      ------  plots h at quad nodes
   if ( auxv(1,ipoig).le.1.e-6) CYCLE
   write(gid_res,*) ipoig, ' 0.   0.   ', auxv(1,ipoig)
enddo
write(gid_res,*) ' End Values ' 


write(gid_res,*) 'Result "dis" "Disp" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values  '  

Call Check_Out_Domain_SW

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
      write(gid_res,*) i + npoig, xdisp, ydisp, zdisp
ENDDO
write(gid_res,*) ' End Values ' 

!      ------  Velocities  -----------------------------

write(gid_res,*) 'Result "vel" "veloc" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoig
   if ( auxv(1,ipoig).le.1.e-6) CYCLE
   write(gid_res,*) ipoig, auxv(2,ipoig), auxv(3,ipoig), ' 0.0 '
enddo
do i = 1,npoin 
   if ( if_Out_Domain(i).eq.1) CYCLE
   vxx = vx(1,i)
   if (ndimn.eq.1) then
       vyy = 0.0
   else
       vyy = vx(2,i)
   endif
   write(gid_res,*) i + npoig, vxx, vyy, ' 0.0 ' 
enddo
write(gid_res,*) ' End Values '

!      ------  And PWPs  -----------------------------

if (icpwp.eq.1) then
   write(gid_res,*) 'Result "Pwp" "pwpress" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig
      if ( auxv(1,ipoig).le.1.e-6) CYCLE
      write(gid_res,*) ipoig, '  0.0   0.0   ', auxv(4,ipoig)
   enddo
!   do i = 1,npoin 
!      if ( if_Out_Domain(i).eq.1) CYCLE
!      write(gid_res,*) i + npoig, '  0.0   0.0   ', u(i)
!   enddo
   write(gid_res,*) ' End Values '
endif   

!      ------  Slopes  -----------------------------
!
!write(gid_res,*) 'Result "slp" "slope" ',time,' Vector OnNodes "where" '
!write(gid_res,*) ' Values '
!
!do ipoig = 1,npoig
!   if (poig_in_model(ipoig).ne.0) then
!      vxx = topol(2,ipoig)
!      vyy = topol(3,ipoig)
!   endif
!   write(gid_res,*) ipoig, vxx, vyy, ' 0.0 ' 
!enddo
!
!write(gid_res,*) ' End Values '
!


deallocate ( disp )
deallocate ( auxv )

 
End Subroutine OutputRes_Classic 


!-------------------------------------------------------------------

       Subroutine OutputRes_WIR_BAK23 
       
!-------------------------------------------------------------------


!      ------  OUT for waves in reservoirs --------



!      Filter for GID
!
!   We will use a mask to plot variables, nGidMaskSw = 12
!   Gid_Mask_SW (nGidMask) can be zero (no plot) or one (plot)
!   SW case:
!          
!       1 ... h soil       
!       2 ... displacement          
!       3 ... velocity          
!       4 ... Pwp          
!       5 ... erosion          
!       6 ... Z         
!       7 ... h soil relative  (used in Lausanne)         
!       8 ... h water          
!       9 ... eta  (when (ic_ws_Interact.eq.2) is porosity n        
!      10 ... hs + hw
!      11 ... hsml OR hrelSat_DF
!      12 .   indicates full pwp output in 1d problems

implicit none

integer  (ink)  i, ipoin, ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer  (ink)  inx0, inx1, inx, iny0, iny1, iny, npoing
integer  (ink)  idest, iaux
real     (irk)  xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp, dhh
integer  (ink)  dinx, diny     ! MP 01 may 2019  
real     (irk)  xdisp, ydisp, zdisp, vxx, vyy, vvv, hhh, hh0, cgra, Froude
real     (irk), allocatable::disp(:)
real     (irk), allocatable::auxv(:,:)
real     (irk), allocatable, SAVE:: etag   (:)       
integer  (ink), allocatable, SAVE:: ictrack(:)      
real     (irk)  distx2, disty2, dist, distg 
real     (irk)  dArea, Mtotal_out
real     (irk)  xcg, ycg, zcg       ! centre of gravity avalanche MP 1st dec 2017
real(irk), allocatable:: height_g(:)

real    (irk) dens, denss, densw        !   ic_fct_dens = C20  0 no correction  1.. correct 
real    (irk) porosity, Sr              !   Porosity and degree of saturation
real    (irk) dens1, dens2              !   mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    !   average and effective average dens of 2 layers
real    (irk) alphaS                    !   ratio hsat/h, from hrelpw. If below 0.1 -> 0

real    (irk) deng, gravz
real    (irk) height, pwp_total, Zpx, Zpy

real   (irk) hh,tt,xt
real   (irk) vv, ee, Te, Tc, Pe, Cvv

integer (ink) ic_fct_dens, ic_cosTheta    

if(.NOT.allocated(disp) ) allocate ( disp(ndimn) ) 
if(.NOT.allocated(auxv) ) allocate ( auxv(14,npoig))

cgra        = const  (1)                
dens        = const  (2)            
denss       = const (17)             
densw       = const (18)
Sr          = const (19)
alphaS      = const (14)
ic_fct_dens = const (20)
ic_cosTheta = const (16)
if (( Gid_Mask_SW (4).eq.-1 ).and. ic_fct_dens.eq.0) then
    write (*,*) ' take care C20 ic_fct_dens should be 1 or 2 '
    write (*,*) ' hit a numerical key ' 
    read  (*,*) i
    stop
endif

Call Check_Out_Domain_SW

!      ------   First of all we get values at topo grid. 
!               Add for each cell value*(1/dist) and divide by sum (1/dist)

auxv    = 0.0

if(.not.allocated(ictrack)) then                   
   allocate (ictrack(npoig))
   ictrack = 0
endif
                                                   
distg   = (deltxg*deltxg + deltyg*deltyg)**0.5

DO ipoin = 1,npoin    
   if (allocated (if_correct)) then             ! CHGD MP  april 2020. if_correct is not always allocated
      if (if_correct(ipoin).eq.1) CYCLE         ! CHGD MP 17th March 2020  skip node to avoid high porosities
   endif
   pwp_total = 0.0                                                  ! ** mp 20th Jan 2019
   if ( (Gid_Mask_SW (4).eq.-1) .AND. (ipoin.LE.npoin_s) ) then     ! ** mp 20th Jan 2019 avoid w points
       xx  = x(1,ipoin)            !  fixing bug in Hungr: before the rk if 
       yy  = (yming+ymaxg)/2     
       if (ndimn.eq.2) yy = x(2,ipoin)                                                      
       if ( ic_cosTheta.eq.1) then                          !  if we use the theta formulation   
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)     !  we will define them again              
          deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
          gravz =  cgra/deng                              ! *** think of joining with nfrict 7 
       elseif ( ic_cosTheta.eq.0) then 
!          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)                 
!          deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
          gravz =  cgra 
       elseif ( ic_cosTheta.eq.2) then                    !  do nothing, we have gx gy gz 
          deng  = deng
       endif
          
       if (ic_ws_Interact.eq.2) then
          height   = h_DF(ipoin)
          porosity = n_DF(ipoin)
       else
          height   = rho(ipoin)  
          porosity = (denss-dens)/(denss-densw)
       endif 
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !     MP 19 January 2019 
     
       pwp_total = gravz * (densw * height * alphaS + dens_bar_eff * u(ipoin) )
       pwp_total = pwp_total/1000.   ! KPa
     
   endif
                        
   if ( if_Out_Domain(ipoin).eq.1 )  CYCLE  ! ** mp 20th Jan 2019 avois water pts
   idest = 5*(iabs(itype(ipoin))-1)                                 !  = 0 for soil and 5 for water
   xx = x(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
!   inx0      = max(1,ipoigx-1)
!   iny0      = max(1,ipoigy-1)
!   inx1      = min(npoigx,ipoigx+2)
!   iny1      = min(npoigy,ipoigy+2)
   dinx      = hsml(ipoin)/deltxg 
   diny      = hsml(ipoin)/deltyg 
   if (ndimn.eq.1) diny = 0 
   inx0      = max(1,ipoigx-dinx)
   iny0      = max(1,ipoigy-diny)
   inx1      = min(npoigx,ipoigx + dinx+1)  !  min(npoigx,ipoigx+1) ! min(npoigx,ipoigx )  MP 1 May 19
   iny1      = min(npoigy,ipoigy + diny+1)  ! min(npoigy,ipoigy+1)  ! min(npoigy,ipoigy )
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      distx2 = ( x(1,ipoin) - coorg(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x(2,ipoin) - coorg(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      dist   = dist/distg                                                         ! distance is normalized
      auxv(1+idest,ipoig) = auxv(1+idest,ipoig) + rho(ipoin)/dist                 ! h  s,w
      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + vx(1,ipoin)/dist                ! vx        
      if(ndimn.eq.2) auxv(3+idest,ipoig) = auxv(3+idest,ipoig) + vx(2,ipoin)/dist ! vy
!      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + hsml(ipoin)/dist               ! tmp hsml s and w
      if ( Gid_Mask_SW (4).eq.1 ) then
          if (ic_ws_interact.eq.2) then
              if (ipoin.le.npoin_s) then
                  auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/h_df(ipoin))/dist 
              endif
          else
              auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/rho(ipoin))/dist      ! pwp
          endif
      elseif ( Gid_Mask_SW (4).eq.-1 ) then
          auxv(4+idest,ipoig) = auxv(4+idest,ipoig) +  pwp_total/dist             ! changed mp Total press
      endif
      auxv(5+idest,ipoig) = auxv(5+idest,ipoig) + 1.00/dist                       ! dist  
      auxv(11,ipoig)      = auxv(11,ipoig) +  rho(ipoin)/dist
      auxv(12,ipoig)      = auxv(12,ipoig) + 1.00/dist                  ! dist. we are adding Hsoil and water 
      if (ic_hrelSat.eq.1.AND.ipoin.le.npoin_s) then
           auxv(13,ipoig)      = auxv(13,ipoig) + hrelSat_DF(ipoin)/dist
      elseif (ic_hrelSat.ne.1) then
           auxv(13,ipoig)      = auxv(13,ipoig) + hsml(ipoin)/dist
      endif

     ! --------   to Plot peclet number when a Stefan problem July 2021 MP

      if (ic_Stef_eros.eq.1.AND.Gid_Mask_SW (7).eq.1) then  
          hh = rho (ipoin)
          vv = vx(1,ipoin)*vx(1,ipoin)
          if (ndimn.eq.2) then
             vv = vv + vx(2,ipoin)*vx(2,ipoin)  
          endif
          vv = vv**0.5
          ee = 0.0
          if ( allocated (activ_Pwp_Stef_eros_Pt ) ) then
             ee = erosion_rate_V(ipoin)
          endif
          Te = 1.e6
          if (ee.gt.1.e-6) then
             Te = hh/ee
          endif
          Cvv = 4.0*const(13)/pi**2
          Tc  = hh*hh/Cvv
          Pe  = min(333.0, Tc/Te)
          auxv(14,ipoig)      = auxv(14,ipoig) + Pe/dist
      endif
   enddo
   enddo
ENDDO

DO ipoig = 1, npoig
   if (auxv(5,ipoig).ge.1.e-6) then
       auxv(1:4,ipoig) = auxv(1:4,ipoig)/auxv(5,ipoig)
       if (ic_hrelSat.eq.1) then
           auxv(13,ipoig)      = auxv(13,ipoig)/auxv(5,ipoig)
       endif
       if ( allocated (activ_Pwp_Stef_eros_Pt)) then  ! Peclet computation at nodes
           auxv(14,ipoig)      = auxv(14,ipoig)/auxv(5,ipoig)
       endif    
   endif
   if (auxv(10,ipoig).ge.1.e-6) then
       auxv(6:9,ipoig) = auxv(6:9,ipoig)/auxv(10,ipoig) 
   endif
   if (auxv(12,ipoig).ge.1.e-6) then
       auxv(11, ipoig) = auxv(11,ipoig)/auxv(12,ipoig)
   endif      

   if (auxv(1,ipoig).ge.1.e-6) ictrack(ipoig)=1      !   ** changed to 0.001 <- 0.1
enddo

auxv(11,:)=auxv(1,:)+auxv(6,:)                      ! This is hs(1) + hw(6)

!   ************   changed 21st July 2014 : avoid white walls

DO ipoig  = 1, npoig
   IF (ictrack (ipoig).eq.1) THEN
       xx     = coorg(1,ipoig)
       yy     = coorg(2,ipoig)
       ipoigx = (xx-xming)/deltxg + 1
       ipoigy = (yy-yming)/deltyg + 1
       inx0   = max(1,ipoigx-1)
       iny0   = max(1,ipoigy-1)
       inx1   = min(npoigx,ipoigx+1)
       iny1   = min(npoigy,ipoigy+1)
       DO inx = inx0, inx1
       DO iny   = iny0, iny1
          idest = (iny-1)*npoigx + inx
          if (ictrack(idest).eq.0) then
              ictrack(idest) = 2
          endif
       ENDDO
       ENDDO
   ENDIF
ENDDO

! where(ictrack.eq.2) ictrack = 1

!   ************   changed 21st July 2014 : avoid white walls


!      ------   And now plot:
 

if(.not.allocated(etag)) then    
   allocate (etag(npoig))
   etag(:)=auxv(11,:) + topol(1,:)
endif
 
   
!      ------ 1:  h soil  -----------------------------

if ( Gid_Mask_SW (1).eq.1 ) then  
   write(gid_res,*) 'Result "height soil" "Height soil" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots h SOIL at quad nodes  ADDED
      if ( ictrack(ipoig).eq.0) CYCLE   
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(1,ipoig)             ! max (0.1,auxv(1,ipoig))
   enddo
   write(gid_res,*) ' End Values ' 
endif 

!      ------  2 Displacements at particles   -------------

if ( Gid_Mask_SW (2).eq.1 ) then    
   write(gid_res,*) 'Result "dis" "Disp" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values  '  
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
      write(gid_res,*) i + 2*npoig, xdisp, ydisp, zdisp
   ENDDO
   write(gid_res,*) ' End Values ' 
endif    
  

!      ------  Velocities  -----------------------------

if ( Gid_Mask_SW (3).eq.1 ) then

    write(gid_res,*) 'Result "vel" "veloc" ',time,' Vector OnNodes "where" '
    write(gid_res,*) ' Values '
    do ipoig = 1, npoig
       if ( auxv(1,ipoig).le.1.e-6) CYCLE
       write(gid_res,*) ipoig, auxv(2,ipoig), auxv(3,ipoig), ' 0.0 '
    enddo
    do i = 1,npoin 
       if ( if_Out_Domain(i).eq.1) CYCLE
       vxx = vx(1,i)
       if (ndimn.eq.1) then
           vyy = 0.0
       else
           vyy = vx(2,i)
       endif
       write(gid_res,*) i + 2*npoig, vxx, vyy, ' 0.0 ' 
    enddo
    write(gid_res,*) ' End Values '    
 endif   
   
!      ------ 4 And PWPs  at the base -----------------------

if ( iabs(Gid_Mask_SW (4)).eq.1 ) then   
   if (icpwp.eq.1.or.icpwp.eq.3.or.icpwp.eq.2) then
      write(gid_res,*) 'Result "Pwp" "pwpress" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig 
         if ( auxv(1,ipoig).le.1.e-6) CYCLE
         write(gid_res,*) ipoig, '   0.0   0.0   ', auxv(4,ipoig)
      enddo
   !   do i = 1,npoin 
   !      if ( if_Out_Domain(i).eq.1) CYCLE
   !      write(gid_res,*) i + 2*npoig, '  0.0   0.0   ', u(i)
   !   enddo
      write(gid_res,*) ' End Values '
   endif   
endif    

!      ------ 5 ErodedTopo  -----------------------------

if ( Gid_Mask_SW (5).eq.1 ) then  
   write(gid_res,*) 'Result "erodedTopo" "-erodedTopo-" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      write(gid_res,*) ipoig, ' 0.   0.   ',   topol(14,ipoig)      
   enddo
   write(gid_res,*) ' End Values ' 

   write(gid_res,*) 'Result "AvailErosD" "ErosDepth" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      write(gid_res,*) ipoig, ' 0.   0.   ',   topol(15,ipoig)      
   enddo
   write(gid_res,*) ' End Values ' 
    
   write(gid_res,*) 'Result "hs over topo" "-hs over topo-" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      if (ic_basal_erosion.eq.1) then
         write(gid_res,*) ipoig, ' 0.   0.   ',   topol(1,ipoig)  - topol(14,ipoig)  
      else
         write(gid_res,*) ipoig, ' 0.   0.   ',   topol(1,ipoig) ! - topol(14,ipoig)  
      endif    
   enddo
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      i = npoig + ipoig
      if (ic_basal_erosion.eq.1) then
         if (ic_ws_interact.eq.2) then
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(11,ipoig)   - topol(14,ipoig) 
         else 
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(1,ipoig)     - topol(14,ipoig)
         endif 
      else
         if (ic_ws_interact.eq.2) then
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(11,ipoig)   
         else 
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(1,ipoig)     
         endif 
      endif            
   enddo                                   !      ------  plots TWO MESHES: S                                 

   write(gid_res,*) ' End Values ' 
   
endif   

!      ------ 6  Topo Z  -----------------------------

if ( Gid_Mask_SW (6).eq.1 ) then
   write(gid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig 
      if (allocated(hs_plus_Z_topo)) then 
          write(gid_res,*) ipoig, ' 0.   0.   ', hs_plus_Z_topo (ipoig)  
      else
          write(gid_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
      endif
   enddo
   write(gid_res,*) ' End Values ' 
endif

!      ------7   h REL soil  -----------------------------  
!              This is an alternative evaluation and plotting of height

if ( Gid_Mask_SW (7).eq.999 ) then
   
   write(gid_res,*) 'Result "hs" "hs" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '

   dArea    = deltxg*deltyg
   if(.NOT.allocated(height_g) ) allocate ( height_g(npoig) )
   height_g = 0.0
   do ipoin = 1, npoin                                   
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      xx = x(1,ipoin)  
      yy = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,ipoin)
      ipoigx = (xx-xming)/deltxg  + 1
      ipoigy = (yy-yming)/deltyg  + 1
      if (ipoigx.eq.npoigx) then
         if (abs(xmaxg-xx).le.1.e-3) then
            ipoigx = ipoigx -1
         else
            write(*,*) ' point x outside topo mesh', xx, xmaxg
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
      xref   = xming + (ipoigx-1)*deltxg
      yref   = yming + (ipoigy-1)*deltyg
      xi     = (xx-xref)/deltxg
      eta    = (yy-yref)/deltyg
      sh1    = (1.-xi) * (1.-eta)
      sh2    =     xi  * (1.-eta)
      sh3    =     xi  * eta
      sh4    = (1.-xi) * eta
      height_g(n1) = height_g(n1) + mass(ipoin) 
      height_g(n2) = height_g(n2) + mass(ipoin) 
      height_g(n3) = height_g(n3) + mass(ipoin) 
      height_g(n4) = height_g(n4) + mass(ipoin) 
   enddo
   height_g = height_g/(4.*dArea)
   do ipoig = 1, npoig                                     
      if ( height_g(ipoig).le.0.005) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', height_g(ipoig)
   enddo
   deallocate ( height_g )

   write(gid_res,*) ' End Values ' 

elseif ( Gid_Mask_SW (7).eq.1 ) then   ! Peclet
   if ( allocated (activ_Pwp_Stef_eros_Pt) ) then 
     write(gid_res,*) 'Result "Peclet nr" "Peclet" ',time,' Vector OnNodes "where" '
     write(gid_res,*) ' Values '
     do ipoig = 1, npoig                                    !      ------  plots Peclet  ADDED
        if ( ictrack(ipoig).eq.0) CYCLE   
        write(gid_res,*) ipoig, ' 0.   0.   ', auxv(14,ipoig)          ! max (0.1,auxv(1,ipoig))
     enddo
     write(gid_res,*) ' End Values '  
  endif     
endif

!      ------ 8   h water -----------------------------

if ( Gid_Mask_SW (8).eq.1 ) then
   write(gid_res,*) 'Result "height water" "Height water" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots h WATER at quad nodes
!      if ( auxv(6,ipoig).le.0.3) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(6,ipoig)
   enddo
   write(gid_res,*) ' End Values ' 
endif   

!      ------ 9  Eta (wave height)  -------------------- 

if ( Gid_Mask_SW (9).eq.1 ) then

   if (ic_ws_Interact.eq.2) then
      write(gid_res,*) 'Result " n " " porosity" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
         do ipoig = 1, npoig                                     
            if ( (auxv(1,ipoig)+auxv(6,ipoig)).le.0.01) CYCLE
            xi = auxv(6,ipoig)/(auxv(1,ipoig)+auxv(6,ipoig))
            if (ic_hrelSat.eq.1.AND. xi.LT.const(35)) xi = const(35)           !   porosity limited by C35   
            write(gid_res,*) ipoig, ' 0.   0.   ', xi                          !      ------  plots h WATER at quad nodes ONLY Z>SWL=0
         enddo 
   else 
      write(gid_res,*) 'Result "eta" " sobreelev" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig
         if ( auxv(6,ipoig).le.0.3) CYCLE
         write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig) + topol(1,ipoig)-etag(ipoig)!   ------  CHANGED
      enddo
   endif
   
   write(gid_res,*) ' End Values ' 
   
endif   

!      ------ 10  Joint h soil + water  ----------------

if ( Gid_Mask_SW (10).eq.1 ) then 
   write(gid_res,*) 'Result "jointh" "JOINTh" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S+W
!      if ( auxv(11,ipoig).le.0.3.or.auxv(6,ipoig).le.0.3) CYCLE  *** 23 July 14
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)
   enddo
   do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S
      i = npoig + ipoig
!      if ( auxv(1,ipoig).le.0.3) CYCLE *** 23 July 14
      write(gid_res,*) i, ' 0.   0.   ', auxv(1,ipoig)
   enddo
   write(gid_res,*) ' End Values ' 
endif   
    
!      ------11   hsml  OR hrelSat_DF  -----------------------------


if ( Gid_Mask_SW (11).eq.1 ) then  
      if (ic_hrelSat.eq.1) then
          write(gid_res,*) 'Result "hrelSat" "-hrelSat-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)     
          enddo
          do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S
             i = npoig + ipoig
             write(gid_res,*)      i, ' 0.   0.   ', auxv(13,ipoig)*auxv(11,ipoig)
          enddo
          write(gid_res,*) ' End Values ' 
      elseif (ic_hrelSat.ne.1) then
          write(gid_res,*) 'Result "hsml" "-HSML-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(13,ipoig)     
          enddo
          write(gid_res,*) ' End Values ' 
      endif
      
      if (ic_hrelSat.eq.1) then
          write(gid_res,*) 'Result "hrSat_only" "-hrSat-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(13,ipoig)     
          enddo
          write(gid_res,*) ' End Values ' 
      endif
endif  

if ( Gid_Mask_SW (11).eq.-1 ) then  
   write(gid_res,*) 'Result "F or I" "-F or I -" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  
      vxx = auxv(2,ipoig) * auxv(2,ipoig)
      vyy = 0.0
      if (ndimn.eq.2) vyy = auxv(3,ipoig) * auxv(3,ipoig)
      vvv = (vxx + vyy)**0.5
      hhh = auxv(1,ipoig)
      cgra = const(1)
      if (hhh.gt.const(10)) then
         Froude = vvv/((hhh*gravz)**0.5)
         if (const(5).eq.8) then
             if (hhh.lt.const(6)) hhh=const(6)
             Froude = Froude * (const(6)/hhh)
         endif
      else
         Froude = 0.0
      endif
      write(gid_res,*) ipoig, ' 0.   0.   ', Froude     
   enddo
   write(gid_res,*) ' End Values ' 
endif 


!      ------  Pwp (1d full)  -------------------- 

if ( Gid_Mask_SW (12).eq.1.and.ndimn.eq.1.and.nUaux.gt.1 ) then

   write(gid_res,*) 'Result "Disp pwp" "disp pore press" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   ipoig = 0
   do ipoin = 1, npoin_s                                         !      ------  plots disp at particles

      !IF ( if_Out_Domain(ipoin).eq.1) then
      !        ipoig = ipoig + nUaux
      !endif

      IF ( if_Out_Domain(ipoin).eq.1) CYCLE

      disp(:) = x(:,ipoin) - x00(:,ipoin)
      xdisp   = disp(1)
      if (ndimn.eq.1) then
          ydisp = 0.0
      else
         ydisp = disp(2)
      endif
      xx     =  x(1,ipoin)
      yy     = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,ipoin) 
         
      call Get_Z_Topo (ipoin, xx,yy,Zp)
   
      zdisp  = Zp - Z00(ipoin)
      
      hhh    = rho(ipoin)   
      hh0    = rho00(ipoin)                         !  ch mp 20th Jan 2019 pwp plot
      if (allocated (h_df)) then
            hhh    = h_df (ipoin)
            hh0    = rho00(ipoin) + rho00(ipoin+npoin_s)
      endif
      dhh     = (hhh-hh0)/(nUaux-1)                 !******* 
      
      do iaux  = 1, nUaux
         ipoig = ipoig + 1
         idest = 2*npoig + npoin + ipoig
         write(gid_res,*) idest, xdisp, ydisp, zdisp + (iaux-1)*dhh  
      enddo
   enddo
   write(gid_res,*) ' End Values ' 

   write(gid_res,*) 'Result "Full  pwp" "full pwp" ',time,' Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   ipoig = 0
   do ipoin = 1, npoin_s  

      IF ( if_Out_Domain(ipoin).eq.1) CYCLE

      do iaux  = 1, nUaux
         ipoig = ipoig + 1
         idest = 2*npoig + npoin + ipoig
         hhh    = rho(ipoin)                            !  ch mp 02 july 17 for df's
         if (allocated (h_df)) hhh    = h_df (ipoin)
         write(gid_res,*) idest,   Uaux (iaux,ipoin)/hhh 
      enddo
   enddo
   write(gid_res,*) ' End Values ' 
   
   
endif   

!      ------  Output for file 3  GIS nodes

xcg = 0.0; ycg = 0.0; zcg = 0.0; Mtotal_out = 0.0

npoing = 0
do ipoin = 1, npoin
   if (if_out_domain(ipoin).eq.0)  then  !   *** only for pts insid
      npoing = npoing +1
      xx = x(1,ipoin)
      yy = (ymaxg+yming)/2.
      if (ndimn.eq.2) yy = x(2,ipoin)
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      xcg = xcg + mass(ipoin)*xx
      ycg = ycg + mass(ipoin)*yy
      zcg = zcg + mass(ipoin)*Zp
      Mtotal_out = Mtotal_out + mass(ipoin)
   endif
enddo   
xcg = xcg / Mtotal_out
ycg = ycg / Mtotal_out
zcg = zcg / Mtotal_out
write(res_file3,9998) time, xcg, ycg, zcg, Mtotal_out   ! mp 1st Dec 2017
9998 format ( ' t xcg ycg zcg Mtot  ', 5(2x,g10.4))


! *******  changed 3 march

!      ------  Output for file 4  GIS nodes

write(res_file4,*) ' ipoin   x   y   h  vx  vy  pwp  r w  t = ', time, ' M ' , Mtotal_out
do ipoin=1,npoin
   if (ic_crush.eq.1) then
      write(res_file4,9999) ipoin, (x(i,ipoin),i=1,ndimn), rho(ipoin), (vx(i,ipoin),i=1,ndimn), u(ipoin), r_crush(ipoin), w_crush(ipoin)
   else
      write(res_file4,9999) ipoin, (x(i,ipoin),i=1,ndimn), rho(ipoin), (vx(i,ipoin),i=1,ndimn), u(ipoin) 
   endif
enddo 
9999 format (2x, i6, 2x,8(2x,f10.4))

!      ------  Output for file 5  GIS erosion

if (gid_mask_SW(5).eq.1) then
    write(res_file5,*) ' ipoig   x   y   z  eroded  ', '  time = ', time
    do ipoig=1,npoig
       write(res_file5,9999) ipoig,  coorg(1,ipoig), coorg(2,ipoig), topol(14,ipoig)
    enddo
endif      

!      ------  Output for file 6   pwp profiles   **** MP 11-Jan 2019

if (gid_mask_SW(12).eq.1.AND.nUAux.GT.1) then
    write(res_file6,*) ' ..... time = ',  time, ' pwp bottom -> pwp top relatives '
    do ipoin=1,npoin_s
       if (ic_ws_interact.eq.2) then
           hhh = h_DF(ipoin)
       else
           hhh = rho(ipoin)
       endif
       write (res_file6,9997) ipoin, (( Uaux (iaux,ipoin)/hhh ),iaux=1,nUaux) 
    enddo
endif  
9997 format ( 2x, i6, 20(1x,f8.4))    

! *******  

deallocate ( disp )
deallocate ( auxv )

 
End Subroutine OutputRes_WIR_BAK23


!-------------------------------------------------------------------

       Subroutine OutputRes_WIR
       
!-------------------------------------------------------------------


!      ------  OUT for waves in reservoirs --------



!      Filter for GID
!
!   We will use a mask to plot variables, nGidMaskSw = 12
!   Gid_Mask_SW (nGidMask) can be zero (no plot) or one (plot)
!   SW case:
!          
!       1 ... h soil       
!       2 ... displacement          
!       3 ... velocity          
!       4 ... Pwp ( -1 -> TOTAL kpA 1 PWP REL)         
!       5 ... erosion          
!       6 ... Z         
!       7 ... h soil relative  (used in Lausanne)// Peclet         
!       8 ... h water          
!       9 ... eta  (when (ic_ws_Interact.eq.2) is porosity n        
!      10 ... hs + hw
!      11 ... hsml OR hrelSat_DF (-1 IS FROUDE)
!      12 .   indicates full pwp output in 1d problems
!
!      mask 9 -> do not plot void points
implicit none

integer  (ink)  i, ipoin, ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer  (ink)  im, nm, index_topo_zone        ! Saeidmt 12Nov2023      im is the index of mass in zones    index_topo_zone is added                              
integer  (ink)  inx0, inx1, inx, iny0, iny1, iny, npoing
integer  (ink)  idest, iaux
real     (irk)  xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp, dhh
integer  (ink)  dinx, diny     ! MP 01 may 2019  
real     (irk)  xdisp, ydisp, zdisp, vxx, vyy, vvv, hhh, hh0, cgra, Froude
real     (irk), allocatable::disp(:)
real     (irk), allocatable::auxv(:,:)
real     (irk), allocatable, SAVE:: etag   (:)       
integer  (ink), allocatable, SAVE:: ictrack(:)      
real     (irk)  distx2, disty2, dist, distg 
real     (irk)  dArea, Mtotal_out
real     (irk)  xcg, ycg, zcg       ! centre of gravity avalanche MP 1st dec 2017
real(irk), allocatable:: height_g(:)
real(irk), allocatable:: mzones(:) ! Saeidmt 12Nov2023 mzones is the total mass in each zone
   
integer  (ink), allocatable::  Tmp_MASk(:)   ! we will modify Mask.   June 2023 MP
                                             ! Keep orig here, recover at the end

real    (irk) dens, denss, densw        !   ic_fct_dens = C20  0 no correction  1.. correct 
real    (irk) porosity, Sr              !   Porosity and degree of saturation
real    (irk) dens1, dens2              !   mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    !   average and effective average dens of 2 layers
real    (irk) alphaS                    !   ratio hsat/h, from hrelpw. If below 0.1 -> 0

real    (irk) deng, gravz
real    (irk) height, pwp_total, Zpx, Zpy

real   (irk) hh,tt,xt
real   (irk) vv, ee, Te, Tc, Pe, Cvv
real   (irk) Voellmy_Gid
real   (irk) total_mass, percentage_mass              ! Saeidmt 12Nov2023
integer (ink) ic_fct_dens, ic_cosTheta
integer (ink) nfric_gid                 !   MP July 2023 to plot 202 frict case    

if(.NOT.allocated(disp) )    allocate ( disp(ndimn) ) 
if(.NOT.allocated(auxv) )    allocate ( auxv(15,npoig))    !  MP July 2023 15 for porosity 202 frict
if(.NOT.allocated(TMP_MAsk)) allocate ( Tmp_Mask (nGidMaskSw) )
auxv = 0.0                  ! MP July 2023 
Tmp_Mask = Gid_Mask_SW     !  MP June 2023
where (Gid_Mask_SW.eq.-9) Gid_Mask_SW = -1
where (Gid_Mask_SW.eq. 9) Gid_Mask_SW =  1

cgra        = const  (1)                
dens        = const  (2)            
denss       = const (17)             
densw       = const (18)
Sr          = const (19)
alphaS      = const (14)
ic_fct_dens = const (20)
ic_cosTheta = const (16)
if (( Gid_Mask_SW (4).eq.-1 ).and. ic_fct_dens.eq.0) then
    write (*,*) ' take care C20 ic_fct_dens should be 1 or 2 '
    write (*,*) ' hit a numerical key ' 
    read  (*,*) i
    stop
endif

Call Check_Out_Domain_SW

!      ------   First of all we get values at topo grid. 
!               Add for each cell value*(1/dist) and divide by sum (1/dist)

auxv    = 0.0

if(.not.allocated(ictrack)) then                   
   allocate (ictrack(npoig))
   ictrack = 0
endif
                                                   
distg   = (deltxg*deltxg + deltyg*deltyg)**0.5

DO ipoin = 1,npoin    
   if (allocated (if_correct)) then             ! CHGD MP  april 2020. if_correct is not always allocated
      if (if_correct(ipoin).eq.1) CYCLE         ! CHGD MP 17th March 2020  skip node to avoid high porosities
   endif
   pwp_total = 0.0                                                  ! ** mp 20th Jan 2019
   if ( (Gid_Mask_SW (4).eq.-1) .AND. (ipoin.LE.npoin_s) ) then     ! ** mp 20th Jan 2019 avoid w points
       xx  = x(1,ipoin)             
       yy  = (yming+ymaxg)/2     
       if (ndimn.eq.2) yy = x(2,ipoin)                                                      
       if ( ic_cosTheta.eq.1) then                          !  if we use the theta formulation   
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)     !  we will define them again              
          deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
          gravz =  cgra/deng                              ! *** think of joining with nfrict 7 
       elseif ( ic_cosTheta.eq.0) then 
!          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)                 
!          deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
          gravz =  cgra 
       elseif ( ic_cosTheta.eq.2) then                    !  do nothing, we have gx gy gz 
          deng  = deng
       endif
          
       if (ic_ws_Interact.eq.2) then
          height   = h_DF(ipoin)
          porosity = n_DF(ipoin)
       else
          height   = rho(ipoin)  
          porosity = (denss-dens)/(denss-densw)
       endif 
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !     MP 19 January 2019 
     
       pwp_total = gravz * (densw * height * alphaS + dens_bar_eff * u(ipoin) )
       pwp_total = pwp_total/1000.   ! KPa
     
   endif
                        
   if ( if_Out_Domain(ipoin).eq.1 )  CYCLE  ! ** mp 20th Jan 2019 avois water pts
   idest = 5*(iabs(itype(ipoin))-1)                                 !  = 0 for soil and 5 for water
   xx = x(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
    inx0      = max(1,ipoigx-1)            ! activated mp July 2023
    iny0      = max(1,ipoigy-1)
    inx1      = min(npoigx,ipoigx+2)
    iny1      = min(npoigy,ipoigy+2)
   dinx      = hsml(ipoin)/deltxg 
   diny      = hsml(ipoin)/deltyg 
   if (ndimn.eq.1) diny = 0 
 !  inx0      = max(1,ipoigx-dinx)
 !  iny0      = max(1,ipoigy-diny)
 !  inx1      = min(npoigx,ipoigx + dinx+1)  !  min(npoigx,ipoigx+1) ! min(npoigx,ipoigx )  MP 1 May 19
 !  iny1      = min(npoigy,ipoigy + diny+1)  ! min(npoigy,ipoigy+1)  ! min(npoigy,ipoigy )
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      distx2 = ( x(1,ipoin) - coorg(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x(2,ipoin) - coorg(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      dist   = dist/distg                                                         ! distance is normalized
      auxv(1+idest,ipoig) = auxv(1+idest,ipoig) + rho(ipoin)/dist                 ! h  s,w
      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + vx(1,ipoin)/dist                ! vx        
      if(ndimn.eq.2) auxv(3+idest,ipoig) = auxv(3+idest,ipoig) + vx(2,ipoin)/dist ! vy
!      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + hsml(ipoin)/dist               ! tmp hsml s and w
      if ( Gid_Mask_SW (4).eq.1 ) then
          if (ic_ws_interact.eq.2) then
              if (ipoin.le.npoin_s) then
                  auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/h_df(ipoin))/dist 
              endif
          else
              auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/rho(ipoin))/dist      ! pwp
          endif
      elseif ( Gid_Mask_SW (4).eq.-1 ) then
          auxv(4+idest,ipoig) = auxv(4+idest,ipoig) +  pwp_total/dist             ! changed mp Total press
      endif
      auxv(5+idest,ipoig) = auxv(5+idest,ipoig) + 1.00/dist                       ! dist  
      auxv(11,ipoig)      = auxv(11,ipoig) +  rho(ipoin)/dist
      auxv(12,ipoig)      = auxv(12,ipoig) + 1.00/dist                  ! dist. we are adding Hsoil and water 
      if (ic_Fs.ne.1.AND.ic_hrelSat.eq.1.AND.ipoin.le.npoin_s) then
           auxv(13,ipoig)      = auxv(13,ipoig) + hrelSat_DF(ipoin)/dist
      elseif (ic_hrelSat.ne.1) then
           auxv(13,ipoig)      = auxv(13,ipoig) + hsml(ipoin)/dist
      endif

     ! --------   to Plot peclet number when a Stefan problem July 2021 MP

      if (ic_Stef_eros.eq.1.AND.Gid_Mask_SW (7).eq.1) then  
          hh = rho (ipoin)
          vv = vx(1,ipoin)*vx(1,ipoin)
          if (ndimn.eq.2) then
             vv = vv + vx(2,ipoin)*vx(2,ipoin)  
          endif
          vv = vv**0.5
          ee = 0.0
          if ( allocated (activ_Pwp_Stef_eros_Pt ) ) then
             ee = erosion_rate_V(ipoin)
          endif
          Te = 1.e6
          if (ee.gt.1.e-6) then
             Te = hh/ee
          endif
          Cvv = 4.0*const(13)/pi**2
          Tc  = hh*hh/Cvv
          Pe  = min(333.0, Tc/Te)
          auxv(14,ipoig)      = auxv(14,ipoig) + Pe/dist
      endif

      nfric_Gid = 0 
      if (Gid_Mask_SW (7).eq.1) then    ! Voellmy      MP July 2023 
          nfric_Gid = const(5) + 0.001 
          if (nfric_Gid.eq.202.OR.nfric_Gid.eq.203) then   ! MP Dec 2023
             if ( allocated (ns_dep_props ) ) then
                Voellmy_Gid = ns_dep_props(1,ipoin)
             else
                Voellmy_Gid = 0.0
             endif
             auxv(15,ipoig) = auxv(15,ipoig) + n_DF(ipoin)/dist
             auxv(14,ipoig) = auxv(14,ipoig) + Voellmy_Gid/dist
          endif 
      endif


   enddo
   enddo
ENDDO

DO ipoig = 1, npoig
   if (auxv(5,ipoig).ge.1.e-6) then
       auxv(1:4,ipoig) = auxv(1:4,ipoig)/auxv(5,ipoig)
       if (ic_hrelSat.eq.1.and.nfric_GID.ne.202) then          !  MP July 2023
           auxv(13,ipoig)      = auxv(13,ipoig)/auxv(5,ipoig)
       endif
       if ( allocated (activ_Pwp_Stef_eros_Pt)) then  ! Peclet computation at nodes
           auxv(14,ipoig)      = auxv(14,ipoig)/auxv(5,ipoig)
       endif    
   endif
   if (auxv(10,ipoig).ge.1.e-6) then
       auxv(6:9,ipoig) = auxv(6:9,ipoig)/auxv(10,ipoig) 
   endif
   if (auxv(12,ipoig).ge.1.e-6) then
       auxv(11, ipoig) = auxv(11,ipoig)/auxv(12,ipoig)
   endif  
   
   if (Gid_Mask_SW (7).eq.1.AND.nfric_Gid.eq.202) then    ! Voellmy      MP July 2023 
       auxv(14,ipoig)      = auxv(14,ipoig)/auxv(5,ipoig)  
       auxv(15,ipoig)      = auxv(15,ipoig)/auxv(5,ipoig)                
   endif

   if (auxv(1,ipoig).ge.1.e-6) ictrack(ipoig)=1      !   ** changed to 0.001 <- 0.1
enddo

auxv(11,:)=auxv(1,:)+auxv(6,:)                      ! This is hs(1) + hw(6)

!   ************   changed 21st July 2014 : avoid white walls

DO ipoig  = 1, npoig
   IF (ictrack (ipoig).eq.1) THEN
       xx     = coorg(1,ipoig)
       yy     = coorg(2,ipoig)
       ipoigx = (xx-xming)/deltxg + 1
       ipoigy = (yy-yming)/deltyg + 1
       inx0   = max(1,ipoigx-1)
       iny0   = max(1,ipoigy-1)
       inx1   = min(npoigx,ipoigx+1)
       iny1   = min(npoigy,ipoigy+1)
       DO inx = inx0, inx1
       DO iny   = iny0, iny1
          idest = (iny-1)*npoigx + inx
          if (ictrack(idest).eq.0) then
              ictrack(idest) = 2
          endif
       ENDDO
       ENDDO
   ENDIF
ENDDO

! where(ictrack.eq.2) ictrack = 1

!   ************   changed 21st July 2014 : avoid white walls


!      ------   And now plot:
 

if(.not.allocated(etag)) then    
   allocate (etag(npoig))
   etag(:)=auxv(11,:) + topol(1,:)
endif
 
   
!      ------ 1:  h soil  -----------------------------

if ( Gid_Mask_SW (1).eq.1 ) then  
   write(gid_res,*) 'Result "height soil" "Height soil" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots h SOIL at quad nodes  ADDED
      !if ( ictrack(ipoig).eq.0) CYCLE                                  !  Saeidmt 31March2022   
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(1,ipoig)             ! max (0.1,auxv(1,ipoig))
   enddo
   write(gid_res,*) ' End Values ' 
endif 

!      ------  2 Displacements at particles   -------------

if ( Gid_Mask_SW (2).eq.1 ) then    
   write(gid_res,*) 'Result "dis" "Disp" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values  '  
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
      write(gid_res,*) i + 2*npoig, xdisp, ydisp, zdisp
   ENDDO
   write(gid_res,*) ' End Values ' 
endif    
  

!      ------  Velocities  -----------------------------

if ( Gid_Mask_SW (3).eq.1 ) then

    write(gid_res,*) 'Result "vel" "veloc" ',time,' Vector OnNodes "where" '
    write(gid_res,*) ' Values '
    do ipoig = 1, npoig
       ! if ( auxv(1,ipoig).le.1.e-6) CYCLE
       if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(3).eq.9) CYCLE     !  MP June 2023
       write(gid_res,*) ipoig, auxv(2,ipoig), auxv(3,ipoig), ' 0.0 '
    enddo
    do i = 1,npoin 
       if ( if_Out_Domain(i).eq.1) CYCLE
       vxx = vx(1,i)
       if (ndimn.eq.1) then
           vyy = 0.0
       else
           vyy = vx(2,i)
       endif
       write(gid_res,*) i + 2*npoig, vxx, vyy, ' 0.0 ' 
    enddo
    write(gid_res,*) ' End Values '    
 endif   
   
!      ------ 4 And PWPs  at the base -----------------------

if ( iabs(Gid_Mask_SW (4)).eq.1 ) then   
   if (icpwp.eq.1.or.icpwp.eq.3.or.icpwp.eq.2) then
      write(gid_res,*) 'Result "Pwp" "pwpress" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig 
   !      if ( auxv(1,ipoig).le.1.e-6) CYCLE
           if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(4).eq. 9) CYCLE     !  MP June 2023
         write(gid_res,*) ipoig, '   0.0   0.0   ', auxv(4,ipoig)
      enddo
   !   do i = 1,npoin 
   !      if ( if_Out_Domain(i).eq.1) CYCLE
   !      write(gid_res,*) i + 2*npoig, '  0.0   0.0   ', u(i)
   !   enddo
      write(gid_res,*) ' End Values '
   endif   
endif    

!      ------ 5 ErodedTopo  -----------------------------

if ( Gid_Mask_SW (5).eq.1 ) then  
   write(gid_res,*) 'Result "erodedTopo" "-erodedTopo-" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      write(gid_res,*) ipoig, ' 0.   0.   ',   topol(14,ipoig)      
   enddo
   write(gid_res,*) ' End Values ' 

   write(gid_res,*) 'Result "AvailErosD" "ErosDepth" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      write(gid_res,*) ipoig, ' 0.   0.   ',   topol(15,ipoig)      
   enddo
   write(gid_res,*) ' End Values ' 
    
   write(gid_res,*) 'Result "hs over topo" "-hs over topo-" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      if (ic_basal_erosion.eq.1) then
         write(gid_res,*) ipoig, ' 0.   0.   ',   topol(1,ipoig)  - topol(14,ipoig)  
      else
         write(gid_res,*) ipoig, ' 0.   0.   ',   topol(1,ipoig) ! - topol(14,ipoig)  
      endif    
   enddo
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      i = npoig + ipoig
      if (ic_basal_erosion.eq.1) then
         if (ic_ws_interact.eq.2) then
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(11,ipoig)   - topol(14,ipoig) 
         else 
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(1,ipoig)     - topol(14,ipoig)
         endif 
      else
         if (ic_ws_interact.eq.2) then
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(11,ipoig)   
         else 
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(1,ipoig)     
         endif 
      endif            
   enddo                                   !      ------  plots TWO MESHES: S                                 

   write(gid_res,*) ' End Values ' 
   
endif   

!      ------ 6  Topo Z  -----------------------------

if (.NOT.allocated(hs_plus_Z_topo)) then      ! MP 21st June 2023
   if ( Gid_Mask_SW (6).eq.1.AND.itimestep.EQ.0) then
      write(gid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig
          write(gid_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
      enddo
      write(gid_res,*) ' End Values '
   endif
elseif (allocated(hs_plus_Z_topo)) then    
   if ( Gid_Mask_SW (6).eq.1 ) then
      write(gid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig
          write(gid_res,*) ipoig, ' 0.   0.   ', hs_plus_Z_topo (ipoig) 
      enddo
      write(gid_res,*) ' End Values '
   endif    
endif 

!      ------7   h REL soil  -----------------------------  
!              This is an alternative evaluation and plotting of height

if ( Gid_Mask_SW (7).eq.999 ) then
   
   write(gid_res,*) 'Result "hs" "hs" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '

   dArea    = deltxg*deltyg
   if(.NOT.allocated(height_g) ) allocate ( height_g(npoig) )
   height_g = 0.0
   do ipoin = 1, npoin                                   
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      xx = x(1,ipoin)  
      yy = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,ipoin)
      ipoigx = (xx-xming)/deltxg  + 1
      ipoigy = (yy-yming)/deltyg  + 1
      if (ipoigx.eq.npoigx) then
         if (abs(xmaxg-xx).le.1.e-3) then
            ipoigx = ipoigx -1
         else
            write(*,*) ' point x outside topo mesh', xx, xmaxg
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
      xref   = xming + (ipoigx-1)*deltxg
      yref   = yming + (ipoigy-1)*deltyg
      xi     = (xx-xref)/deltxg
      eta    = (yy-yref)/deltyg
      sh1    = (1.-xi) * (1.-eta)
      sh2    =     xi  * (1.-eta)
      sh3    =     xi  * eta
      sh4    = (1.-xi) * eta
      height_g(n1) = height_g(n1) + mass(ipoin) 
      height_g(n2) = height_g(n2) + mass(ipoin) 
      height_g(n3) = height_g(n3) + mass(ipoin) 
      height_g(n4) = height_g(n4) + mass(ipoin) 
   enddo
   height_g = height_g/(4.*dArea)
   do ipoig = 1, npoig                                     
      if ( height_g(ipoig).le.0.005) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', height_g(ipoig)
   enddo
   deallocate ( height_g )

   write(gid_res,*) ' End Values ' 

elseif ( Gid_Mask_SW (7).eq.1 ) then   ! Peclet
   if ( allocated (activ_Pwp_Stef_eros_Pt) ) then 
     write(gid_res,*) 'Result "Peclet nr" "Peclet" ',time,' Vector OnNodes "where" '
     write(gid_res,*) ' Values '
     do ipoig = 1, npoig                                    !      ------  plots Peclet  ADDED
!        if ( ictrack(ipoig).eq.0) CYCLE 
         if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(7).eq.9) CYCLE     !  MP June 2023  
        write(gid_res,*) ipoig, ' 0.   0.   ', auxv(14,ipoig)          ! max (0.1,auxv(1,ipoig))
     enddo
     write(gid_res,*) ' End Values '  
  elseif (nfric_Gid.EQ.202) then  !   plot voellmy and pososity  MP July 2023
     write(gid_res,*) 'Result "porosity" "poros" ',time,' Vector OnNodes "where" '
     write(gid_res,*) ' Values '
     do ipoig = 1, npoig    
        if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(7).eq.9) CYCLE     !  MP June 2023  
        write(gid_res,*) ipoig, ' 0.   0.   ', auxv(15,ipoig)          !  
     enddo
     write(gid_res,*) ' End Values '  
     
     write(gid_res,*) 'Result "Voellmy" "Voellm" ',time,' Vector OnNodes "where" '
     write(gid_res,*) ' Values '
     do ipoig = 1, npoig    
        if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(7).eq.9) CYCLE     !  MP June 2023  
        write(gid_res,*) ipoig, ' 0.   0.   ', auxv(14,ipoig)          !  
     enddo
     write(gid_res,*) ' End Values '  
  
  endif     
endif

!      ------ 8   h water -----------------------------

if ( Gid_Mask_SW (8).eq.1 ) then
   write(gid_res,*) 'Result "height water" "Height water" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots h WATER at quad nodes
!      if ( auxv(6,ipoig).le.0.3) CYCLE
      if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(8).eq.9) CYCLE 
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(6,ipoig)
   enddo
   write(gid_res,*) ' End Values ' 
endif   

!      ------ 9  Eta (wave height)  -------------------- 

if ( Gid_Mask_SW (9).eq.1 ) then

   if (ic_ws_Interact.eq.2) then
      write(gid_res,*) 'Result " n " " porosity" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
         do ipoig = 1, npoig                                     
            if ( (auxv(1,ipoig)+auxv(6,ipoig)).le.0.01) CYCLE
            xi = auxv(6,ipoig)/(auxv(1,ipoig)+auxv(6,ipoig))
            if (ic_hrelSat.eq.1.AND. xi.LT.const(35)) xi = const(35)      !   porosity limited by C35  
            if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(9).eq.9) CYCLE     !  MP June 2023 
            write(gid_res,*) ipoig, ' 0.   0.   ', xi                     !   plots h WATER at quad nodes ONLY Z>SWL=0
         enddo 
   else 
      write(gid_res,*) 'Result "eta" " sobreelev" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig
        ! if ( auxv(6,ipoig).le.0.3) CYCLE 
         if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(9).eq.9) CYCLE     !  MP June 2023 
         write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig) + topol(1,ipoig)-etag(ipoig)!   ------  CHANGED
      enddo
   endif
   
   write(gid_res,*) ' End Values ' 
   
endif   

!      ------ 10  Joint h soil + water  ----------------

if ( Gid_Mask_SW (10).eq.1 ) then 
   write(gid_res,*) 'Result "jointh" "JOINTh" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S+W
!      if ( auxv(11,ipoig).le.0.3.or.auxv(6,ipoig).le.0.3) CYCLE  *** 23 July 14
      if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(10).eq.9) CYCLE     !  MP June 2023                                                                                                                      
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)
   enddo
   do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S
      i = npoig + ipoig
!      if ( auxv(1,ipoig).le.0.3) CYCLE *** 23 July 14
      if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(10).eq.9) CYCLE     !  MP June 2023
      write(gid_res,*) i, ' 0.   0.   ', auxv(1,ipoig)
   enddo
   write(gid_res,*) ' End Values ' 
endif   
    
!      ------11   hsml  OR hrelSat_DF  -----------------------------


if ( Gid_Mask_SW (11).eq.1 ) then  
      if (ic_hrelSat.eq.1) then
          write(gid_res,*) 'Result "hrelSat" "-hrelSat-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                               !   hsml at quad nodes   
              if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(11).eq.9) CYCLE     !  MP June 2023
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)     
          enddo
          do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S
             i = npoig + ipoig
             if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(11).eq.9) CYCLE
             write(gid_res,*)      i, ' 0.   0.   ', auxv(13,ipoig)*auxv(11,ipoig)
          enddo
          write(gid_res,*) ' End Values ' 
      elseif (ic_hrelSat.ne.1) then
          write(gid_res,*) 'Result "hsml" "-HSML-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig  
              if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(11).eq.9) CYCLE  !  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(13,ipoig)     
          enddo
          write(gid_res,*) ' End Values ' 
      endif
      
      if (ic_hrelSat.eq.1) then
          write(gid_res,*) 'Result "hrSat_only" "-hrSat-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(11).eq.9) CYCLE
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(13,ipoig)     
          enddo
          write(gid_res,*) ' End Values ' 
      endif
endif  

if ( Gid_Mask_SW (11).eq.-1 ) then  
   write(gid_res,*) 'Result "F or I" "-F or I -" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  
      vxx = auxv(2,ipoig) * auxv(2,ipoig)
      vyy = 0.0
      if (ndimn.eq.2) vyy = auxv(3,ipoig) * auxv(3,ipoig)
      vvv = (vxx + vyy)**0.5
      hhh = auxv(1,ipoig)
      cgra = const(1)
      if (hhh.gt.const(10)) then
         Froude = vvv/((hhh*gravz)**0.5)
         if (const(5).eq.8) then
             if (hhh.lt.const(6)) hhh=const(6)
             Froude = Froude * (const(6)/hhh)
         endif
      else
         Froude = 0.0
      endif
      if ( auxv(12,ipoig).le.1.e-6.AND. tmp_Mask(11).eq.-9) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', Froude     
   enddo
   write(gid_res,*) ' End Values ' 
endif 


!      ------  Pwp (1d full)  -------------------- 

if ( Gid_Mask_SW (12).eq.1.and.ndimn.eq.1.and.nUaux.gt.1 ) then

   write(gid_res,*) 'Result "Disp pwp" "disp pore press" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   ipoig = 0
   do ipoin = 1, npoin_s                                         !      ------  plots disp at particles

      !IF ( if_Out_Domain(ipoin).eq.1) then
      !        ipoig = ipoig + nUaux
      !endif

      IF ( if_Out_Domain(ipoin).eq.1) CYCLE

      disp(:) = x(:,ipoin) - x00(:,ipoin)
      xdisp   = disp(1)
      if (ndimn.eq.1) then
          ydisp = 0.0
      else
         ydisp = disp(2)
      endif
      xx     =  x(1,ipoin)
      yy     = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,ipoin) 
         
      call Get_Z_Topo (ipoin, xx,yy,Zp)
   
      zdisp  = Zp - Z00(ipoin)
      
      hhh    = rho(ipoin)   
      hh0    = rho00(ipoin)                         !  ch mp 20th Jan 2019 pwp plot
      if (allocated (h_df)) then
            hhh    = h_df (ipoin)
            hh0    = rho00(ipoin) + rho00(ipoin+npoin_s)
      endif
      dhh     = (hhh-hh0)/(nUaux-1)                 !******* 
      
      do iaux  = 1, nUaux
         ipoig = ipoig + 1
         idest = 2*npoig + npoin + ipoig
         write(gid_res,*) idest, xdisp, ydisp, zdisp + (iaux-1)*dhh  
      enddo
   enddo
   write(gid_res,*) ' End Values ' 

   write(gid_res,*) 'Result "Full  pwp" "full pwp" ',time,' Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   ipoig = 0
   do ipoin = 1, npoin_s  

      IF ( if_Out_Domain(ipoin).eq.1) CYCLE

      do iaux  = 1, nUaux
         ipoig = ipoig + 1
         idest = 2*npoig + npoin + ipoig
         hhh    = rho(ipoin)                            !  ch mp 02 july 17 for df's
         if (allocated (h_df)) hhh    = h_df (ipoin)
         write(gid_res,*) idest,   Uaux (iaux,ipoin)/hhh 
      enddo
   enddo
   write(gid_res,*) ' End Values ' 
   
   
endif   

!      ------  Output for file 3  GIS nodes

xcg = 0.0; ycg = 0.0; zcg = 0.0; Mtotal_out = 0.0

npoing = 0
do ipoin = 1, npoin
   if (if_out_domain(ipoin).eq.0)  then  !   *** only for pts insid
      npoing = npoing +1
      xx = x(1,ipoin)
      yy = (ymaxg+yming)/2.
      if (ndimn.eq.2) yy = x(2,ipoin)
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      xcg = xcg + mass(ipoin)*xx
      ycg = ycg + mass(ipoin)*yy
      zcg = zcg + mass(ipoin)*Zp
      Mtotal_out = Mtotal_out + mass(ipoin)
   endif
enddo   
xcg = xcg / Mtotal_out
ycg = ycg / Mtotal_out
zcg = zcg / Mtotal_out
write(res_file3,9998) time, xcg, ycg, zcg, Mtotal_out   ! mp 1st Dec 2017
9998 format ( ' t xcg ycg zcg Mtot  ', 5(2x,g10.4))


! *******  changed 3 march

!      ------  Output for file 4  GIS nodes

write(res_file4,*) ' ipoin   x   y   h  vx  vy  pwp  r w  t = ', time, ' M ' , Mtotal_out
do ipoin=1,npoin
   if (ic_crush.eq.1) then
      write(res_file4,9999) ipoin, (x(i,ipoin),i=1,ndimn), rho(ipoin), (vx(i,ipoin),i=1,ndimn), u(ipoin), r_crush(ipoin), w_crush(ipoin)
   else
      write(res_file4,9999) ipoin, (x(i,ipoin),i=1,ndimn), rho(ipoin), (vx(i,ipoin),i=1,ndimn), u(ipoin) 
   endif
enddo 
9999 format (2x, i6, 2x,8(2x,g10.4))

!      ------  Output for file 5  GIS erosion

if (gid_mask_SW(5).eq.1) then
    write(res_file5,*) ' ipoig   x   y   z  eroded  ', '  time = ', time
    do ipoig=1,npoig
       write(res_file5,9999) ipoig,  coorg(1,ipoig), coorg(2,ipoig), topol(14,ipoig)
    enddo
endif      

!      ------  Output for file 6   pwp profiles   **** MP 11-Jan 2019

if (gid_mask_SW(12).eq.1.AND.nUAux.GT.1) then
    write(res_file6,*) ' ..... time = ',  time, ' pwp bottom -> pwp top relatives '
    do ipoin=1,npoin_s
       if (ic_ws_interact.eq.2) then
           hhh = h_DF(ipoin)
       else
           hhh = rho(ipoin)
       endif
       write (res_file6,9997) ipoin, (( Uaux (iaux,ipoin)/hhh ),iaux=1,nUaux) 
    enddo
endif  
9997 format ( 2x, i6, 20(1x,f8.4))    

! *******  

Gid_Mask_SW = Tmp_Mask   ! recover values of Gid_MAsk_SW

!      ------  Output for file 8   Volume of the zones   ****   ! Saeidmt 12Nov2023

nm = 0
if (mass_zones.eq.1.OR.mass_zones.eq.2) then
    total_mass=0.0
    do ipoin=1,npoin
        total_mass = total_mass + mass(ipoin)
    end do
    do ipoin=1,npoin
        im = 0
        xx     =  x(1,ipoin)
        yy     = (yming+ymaxg)/2
        if (ndimn.eq.2) yy = x(2,ipoin)
        if (mass_zones.eq.1) call Which_Topo_Zone(xx, yy, 5, index_topo_zone)
        if (mass_zones.eq.2) call Which_Topo_Zone(xx, yy, 55, index_topo_zone)
        if (index_topo_zone.gt.0) then
            nm = basal_law (index_topo_zone)
            if(.NOT.allocated(mzones) ) then 
                allocate ( mzones(nm))
                 mzones = 0
            end if
            do im = 1, nm
                if (im.ne.index_topo_zone) cycle
                mzones(im) = mzones(im) + mass(ipoin)
            end do
        end if
    end do
    
    do im = 1, nm
        percentage_mass = (mzones(im)/total_mass)*100
        if (mzones(im).LT.0.0001) cycle
        write(res_file7,*) time, im, mzones(im), total_mass, percentage_mass
        !write(res_file8,9996) ' ..... time = ',  time, ' zone = ', im, mzones(im), total_mass, percentage_mass
    end do
    if (index_topo_zone.gt.0) then
        deallocate ( mzones )
    endif
endif 

deallocate ( disp )
deallocate ( auxv )
deallocate (Tmp_Mask )   !  MP June 2023

 
End Subroutine OutputRes_WIR



!-------------------------------------------------------------------

       Subroutine OutputRes_WIR_SAEID_bug 
       
!-------------------------------------------------------------------


!      ------  OUT for waves in reservoirs --------



!      Filter for GID
!
!   We will use a mask to plot variables, nGidMaskSw = 12
!   Gid_Mask_SW (nGidMask) can be zero (no plot) or one (plot)
!   SW case:
!          
!       1 ... h soil       
!       2 ... displacement          
!       3 ... velocity          
!       4 ... Pwp          
!       5 ... erosion          
!       6 ... Z         
!       7 ... h soil relative  (used in Lausanne)         
!       8 ... h water          
!       9 ... eta  (when (ic_ws_Interact.eq.2) is porosity n        
!      10 ... hs + hw
!      11 ... hsml OR hrelSat_DF
!      12 .   indicates full pwp output in 1d problems

implicit none

integer  (ink)  i, ipoin, ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer  (ink)  inx0, inx1, inx, iny0, iny1, iny, npoing
integer  (ink)  idest, iaux
real     (irk)  xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp, dhh
integer  (ink)  dinx, diny     ! MP 01 may 2019  
real     (irk)  xdisp, ydisp, zdisp, vxx, vyy, vvv, hhh, hh0, cgra, Froude
real     (irk), allocatable::disp(:)
real     (irk), allocatable::auxv(:,:)
real     (irk), allocatable, SAVE:: etag   (:)       
integer  (ink), allocatable, SAVE:: ictrack(:)      
real     (irk)  distx2, disty2, dist, distg 
real     (irk)  dArea, Mtotal_out
real     (irk)  xcg, ycg, zcg       ! centre of gravity avalanche MP 1st dec 2017
real(irk), allocatable:: height_g(:)

real    (irk) dens, denss, densw        !   ic_fct_dens = C20  0 no correction  1.. correct 
real    (irk) porosity, Sr              !   Porosity and degree of saturation
real    (irk) dens1, dens2              !   mixture densities of layer 1 (bottom, sat) and 2 top
real    (irk) dens_bar, dens_bar_eff    !   average and effective average dens of 2 layers
real    (irk) alphaS                    !   ratio hsat/h, from hrelpw. If below 0.1 -> 0

real    (irk) deng, gravz
real    (irk) height, pwp_total, Zpx, Zpy

real   (irk) hh,tt,xt
real   (irk) vv, ee, Te, Tc, Pe, Cvv

integer (ink) ic_fct_dens, ic_cosTheta    

if(.NOT.allocated(disp) ) allocate ( disp(ndimn) ) 
if(.NOT.allocated(auxv) ) allocate ( auxv(16,npoig))  !   Saeidmt 26May2021  old**  if(.NOT.allocated(auxv) ) allocate ( auxv(14,npoig))

cgra        = const  (1)                
dens        = const  (2)            
denss       = const (17)             
densw       = const (18)
Sr          = const (19)
alphaS      = const (14)
ic_fct_dens = const (20)
ic_cosTheta = const (16)
if (( Gid_Mask_SW (4).eq.-1 ).and. ic_fct_dens.eq.0) then
    write (*,*) ' take care C20 ic_fct_dens should be 1 or 2 '
    write (*,*) ' hit a numerical key ' 
    read  (*,*) i
    stop
endif

Call Check_Out_Domain_SW

!      ------   First of all we get values at topo grid. 
!               Add for each cell value*(1/dist) and divide by sum (1/dist)

auxv    = 0.0

if(.not.allocated(ictrack)) then                   
   allocate (ictrack(npoig))
   ictrack = 0
endif
                                                   
distg   = (deltxg*deltxg + deltyg*deltyg)**0.5

DO ipoin = 1,npoin    
   if (allocated (if_correct)) then             ! CHGD MP  april 2020. if_correct is not always allocated
      if (if_correct(ipoin).eq.1) CYCLE         ! CHGD MP 17th March 2020  skip node to avoid high porosities
   endif
   pwp_total = 0.0                                                  ! ** mp 20th Jan 2019
   if ( (Gid_Mask_SW (4).eq.-1) .AND. (ipoin.LE.npoin_s) ) then     ! ** mp 20th Jan 2019 avoid w points
       xx  = x(1,ipoin)            !  fixing bug in Hungr: before the rk if 
       yy  = (yming+ymaxg)/2     
       if (ndimn.eq.2) yy = x(2,ipoin)                                                      
       if ( ic_cosTheta.eq.1) then                          !  if we use the theta formulation   
          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)     !  we will define them again              
          deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
          gravz =  cgra/deng                              ! *** think of joining with nfrict 7 
       elseif ( ic_cosTheta.eq.0) then 
!          call Get_Z_Topo (ipoin, xx, yy, Zp, Zpx, Zpy)                 
!          deng  = (1. +Zpx*Zpx + Zpy*Zpy)**0.5
          gravz =  cgra 
       elseif ( ic_cosTheta.eq.2) then                    !  do nothing, we have gx gy gz 
          deng  = deng
       endif
          
       if (ic_ws_Interact.eq.2) then
          height   = h_DF(ipoin)
          porosity = n_DF(ipoin)
       else
          height   = rho(ipoin)  
          porosity = (denss-dens)/(denss-densw)
       endif 
       dens2        = (1-porosity)*denss + porosity*densw*Sr  !  -- dens of upper layer (can be unsat)
       dens1        = (1-porosity)*denss + porosity*densw     !     lower layer (sat)   Sr=alphaS
       dens_bar     = (1-alphaS)*dens2 + alphaS*dens1         !     average dens
       dens_bar_eff = dens_bar - densw*alphaS                 !     MP 19 January 2019 
     
       pwp_total = gravz * (densw * height * alphaS + dens_bar_eff * u(ipoin) )
       pwp_total = pwp_total/1000.   ! KPa
     
   endif
                        
   if ( if_Out_Domain(ipoin).eq.1 )  CYCLE  ! ** mp 20th Jan 2019 avois water pts
   idest = 5*(iabs(itype(ipoin))-1)                                 !  = 0 for soil and 5 for water
   xx = x(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
!   inx0      = max(1,ipoigx-1)
!   iny0      = max(1,ipoigy-1)
!   inx1      = min(npoigx,ipoigx+2)
!   iny1      = min(npoigy,ipoigy+2)
   dinx      = hsml(ipoin)/deltxg 
   diny      = hsml(ipoin)/deltyg 
   if (ndimn.eq.1) diny = 0 
   inx0      = max(1,ipoigx-dinx)
   iny0      = max(1,ipoigy-diny)
   inx1      = min(npoigx,ipoigx + dinx+1)  !  min(npoigx,ipoigx+1) ! min(npoigx,ipoigx )  MP 1 May 19
   iny1      = min(npoigy,ipoigy + diny+1)  ! min(npoigy,ipoigy+1)  ! min(npoigy,ipoigy )
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoigx + inx
      distx2 = ( x(1,ipoin) - coorg(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x(2,ipoin) - coorg(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      dist   = dist/distg                                                         ! distance is normalized
      auxv(1+idest,ipoig) = auxv(1+idest,ipoig) + rho(ipoin)/dist                 ! h  s,w
      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + vx(1,ipoin)/dist                ! vx        
      if(ndimn.eq.2) auxv(3+idest,ipoig) = auxv(3+idest,ipoig) + vx(2,ipoin)/dist ! vy
!      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + hsml(ipoin)/dist               ! tmp hsml s and w
      if ( Gid_Mask_SW (4).eq.1 ) then
          if (ic_ws_interact.eq.2) then
              if (ipoin.le.npoin_s) then
                  auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/h_df(ipoin))/dist 
              endif
          else
              auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/rho(ipoin))/dist      ! pwp
          endif
      elseif ( Gid_Mask_SW (4).eq.-1 ) then
          auxv(4+idest,ipoig) = auxv(4+idest,ipoig) +  pwp_total/dist             ! changed mp Total press
      endif
      auxv(5+idest,ipoig) = auxv(5+idest,ipoig) + 1.00/dist                       ! dist  
      auxv(11,ipoig)      = auxv(11,ipoig) +  rho(ipoin)/dist
      auxv(12,ipoig)      = auxv(12,ipoig) + 1.00/dist                  ! dist. we are adding Hsoil and water 
      if (ic_hrelSat.eq.1.AND.ipoin.le.npoin_s) then
           auxv(13,ipoig)      = auxv(13,ipoig) + hrelSat_DF(ipoin)/dist
      elseif (ic_hrelSat.ne.1) then
           auxv(13,ipoig)      = auxv(13,ipoig) + hsml(ipoin)/dist
      endif
      
     ! --------   to Plot peclet number when a Stefan problem July 2021 MP

      if (ic_Stef_eros.eq.1.AND.Gid_Mask_SW (7).eq.1) then  
          hh = rho (ipoin)
          vv = vx(1,ipoin)*vx(1,ipoin)
          if (ndimn.eq.2) then
             vv = vv + vx(2,ipoin)*vx(2,ipoin)  
          endif
          vv = vv**0.5
          ee = 0.0
          if ( allocated (activ_Pwp_Stef_eros_Pt ) ) then
             ee = erosion_rate_V(ipoin)
          endif
          Te = 1.e6
          if (ee.gt.1.e-6) then
             Te = hh/ee
          endif
          Cvv = 4.0*const(13)/pi**2
          Tc  = hh*hh/Cvv
          Pe  = min(333.0, Tc/Te)
          auxv(14,ipoig)      = auxv(14,ipoig) + Pe/dist
      endif
   enddo
   enddo
ENDDO

DO ipoig = 1, npoig
   if (auxv(5,ipoig).ge.1.e-6) then
       auxv(1:4,ipoig) = auxv(1:4,ipoig)/auxv(5,ipoig)
       if (ic_hrelSat.eq.1) then
           auxv(13,ipoig)      = auxv(13,ipoig)/auxv(5,ipoig)
       endif
       if ( allocated (activ_Pwp_Stef_eros_Pt)) then  ! Peclet computation at nodes
           auxv(14,ipoig)      = auxv(14,ipoig)/auxv(5,ipoig)
       endif    
   endif
   if (auxv(10,ipoig).ge.1.e-6) then
       auxv(6:9,ipoig) = auxv(6:9,ipoig)/auxv(10,ipoig) 
   endif
   if (auxv(12,ipoig).ge.1.e-6) then
       auxv(11, ipoig) = auxv(11,ipoig)/auxv(12,ipoig)
   endif      

   if (auxv(1,ipoig).ge.1.e-6) ictrack(ipoig)=1      !   ** changed to 0.001 <- 0.1
enddo

do ipoig = 1, npoig                                                                                                  !   Saeidmt 26May2021   ** Add  porosity
    if (auxv(6,ipoig)/(auxv(1,ipoig)+auxv(6,ipoig)).GE.const(35)+0.001.AND.ic_hrelSat.eq.1) then !if (((auxv(13,ipoig)+(Sr*(1-auxv(13,ipoig))))*(auxv(1,ipoig)))+auxv(6,ipoig).le.0.01) CYCLE
        auxv(13,ipoig) = 1.0
        auxv(16,ipoig) = auxv(6,ipoig)/(auxv(1,ipoig)+auxv(6,ipoig))        ! Saeidmt 21July2021 
    elseif (ic_hrelSat.eq.1.AND.auxv(13,ipoig).LT.0.99) then
        auxv(16,ipoig) = const(35)           !   porosity limited by C35 ! Saeidmt 08Aug2021
        auxv(13,ipoig) = (1./(1.-Sr)) * (((1.-auxv(16,ipoig))*auxv(6,ipoig)/(auxv(16,ipoig)*auxv(1,ipoig))) - Sr)
    endif
enddo

auxv(15,:)=auxv(1,:)+auxv(6,:)      !   Saeidmt 26May2021  **Add hs(1) + hw(6)
auxv(11,:)=auxv(1,:)+auxv(6,:) ! Saeidmt 26May2021  **old auxv(11,:)=(auxv(1,:)+auxv(6,:) This is hs(1) + hw(6)

do ipoig = 1, npoig
    if (ic_hrelSat.eq.1) then
        if (1.0-((auxv(16,ipoig)*(1.0-Sr)*(1.0-auxv(13,ipoig)))).le.0.001) cycle    ! Saeidmt 21July2021
        auxv(11,ipoig)=(auxv(1,ipoig)+auxv(6,ipoig))/(1.0-((auxv(16,ipoig)*(1.0-Sr)*(1.0-auxv(13,ipoig)))))               ! Saeidmt 26May2021  hs+hw+ha **old auxv(11,:)=(auxv(1,:)+auxv(6,:) This is hs(1) + hw(6)
        if (auxv(1,ipoig).le.0.001) auxv(11,ipoig) = 0.0
    endif
end do

!   ************   changed 21st July 2014 : avoid white walls

DO ipoig  = 1, npoig
   IF (ictrack (ipoig).eq.1) THEN
       xx     = coorg(1,ipoig)
       yy     = coorg(2,ipoig)
       ipoigx = (xx-xming)/deltxg + 1
       ipoigy = (yy-yming)/deltyg + 1
       inx0   = max(1,ipoigx-1)
       iny0   = max(1,ipoigy-1)
       inx1   = min(npoigx,ipoigx+1)
       iny1   = min(npoigy,ipoigy+1)
       DO inx = inx0, inx1
       DO iny   = iny0, iny1
          idest = (iny-1)*npoigx + inx
          if (ictrack(idest).eq.0) then
              ictrack(idest) = 2
          endif
       ENDDO
       ENDDO
   ENDIF
ENDDO

! where(ictrack.eq.2) ictrack = 1

!   ************   changed 21st July 2014 : avoid white walls


!      ------   And now plot:
 

if(.not.allocated(etag)) then    
   allocate (etag(npoig))
   etag(:)=auxv(11,:) + topol(1,:)
endif
 
   
!      ------ 1:  h soil  -----------------------------

if ( Gid_Mask_SW (1).eq.1 ) then  
   write(gid_res,*) 'Result "height soil" "Height soil" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots h SOIL at quad nodes  ADDED
      !if ( ictrack(ipoig).eq.0) CYCLE                                  ! Saeidmt 12Nov2023   
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(1,ipoig)             ! max (0.1,auxv(1,ipoig))
   enddo
   write(gid_res,*) ' End Values ' 
endif 

!      ------  2 Displacements at particles   -------------

if ( Gid_Mask_SW (2).eq.1 ) then    
   write(gid_res,*) 'Result "dis" "Disp" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values  '  
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
      write(gid_res,*) i + 3*npoig, xdisp, ydisp, zdisp                      !   Saeidmt 26May2021  ** old write(gid_res,*) i + 2*npoig, xdisp, ydisp, zdisp
   ENDDO
   write(gid_res,*) ' End Values ' 
endif    
  

!      ------  Velocities  -----------------------------

if ( Gid_Mask_SW (3).eq.1 ) then

    write(gid_res,*) 'Result "vel" "veloc" ',time,' Vector OnNodes "where" '
    write(gid_res,*) ' Values '
    do ipoig = 1, npoig
       if ( auxv(1,ipoig).le.1.e-6) CYCLE
       write(gid_res,*) ipoig, auxv(2,ipoig), auxv(3,ipoig), ' 0.0 '
    enddo
    do i = 1,npoin 
       if ( if_Out_Domain(i).eq.1) CYCLE
       vxx = vx(1,i)
       if (ndimn.eq.1) then
           vyy = 0.0
       else
           vyy = vx(2,i)
       endif
       write(gid_res,*) i + 3*npoig, vxx, vyy, ' 0.0 '            !   Saeidmt 26May2021  ** old    write(gid_res,*) i + 2*npoig, vxx, vyy, ' 0.0 '
    enddo
    write(gid_res,*) ' End Values '    
 endif   
   
!      ------ 4 And PWPs  at the base -----------------------

if ( iabs(Gid_Mask_SW (4)).eq.1 ) then   
   if (icpwp.eq.1.or.icpwp.eq.3.or.icpwp.eq.2) then
      write(gid_res,*) 'Result "Pwp" "pwpress" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig 
         if ( auxv(1,ipoig).le.1.e-6) CYCLE
         write(gid_res,*) ipoig, '   0.0   0.0   ', auxv(4,ipoig)
      enddo
   !   do i = 1,npoin 
   !      if ( if_Out_Domain(i).eq.1) CYCLE
   !      write(gid_res,*) i + 2*npoig, '  0.0   0.0   ', u(i)
   !   enddo
      write(gid_res,*) ' End Values '
   endif   
endif    

!      ------ 5 ErodedTopo  -----------------------------

if ( Gid_Mask_SW (5).eq.1 ) then  
   write(gid_res,*) 'Result "erodedTopo" "-erodedTopo-" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      write(gid_res,*) ipoig, ' 0.   0.   ',   topol(14,ipoig)      
   enddo
   write(gid_res,*) ' End Values ' 

   write(gid_res,*) 'Result "AvailErosD" "ErosDepth" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      write(gid_res,*) ipoig, ' 0.   0.   ',   topol(15,ipoig)      
   enddo
   write(gid_res,*) ' End Values ' 
    
   write(gid_res,*) 'Result "hs over topo" "-hs over topo-" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      if (ic_basal_erosion.eq.1) then
         write(gid_res,*) ipoig, ' 0.   0.   ',   topol(1,ipoig)  - topol(14,ipoig)  
      else
         write(gid_res,*) ipoig, ' 0.   0.   ',   topol(1,ipoig) ! - topol(14,ipoig)  
      endif    
   enddo
   do ipoig = 1, npoig          !      ------  plots hsml at quad nodes  ADDED
      i = npoig + ipoig
      if (ic_basal_erosion.eq.1) then
         if (ic_ws_interact.eq.2) then
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(11,ipoig)   - topol(14,ipoig) 
         else 
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(1,ipoig)     - topol(14,ipoig)
         endif 
      else
         if (ic_ws_interact.eq.2) then
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(11,ipoig)   
         else 
             write(gid_res,*) i, ' 0.   0.   ',   topol(1,ipoig) + auxv(1,ipoig)     
         endif 
      endif            
   enddo                                   !      ------  plots TWO MESHES: S                                 

   write(gid_res,*) ' End Values ' 
   
endif   

!      ------ 6  Topo Z  -----------------------------

if ( Gid_Mask_SW (6).eq.1 ) then
   write(gid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig 
      if (allocated(hs_plus_Z_topo)) then 
          write(gid_res,*) ipoig, ' 0.   0.   ', hs_plus_Z_topo (ipoig)  
      else
          write(gid_res,*) ipoig, ' 0.   0.   ', topol(1,ipoig)
      endif
   enddo
   write(gid_res,*) ' End Values ' 
endif

!      ------7   h REL soil  -----------------------------  
!              This is an alternative evaluation and plotting of height

if ( Gid_Mask_SW (7).eq.999 ) then
   
   write(gid_res,*) 'Result "hs" "hs" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '

   dArea    = deltxg*deltyg
   if(.NOT.allocated(height_g) ) allocate ( height_g(npoig) )
   height_g = 0.0
   do ipoin = 1, npoin                                   
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      xx = x(1,ipoin)  
      yy = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,ipoin)
      ipoigx = (xx-xming)/deltxg  + 1
      ipoigy = (yy-yming)/deltyg  + 1
      if (ipoigx.eq.npoigx) then
         if (abs(xmaxg-xx).le.1.e-3) then
            ipoigx = ipoigx -1
         else
            write(*,*) ' point x outside topo mesh', xx, xmaxg
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
      xref   = xming + (ipoigx-1)*deltxg
      yref   = yming + (ipoigy-1)*deltyg
      xi     = (xx-xref)/deltxg
      eta    = (yy-yref)/deltyg
      sh1    = (1.-xi) * (1.-eta)
      sh2    =     xi  * (1.-eta)
      sh3    =     xi  * eta
      sh4    = (1.-xi) * eta
      height_g(n1) = height_g(n1) + mass(ipoin) 
      height_g(n2) = height_g(n2) + mass(ipoin) 
      height_g(n3) = height_g(n3) + mass(ipoin) 
      height_g(n4) = height_g(n4) + mass(ipoin) 
   enddo
   height_g = height_g/(4.*dArea)
   do ipoig = 1, npoig                                     
      if ( height_g(ipoig).le.0.005) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', height_g(ipoig)
   enddo
   deallocate ( height_g )

   write(gid_res,*) ' End Values ' 

elseif ( Gid_Mask_SW (7).eq.1 ) then   ! Peclet
   if ( allocated (activ_Pwp_Stef_eros_Pt) ) then 
     write(gid_res,*) 'Result "Peclet nr" "Peclet" ',time,' Vector OnNodes "where" '
     write(gid_res,*) ' Values '
     do ipoig = 1, npoig                                    !      ------  plots Peclet  ADDED
        if ( ictrack(ipoig).eq.0) CYCLE   
        write(gid_res,*) ipoig, ' 0.   0.   ', auxv(14,ipoig)          ! max (0.1,auxv(1,ipoig))
     enddo
     write(gid_res,*) ' End Values '  
  endif     
endif

!      ------ 8   h water -----------------------------

if ( Gid_Mask_SW (8).eq.1 ) then
   write(gid_res,*) 'Result "height water" "Height water" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  plots h WATER at quad nodes
!      if ( auxv(6,ipoig).le.0.3) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(6,ipoig)
   enddo
   write(gid_res,*) ' End Values ' 
endif   

!      ------ 9  Eta (wave height)  -------------------- 

if ( Gid_Mask_SW (9).eq.1 ) then

   if (ic_ws_Interact.eq.2) then
      write(gid_res,*) 'Result " n " " porosity" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
         do ipoig = 1, npoig                                     
            if ( (auxv(1,ipoig)+auxv(6,ipoig)).le.0.01) CYCLE
            xi = auxv(6,ipoig)/(auxv(6,ipoig)+((auxv(13,ipoig)+Sr*(1.0-auxv(13,ipoig)))*auxv(1,ipoig)))                  ! Saeidmt 21July2021
            if (ic_hrelSat.eq.1.AND.auxv(13,ipoig).LT.0.99) xi = const(35)           !   porosity limited by C35         ! Saeidmt 08Aug2021 
            write(gid_res,*) ipoig, ' 0.   0.   ', xi                          !      ------  plots h WATER at quad nodes ONLY Z>SWL=0
         enddo 
   else 
      write(gid_res,*) 'Result "eta" " sobreelev" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
      do ipoig = 1, npoig
         if ( auxv(6,ipoig).le.0.3) CYCLE
         write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig) + topol(1,ipoig)-etag(ipoig)!   ------  CHANGED
      enddo
   endif
   
   if (ic_ws_Interact.eq.2) then                                                          !      Saeidmt 07July2021 Showing the real porosity
      write(gid_res,*) 'Result " n_sw " " n_sw" ',time,' Vector OnNodes "where" '
      write(gid_res,*) ' Values '
         do ipoig = 1, npoig                                     
            if ( (auxv(1,ipoig)+auxv(6,ipoig)).le.0.01) CYCLE
            xi = auxv(6,ipoig)/(auxv(6,ipoig)+((auxv(13,ipoig)+Sr*(1.0-auxv(13,ipoig)))*auxv(1,ipoig)))                 ! Saeidmt 21July2021
            write(gid_res,*) ipoig, ' 0.   0.   ', xi                          !      ------  plots h WATER at quad nodes ONLY Z>SWL=0
         enddo 
   endif    
        
   write(gid_res,*) ' End Values '
   
endif   

!      ------ 10  Joint h soil + water  ----------------

if ( Gid_Mask_SW (10).eq.1 ) then 
   write(gid_res,*) 'Result "jointh" "JOINTh" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !   Saeidmt 26May2021      ------  plots THREE MESHES: S+W+a
      write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)
   enddo
   do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S+W
       i = npoig + ipoig                                  !   ** Add  Saeidmt 26May2021 
!      if ( auxv(11,ipoig).le.0.3.or.auxv(6,ipoig).le.0.3) CYCLE  *** 23 July 14
      write(gid_res,*) i, ' 0.   0.   ', auxv(15,ipoig) !   Saeidmt 26May2021   **old   write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)
   enddo
   do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S
      i = 2*npoig + ipoig                                 !   ** Add  Saeidmt 26May2021 
!      if ( auxv(1,ipoig).le.0.3) CYCLE *** 23 July 14
      write(gid_res,*) i, ' 0.   0.   ', auxv(1,ipoig)
   enddo
   write(gid_res,*) ' End Values ' 
endif   
    
!      ------11   hsml  OR hrelSat_DF  -----------------------------


if ( Gid_Mask_SW (11).eq.1 ) then  
      if (ic_hrelSat.eq.1) then
          write(gid_res,*) 'Result "hrelSat" "-hrelSat-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)     
          enddo
          do ipoig = 1, npoig                                    !      ------  plots TWO MESHES: S
             i = npoig + ipoig
             write(gid_res,*)      i, ' 0.   0.   ', auxv(13,ipoig)*auxv(11,ipoig)
          enddo
          write(gid_res,*) ' End Values ' 
      elseif (ic_hrelSat.ne.1) then
          write(gid_res,*) 'Result "hsml" "-HSML-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(13,ipoig)     
          enddo
          write(gid_res,*) ' End Values ' 
      endif
      
      if (ic_hrelSat.eq.1) then
          write(gid_res,*) 'Result "hrSat_only" "-hrSat-" ',time,' Vector OnNodes "where" '
          write(gid_res,*) ' Values '
          do ipoig = 1, npoig                                    !      ------  plots hsml at quad nodes  ADDED
              write(gid_res,*) ipoig, ' 0.   0.   ', auxv(13,ipoig)     
          enddo
          write(gid_res,*) ' End Values ' 
      endif
endif  

if ( Gid_Mask_SW (11).eq.-1 ) then  
   write(gid_res,*) 'Result "F or I" "-F or I -" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoig                                    !      ------  
      vxx = auxv(2,ipoig) * auxv(2,ipoig)
      vyy = 0.0
      if (ndimn.eq.2) vyy = auxv(3,ipoig) * auxv(3,ipoig)
      vvv = (vxx + vyy)**0.5
      hhh = auxv(1,ipoig)
      cgra = const(1)
      if (hhh.gt.const(10)) then
         Froude = vvv/((hhh*gravz)**0.5)
         if (const(5).eq.8) then
             if (hhh.lt.const(6)) hhh=const(6)
             Froude = Froude * (const(6)/hhh)
         endif
      else
         Froude = 0.0
      endif
      write(gid_res,*) ipoig, ' 0.   0.   ', Froude     
   enddo
   write(gid_res,*) ' End Values ' 
endif 


!      ------  Pwp (1d full)  -------------------- 

if ( Gid_Mask_SW (12).eq.1.and.ndimn.eq.1.and.nUaux.gt.1 ) then

   write(gid_res,*) 'Result "Disp pwp" "disp pore press" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   ipoig = 0
   do ipoin = 1, npoin_s                                         !      ------  plots disp at particles

      !IF ( if_Out_Domain(ipoin).eq.1) then
      !        ipoig = ipoig + nUaux
      !endif

      IF ( if_Out_Domain(ipoin).eq.1) CYCLE

      disp(:) = x(:,ipoin) - x00(:,ipoin)
      xdisp   = disp(1)
      if (ndimn.eq.1) then
          ydisp = 0.0
      else
         ydisp = disp(2)
      endif
      xx     =  x(1,ipoin)
      yy     = (yming+ymaxg)/2
      if (ndimn.eq.2) yy = x(2,ipoin) 
         
      call Get_Z_Topo (ipoin, xx,yy,Zp)
   
      zdisp  = Zp - Z00(ipoin)
      
      hhh    = rho(ipoin)   
      hh0    = rho00(ipoin)                         !  ch mp 20th Jan 2019 pwp plot
      if (allocated (h_df)) then
            hhh    = h_df (ipoin)
            hh0    = rho00(ipoin) + rho00(ipoin+npoin_s)
      endif
      dhh     = (hhh-hh0)/(nUaux-1)                 !******* 
      
      do iaux  = 1, nUaux
         ipoig = ipoig + 1
         idest = 3*npoig + npoin + ipoig                                  !   Saeidmt 26May2021  ** old    idest = 2*npoig + npoin + ipoig 
         write(gid_res,*) idest, xdisp, ydisp, zdisp + (iaux-1)*dhh  
      enddo
   enddo
   write(gid_res,*) ' End Values ' 

   write(gid_res,*) 'Result "Full  pwp" "full pwp" ',time,' Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   ipoig = 0
   do ipoin = 1, npoin_s  

      IF ( if_Out_Domain(ipoin).eq.1) CYCLE

      do iaux  = 1, nUaux
         ipoig = ipoig + 1
         idest = 3*npoig + npoin + ipoig                !   Saeidmt 26May2021  ** old    idest = 2*npoig + npoin + ipoig 
         hhh    = rho(ipoin)                            !  ch mp 02 july 17 for df's
         if (allocated (h_df)) hhh    = h_df (ipoin)
         write(gid_res,*) idest,   Uaux (iaux,ipoin)/hhh 
      enddo
   enddo
   write(gid_res,*) ' End Values ' 
   
   
endif   

!      ------  Output for file 3  GIS nodes

xcg = 0.0; ycg = 0.0; zcg = 0.0; Mtotal_out = 0.0

npoing = 0
do ipoin = 1, npoin
   if (if_out_domain(ipoin).eq.0)  then  !   *** only for pts insid
      npoing = npoing +1
      xx = x(1,ipoin)
      yy = (ymaxg+yming)/2.
      if (ndimn.eq.2) yy = x(2,ipoin)
      call Get_Z_Topo ( ipoin, xx, yy, Zp)
      xcg = xcg + mass(ipoin)*xx
      ycg = ycg + mass(ipoin)*yy
      zcg = zcg + mass(ipoin)*Zp
      Mtotal_out = Mtotal_out + mass(ipoin)
   endif
enddo   
xcg = xcg / Mtotal_out
ycg = ycg / Mtotal_out
zcg = zcg / Mtotal_out
write(res_file3,9998) time, xcg, ycg, zcg, Mtotal_out   ! mp 1st Dec 2017
9998 format ( ' t xcg ycg zcg Mtot  ', 5(2x,g10.4))


! *******  changed 3 march

!      ------  Output for file 4  GIS nodes

write(res_file4,*) ' ipoin   x   y   h  vx  vy  pwp  r w  t = ', time, ' M ' , Mtotal_out
do ipoin=1,npoin
   if (ic_crush.eq.1) then
      write(res_file4,9999) ipoin, (x(i,ipoin),i=1,ndimn), rho(ipoin), (vx(i,ipoin),i=1,ndimn), u(ipoin), r_crush(ipoin), w_crush(ipoin)
   else
      write(res_file4,9999) ipoin, (x(i,ipoin),i=1,ndimn), rho(ipoin), (vx(i,ipoin),i=1,ndimn), u(ipoin) 
   endif
enddo 
9999 format (2x, i6, 2x,8(2x,f10.4))

!      ------  Output for file 5  GIS erosion

if (gid_mask_SW(5).eq.1) then
    write(res_file5,*) ' ipoig   x   y   z  eroded  ', '  time = ', time
    do ipoig=1,npoig
       write(res_file5,9999) ipoig,  coorg(1,ipoig), coorg(2,ipoig), topol(14,ipoig)
    enddo
endif      

!      ------  Output for file 6   pwp profiles   **** MP 11-Jan 2019

if (gid_mask_SW(12).eq.1.AND.nUAux.GT.1) then
    write(res_file6,*) ' ..... time = ',  time, ' pwp bottom -> pwp top relatives '
    do ipoin=1,npoin_s
       if (ic_ws_interact.eq.2) then
  !         if (ic_hrelSat.eq.1) then                ! Saeidmt 31May2021   **Add hrelSat_DF(ipoin)
  !             hhh = h_DF(ipoin)*hrelSat_DF(ipoin)
  !         else
               hhh = h_DF(ipoin)
   !        endif
       else
           hhh = rho(ipoin)
       endif
       write (res_file6,9997) ipoin, (( Uaux (iaux,ipoin)/hhh ),iaux=1,nUaux) 
    enddo
endif  
9997 format ( 2x, i6, 10(2x,f10.4))    

! *******  
 

deallocate ( disp )
deallocate ( auxv )

 
End Subroutine OutputRes_WIR_SAEID_bug


!-------------------------------------------------------------------

       Subroutine OutputRes_coarse_SW 
       
!-------------------------------------------------------------------


!      ------  OUT for waves in reservoirs --------

implicit none

integer  (ink)  i, ipoin, ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer  (ink)  inx0, inx1, inx, iny0, iny1, iny
integer  (ink)  idest
real     (irk)  xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp 
real     (irk)  xdisp, ydisp, zdisp, vxx, vyy
real     (irk), allocatable::disp(:)
real     (irk), allocatable::auxv(:,:)
real     (irk), allocatable, SAVE:: etag   (:)                  !      ------  Added
integer  (ink), allocatable, SAVE:: ictrack(:)                  !      ------  Added
real     (irk)  distx2, disty2, dist, distg 
real     (irk)  dArea
real(irk), allocatable:: height_g(:)

if(.NOT.allocated(disp) ) allocate ( disp(ndimn) ) 
if(.NOT.allocated(auxv) ) allocate ( auxv(12,npoi_plt))

Call Check_Out_Domain_SW


!      ------   First of all we get values at topo grid. Add for each cell value*(1/dist) and divide by sum (1/dist)

auxv    = 0.0

if(.not.allocated(ictrack)) then                         
   allocate (ictrack(npoi_plt))
   ictrack = 0
endif
                                                         
distg   = (deltx_plt*deltx_plt + delty_plt*delty_plt)**0.5

DO ipoin = 1,npoin                                       
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   idest = 5*(iabs(itype(ipoin))-1)                     !     = 0 for soil and 5 for water
   xx = x(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)  
   ipoigx    = (xx-xming)/deltx_plt + 1
   ipoigy    = (yy-yming)/delty_plt + 1
   inx0      = max(1,ipoigx-1)
   iny0      = max(1,ipoigy-1)
   inx1      = min(npoix_plt,ipoigx+2)
   iny1      = min(npoiy_plt,ipoigy+2)
!   inx0      = max(1,ipoigx)
!   iny0      = max(1,ipoigy)
!   inx1      = min(npoix_plt,ipoigx+1)
!   iny1      = min(npoiy_plt,ipoigy+1)
   DO inx    = inx0,inx1
   DO iny    = iny0,iny1
      ipoig  = (iny-1)*npoix_plt + inx
      distx2 = ( x(1,ipoin) - coor_plt(1,ipoig) )
      distx2 = distx2*distx2
      disty2 = 0.0; if (ndimn.eq.2) disty2 = ( x(2,ipoin) - coor_plt(2,ipoig) )
      disty2 = disty2*disty2
      dist   = max( 1.e-6, (distx2+disty2)**0.5 )
      dist   = dist/distg
      auxv(1+idest,ipoig) = auxv(1+idest,ipoig) + rho(ipoin)/dist               ! h        soil  and water
      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + hsml(ipoin)/dist              ! temporal hsml s and w
      auxv(4+idest,ipoig) = auxv(4+idest,ipoig) + (u(ipoin)/rho(ipoin))/dist    ! pwp
      auxv(5+idest,ipoig) = auxv(5+idest,ipoig) + 1.00/dist                     ! dist for soil and water
!      auxv(2+idest,ipoig) = auxv(2+idest,ipoig) + vx(1,ipoin)/dist             ! vx       soil and water
!      if(ndimn.eq.2) auxv(3+idest,ipoig) = auxv(3+idest,ipoig) + vx(2,ipoin)/dist
!      auxv(11,ipoig)      = auxv(11,ipoig) + + rho(ipoin)/dist
!      auxv(12,ipoig)      = auxv(12,ipoig) + 1.00/dist                         ! dist. we are adding Hsoil and water  
   enddo
   enddo
ENDDO



DO ipoig = 1, npoi_plt
   if (auxv(5,ipoig).ge.1.e-6) then
       auxv(1:4,ipoig) = auxv(1:4,ipoig)/auxv(5,ipoig) 
   endif
   if (auxv(10,ipoig).ge.1.e-6) then
       auxv(6:9,ipoig) = auxv(6:9,ipoig)/auxv(10,ipoig) 
   endif
   if (auxv(1,ipoig).ge.1.e-6) ictrack(ipoig)=1                 !    ADDED and changed to 0.001 from 0.1
enddo


auxv(11,:)=auxv(1,:)+auxv(6,:)                          ! This is hs(1) + hw(6)


!      ------   And now plot:
 

if(.not.allocated(etag)) then   !   ------  CHANGED
   allocate (etag(npoi_plt))
   etag(:)=auxv(11,:) + topol(1,:)
endif

!      ------  Erosion  -----------------------------


If (const(4).gt.1.e-6) then
   write(gid_res,*) 'Result "erosion" " erosion " ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoi_plt                                    !      ------  plots erosioned soil
      if ( topol(14,ipoig).le.1.e-6) CYCLE
      write(gid_res,*) ipoig, ' 0.   0.   ', -topol(14,ipoig) 
   enddo
   write(gid_res,*) ' End Values ' 
endif   

!      ------  Eta (wave height)  -------------------- 


write(gid_res,*) 'Result "eta" " sobreelev" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoi_plt                                    !      ------  plots h WATER at quad nodes ONLY Z>SWL=0
   if ( auxv(6,ipoig).le.0.3) CYCLE
!   if (( auxv(6,ipoig).le.1.e-6).or.(topol(1,ipoig).gt.0.0).or.(auxv(6,ipoig).lt.50) ) CYCLE
   write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig) + topol(1,ipoig)-etag(ipoig)!   ------  CHANGED
enddo
write(gid_res,*) ' End Values ' 

!      ------  h water -----------------------------

write(gid_res,*) 'Result "height water" "Height water" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoi_plt                                    !      ------  plots h WATER at quad nodes
   if ( auxv(6,ipoig).le.0.3) CYCLE
   write(gid_res,*) ipoig, ' 0.   0.   ', auxv(6,ipoig)
enddo
write(gid_res,*) ' End Values ' 

!      ------  h soil  -----------------------------

write(gid_res,*) 'Result "height soil" "Height soil" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoi_plt                                    !      ------  plots h SOIL at quad nodes  ADDED
   if ( ictrack(ipoig).ne.1) CYCLE
   write(gid_res,*) ipoig, ' 0.   0.   ', auxv(1,ipoig)             ! max (0.1,auxv(1,ipoig))
enddo
write(gid_res,*) ' End Values ' 

!      ------  h REL soil  -----------------------------  
!              This is an alternative evaluation and plotting of height

write(gid_res,*) 'Result "hs" "hs" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '

dArea    = deltx_plt*delty_plt
if(.NOT.allocated(height_g) ) allocate ( height_g(npoi_plt) )
height_g = 0.0
do ipoin = 1, npoin                                      
   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   xx = x(1,ipoin)  
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)
   ipoigx = (xx-xming)/deltx_plt  + 1
   ipoigy = (yy-yming)/delty_plt  + 1
   if (ipoigx.eq.npoix_plt) then
      if (abs(xmaxg-xx).le.1.e-3) then
         ipoigx = ipoigx -1
      else
         write(*,*) ' point x outside topo mesh', xx, xmaxg
      endif
   endif   
   if (ipoigy.eq.npoiy_plt) then
      if (abs(ymaxg-yy).le.1.e-3) then
         ipoigy = ipoigy -1
      else
         write(*,*) ' point y outside topo mesh', xx, xmaxg
      endif
   endif
   ipoig  = (ipoigy-1)*npoix_plt + ipoigx
   n1     = ipoig
   n2     = ipoig + 1
   n3     = ipoig + 1 + npoix_plt
   n4     = n3 -1 
   xref   = xming + (ipoigx-1)*deltx_plt
   yref   = yming + (ipoigy-1)*delty_plt
   xi     = (xx-xref)/deltx_plt
   eta    = (yy-yref)/delty_plt
   sh1    = (1.-xi) * (1.-eta)
   sh2    =     xi  * (1.-eta)
   sh3    =     xi  * eta
   sh4    = (1.-xi) * eta
   height_g(n1) = height_g(n1) + mass(ipoin) 
   height_g(n2) = height_g(n2) + mass(ipoin) 
   height_g(n3) = height_g(n3) + mass(ipoin) 
   height_g(n4) = height_g(n4) + mass(ipoin) 
enddo
height_g = height_g/(4.*dArea)
do ipoig = 1, npoi_plt                                     
   if ( height_g(ipoig).le.0.005) CYCLE
   write(gid_res,*) ipoig, ' 0.   0.   ', height_g(ipoig)
enddo
deallocate ( height_g )

write(gid_res,*) ' End Values '
!      ------  Joint h soil + water  ----------------
 
write(gid_res,*) 'Result "jointh" "JOINTh" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoi_plt                                    !      ------  plots TWO MESHES: S+W
   if ( auxv(11,ipoig).le.0.3.or.auxv(6,ipoig).le.0.3) CYCLE
   write(gid_res,*) ipoig, ' 0.   0.   ', auxv(11,ipoig)
enddo
do ipoig = 1, npoi_plt                                    !      ------  plots TWO MESHES: S
   i = npoi_plt + ipoig
   if ( auxv(1,ipoig).le.0.3) CYCLE
   write(gid_res,*) i, ' 0.   0.   ', auxv(1,ipoig)
enddo
write(gid_res,*) ' End Values ' 

!      ------  Topo Z  -----------------------------

write(gid_res,*) 'Result " Z topo " "   Z   " ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values '
do ipoig = 1, npoi_plt 
   write(gid_res,*) ipoig, ' 0.   0.   ', coor_plt(3,ipoig)
enddo
write(gid_res,*) ' End Values ' 

!      ------  Displacements at particles   -------------

write(gid_res,*) 'Result "dis" "Disp" ',time,' Vector OnNodes "where" '
write(gid_res,*) ' Values  '  
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
      write(gid_res,*) i + 2*npoi_plt, xdisp, ydisp, zdisp
ENDDO
write(gid_res,*) ' End Values ' 

!      ------  Velocities  -----------------------------

!write(gid_res,*) 'Result "vel" "veloc" ',time,' Vector OnNodes "where" '
!write(gid_res,*) ' Values '
!do ipoig = 1, npoi_plt
!   if ( auxv(1,ipoig).le.1.e-6) CYCLE
!   write(gid_res,*) ipoig, auxv(2,ipoig), auxv(3,ipoig), ' 0.0 '
!enddo
!do i = 1,npoin 
!   if ( if_Out_Domain(i).eq.1) CYCLE
!   vxx = vx(1,i)
!   if (ndimn.eq.1) then
!       vyy = 0.0
!   else
!       vyy = vx(2,i)
!   endif
!   write(gid_res,*) i + 2*npoi_plt, vxx, vyy, ' 0.0 ' 
!enddo
!write(gid_res,*) ' End Values '

!      ------  And PWPs  -----------------------------

if (icpwp.eq.1) then
   write(gid_res,*) 'Result "Pwp" "pwpress" ',time,' Vector OnNodes "where" '
   write(gid_res,*) ' Values '
   do ipoig = 1, npoi_plt
      if ( auxv(1,ipoig).le.1.e-6) CYCLE
      write(gid_res,*) ipoig, '  0.0   0.0   ', auxv(4,ipoig)
   enddo
!   do i = 1,npoin 
!      if ( if_Out_Domain(i).eq.1) CYCLE
!      write(gid_res,*) i + 2*npoi_plt, '  0.0   0.0  ', u(i)
!   enddo
   write(gid_res,*) ' End Values '
endif   

deallocate ( disp )
deallocate ( auxv )

 
End Subroutine OutputRes_coarse_SW


!-------------------------------------------------------------------

       Subroutine Moni_SW 
            
!-------------------------------------------------------------------


!      ------  Control output for a given set of check points  --------

implicit none

integer(ink) ipoin   
real   (irk) vtemp, vtempx, vtempy

integer(ink) ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4
integer(ink) inx0, iny0, inx1, iny1, inx, iny

real   (irk) xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4 
real   (irk) distg, distx2, disty2, dist
real   (irk), allocatable:: auxv(:,:), h_chk(:)


1000 format  (5(2x,g12.5))

!      ------  Obtain vtemp 

vtemp = 0.0
do ipoin = 1, npoin
   if (If_Out_Domain(ipoin).eq.1) cycle
   vtempx = vx(1,ipoin)
   vtempy = 0.0
   if (ndimn.eq.2) vtempy = vx(2,ipoin)
   vtemp  = vtemp + (vtempx*vtempx + vtempy*vtempy)**0.5
enddo
vtemp = vtemp/npoin
if (ic_TRIGGER.NE.0) write(*,*) '   vtemp =  ', vtemp   !  MP July 2022

!      ------  Obtain variables at topo grid nodes 

IF (nchk_pts.GT.0) THEN       ! MP 14th June 2022
   distg   = (deltxg*deltxg + deltyg*deltyg)**0.5
   if(.NOT.allocated(auxv) ) allocate ( auxv(2,npoig) )
   auxv = 0.0

   DO ipoin = 1,npoin                                       
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
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
         auxv(1,ipoig) = auxv(1,ipoig) + rho(ipoin)/dist           ! h        soil   
         auxv(2,ipoig) = auxv(2,ipoig) + 1.00/dist         
      enddo
      enddo
   ENDDO

   DO ipoig = 1, npoig
      if (auxv(2,ipoig).ge.1.e-6) then
          auxv(1,ipoig) = auxv(1,ipoig)/auxv(2,ipoig) 
      endif
   enddo

!      ------  Now we know h at grip nodes, interpolate at chk pts

   if(.NOT.allocated(h_chk) ) allocate ( h_chk(nchk_pts) )

   DO ipoin = 1, nchk_pts

      xx  = coor_chk_pts (1,ipoin)
      yy  = coor_chk_pts (2,ipoin)
      ipoigx = (xx-xming)/deltxg  + 1
      ipoigy = (yy-yming)/deltyg  + 1
   
      if (ipoigx.eq.npoigx) then
         if (abs(xmaxg-xx).le.1.e-3) then
            ipoigx = ipoigx -1
         else
            write(*,*) ' point x outside topo mesh', xx, xmaxg
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
   
      h_chk(ipoin)  = sh1*auxv(1,n1) + sh2*auxv(1,n2) + sh3*auxv(1,n3) + sh4*auxv(1,n4)
   
   
   ENDDO   

   write(res_file2,3000) time, (h_chk(ipoin),ipoin=1, nchk_pts)
   deallocate ( auxv,h_chk ) 
   3000 format (8(f10.4,2x))
   
ENDIF   


End Subroutine Moni_SW

!-------------------------------------------------------------------

       Subroutine output_points_SW 
            
!-------------------------------------------------------------------


!      ------  Output x h v  h_exact DBK dry   v_exact DBK_dry  --------

implicit none

integer(ink) ipoin, id 
real   (irk) hh,tt,xt
real   (irk) vv, ee, Te, Tc, Pe, Cvv
real   (irk) pi, deg2rads, theta, phi, h0, grav, c0, m, x_l, x_r
real   (irk) tanPhi

1000 format  (7(2x,g12.5))
1001 format  (i6, 2x, 7(2x,g12.5))

tt =time
if (time.lt.1.e-6) tt=1.e-6
pi       = 4.*atan(1.0)
deg2rads = pi/180
                                      !!!!  theta    = 0.0                        
                                      !!!!! 30.0*deg2rads
                                      !!!!!  phi      = 0.0    ! 25.0*deg2rads
if ( ic_Stef_eros.eq.1) then
   write(res_file1,*) ' ipoin  hh                   vv         ee            Te            Tc         Peclet at base  time =' , time
   do ipoin = 1, npoin_s
      hh = rho (ipoin)
      vv = vx(1,ipoin)*vx(1,ipoin)
      if (ndimn.eq.2) then
         vv = vv + vx(2,ipoin)*vx(2,ipoin)  
      endif
      vv = vv**0.5
      ee = 0.0
      if (allocated (activ_Pwp_Stef_eros_Pt) ) then
         ee = erosion_rate_V(ipoin)
      endif
      Te = 1.e6
      if (ee.gt.1.e-6) then
         Te = hh/ee
      endif
      Cvv = 4.0*const(13)/pi**2
      Tc  = hh*hh/Cvv
      Pe  = Tc/Te
      write(res_file1,1001) ipoin, hh, vv, ee, Te, Tc, Pe
   enddo 
elseif (ic_hrelSat.eq.1) then
   if (nfric.eq.202) then
      write(res_file1,*) ' ipoin   x    y    h   n_DF  Voellmy  mass  hsml time =' , time  
      do ipoin = 1, npoin_s
         if (ndimn.eq.2) then 
             write(res_file1,1001) ipoin, x(1,ipoin),x(2,ipoin), & 
                                   rho(ipoin),n_DF(ipoin), ns_dep_props(1,ipoin) , &
                                   hsml(ipoin), mass(ipoin) 
         endif                           
      enddo
   else
      write(res_file1,*) ' ipoin  hs  hw  h_DF  HrelSat  n  vx    time =' , time  
      do ipoin = 1, npoin_s 
         if (ic_DF.eq.1) then                    !  MP March 2023 FS can have icherelsat = 1
            write(res_file1,1001) ipoin, rho(ipoin),rho(ipoin+npoin_s) , & 
                               h_DF(ipoin), hrelSat_DF(ipoin), n_DF(ipoin), vx(1,ipoin) 
         elseif (ic_FS.EQ.1) then
            write(res_file1,1001) ipoin, rho(ipoin),rho(ipoin) , & 
                            rho(ipoin), hrelSat_DF(ipoin), n_DF(ipoin), vx(1,ipoin) 
         endif                               
     enddo
   endif

else

   h0       = 10.
   tanPhi   = const(9)
   grav     = const(1)
   if (gravz.le.1.e-6) gravz = const(1)
   c0       = sqrt (gravz*h0)            !!!!! sqrt(grav*h0*cos(Theta))
   m        = gravz*tanPhi - gravx       !!!!! grav*(cos(Theta)*tanPhi-sin(Theta))
   x_r      = 2*c0*tt - 0.5*m*tt*tt
   x_l      =  -c0*tt - 0.5*m*tt*tt

   write(res_file1,*) ' x   h_sph h_analyt at time=', time
   do ipoin = 1, npoin_s   !10                               ! Saeidmt 19Oct2021    ** do ipoin = 1, npoin
      if (x(1,ipoin).le.x_l) then
          hh = h0
      elseif ( x(1,ipoin).le.x_r) then
          hh = (1./(9.*gravz))*(2*c0-(x(1,ipoin)/tt)-0.5*m*tt)**2.
      else
          hh = 0.
      endif
      !   xt = 2*(10.0*const(1))**0.5 - x(1,ipoin)/tt  
      !  hh = xt**2./(9.0*const(1))
      if (hh.lt.0.0)  hh = 0.0
      if (hh.gt. h0)  hh = h0 
      if (allocated (h_DF) ) then   
         write(res_file1,1000) x(1,ipoin), h_DF(ipoin), hh, vx(1,ipoin), hsml(ipoin)   ! Saeidmt 19Oct2021    ** 
      else
         write(res_file1,1000) x(1,ipoin), rho(ipoin), hh, vx(1,ipoin), hsml(ipoin)    ! MP Nov 2021
      endif
   enddo

endif

End Subroutine output_points_SW


! ----------------------------------------------------------------------   

      subroutine  Read_vn0Bcs

! ----------------------------------------------------------------------   

!      This subroutine provides following output:
!               nvn0                    number of Abs Bcs
!               lvn0        (nvn0)      list of vn=0 nodes
!               x_vn0 (ndimn,nvn0)      coordinates x and y(if ndimn=2) 
!               vn_vn0(ndimn,nvn0)      unit normal vector, pointing outwards
!      All those variables defined as allocatable in SW and dimensioned here

!      We assume that nvn0 will be smaller than npoig, and dimension temp arrays so



implicit none

integer(ink) Narc, iarc, ndiv_arc,  idiv_arc,  Nseg, iseg, idimn, iivn0, jvn0 
integer(ink) ipoin, idest
integer(ink) ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4

real   (irk) R_arc, s_arc, ds_arc, Zext_arc,d_arc, d_arc2, h_arc, h_arc2
integer(ink) Type_arc 
real   (irk) vn_mod, theta0, thetaf, theta, dtheta
real   (irk) pi
real   (irk) xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp
real   (irk) tmpaux1, tmpaux2
real   (irk) fact_hsml_vn0

real   (irk), allocatable:: x_temp(:,:), vn_temp(:,:), Z_temp(:), hsml_temp(:) 
integer(ink), allocatable:: type_temp(:)                                       ! MP 28 Nov 2022
real   (irk), allocatable:: x0_arc(:),   xf_arc(:),    xm_arc(:), xc_arc(:)
real   (irk), allocatable:: vt_arc(:),   vn_arc(:),    r0_arc(:), rf_arc(:)

!      ------  

nvn0 = 0
pi   = 4.*atan(1.0)

read (dat_file,1005) text         !   --- chgd MP 28 Nov 2022
write(chk_file,1005) text         ! max nr of vn0 pts
read (dat_file,*)    max_nvn0     ! instead of npoig (very large)
write(chk_file,*)    max_nvn0

IF (ndimn.eq.2) THEN
   if(.NOT.allocated(x_temp) )  then
     allocate (x_temp(ndimn,max_nvn0), vn_temp(ndimn,max_nvn0))
     allocate (type_temp (max_nvn0),   Z_temp(max_nvn0), hsml_temp(max_nvn0))      ! MP 2022 type 1 not perm 2 perm
     allocate (x0_arc(ndimn),     xf_arc(ndimn),   xm_arc(ndimn), xc_arc(ndimn) )
     allocate (vt_arc(ndimn),     vn_arc(ndimn),   r0_arc(ndimn), rf_arc(ndimn))
   endif
   x_temp = 0; vn_temp=0.; Z_temp=0.; hsml_temp=0

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    Narc, Nseg    
   write(chk_file,*)    Narc, Nseg

   DO iarc = 1, Narc            !      ------  Begin loop around circular arcs

      read (dat_file,1005) text
      write(chk_file,1005) text
      read (dat_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn) !  --  Initial and final points
      write(chk_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)
      read (dat_file,*)    R_arc, ds_arc, Zext_arc, Type_arc                            !  --  Radius, sep between points, Depth
      write(chk_file,*)    R_arc, ds_arc, Zext_arc, Type_arc                            !  Type-arc 1 not permeable 2 permeable
                                                                                        !  MP 28 Nov 2022
      xm_arc    = (xf_arc + x0_arc)/2.  
      vt_arc    = (xf_arc - x0_arc)/2.
      d_arc2    = ( vt_arc(1)*vt_arc(1) + vt_arc(2)*vt_arc(2) ) 
      h_arc     = sqrt(R_arc*R_arc - d_arc2) 
      vn_arc(1) = - vt_arc(2) 
      vn_arc(2) =   vt_arc(1) 
      vn_mod    = sqrt(vn_arc(1)**2. + vn_arc(2)**2.) 
      xc_arc    = xm_arc + ( h_arc/vn_mod)*vn_arc
      R0_arc    = x0_arc - xc_arc
      Rf_arc    = xf_arc - xc_arc
      Theta0    = atan2 (r0_arc(2),r0_arc(1))    
      Thetaf    = atan2 (rf_arc(2),rf_arc(1))
      if (Thetaf.le.Theta0) Thetaf = Thetaf + 2.*pi
      s_arc     = (Thetaf-Theta0)*R_arc
      ndiv_arc  = (s_arc/ds_arc)
      ds_arc    = s_arc/ndiv_arc
      dtheta    = ds_arc/R_arc
      
      !      -------  Now we generate all points, centered at segments dtheta
     
      Do idiv_arc = 1, ndiv_arc
         idest    = nvn0 + idiv_arc
         if (nvn0+idiv_arc.GT.max_nvn0) then      !  MP 28 Nov 2020
            write(*,*) ' in routine readvn0 max nr of pts exceeded...modify code to 2*npoig '
            write(*,*) npoig, nvn0+idiv_arc
            stop
         endif            
         theta                       = Theta0 + dtheta*(idiv_arc - 1) + dtheta/2.
         x_temp (1, idest)   = R_arc*cos(theta)
         x_temp (2, idest)   = R_arc*sin(theta)
         vn_temp(1, idest)   = cos(theta)
         vn_temp(2, idest)   = sin(theta)
         Z_temp    (idest)   = Zext_arc
         hsml_temp (idest)   = 2.*ds_arc 
         type_temp (idest)   = Type_arc      !  MP 28 Nov 2020
      Enddo
      
      nvn0  = nvn0 + ndiv_arc
      idest = nvn0      
      
   ENDDO                        !      ------  End of loop around circular arcs
   
   DO iseg = 1, Nseg            !      ------  Begin loop around Rectilinear segments (a-clockwise)
   
      read (dat_file,1005) text
      write(chk_file,1005) text
      read (dat_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)  !  --  Initial and final points
      write(chk_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)    
      read (dat_file,*)     ds_arc, fact_hsml_vn0, Type_arc                              !  --  Sep between points, Depth
      write(chk_file,*)     ds_arc, fact_hsml_vn0, Type_arc                              !  MP 28 Nov 2022 
      vt_arc    = (xf_arc - x0_arc)
      d_arc2    = ( vt_arc(1)*vt_arc(1) + vt_arc(2)*vt_arc(2) )**0.5
      ndiv_arc  = d_arc2/ds_arc
      ds_arc    = d_arc2/ndiv_arc
      vt_arc    = vt_arc / d_arc2
      vn_arc(1) =   vt_arc(2)
      vn_arc(2) =  -vt_arc(1)
      
      DO idiv_arc = 1, ndiv_arc
         idest    = nvn0 + idiv_arc
         if (nvn0+idiv_arc.GT.max_nvn0) then      !  MP 28 Nov 2020
            write(*,*) ' in routine readvn0 max nr of pts exceeded...modify code to 2*npoig '
            write(*,*) npoig, nvn0+idiv_arc
            stop
         endif            
         x_temp (:,idest) = x0_arc(:) + (ds_arc/2.)*vt_arc(:) + (idiv_arc-1)*ds_arc*vt_arc(:)
         vn_temp(:,idest) = vn_arc(:)
         Z_temp   (idest) = Zext_arc
         hsml_temp(idest) = fact_hsml_vn0*ds_arc 
         type_temp (idest)   = Type_arc      !  MP 28 Nov 2020
      Enddo
      
      nvn0  = nvn0 + ndiv_arc
      idest = nvn0                      ! just to see it through debugger  
      
   ENDDO 
   
   if (.NOT.allocated (lvn0)) then
      allocate (lvn0(nvn0),  x_vn0(ndimn,nvn0), vn_vn0(ndimn,nvn0), hsml_vn0(nvn0), type_vn0(nvn0) )   !  MP 28 nov 2022  
   endif   
   
   idest = 0
   DO iivn0 = 1, nvn0
      xx    = x_temp (1,iivn0)          !      ------  Zp is bottom Z, Z-temp is surface. Then, h = Z-temp - Zp
      yy    = ( ymaxg + yming )/2.
      if (ndimn.eq.2) yy = x_temp (2,iivn0)
      idest = idest +1
      hsml_vn0(idest)  = hsml_temp(iivn0)
      x_vn0 (:,idest)  = x_temp (:,iivn0)
      vn_vn0(:,idest)  = vn_temp(:,iivn0)
      type_vn0(idest)  = type_temp(iivn0)      !  MP 28 nov 2022  
   enddo  
   nvn0 = idest                                         !     ------  we only add nodes below SWL    
      
   deallocate (x_temp, vn_temp, Z_temp, hsml_temp, type_temp )
   deallocate (x0_arc,  xf_arc, xm_arc, xc_arc)
   deallocate (vt_arc,  vn_arc, r0_arc, rf_arc)
 
   
ELSEIF (ndimn.eq.1) THEN                !      ------  This is for 1D cases
 
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    nvn0
   write(chk_file,*)    nvn0   
  
   if (.NOT.allocated (lvn0)) then
      allocate (lvn0(nvn0),  x_vn0(ndimn,nvn0), vn_vn0(ndimn,nvn0), hsml_vn0(nvn0), type_vn0(nvn0) )   ! MP 28 Nov 2022 
   endif   
   
   read (dat_file,1005) text
   write(chk_file,1005) text
   DO jvn0 = 1,nvn0
      read (dat_file,*) lvn0(jvn0), x_vn0(1,jvn0), vn_vn0(1,jvn0), hsml_vn0(jvn0), type_vn0(jvn0)
      write(chk_file,*) lvn0(jvn0), x_vn0(1,jvn0), vn_vn0(1,jvn0), hsml_vn0(jvn0), type_vn0(jvn0)
   enddo
   
ENDIF

!      ------   

1005     format(a60) 


END SUBROUTINE Read_vn0Bcs


! ----------------------------------------------------------------------   

      subroutine  Read_AbsBcs

! ----------------------------------------------------------------------   

!      This subroutine provides following output:
!               nabs                    number of Abs Bcs
!               labs        (nabs)      list of absorbing nodes
!               Zext_abs    (nabs)      Depth of water otside
!               x_abs (ndimn,nabs)      coordinates x and y(if ndimn=2) 
!               vn_abs(ndimn,nabs)      unit normal vector, pointing outwards
!      All those variables defined as allocatable in SW and dimensioned here

!      We assume that nabs will be smaller than npoig, and dimension temp arrays so



implicit none

integer(ink) Narc, iarc, ndiv_arc,  idiv_arc,  Nseg, iseg, idimn, iiabs, jabs 
integer(ink) ipoin, idest
integer(ink) ipoigx, ipoigy, ieleg, ipoig, n1, n2, n3, n4

real   (irk) R_arc, s_arc, ds_arc, Zext_arc,d_arc, d_arc2, h_arc, h_arc2  
real   (irk) vn_mod, theta0, thetaf, theta, dtheta
real   (irk) pi
real   (irk) xx, yy, xref, yref, xi, eta, sh1, sh2, sh3, sh4, Zp
real   (irk) tmpaux1, tmpaux2

real   (irk), allocatable:: x_temp(:,:), vn_temp(:,:), Z_temp(:), hsml_temp(:) 
real   (irk), allocatable:: x0_arc(:),   xf_arc(:),    xm_arc(:), xc_arc(:)
real   (irk), allocatable:: vt_arc(:),   vn_arc(:),    r0_arc(:), rf_arc(:)

!      ------  

nabs = 0
pi   = 4.*atan(1.0)

IF (ndimn.eq.2) THEN
   if(.NOT.allocated(x_temp) ) then
     allocate (x_temp(ndimn,npoig), vn_temp(ndimn,npoig),   Z_temp(npoig), hsml_temp(npoig))
     allocate (x0_arc(ndimn),     xf_arc(ndimn),   xm_arc(ndimn), xc_arc(ndimn) )
     allocate (vt_arc(ndimn),     vn_arc(ndimn),   r0_arc(ndimn), rf_arc(ndimn))
   endif
   x_temp = 0; vn_temp=0.; Z_temp=0.; hsml_temp=0

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    Narc, Nseg    
   write(chk_file,*)    Narc, Nseg

   DO iarc = 1, Narc            !      ------  Begin loop around circular arcs

      read (dat_file,1005) text
      write(chk_file,1005) text
      read (dat_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn) !  --  Initial and final points
      write(chk_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)
      read (dat_file,*)    R_arc, ds_arc, Zext_arc                                  !  --  Radius, sep between points, Depth
      write(chk_file,*)    R_arc, ds_arc, Zext_arc  
   
      xm_arc    = (xf_arc + x0_arc)/2.  
      vt_arc    = (xf_arc - x0_arc)/2.
      d_arc2    = ( vt_arc(1)*vt_arc(1) + vt_arc(2)*vt_arc(2) ) 
      h_arc     = sqrt(R_arc*R_arc - d_arc2) 
      vn_arc(1) = - vt_arc(2) 
      vn_arc(2) =   vt_arc(1) 
      vn_mod    = sqrt(vn_arc(1)**2. + vn_arc(2)**2.) 
      xc_arc    = xm_arc + ( h_arc/vn_mod)*vn_arc
      R0_arc    = x0_arc - xc_arc
      Rf_arc    = xf_arc - xc_arc
      Theta0    = atan2 (r0_arc(2),r0_arc(1))    
      Thetaf    = atan2 (rf_arc(2),rf_arc(1))
      if (Thetaf.le.Theta0) Thetaf = Thetaf + 2.*pi
      s_arc     = (Thetaf-Theta0)*R_arc
      ndiv_arc  = (s_arc/ds_arc)
      ds_arc    = s_arc/ndiv_arc
      dtheta    = ds_arc/R_arc
      
      !      -------  Now we generate all points, centered at segments dtheta
     
      Do idiv_arc = 1, ndiv_arc
         idest    = nabs + idiv_arc
         if (nabs+idiv_arc.GT.npoig) then
            write(*,*) ' in routine readabs max nr of pts exceeded...modify code to 2*npoig '
            write(*,*) npoig, nabs+idiv_arc
            stop
         endif            
         theta                       = Theta0 + dtheta*(idiv_arc - 1) + dtheta/2.
         x_temp (1, idest)   = R_arc*cos(theta)
         x_temp (2, idest)   = R_arc*sin(theta)
         vn_temp(1, idest)   = cos(theta)
         vn_temp(2, idest)   = sin(theta)
         Z_temp    (idest)   = Zext_arc
         hsml_temp (idest) = 2.*ds_arc 
      Enddo
      
      nabs  = nabs + ndiv_arc
      idest = nabs      
      
   ENDDO                        !      ------  End of loop around circular arcs
   
   DO iseg = 1, Nseg            !      ------  Begin loop around Rectilinear segments (a-clockwise)
   
      read (dat_file,1005) text
      write(chk_file,1005) text
      read (dat_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)  !  --  Initial and final points
      write(chk_file,*)    (x0_arc(idimn),idimn=1,ndimn), (xf_arc(idimn),idimn=1,ndimn)    
      read (dat_file,*)     ds_arc, Zext_arc                                            !  --  Sep between points, Depth
      write(chk_file,*)     ds_arc, Zext_arc  
      vt_arc    = (xf_arc - x0_arc)
      d_arc2    = ( vt_arc(1)*vt_arc(1) + vt_arc(2)*vt_arc(2) )**0.5
      ndiv_arc  = d_arc2/ds_arc
      ds_arc    = d_arc2/ndiv_arc
      vt_arc    = vt_arc / d_arc2
      vn_arc(1) =   vt_arc(2)
      vn_arc(2) =  -vt_arc(1)
      
      DO idiv_arc = 1, ndiv_arc
         idest    = nabs + idiv_arc
         if (nabs+idiv_arc.GT.npoig) then
            write(*,*) ' in routine readabs max nr of pts exceeded...modify code to 2*npoig '
            write(*,*) npoig, nabs+idiv_arc
            stop
         endif            
         x_temp (:,idest) = x0_arc(:) + (ds_arc/2.)*vt_arc(:) + (idiv_arc-1)*ds_arc*vt_arc(:)
         vn_temp(:,idest) = vn_arc(:)
         Z_temp   (idest) = Zext_arc
         hsml_temp(idest) = 2.*ds_arc 
      Enddo
      
      nabs  = nabs + ndiv_arc
      idest = nabs                      ! just to see it through debugger  
      
   ENDDO 

   if(.NOT.allocated(labs) ) then  
     allocate (labs(nabs), Zext_abs(nabs), x_abs(ndimn,nabs), vn_abs(ndimn,nabs), hsml_abs(nabs) )
   endif
   idest = 0
   DO iiabs = 1, nabs
      xx    = x_temp (1,iiabs)          !      ------  Zp is bottom Z, Z-temp is surface. Then, h = Z-temp - Zp
      yy    = ( ymaxg + yming )/2.
      if (ndimn.eq.2) yy = x_temp (2,iiabs)
      call Get_Z_Topo (iiabs, xx,yy,Zp) 
      if ( Zp.gt.Z_temp(iiabs)) CYCLE
      idest = idest +1
      labs    (idest) = npoin + nvirt + iiabs   !       WRONG  npoin not defined yet Go to SetUpGlobal routine. Not used anyway
      Zext_abs(idest) = Z_temp(iiabs) - Zp 
      hsml_abs(idest) = hsml_temp(iiabs)
      x_abs (:,idest) = x_temp (:,iiabs)
      vn_abs(:,idest) = vn_temp(:,iiabs)
      
      tmpaux1 = Zext_abs(idest)
      tmpaux2 = labs(idest)
      
   enddo  
   nabs = idest                                         !     ------  we only add nodes below SWL    
      
   deallocate (x_temp, vn_temp, Z_temp, hsml_temp )
   deallocate (x0_arc,  xf_arc, xm_arc, xc_arc)
   deallocate (vt_arc,  vn_arc, r0_arc, rf_arc)
 
   
ELSEIF (ndimn.eq.1) THEN                !      ------  This is for 1D cases
 
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    nabs
   write(chk_file,*)    nabs   
   if(.NOT.allocated(labs) )  then
     allocate (labs(nabs), Zext_abs(nabs), x_abs(ndimn,nabs), vn_abs(ndimn,nabs), hsml_abs(nabs) )   
   endif   
   read (dat_file,1005) text
   write(chk_file,1005) text
   
   DO jabs = 1,nabs
      read (dat_file,*) labs(jabs), x_abs(1,jabs), vn_abs(1,jabs), Zext_abs(jabs), hsml_abs(jabs)
      write(chk_file,*) labs(jabs), x_abs(1,jabs), vn_abs(1,jabs), Zext_abs(jabs), hsml_abs(jabs)
   enddo
   
ENDIF

!      ------   

1005     format(a60) 

Call Check_Abs_BCs_SW (nabs, labs, Zext_abs, hsml_abs, x_abs, vn_abs )

END SUBROUTINE Read_AbsBcs



! ----------------------------------------------------------------------   

      subroutine  Abs_BCs_SW 

! ----------------------------------------------------------------------   

!      This subroutine  applies BCs of absorbing type 
!

implicit none

integer(ink) i, j, jabs, idimn, ipoin

real   (irk) distij, anx, any, Zext, vn, vtx, vty, h, cc, c0
real   (irk) r10, r20, hnew, vnew, cgrav, vmod, anmod
real   (irk)  xx, yy, zp

real   (irk), allocatable:: Aux_Bcs (:,:)

if (nabs.ge.1) then   

   Call Check_Abs_BCs_SW (nabs, labs, Zext_abs, hsml_abs, x_abs, vn_abs )

   if(.NOT.allocated(Aux_Bcs) ) allocate (Aux_Bcs(5,npoin))
   Aux_Bcs = 0.0
   cgrav   = const(1)

   current=> last
   Do while (associated(current))
      i = current%pair_i 
      j = current%pair_j
      if (current%Pint_type.eq.-11) then                !  -11 means ABS BCs
          jabs   = j-npoin-nvirt
          distij = 0.0
          do idimn=1,ndimn
             distij = distij + ( x(idimn,i)-x(idimn,j) )**2.
          enddo
             distij = distij**0.5
             if (distij.lt.1.e-6) distij=1.e-6
             aux_bcs(1,i) = aux_bcs(1,i) + 1.
             aux_bcs(2,i) = aux_bcs(2,i) + (1./distij)
             aux_bcs(3,i) = aux_bcs(3,i) + Zext_abs(jabs)/distij
             aux_bcs(4,i) = aux_bcs(4,i) + vn_abs(1,jabs)/distij
             if (ndimn.eq.2) aux_bcs(5,i) = aux_bcs(5,i) + vn_abs(2,jabs)/distij
      endif   
      current=>current%next
   enddo

   DO ipoin = 1,npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      if (aux_bcs(1,ipoin).gt.1.e-6) then
          anx    = aux_bcs(4,ipoin)/aux_bcs(2,ipoin)
          any    = 0.0
          if (ndimn.eq.2)  any   = aux_bcs(5,ipoin)/aux_bcs(2,ipoin)
          anmod  = sqrt(anx*anx + any*any)
          anx    = anx/anmod
          any    = any/anmod
          Zext   = aux_bcs(3,ipoin)/aux_bcs(2,ipoin)      
!         -----------------------
           xx     = x(1,ipoin)                  ! ***** changed. Apply them at point
           yy     = (yming + ymaxg)/2.
           if (ndimn.eq.2) yy = x(2,ipoin)
           call Get_Z_Topo (ipoin, xx,yy,Zp) 
           Zext   = 0.0 - Zp                                ! *******  SWL is ZERO
                  if (Zext.le.0.1) CYCLE
!         ---------------------          
          vmod   = (vx(1,ipoin)*vx(1,ipoin))
          if (ndimn.eq.2) vmod = vmod + (vx(2,ipoin)*vx(2,ipoin))
          vmod   = vmod**0.5
          vn     = vx(1,ipoin) * anx
          if (ndimn.eq.2) vn = vn + vx(2,ipoin) * any
          vtx    = vx(1,ipoin) - anx*vn
          vty  = 0.0
          if (ndimn.eq.2) vty = vx(2,ipoin) - any*vn  
          cc     = (cgrav*rho(ipoin))**0.5
          c0     = (cgrav*Zext)**0.5
          R10    = 2.*cc + vn
          R20    = 2.*c0
          hnew   = (1./(16.*cgrav))*(r10 + r20)**2.         
          vnew   = 0.5*(R10-R20)                            
          vx(1,ipoin) = vtx + vnew*anx
          if (ndimn.eq.2) vx(2,ipoin)= vty + vnew*any
          rho (ipoin) = hnew
          x(1,ipoin)   = x0(1,ipoin) + dt_sph*vx(1,ipoin)          
          if (ndimn.eq.2) x(2,ipoin)   = x0(2,ipoin) + dt_sph*vx(2,ipoin)  
       endif
   enddo

   deallocate (Aux_bcs)

endif

END Subroutine Abs_BCs_SW 


! ----------------------------------------------------------------------   

      subroutine  Abs_BCs_SW_TRY

! ----------------------------------------------------------------------   

!      This subroutine  applies BCs of absorbing type 
!

implicit none

integer(ink) i, j, jabs, idimn, ipoin

real   (irk) distij, anx, any, Zext, vn, vtx, vty, h, cc, c0
real   (irk) r10, r20, hnew, vnew, cgrav, vmod, anmod
real   (irk)  xx, yy, zp, rr

real   (irk), allocatable:: Aux_Bcs (:,:)

if (nabs.ge.1) then   

   Call Check_Abs_BCs_SW (nabs, labs, Zext_abs, hsml_abs, x_abs, vn_abs )

   cgrav   = const(1)

   DO ipoin = 1,npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
          xx     = x(1,ipoin)                   
          yy     = (yming + ymaxg)/2.
          if (ndimn.eq.2) yy = x(2,ipoin)
          rr     = ( xx*xx + yy*yy)**0.5
      if  ( abs(rr-6000).le.2.*hsml(ipoin)) then
          anx    = xx
          any    = 0.0
          if (ndimn.eq.2)  any   = yy
          anmod  = sqrt(anx*anx + any*any)
          anx    = anx/anmod
          any    = any/anmod
          Zext   = 800. 
          vmod   = (vx(1,ipoin)*vx(1,ipoin))
          if (ndimn.eq.2) vmod = vmod + (vx(2,ipoin)*vx(2,ipoin))
          vmod   = vmod**0.5
          vn     = vx(1,ipoin) * anx
          if (ndimn.eq.2) vn = vn + vx(2,ipoin) * any
          vtx    = vx(1,ipoin) - anx*vn
          vty  = 0.0
          if (ndimn.eq.2) vty = vx(2,ipoin) - any*vn  
          cc     = (cgrav*rho(ipoin))**0.5
          c0     = (cgrav*Zext)**0.5
          R10    = 2.*cc + vn
          R20    = 2.*c0
          hnew   = (1./(16.*cgrav))*(r10 + r20)**2.         
          vnew   = 0.5*(R10-R20)                            
          vx(1,ipoin) = vtx + vnew*anx
          if (ndimn.eq.2) vx(2,ipoin)= vty + vnew*any
          rho (ipoin) = hnew
          x(1,ipoin)   = x0(1,ipoin) + dt_sph*vx(1,ipoin)          
          if (ndimn.eq.2) x(2,ipoin)   = x0(2,ipoin) + dt_sph*vx(2,ipoin)  
       endif
   enddo


endif

END Subroutine Abs_BCs_SW_TRY


! ----------------------------------------------------------------------   

      subroutine  Vn0_BCs_SW (dtvn0) 

! ----------------------------------------------------------------------   

!      This subroutine  applies BCs of Vn=0 type 
!

implicit none

integer(ink) i, j, jvn0, idimn, ipoin

real   (irk) distij, anx, any,  vn, vtx, vty, h, cc, c0
real   (irk) r10, r20, hnew, vnew,  vmod, anmod
real   (irk)  xx, yy, zp
real   (irk)  dtvn0

real   (irk), allocatable:: Aux_Bcs (:,:)

if (nvn0.ge.1) then   
   
   if(.NOT.allocated(Aux_Bcs)) allocate (Aux_Bcs(5,npoin))
   Aux_Bcs = 0.0

   current=> last
   Do while (associated(current))
      i = current%pair_i 
      j = current%pair_j
      IF (current%Pint_type.eq.-21) then         
          jvn0   = j-npoin-nvirt-nabs      !  bug (a3) 21 nov 2010
          IF (itype(i).eq.2.AND.type_vn0(jvn0).eq.2) then
             idimn = idimn                 ! do nothing perm wall  MP 28  nov 2022
             xx = 0.0
          ELSE
             distij = 0.0
             do idimn=1,ndimn
                distij = distij + ( x(idimn,i)-x(idimn,j) )**2.
             enddo
             distij = distij**0.5
             if (distij.lt.1.e-6) distij=1.e-6
             aux_bcs(1,i) = aux_bcs(1,i) + 1.
             aux_bcs(2,i) = aux_bcs(2,i) + (1./distij)
             aux_bcs(4,i) = aux_bcs(4,i) + vn_vn0(1,jvn0)/distij
             if (ndimn.eq.2) aux_bcs(5,i) = aux_bcs(5,i) + vn_vn0(2,jvn0)/distij
          ENDIF
      ENDIF   
      current=>current%next
   enddo

   if_correct = 0
   
   DO ipoin = 1,npoin
      if ( if_Out_Domain(ipoin).eq.1) CYCLE
      if (aux_bcs(1,ipoin).gt.1.e-6) then
      
          if_correct (ipoin) = 1
      
          anx    = aux_bcs(4,ipoin)/aux_bcs(2,ipoin)
          any    = 0
          if (ndimn.eq.2) any    = aux_bcs(5,ipoin)/aux_bcs(2,ipoin)
          anmod  = sqrt(anx*anx + any*any)
          anx    = anx/anmod
          any    = any/anmod 
          
          vn     =                  vx (1,ipoin) * anx
          if (ndimn.eq.2) vn = vn + vx (2,ipoin) * any
          
          if (vn.lt.0.0) then 
              CYCLE   ! *** MP 27th June BACK 7march
          endif

         ! if (itype(ipoin).eq.2.AND.ic_vn0.eq.2) then  ! deactivated MP 28 Nov 2022
         !     CYCLE   ! *** MP 20th december. Skip wall-water interactions, leaving pics to pass trough
         ! endif
          
          vtx    = vx(1,ipoin) - anx*vn         
          vx(1,ipoin) = vtx  
          x(1,ipoin)  = x0(1,ipoin) + dtvn0*vx(1,ipoin)     !  ***** CHANGED 6 march 16: dt*f1rk
          
          if (ndimn.eq.2) then
              vty    = vx(2,ipoin) - any*vn
              vx(2,ipoin) = vty               
              x(2,ipoin)  = x (2,ipoin) + dtvn0*vx(2,ipoin)  !  ***** CHANGED 6 march 16: commented
        !      x(2,ipoin)  = x0(2,ipoin) + dtvn0*vx(2,ipoin)  !  ***** CHANGED 6 march 16: commented
          endif
       endif
   enddo

   deallocate (Aux_bcs)

endif

END Subroutine vn0_BCs_SW 


! ----------------------------------------------------------------------   

      subroutine  Check_Abs_BCs_SW (nabs, labs, Zext_abs, hsml_abs, x_abs, vn_abs )

! ----------------------------------------------------------------------   

!      This subroutine  checks abs arrays
!

implicit none
    
integer(ink) iabs, nabs, labs(nabs)
real   (irk) Zext_abs(nabs), hsml_abs(nabs), x_abs(ndimn,nabs), vn_abs(ndimn,nabs)
real   (irk) z, xx, yy, hh

do iabs=1,nabs
   z   = zext_abs(iabs)
   xx  = x_abs(1,iabs)
   if (ndimn.eq.2)  yy  = x_abs(2,iabs)
   hh  = hsml_abs(iabs)
enddo


END Subroutine Check_Abs_BCs_SW


! ----------------------------------------------------------------------   

      subroutine  Check_Out_Domain_SW  

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

END  SUBROUTINE Check_Out_Domain_SW


! ----------------------------------------------------------------------   

      subroutine  change_to_w_SW  

! ----------------------------------------------------------------------   

!   Subroutine to check wheter a given point is outside the computational domain  
 
implicit none

integer(ink) ipoin, idimn
integer(ink) ipoigx, ipoigy, ipoig, inx0, inx1, inx, iny0,iny1, iny 

real   (irk) xx, yy, distx2, disty2, dist
real   (irk), allocatable:: auxv(:,:)

if(.NOT.allocated(auxv) ) allocate (auxv(2,npoig))
auxv = 0.0

DO ipoin=1,npoin

   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   if ( itype(ipoin).ne.1) CYCLE
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
      auxv(1,ipoig) = auxv(1,ipoig) + rho(ipoin)/dist            
      auxv(2,ipoig) = auxv(2,ipoig) + 1.0/dist 
   enddo
   enddo
   If_Out_Domain(ipoin) = 1
   rho (ipoin)          = 0
   
ENDDO      

Do ipoig =1,npoig
   if (auxv(2,ipoig).gt.1.e-6) then
      topol(1,ipoig) = topol(1,ipoig) +  auxv(1,ipoig)/auxv(2,ipoig)
   endif
enddo

Call Update_topol

deallocate (auxv)


END  SUBROUTINE change_to_w_SW

! ----------------------------------------------------------------------   

      subroutine Read_histogram_SW

! ---------------------------------------------------------------------- 

!     This subroutine is used to read histogram

implicit none  
 
integer(ink)  idimn, ip_gate, il_gate, ipoin, ip_his, i_hist 
real   (irk), allocatable:: x0_gate(:,:), xf_gate(:,:), tv_gate(:,:) 
real   (irk) cgra05
real   (irk) dx_gate

real   (irk) xx, yy     !  debugger
integer(ink) nn

pi      = 4.*atan(1.0)

read (dat_file,1005) text
read (dat_file,*)    type_of_histogram
write(chk_file,1005) text
write(chk_file,*) 
write(chk_file,*)  '   Type of histogram: 1 for h (m) ,  2 for q (m3/m.s) '
write(chk_file,*)  type_of_histogram
write(chk_file,*)

read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*) mp_hist, mp_gate, ml_gate    
write(chk_file,*) mp_hist, mp_gate, ml_gate 

if(.NOT.allocated(h_his) )  then
  allocate ( h_his   (n_hist,mp_hist),         v_his(n_hist,mp_hist), t_his(n_hist,mp_hist) )
  allocate ( x_gate  (n_hist,ndimn,mp_gate), nv_gate(n_hist,ndimn), tv_gate(n_hist,ndimn) )
  allocate ( x0_gate (n_hist,ndimn),         xf_gate(n_hist,ndimn) )
  allocate ( xl_gate (n_hist) )
  allocate ( np_gate (n_hist), nl_gate   (n_hist), np_hist(n_hist) )
  allocate ( np_track(n_hist), ndof_track(n_hist,ml_gate) )
  allocate ( np_bakG (n_hist) )
endif

np_track = 0   ; ndof_track = 0
np_gate  = 0   ; nl_gate    = 0   ; np_hist  = 0 
h_his    = 0.0 ; v_his      = 0.0 ; t_his    = 0.0
x_gate   = 0.0 ; nv_gate    = 0.0 ; tv_gate  = 0.0
x0_gate  = 0.0 ; xf_gate    = 0.0 ; xl_gate  = 0.0
np_bakG  = 0.0

!      ------  We read 3 sets of np_his points: time, h and v (directed along normal to gate)

DO i_hist = 1, n_hist


   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)    np_hist(i_hist), np_gate(i_hist), nl_gate(i_hist)
   write(chk_file,*)    np_hist(i_hist), np_gate(i_hist), nl_gate(i_hist)

   if ( nl_gate(i_hist).ne.ml_gate) then
      write(*,*) ' we are reading nl_gate fo r histogram ', i_hist
          write(*,*) ' and the value read ', nl_gate(i_hist), ' is different from ml_gate'
          write(*,*) ' ml_gate = ', ml_gate
          write(*,*) ' please re run code after modifying file.dat '
          write(*,*) ' All values of nl_gate should be equal to ml_gate '
          Pause 
          STOP
   endif

   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) ( t_his(i_hist,ip_his),ip_his = 1,np_hist(i_hist)) 
   write(chk_file,*) ( t_his(i_hist,ip_his),ip_his = 1,np_hist(i_hist))
 
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) ( h_his(i_hist,ip_his),ip_his = 1,np_hist(i_hist)) 
   write(chk_file,*) ( h_his(i_hist,ip_his),ip_his = 1,np_hist(i_hist))
 
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*) ( v_his(i_hist,ip_his),ip_his = 1,np_hist(i_hist)) 
   write(chk_file,*) ( v_his(i_hist,ip_his),ip_his = 1,np_hist(i_hist))

   !      ------  Gate information
   
   read (dat_file,1005) text
   write(chk_file,1005) text
   read (dat_file,*)  (x0_gate(i_hist,idimn), idimn=1,ndimn),(xf_gate(i_hist,idimn), idimn=1,ndimn)  
   write(chk_file,*)  (x0_gate(i_hist,idimn), idimn=1,ndimn),(xf_gate(i_hist,idimn), idimn=1,ndimn)  

   xl_gate(i_hist) = 0.0
   do idimn = 1, ndimn
      xl_gate(i_hist) = xl_gate(i_hist) + ( xf_gate(i_hist,idimn)-x0_gate(i_hist,idimn) )**2.
      tv_gate(i_hist,idimn) = (xf_gate(i_hist,idimn)-x0_gate(i_hist,idimn))
   enddo
   xl_gate(i_hist) = sqrt(xl_gate(i_hist))

   If (np_gate(i_hist).ne.1) THEN                  !      -   This is 2D case
      tv_gate(i_hist,:) = tv_gate(i_hist,:)/(np_gate(i_hist) - 1)       !      -   incremental vector along gate from pt to next
      nv_gate(i_hist,1) = - (xf_gate(i_hist,2) -x0_gate(i_hist,2))/xl_gate(i_hist)
      nv_gate(i_hist,2) =   (xf_gate(i_hist,1) -x0_gate(i_hist,1))/xl_gate(i_hist)
      DO ip_gate = 1, np_gate(i_hist)
         do idimn =1,ndimn
            x_gate(i_hist,idimn,ip_gate) = x0_gate(i_hist,idimn) + (ip_gate-1)*tv_gate(i_hist,idimn)
         enddo
      enddo
   else                                 !      -  This is the 1D case. Input of normal using XF
      x_gate (i_hist,:,1) = x0_gate(i_hist,:)
      nv_gate(i_hist,:)   = xf_gate(i_hist,:)
   endif  
   
ENDDO

write(chk_file,*) ' give factor for smoothing length. Total is 2*facthsml'      !      ------  used at injection time
read (dat_file,1005) text
write(chk_file,1005) text
read (dat_file,*)  hsml_gfact
write(chk_file,*)  hsml_gfact   

!      ------  Generate points

npoin_s = 0
DO i_hist = 1, n_hist
   npoin_s = npoin_s + np_gate(i_hist)*nl_gate(i_hist)
ENDDO

if(.NOT.allocated(x_s)) then
  allocate ( x_s(ndimn,npoin_s),  rho_s (npoin_s), mass_s (npoin_s), hsml_s (npoin_s) )
  allocate (vx_s(ndimn,npoin_s),  u_s   (npoin_s), itype_s(npoin_s) )
endif
                                                                        !  deleted p_s
ipoin = 0

DO il_gate = 1, ml_gate             !  LIMITATION OF CODE: all gates same nl_gate
   DO i_hist = 1, n_hist
      DO ip_gate =1, np_gate(i_hist)
         ipoin = ipoin +1
         if (ipoin.gt.npoin_s) then
             write(*,*) ' npoin is larger than mpoin', npoin, npoin_s
             write(*,*) ' program is stopped '
             write(chk_file,*) ' npoin is larger than mpoin', npoin, npoin_s
             write(chk_file,*) ' program is stopped '
             stop
         endif
         x_s   (:,ipoin) = x_gate(i_hist,:,ip_gate)
         itype_s (ipoin) = 1
      enddo
   enddo
ENDDO   

!      ------  all other vars set to zero because they will be read later on

rho_s  = 0.
mass_s = 0.
hsml_s = 0.
vx_s   = 0.
u_s    = 0.

!      ------  We add Back points in 4 lines at gates

if (ndimn.eq.2) then
   ipoin         = 0
   np_bakG_Total = 0
   DO i_hist = 1, n_hist
      np_bakG (i_hist) =  4*np_gate(i_hist)
      nn = np_bakG (i_hist)
      np_bakG_Total    = np_bakG_Total + np_bakG (i_hist)
      dx_gate          = xl_gate(i_hist)/(np_gate(i_hist)-1)
      DO il_gate = 1, 4
      DO ip_gate = 1, np_gate(i_hist)
         ipoin = ipoin +1
          x_s  (:,ipoin) = x_gate(i_hist,:,ip_gate) - ((2*il_gate-1)/2) * dx_gate * nv_gate(i_hist,:)  ! il-1, then il
          xx = x_s(1,ipoin); yy = x_s(1,ipoin) ! Debug only
          rho_s  (ipoin) = 0.99     
          mass_s (ipoin) = 0.00  ! ** changed. They do not contribute 
          hsml_s (ipoin) = 0.99
          u_s    (ipoin) = 0.99
      ENDDO
      ENDDO
   ENDDO
Endif   
  
npoin_in_model = np_bakG_Total
nn =  np_bakG_Total                            ! Debug only

deallocate ( x0_gate, xf_gate, tv_gate )

  
1005     format(a60)  

END subroutine Read_histogram_SW

!-------------------------------------------------------------------

       Subroutine Clean_SW 
            
!-------------------------------------------------------------------


!      ------  Clean internal allocated variables of module SW  --------

implicit none

if (Ic_abs.eq.1) then
   deallocate ( labs, Zext_abs, vn_abs, x_abs, hsml_abs)
endif

if (ic_histogram.eq.1) then 
    deallocate ( h_his    , v_his   , t_his)
    deallocate ( x_gate   , nv_gate)        
    deallocate ( np_gate  , nl_gate , np_hist)
    deallocate ( np_track , ndof_track)
    deallocate ( xl_gate  , np_bakG )
endif 


End Subroutine Clean_SW 

!-------------------------------------------------------------------

       Subroutine Inject_histogram_SW_hh  (time)  
            
!-------------------------------------------------------------------

implicit none
 
integer(ink) nrow_gate, npg2, iaux, itrack, i_hist 
integer(ink) icheck, ip_his, itotv, itotv1, itotv2, idimn, itotw, irow, itotvG

real   (irk) time, t_ref, xi, t0_his, t1_his, hh, hhh, vv
real   (irk) dx_gate, dy_gate, mass_gate, hsml_gate
real   (irk) fact1, rij, wij, w00, sum_wm

IF ( ic_histogram.ne.0) THEN

   DO i_hist = 1, n_hist
   
      npg2 = np_gate(i_hist)/2
      if (np_gate(i_hist).eq.1) then
         write(*,*) ' code stopped at Inject_histogram because np_gate is 1 '
      endif   

      icheck = 1
      t_ref  = time + dt_sph/2
      itotv1 = npoin_in_model + 1               
      itotv2 = npoin_in_model + np_gate(i_hist)
      itotvG = itotv1 + npg2  
       
      if (t_ref.lt.t_his(i_hist,1).or.t_ref.gt.t_his(i_hist,np_hist(i_hist)))  then
         icheck = -1
      elseif (itotv2.gt.npoin) then
         icheck = -1
      endif

      DO WHILE (icheck.eq.1)                    !      ------  do this only if we are going to inject parts

         dx_gate   = xl_gate(i_hist)/(np_gate(i_hist)-1)
         hsml_gate = hsml_gfact*dx_gate
   
         do ip_his = 1, np_hist(i_hist) - 1     !      ------  obtain interval (t0,t1) to which tref belongs
            t0_his = t_his(i_hist,ip_his )
            t1_his = t_his(i_hist, ip_his + 1 )
            if ( t_ref.ge.t0_his.and.t_ref.le.t1_his)  EXIT
         enddo
                                                !      ------  check if distance to last row is too small
         if ( ndof_track(i_hist,1).ne.0 ) THEN
            iaux  = np_track(i_hist)
            itotw = ndof_track(i_hist,iaux)
            dy_gate    = 0.0
            do idimn   = 1,ndimn
               dy_gate = dy_gate + (x(idimn,itotw)-x_gate(i_hist,idimn,npg2+1))**2.
            enddo
            dy_gate    = sqrt(dy_gate)
            if (dy_gate.lt.0.25*dx_gate) icheck = -1
         endif            
         if (icheck.eq.-1) EXIT   
            
        !      ------  Interpolate within time series
                                
         xi = (t_ref - t0_his)/(t1_his - t0_his)
         hh = (1.-xi)*h_his(i_hist,ip_his) + xi*h_his(i_hist,ip_his+1)
         vv = (1.-xi)*v_his(i_hist,ip_his) + xi*v_his(i_hist,ip_his+1)           
   
        !      ------  Obtain mass and height of nodes to be njected

        !fact1    = 1.0/ hsml_gate                                              **** we will use the 2d factor
         fact1    = 15.e0/(7.e0*pi*hsml_gate*hsml_gate)
         sum_wm   = 0.0
         w00      = (2./3.)*fact1
         itrack   = 1
         
         DO WHILE ( ndof_track(i_hist,itrack).ne.0.AND.itrack.le.nl_gate(i_hist) )
            itotw = ndof_track(i_hist,itrack) 
            rij   = 0.0
            do idimn = 1, ndimn
               rij   = rij + (x(idimn,itotw)-x_gate(i_hist,idimn,npg2+1))**2.
            enddo
            rij = sqrt(rij)
            xi  = rij/hsml_gate
            if ( (xi.ge.0.0).and.(xi.le.1.0) ) then
               wij = fact1*(2./3. -xi*xi + xi**3./2. )
            else if ( (xi.ge.1.0).and.(xi.le.2.0) ) then
               wij = fact1*(1./6.)*(2.-xi)**3.
            else
               wij = 0.0
            endif
            itrack = itrack +1
            if(wij.le.1.e-6) CYCLE
            sum_wm = sum_wm + wij*mass(itotw)
                        if (itrack.gt.6) EXIT
         ENDDO
         
         mass_gate = (hh-sum_wm)/w00
   
!      ------  Now INJECT: All particles at gate have same props

         npoin_in_model  = npoin_in_model  + np_gate(i_hist)
         
         np_track(i_hist) = np_track(i_hist) +1
         iaux = np_track(i_hist)
         ndof_track(i_hist, iaux) = itotvG
   
         DO itotv = itotv1, itotv2
            rho  (itotv) = hh
            mass (itotv) = mass_gate
            mass0(itotv) = mass_gate
            do idimn =1, ndimn
               vx(idimn,itotv) = vv * nv_gate(i_hist, idimn)
           enddo      
           itype(itotv) = 1
           hsml (itotv) = hsml_gate 
           p    (itotv) = (9.81/2.)*rho(itotv)*rho(itotv)
           u    (itotv) = 0.0
           if_out_domain(itotv) = 0
           do idimn=1,ndimn
              x0  (idimn,itotv) = x  (idimn,itotv)   
              x00 (idimn,itotv) = x0 (idimn,itotv)
              vx0 (idimn,itotv) = vx (idimn,itotv)
           enddo
         ENDDO 
   
         icheck = -1         
             
      ENDDO  
      
   ENDDO
   
ENDIF

End Subroutine Inject_histogram_SW_hh  

!-------------------------------------------------------------------

       Subroutine Inject_histogram_SW_BK   (time)  
            
!-------------------------------------------------------------------

!   Modified, we use hh for height and vv is now the weir discharge

implicit none
 
integer(ink) nrow_gate, npg2, iaux, itrack, i_hist, ip_gate 
integer(ink) icheck, ip_his, itotv, itotv1, itotv2, idimn, itotw, irow, itotvG, ipoin
integer(ink) naux

real   (irk) time, t_ref, xi, t0_his, t1_his, hh, hhh, vv, qq, dtime_inject
real   (irk) dx_gate, dy_gate, mass_gate, hsml_gate
real   (irk) fact1, rij, wij, w00, sum_wm
real   (irk) auxx, auxy

real   (irk), SAVE:: t_last_inject

IF ( ic_histogram.ne.0) THEN

!      ------  We check for more available nodes (all that have just exited)

    Do ipoin = npoin0, npoin
       if ( (if_out_domain(ipoin).eq.1).and.(if_out_domain_LAST(ipoin).eq.0)) THEN
           if_out_domain_LAST(ipoin) = 1
           n_avail_nodes = n_avail_nodes +1
           list_avail_nodes (n_avail_nodes) = ipoin
       endif
    enddo

   if (time.lt.1.e-6) t_last_inject=0

   DO i_hist = 1, n_hist
   
      npg2 = np_gate(i_hist)/2
      if (np_gate(i_hist).eq.1) then
         write(*,*) ' code stopped at Inject_histogram because np_gate is 1 '
      endif   

      icheck = 1
      t_ref  = time + dt_sph/2
       
      if (t_ref.lt.t_his(i_hist,1).or.t_ref.gt.t_his(i_hist,np_hist(i_hist)))  then
         icheck = -1
      elseif (itotv2.gt.npoin) then
         icheck = -1
      endif

      DO WHILE (icheck.eq.1)                    !      ------  do this only if we are going to inject parts

         dx_gate   = xl_gate(i_hist)/(np_gate(i_hist)-1)
         hsml_gate = hsml_gfact*dx_gate
   
         do ip_his = 1, np_hist(i_hist) - 1     !      ------  obtain interval (t0,t1) to which tref belongs
            t0_his = t_his(i_hist,ip_his )
            t1_his = t_his(i_hist, ip_his + 1 )
            if ( t_ref.ge.t0_his.and.t_ref.le.t1_his)  EXIT
         enddo
                                                !      ------  check if distance to last row is too small
         if ( ndof_track(i_hist,1).ne.0 ) THEN
            iaux  = np_track(i_hist)
            itotw = ndof_track(i_hist,iaux)
            dy_gate    = 0.0
            do idimn   = 1,ndimn
               dy_gate = dy_gate + (x(idimn,itotw)-x_gate(i_hist,idimn,npg2+1))**2.
            enddo
            dy_gate    = sqrt(dy_gate)
            if (dy_gate.lt.0.66*dx_gate) icheck = -1
         endif            
         if (icheck.eq.-1) EXIT   
            
             !      ------  Interpolate within time series
                                
         xi = (t_ref - t0_his)/(t1_his - t0_his)
         hh = (1.-xi)*h_his(i_hist,ip_his) + xi*h_his(i_hist,ip_his+1)
         vv = (1.-xi)*v_his(i_hist,ip_his) + xi*v_his(i_hist,ip_his+1)
         qq = vv
         vv = 0.0
   
         !      ------  Obtain mass and height of nodes to be njected
         
         dtime_inject = time - t_last_inject
         ! if (time.lt.1.e-6) dtime_inject = dt_sph
         mass_gate = 2.* hh * dx_gate * dx_gate 
   
         !      ------  Now INJECT: All particles at gate have same props

         npoin_in_model  = npoin_in_model  + np_gate(i_hist)
         t_last_inject   = time
         
         if (n_avail_nodes.lt.np_gate(i_hist)) THEN
             write(*,*) ' *******  code stopped at injection: no more available points '
             write(*,*) ' n_avail_nodes = ',n_avail_nodes, ' at gate we have: ' ,np_gate(i_hist)
             PAUSE
             STOP
         endif
         
         itotvG = list_avail_nodes(n_avail_nodes - npg2)
         DO ip_gate = 1, np_gate(i_hist)
            itotv    = list_avail_nodes(n_avail_nodes)
            list_avail_nodes(n_avail_nodes) = 0
            n_avail_nodes = n_avail_nodes -1
            rho  (itotv) = hh
            mass (itotv) = mass_gate
            mass0(itotv) = mass_gate
            do idimn =1, ndimn
               x   (idimn,itotv) =       x_gate(i_hist,idimn,ip_gate)
               vx  (idimn,itotv) = vv * nv_gate(i_hist,idimn)
               x0  (idimn,itotv) = x  (idimn,itotv)
               x00 (idimn,itotv) = x0 (idimn,itotv)
               vx0 (idimn,itotv) = vx (idimn,itotv)
           enddo      
           itype(itotv) = 1
           hsml (itotv) = hsml_gate 
           p    (itotv) = (9.81/2.)*rho(itotv)*rho(itotv)
           u    (itotv) = 0.0
           if_out_domain     (itotv) = 0
           if_out_domain_LAST(itotv) = 0
         ENDDO
         
         np_track(i_hist) = np_track(i_hist) +1
         iaux = np_track(i_hist)
         ndof_track(i_hist, iaux) = itotvG

1111     format(12(f10.2,2x))         
   
         icheck = -1         
             
      ENDDO  
      
   ENDDO
   
ENDIF

icheck = 0

End Subroutine Inject_histogram_SW_BK  

!-------------------------------------------------------------------

       Subroutine Inject_histogram_SW   (time)  
            
!-------------------------------------------------------------------

!   Modified, we use hh for height and vv is now the weir discharge

implicit none
 
integer(ink) nrow_gate, npg2, iaux, itrack, i_hist, ip_gate 
integer(ink) icheck, ip_his, itotv, itotv1, itotv2, idimn, itotw, irow, itotvG, ipoin

real   (irk) time, t_ref, xi, t0_his, t1_his, hh, hhh, vv, qq, dtime_inject
real   (irk) dx_gate, dy_gate, mass_gate, hsml_gate
real   (irk) fact1, rij, wij, w00, sum_wm
real   (irk) auxx, auxy

real   (irk) xx, yy
integer(ink) nn

real   (irk), SAVE:: t_last_inject

IF ( ic_histogram.ne.0) THEN

!      ------  We check for more available nodes (all that have just exited)

    Do ipoin = 1,npoin
       if ( (if_out_domain(ipoin).eq.1).and.(if_out_domain_LAST(ipoin).eq.0)) THEN
           if_out_domain_LAST(ipoin) = 1
           n_avail_nodes = n_avail_nodes +1
           list_avail_nodes (n_avail_nodes) = ipoin
       endif
    enddo

   if (time.lt.1.e-6) t_last_inject=0

   DO i_hist = 1, n_hist
   
      npg2 = np_gate(i_hist)/2
      if (np_gate(i_hist).eq.1) then
         write(*,*) ' code stopped at Inject_histogram because np_gate is 1 '
      endif   

      icheck = 1
      t_ref  = time + dt_sph/2
       
      if (t_ref.lt.t_his(i_hist,1).or.t_ref.gt.t_his(i_hist,np_hist(i_hist)))  then
         icheck = -1
      endif

      DO WHILE (icheck.eq.1)                    !      ------  do this only if we are going to inject parts

         dx_gate   = xl_gate(i_hist)/(np_gate(i_hist)-1)
         hsml_gate = hsml_gfact*dx_gate
   
         do ip_his = 1, np_hist(i_hist) - 1     !      ------  obtain interval (t0,t1) to which tref belongs
            t0_his = t_his(i_hist,ip_his )
            t1_his = t_his(i_hist, ip_his + 1 )
            if ( t_ref.ge.t0_his.and.t_ref.le.t1_his)  EXIT
         enddo
                                                !      ------  check if distance to last row is too small
                                                
         dy_gate = dx_gate                      !    only for first timestep                                            
         if ( ndof_track(i_hist,1).ne.0 ) THEN
            iaux  = np_track(i_hist)
            itotw = ndof_track(i_hist,iaux)
            dy_gate    = 0.0
            do idimn   = 1,ndimn
               auxx    = x_gate(i_hist,idimn,npg2+1)
               auxy    = x(idimn,itotw)
               dy_gate = dy_gate + (x(idimn,itotw)-x_gate(i_hist,idimn,npg2+1))**2.
            enddo
            dy_gate    = sqrt(dy_gate)
            if (dy_gate.lt.0.40*dx_gate) icheck = -1
         endif            
         if (icheck.eq.-1) EXIT   
            
         !      ------  Interpolate within time series
                                
         xi = (t_ref - t0_his)/(t1_his - t0_his)
         hh = (1.-xi)*h_his(i_hist,ip_his) + xi*h_his(i_hist,ip_his+1)
         vv = (1.-xi)*v_his(i_hist,ip_his) + xi*v_his(i_hist,ip_his+1)
   
         !      ------  Obtain mass and height of nodes to be injected
         
         dtime_inject = time - t_last_inject
         if (time.lt.1.e-6) dtime_inject = dt_sph 
   
         !      ------  Now INJECT: All particles at gate have same props

         npoin_in_model  = npoin_in_model  + np_gate(i_hist)
         t_last_inject   = time
         
         if (n_avail_nodes.lt.np_gate(i_hist)) THEN
             write(*,*) ' *******  code stopped at injection: no more available points '
             write(*,*) ' n_avail_nodes = ',n_avail_nodes, ' at gate we have: ' ,np_gate(i_hist)
             PAUSE
             STOP
         endif
         
         itotvG = list_avail_nodes(n_avail_nodes - npg2)
         DO ip_gate = 1, np_gate(i_hist)
            itotv    = list_avail_nodes(n_avail_nodes)
            list_avail_nodes(n_avail_nodes) = 0
            n_avail_nodes = n_avail_nodes -1
            
!            if (type_of_histogram.eq.1) THEN
!                mass_gate    = 2.* hh * dx_gate * dx_gate
!                rho  (itotv) = hh
!            else if (type_of_histogram.eq.2) THEN
!                mass_gate    = hh*dtime_inject*dx_gate      !  ( m3/m.s  * m  * s )
!                rho  (itotv) = mass_gate/(dx_gate*dy_gate)   
!            endif

            mass_gate    = hh*(dtime_inject*vv)*dx_gate         ! *** for efficiency, out of loop
            
            mass (itotv) = mass_gate
            mass0(itotv) = mass_gate            
            do idimn =1, ndimn
               x   (idimn,itotv) =       x_gate(i_hist,idimn,ip_gate)
               vx  (idimn,itotv) = vv * nv_gate(i_hist,idimn)
               x0  (idimn,itotv) = x  (idimn,itotv)
               x00 (idimn,itotv) = x0 (idimn,itotv)
               vx0 (idimn,itotv) = vx (idimn,itotv)
           enddo
                   
           rho  (itotv) = hh 
           rho0 (itotv) = rho(itotv)
           mass (itotv) = mass_gate
           mass0(itotv) = mass_gate       
           itype (itotv) = 1
           hsml  (itotv) = hsml_gate 
           p     (itotv) = (9.81/2.)*rho(itotv)*rho(itotv)
           u     (itotv) = 0.0
           if_out_domain     (itotv) = 0
           if_out_domain_LAST(itotv) = 0

         ENDDO
         
         np_track(i_hist) = np_track(i_hist) +1
         iaux = np_track(i_hist)
         ndof_track(i_hist, iaux) = itotvG
         
         DO itotv = 1, np_bakG(i_hist)
            rho  (itotv) = hh  
            rho0 (itotv) = rho(itotv)
            mass (itotv) = mass_gate
            mass0(itotv) = mass_gate 
            vx (:,itotv) = 0.0
            x0  (:,itotv) = x  (:,itotv)
            x00 (:,itotv) = x0 (:,itotv)
            vx0 (:,itotv) = vx (:,itotv)   
            itype(itotv) = 1
            hsml (itotv) = hsml_gate 
            p    (itotv) = (9.81/2.)*rho(itotv)*rho(itotv)
            u    (itotv) = 0.0
         ENDDO            
   
1111     format(12(f10.2,2x))         
   
         icheck = -1         
             
      ENDDO  
      
   ENDDO
   
ENDIF

icheck = 0

End Subroutine Inject_histogram_SW  


! ----------------------------------------------------------------------   

      subroutine  Get_hs_at_w_SW  

! ----------------------------------------------------------------------   

!   Subroutine  
 
implicit none

integer(ink) ipoin, idimn
integer(ink) ipoigx, ipoigy, ipoig
integer(ink) inx0, iny0, inx1,iny1, inx, iny 
integer(ink) ieleg, n1, n2, n3, n4, kpoin

real   (irk) xx, yy, distx2, disty2, dist
real   (irk) xref, yref, xi, eta
real   (irk) sh1, sh2, sh3, sh4, Zp, dZpdx, dZpdy
real   (irk) dsh1dx, dsh2dx, dsh3dx, dsh4dx
real   (irk) dsh1dy, dsh2dy, dsh3dy, dsh4dy
real   (irk) dxginv, dyginv

real   (irk), allocatable:: auxv(:,:)

if(.NOT.allocated(auxv)) allocate (auxv(2,npoig))
auxv = 0.0

!if (.not.allocated( hs_plus_Z_topo ) ) then
!
!
!   allocate (    hs_plus_Z_topo    (npoig))
!   allocate (    hs_plus_Z_w       (npoin))
!   allocate (gradhs_plus_Z_w (ndimn,npoin)) 
!   hs_plus_Z_topo = 0.0
!   hs_plus_Z_w    = 0.0
!   gradhs_plus_Z_w= 0.0
!   
!endif

! chuan to consider the z+hs for water particles

if (.not.allocated( hs_plus_Z_w ) ) then

   allocate (    hs_plus_Z_w       (npoin))
   allocate (gradhs_plus_Z_w (ndimn,npoin)) 

   hs_plus_Z_w    = 0.0
   gradhs_plus_Z_w= 0.0  
endif

!just in case
if (.not.allocated( hs_plus_Z_topo ) ) then
   allocate (    hs_plus_Z_topo    (npoig))
   hs_plus_Z_topo= 0.0  
endif

DO ipoin=1,npoin

   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   if ( itype(ipoin).ne.1) CYCLE
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
      auxv(1,ipoig) = auxv(1,ipoig) + rho(ipoin)/dist            
      auxv(2,ipoig) = auxv(2,ipoig) + 1.0/dist 
   enddo
   enddo
   
ENDDO      

Do ipoig =1,npoig
   if (auxv(2,ipoig).gt.1.e-6) then
      hs_plus_Z_topo (ipoig) = topol(1,ipoig) +  auxv(1,ipoig)/auxv(2,ipoig)
   else
      hs_plus_Z_topo (ipoig) = topol(1,ipoig)
   endif
enddo

DO ipoin = 1, npoin

   if ( if_Out_Domain(ipoin).eq.1) CYCLE
   if ( itype(ipoin).ne.2) CYCLE 

   xx = x(1,ipoin)   
   yy = (yming+ymaxg)/2
   if (ndimn.eq.2) yy = x(2,ipoin)  
   ipoigx    = (xx-xming)/deltxg + 1
   ipoigy    = (yy-yming)/deltyg + 1
   
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
   
   dsh1dx = (-1.0 + eta )/deltxg
   dsh2dx = ( 1.0 - eta )/deltxg
   dsh3dx = (       eta )/deltxg
   dsh4dx = (     - eta )/deltxg
    
   dsh1dy = (-1.0 +  xi )/deltyg
   dsh2dy = (     -  xi )/deltyg
   dsh3dy = (        xi )/deltyg
   dsh4dy = ( 1.0 -  xi )/deltyg
   
   Zp     =      sh1 * hs_plus_Z_topo (n1) + sh2 * hs_plus_Z_topo (n2) &  
               + sh3 * hs_plus_Z_topo (n3) + sh4 * hs_plus_Z_topo (n4)
   hs_plus_Z_w (ipoin) = Zp     
   
   dZpdx  =    hs_plus_Z_topo (n1) * dsh1dx +  &
               hs_plus_Z_topo (n2) * dsh2dx +  &
               hs_plus_Z_topo (n3) * dsh3dx +  &
               hs_plus_Z_topo (n4) * dsh4dx 
   gradhs_plus_Z_w (1,ipoin) = dZpdx             
             
   if (ndimn.eq.2) then
   
      dZpdy  =    hs_plus_Z_topo (n1) * dsh1dy + &
                  hs_plus_Z_topo (n2) * dsh2dy + &
                  hs_plus_Z_topo (n3) * dsh3dy + &
                  hs_plus_Z_topo (n4) * dsh4dy         !CHUAN 2018.4.7
       gradhs_plus_Z_w (2,ipoin) = dZpdy                  
   endif  

Enddo 

deallocate(auxv)   

END  SUBROUTINE Get_hs_at_w_SW 


! ----------------------------------------------------------------------   

      subroutine  Drum_drive_SW  

! ----------------------------------------------------------------------   

!   Subroutine to drive drum cases, PUBLIC 
 
implicit none

integer(ink) ipoin
real   (irk) xx, dxx, sinth, costh  

if (R_drum.gt.1.e-3) then
   DO ipoin = 1,npoin
      if (if_out_domain(ipoin).eq.1) CYCLE
      xx    = x(1,ipoin)
      sinth = xx/R_drum
      costh = sqrt (1-sinth*sinth)
      dxx   = R_drum*costh*dt_sph*w_drum_radps
      x(1,ipoin) = x(1,ipoin) + dxx
   enddo
endif   


END  SUBROUTINE Drum_drive_SW

!---------------------------------------------------------------

     subroutine Get_s_at_w_DF_SW (if_input)

!---------------------------------------------------------------
 
!   Subroutine to calculate
!
!   (1)     (a) h_sw (ntotal): hw at soil pts and hs at water points
!           (b) v_sw (ndimn, ntotal): same for velocities
!           (c) h_DF (ntotal) = hs + hw
!           (d) n_DF (ntotal) = hw / h_DF
!   We use:         
!            mass(   ntotv) : Particle masses 
!            vx(ndimn,ntov) : particle velocity
!            current%Pint_type: 
!                 0 for particles of same kind
!                 1 pic I soil pic J water
!            wi(ntotal) used for normalization
!
!            h(i) = sum_J (Wij.Mj)
!            V(i) = sum_J (Vj.Wij.Mj/RHOj)
!
!            Normalization denom. 
!                 w(i) = sum_J (Wij.Mj/RHOj)
!
!   (2) Obtain h_DF = hs + hf , n_DF 
!
!       We will consider two cases:  
!                                                  
!          (c1) The column is not fully saturated                                                                    
!                hrelSat_DF < 1                                            
!                Then:  
!                     modify HrelSat_DF using hw & hs keeping n_DF constant
!                     and use it to get the new h_DF   
!                 check: new HrelSat_DF should be LE 1. Otherwise
!                     hrelSat = 1
!                     compute n_DF and H_Df    
!                                                  
!          (c2) The column is fully saturated                                                                    
!                hrelSat_DF = 1                                            
!                Then:  
!                     compute n_DF and h_DF  using hw & hs and hrelSat=1
!                 check: new n_DF should be GE nmin Otherwise
!                     n_DF = n min
!                     compute hrelSat and h_Df                                    

! --------------------------------------------------------------

implicit none

integer(ink), intent(IN), OPTIONAL:: if_input          ! if 1 do not normalize h close to walls 
integer(ink)  I, J , int_sum, iw, ic_input

real   (irk) Sr, alphaS, alph, n_min, hs, hw, n0, denf
real   (irk) rij, rij_inv
real   (irk), allocatable:: sum_rij(:), sum_hsw(:),sum_vsw(:,:)
real   (irk), allocatable:: wi(:)

if (ic_ws_Interact.NE.2) then
   RETURN       ! routine also called from MAIN, where interact not known
endif

ic_input = -999
if ( present(if_input) )   ic_input = if_input
  
if(.NOT.allocated(wi)) allocate ( wi(ntotal) )

wi      = 0.0
int_sum = 0
h_sw    = 0.0
v_sw    = 0.0
!  h_DF    = 0.0            MP sept 2021
!  n_DF    = 0.0
                            !   CH MP 4 march 2019  
Sr         = const (19)     !   upper unsaturated layer Sr
alphaS     = const (14)     !   sat. lower layer height = h * alphaS
ic_hrelSat = const (34)     !   alfa can be ct=C14, or if C34=1 be modified
n_min      = const (35)     !   minimum DF porosity below which desaturates

IF (ic_input.eq.0) then     !  when nodes coincide at t=0 make hs = hw

   do i = 1, npoin_s
      iw = i + npoin_s
      h_sw (i) = rho(iw)
      h_sw (iw)= rho(i)
   enddo
   v_sw = 0.0

ELSEIF (ic_input.eq.1) THEN      ! interpolate
   
   if (.not.allocated(sum_rij)) then
        allocate (sum_rij(npoin), sum_hsw(npoin), sum_vsw(ndimn,npoin))
   endif  

   rij     = 0.0
   sum_rij = 0.0
   sum_hsw = 0.0
   sum_vsw = 0.0
   rij_inv = 0.0

   current=> last
   Do while (associated(current))
      if (current%Pint_type.eq.1) then
         i = current%pair_i 
         j = current%pair_j 
         rij = (x(1,i)-x(1,j))*(x(1,i)-x(1,j))
         if (ndimn.eq.2) then
            rij = rij + (x(2,i)-x(2,j))*(x(2,i)-x(2,j))
         endif
         if (rij.lt.1.e-8) then
            rij = 1.e-8
         endif
         rij_inv       = 1.0/rij
         sum_hsw   (i) = sum_hsw(i)    + rho (j)/rij
         sum_hsw   (j) = sum_hsw(j)    + rho (i)/rij
         sum_vsw (:,i) = sum_vsw (:,i) + vx(:,j)/rij
         sum_vsw (:,j) = sum_vsw (:,j) + vx(:,i)/rij
         sum_rij (i) = sum_rij (i) + rij_inv
         sum_rij (j) = sum_rij (j) + rij_inv
      endif 
      current => current%next
   enddo
   
   do i = 1,npoin
      if (sum_rij(i).gt.1.e-5) then
       h_sw   (i) = sum_hsw(i)/sum_rij(i)
       v_sw (:,i) = sum_vsw(:,i)/sum_rij(i)
      endif
   enddo  

ELSE   !      ------  use the linked list.  when (2)
          
   current=> last                      !  used for normalizing
   Do while (associated(current))
      if (current%Pint_type.eq.1) then
          i = current%pair_i 
          j = current%pair_j 
          wi(i) = wi(i) + mass(j)/rho(j)*current%w 
          wi(j) = wi(j) + mass(i)/rho(i)*current%w
          int_sum = int_sum +1 
     endif   
      current=>current%next
   enddo

   if (int_sum.eq.0) wi = 1.0   !  if no w-s interactions, normalize using 1

   !      ------  Calculate SPH sum for h_sw:

   h_sw = 0.0
   v_sw = 0.0

   current=> last
   Do while (associated(current))
       if (current%Pint_type.eq.1) then
          i = current%pair_i 
          j = current%pair_j 
          h_sw (i)   = h_sw(i)   + mass(j)*current%w 
          h_sw (j)   = h_sw(j)   + mass(i)*current%w
          v_sw (:,i) = v_sw(:,i) + mass(j)*current%w * (vx(:,j)/rho(j)) 
          v_sw (:,j) = v_sw(:,j) + mass(i)*current%w * (vx(:,i)/rho(i))
       endif      
       current=>current%next
   enddo

   !      ------  Compute the normalized rho, rho=sum(rho)/sum(w)
               
   do i = 1, npoin
       if (wi(i).LT.1.e-6) wi(i) = 1         ! new change MP for screens perm to water                  
       if (if_out_domain(i).eq.1) cycle 
       
       if (nor_density) then                !  avoids computing twice if ifcorrect and nor_density
            h_sw  (i) = h_sw  (i)/wi(i)     ! MP changed structure 17 June 2019
       elseif (allocated(if_correct)) then
           if (if_correct(i).eq.1) then      ! for nodes close to vn0 boundaries                                         
               h_sw  (i) = h_sw  (i)/wi(i)   ! ** CHANGED MP 6 June 2017  
           endif
       endif

       if (h_sw(i).LT.h_inf_SW) then         ! ** avoid n .ne.1 at water points
           v_sw(:,i) = 0.0
           h_sw(i)   = h_inf_SW            ! changed MP 14th June 2019
       else   
           v_sw(:,i) = v_sw(:,i)/wi(i)
       endif   
   enddo

ENDIF

!      ------  Obtain h_DF = hs + hf  

DO i = 1, npoin   ! --------------------------------------------------------------------

   alph = alphaS
   if (ic_hrelSat.EQ.1) then
       alph = hrelSat_DF (i)
   endif
   if (i.le.npoin_s) then
      hs = rho(i)
      hw = h_sw(i)
   else
      hw = rho(i)
      hs = h_sw(i) 
   endif

   IF (ic_hrelSat.eq.1) THEN   ! ----------------------------------------------
 
       IF (abs (alph -1.0).LT.1.e-4) THEN   ! colunmn is FULLY Saturated, i.e. hrelSat_DF(i) = 1 
                                                            
          !    (case 1) The column is fully saturated:                                                                
          !                hrelSat_DF = 1                                            
          !                Then: :  
          !                     compute n_DF and h_DF  using hw & hs and hrelSat=1
          !                 check: new n_DF should be GE nmin Otherwise
          !                     n_DF = n min
          !                     compute hrelSat and h_Df
                                 
           denf     = hrelSat_DF (i) + Sr * (1.-hrelSat_DF (i)) 
           n_DF (i) = hw / (hw + denf*hs)
           h_DF (i)      =  (hs+hw) / (1.- n_DF (i)*(1.-Sr)*(1.-hrelSat_DF (i)) ) 
           if (n_DF (i).LT.n_min) then         ! only when using hrelSat and n<nmin
              n_DF (i)      = n_min
              hrelSat_DF(i) = (1./(1.-Sr)) * ( ((1.- n_min )*hw/(n_min*hs)) - Sr)
              h_DF (i)      = (hs+hw) / (1.- n_min*(1.-Sr)*(1.-hrelSat_DF(i)) )
           endif          
                 
       ELSE                                          
          !    (case 2) The column is not fully saturated:  hrelSat_DF(i) < 1                                                              
          !                hrelSat_DF < 1                                            
          !                Then:  
          !                     modify HrelSat_DF using hw & hs keeping n_DF constant
          !                     and use it to get the new h_DF   
          !                 check: new HrelSat_DF should be LE 1. Otherwise
          !                     hrelSat = 1
          !                     compute n_DF and H_Df 

          hrelSat_DF(i) = (1./(1.-Sr)) * ( ((1.- n_DF(i) )*hw/(n_DF(i)*hs)) - Sr)
          h_DF (i)      =  (hs+hw) / (1.- n_DF(i)*(1.-Sr)*(1.-hrelSat_DF(i)) )  
          if (hrelSat_DF(i).GT.1.00) then
              hrelSat_DF(i) = 1                          
              denf     = hrelSat_DF(i) + Sr*(1.- hrelSat_DF(i)) 
              n_DF (i) = hw / (hw + denf*hs)
              h_DF (i)      =  (hs+hw) / (1.- n_DF(i)*(1.-Sr)*(1.-hrelSat_DF(i)) )
          endif
       ENDIF

   ELSE     ! ----------------------------- ! DEFAULT   MP 30 april 19   ic_hrelSat_DF = 0   THEN   

       denf = alph + Sr*(1.-alph)                         
       n0   = hw / (hw + denf*hs)
       n_DF (i) = n0                         
       h_DF (i) = (hs+hw) / (1.- n0*(1.-Sr)*(1.-alph) )

   ENDIF

ENDDO    ! -------------------------end loop i points ----------------------------- 
   
if (icpwp.ge.1) then                !  DFs + PWp RK4 needs to know n_dfaç
    auxS1 (1:npoin) = n_DF(1:npoin)
endif   
   
deallocate (wi)


END SUBROUTINE Get_s_at_w_DF_SW

!---------------------------------------------------------------

     subroutine Get_s_at_w_DF_SW_Saeid (if_input)

!---------------------------------------------------------------
 
!   Subroutine to calculate 
!
!           (a) h_sw (ntotal): hw at soil pts and hs at water points
!           (b) v_sw (ndimn, ntotal): same for velocities
!           (c) h_DF (ntotal) = hs + hw
!           (d) n_DF (ntotal) = hw / h_DF
!   We use:         
!            mass(   ntotv) : Particle masses 
!            vx(ndimn,ntov) : particle velocity
!            current%Pint_type: 
!                 0 for particles of same kind
!                 1 pic I soil pic J water
!            wi(ntotal) used for normalization
!
!            h(i) = sum_J (Wij.Mj)
!            V(i) = sum_J (Vj.Wij.Mj/RHOj)
!
!            Normalization denom. 
!                 w(i) = sum_J (Wij.Mj/RHOj)
!
! --------------------------------------------------------------

implicit none

integer(ink), intent(IN), OPTIONAL:: if_input          ! if 1 do not normalize h close to walls 
integer(ink)  I, J , int_sum, iw, ic_input

real   (irk) Sr, alphaS, alph, n_min, hs, hw, n0, denf
real   (irk) rij, rij_inv
real   (irk), allocatable:: sum_rij(:), sum_hsw(:),sum_vsw(:,:)
real   (irk), allocatable:: wi(:)

if (ic_ws_Interact.NE.2) then
   RETURN       ! routine also called from MAIN, where interact not known
endif

ic_input = -999
if ( present(if_input) )   ic_input = if_input
  
if(.NOT.allocated(wi)) allocate ( wi(ntotal) )

wi      = 0.0
int_sum = 0
h_sw    = 0.0
v_sw    = 0.0
h_DF    = 0.0
n_DF    = 0.0
                            !   CH MP 4 march 2019  
Sr         = const (19)     !   upper unsaturated layer Sr
alphaS     = const (14)     !   sat. lower layer height = h * alphaS
ic_hrelSat = const (34)     !   alfa can be ct=C14, or if C34=1 be modified
n_min      = const (35)     !   minimum DF porosity below which desaturates

IF (ic_input.eq.0) then     !  when nodes coincide at t=0 make hs = hw

   do i = 1, npoin_s
      iw = i + npoin_s
      h_sw (i) = rho(iw)
      h_sw (iw)= rho(i)
   enddo
   v_sw = 0.0

ELSEIF (ic_input.eq.1) THEN      ! interpolate
   
   if (.not.allocated(sum_rij)) then
        allocate (sum_rij(npoin), sum_hsw(npoin), sum_vsw(ndimn,npoin))
   endif  

   rij     = 0.0
   sum_rij = 0.0
   sum_hsw = 0.0
   sum_vsw = 0.0
   rij_inv = 0.0

   current=> last
   Do while (associated(current))
      if (current%Pint_type.eq.1) then
         i = current%pair_i 
         j = current%pair_j 
         rij = (x(1,i)-x(1,j))*(x(1,i)-x(1,j))
         if (ndimn.eq.2) then
            rij = rij + (x(2,i)-x(2,j))*(x(2,i)-x(2,j))
         endif
         if (rij.lt.1.e-8) then
            rij = 1.e-8
         endif
         rij_inv       = 1.0/rij
         sum_hsw   (i) = sum_hsw(i)    + rho (j)/rij
         sum_hsw   (j) = sum_hsw(j)    + rho (i)/rij
         sum_vsw (:,i) = sum_vsw (:,i) + vx(:,j)/rij
         sum_vsw (:,j) = sum_vsw (:,j) + vx(:,i)/rij
         sum_rij (i) = sum_rij (i) + rij_inv
         sum_rij (j) = sum_rij (j) + rij_inv
      endif 
      current => current%next
   enddo
   
   do i = 1,npoin
      if (sum_rij(i).gt.1.e-5) then
       h_sw   (i) = sum_hsw(i)/sum_rij(i)
       v_sw (:,i) = sum_vsw(:,i)/sum_rij(i)
      endif
   enddo  

ELSE   !      ------  use the linked list.
          
   current=> last                      !  used for normalizing
   Do while (associated(current))
      if (current%Pint_type.eq.1) then
          i = current%pair_i 
          j = current%pair_j 
          wi(i) = wi(i) + mass(j)/rho(j)*current%w 
          wi(j) = wi(j) + mass(i)/rho(i)*current%w
          int_sum = int_sum +1 
     endif   
      current=>current%next
   enddo

   if (int_sum.eq.0) wi = 1.0   !  if no w-s interactions, normalize using 1

   !      ------  Calculate SPH sum for h_sw:

   h_sw = 0.0
   v_sw = 0.0

   current=> last
   Do while (associated(current))
       if (current%Pint_type.eq.1) then
          i = current%pair_i 
          j = current%pair_j 
          h_sw (i)   = h_sw(i)   + mass(j)*current%w 
          h_sw (j)   = h_sw(j)   + mass(i)*current%w
          v_sw (:,i) = v_sw(:,i) + mass(j)*current%w * (vx(:,j)/rho(j)) 
          v_sw (:,j) = v_sw(:,j) + mass(i)*current%w * (vx(:,i)/rho(i))
       endif      
       current=>current%next
   enddo

   !      ------  Compute the normalized rho, rho=sum(rho)/sum(w)
               
   do i = 1, npoin
       if (wi(i).LT.1.e-6) wi(i) = 1         ! new change MP for screens perm to water                  
       if (if_out_domain(i).eq.1) cycle 
       
       if (nor_density) then                !  avoids computing twice if ifcorrect and nor_density
            h_sw  (i) = h_sw  (i)/wi(i)     ! MP changed structure 17 June 2019
       elseif (allocated(if_correct)) then
           if (if_correct(i).eq.1) then      ! for nodes close to vn0 boundaries                                         
               h_sw  (i) = h_sw  (i)/wi(i)   ! ** CHANGED MP 6 June 2017  
           endif
       endif

       if (h_sw(i).LT.h_inf_SW) then         ! ** avoid n .ne.1 at water points
           v_sw(:,i) = 0.0
           h_sw(i)   = h_inf_SW            ! changed MP 14th June 2019
       else   
           v_sw(:,i) = v_sw(:,i)/wi(i)
       endif   
   enddo

ENDIF

!      ------  Obtain h_DF = hs + hf  

DO i = 1, npoin
   alph = alphaS
   if (ic_hrelSat.EQ.1) then
       alph = hrelSat_DF (i)
   endif
   if (i.le.npoin_s) then
      hs = rho(i)
      hw = h_sw(i)
   else
      hw = rho(i)
      hs = h_sw(i) 
   endif
   denf = alph + Sr*(1.-alph)                       ! DEFAULT   MP 30 april 19
   n0   = hw / (hw + denf*hs)
   n_DF (i) = n0                         
   h_DF (i) = (hs+hw) / (1.- n0*(1.-Sr)*(1.-alph) )
   if (n0.LT.n_min.AND.ic_hrelSat.eq.1) then         ! only when using hrelSat and n<nmin
      n0 = n_min
      alph = (1./(1.-Sr)) * ( ((1.- n0 )*hw/(n0*hs)) - Sr)
      hrelSat_DF(i) = alph
      n_DF (i)      = n0
      h_DF (i)      = (hs+hw) / (1.- n0*(1.-Sr)*(1.-alph) )
   endif

ENDDO 
   
if (icpwp.ge.1) then                !  DFs + PWp RK4 needs to know n_dfaç
    auxS1 (1:npoin) = n_DF(1:npoin)
endif   
   
deallocate (wi)


END SUBROUTINE Get_s_at_w_DF_SW_Saeid

!---------------------------------------------------------------

     subroutine Get_s_at_w_DF_SW_BAK

!---------------------------------------------------------------
 
!   Subroutine to calculate 
!
!           (a) h_sw (ntotal): hw at soil pts and hs at water points
!           (b) v_sw (ndimn, ntotal): same for velocities
!           (c) h_DF (ntotal) = hs + hw
!           (d) n_DF (ntotal) = hw / h_DF
!   We use:         
!            mass(   ntotv) : Particle masses 
!            vx(ndimn,ntov) : particle velocity
!            current%Pint_type: 
!                 0 for particles of same kind
!                 1 pic I soil pic J water
!            wi(ntotal) used for normalization
!
!            h(i) = sum_J (Wij.Mj)
!            V(i) = sum_J (Vj.Wij.Mj/RHOj)
!
!            Normalization denom. 
!                 w(i) = sum_J (Wij.Mj/RHOj)
!
!        
!
! --------------------------------------------------------------

implicit none

integer(ink)  I, J 

real   (irk), allocatable:: wi(:)

if (ic_ws_Interact.eq.2) then  ! routine also called from MAIN, where interact not known
  
   allocate ( wi(ntotal) )

   wi = 0.0
      
!      ------  use the linked list.

   current=> last                      !  used for normalizing
   Do while (associated(current))
      if (current%Pint_type.eq.1) then
         i = current%pair_i 
         j = current%pair_j 
         wi(i) = wi(i) + mass(j)/rho(j)*current%w 
         wi(j) = wi(j) + mass(i)/rho(i)*current%w
      endif   
      current=>current%next
   enddo

!      ------  Calculate SPH sum for h_sw:

   h_sw = 0.0
   v_sw = 0.0

   current=> last
   Do while (associated(current))
      if (current%Pint_type.eq.1) then
         i = current%pair_i 
        j = current%pair_j 
         h_sw (i)   = h_sw(i)   + mass(j)*current%w 
         h_sw (j)   = h_sw(j)   + mass(i)*current%w
         v_sw (:,i) = v_sw(:,i) + mass(j)*current%w * (vx(:,j)/rho(j)) 
         v_sw (:,j) = v_sw(:,j) + mass(i)*current%w * (vx(:,i)/rho(i))
      endif      
      current=>current%next
   enddo

!      ------  Compute the normalized rho, rho=sum(rho)/sum(w)
               
    do i = 1, npoin                     
       if (if_out_domain(i).eq.1) cycle   
       if (allocated(if_correct)) then
           if (if_correct(i).eq.1 ) then     ! for nodes close to vn0 boundaries                                         
               h_sw  (i) = h_sw  (i)/wi(i)   ! ** CHANGED MP 6 June 2017  
           endif
       elseif (nor_density) then            ! not always allocated. EXE works, debugger stops
           h_sw  (i) = h_sw  (i)/wi(i)       ! MP changed 06 Feb 18
       endif
       if (h_sw(i).LT.h_inf_SW) then         ! ** avoid n .ne.1 at water points
           v_sw(:,i) = 0.0
           h_sw(i)   = h_inf_SW
       else   
           v_sw(:,i) = v_sw(:,i)/wi(i)
       endif   
    enddo

!      ------  Obtain h_DF = hs + hf  

   h_DF = rho + h_sw

!      ------  and n = (h-hs)/h

   do i = 1, npoin

      if (i.le.npoin_s) then
          n_DF(i) = (h_DF(i)- rho(i))/h_DF(i)
      else       
          n_DF(i) = (h_DF(i)-h_sw(i))/h_DF(i) 
      endif
   enddo   
   
   if (icpwp.ge.1) then                !  DFs + PWp RK4 needs to know n_dfaç
      auxS1 (1:npoin) = n_DF(1:npoin)
   endif   
   
   deallocate (wi)

ENDIF

END SUBROUTINE Get_s_at_w_DF_SW_BAK

!---------------------------------------------------------------

     subroutine Get_s_at_w_Wir_SW

!---------------------------------------------------------------
 
!   Subroutine to calculate in WIR problems
!
!           (a) h_sw (ntotal): hw at soil pts and hs at water points in 
!           (b) v_sw (ndimn, ntotal): same for velocities
!           (c) h_DF (ntotal) = hs + hw 
!   We use:         
!            mass(   ntotv) : Particle masses 
!            vx(ndimn,ntov) : particle velocity
!            current%Pint_type: 
!                 0 for particles of same kind
!                 1 pic I soil pic J water
!            wi(ntotal) used for normalization
!
!            h(i) = sum_J (Wij.Mj)
!            V(i) = sum_J (Vj.Wij.Mj/RHOj)
!
!            Normalization denom. 
!                 w(i) = sum_J (Wij.Mj/RHOj)
!
!        
!
! --------------------------------------------------------------

implicit none

integer(ink)  I, J 

real   (irk), allocatable:: wi(:)

  
allocate ( wi(ntotal) )

wi = 0.0
      
!      ------  use the linked list.

current=> last                      !  used for normalizing
Do while (associated(current))
   if (current%Pint_type.eq.1) then
       i = current%pair_i 
       j = current%pair_j 
       wi(i) = wi(i) + mass(j)/rho(j)*current%w 
       wi(j) = wi(j) + mass(i)/rho(i)*current%w
   endif   
   current=>current%next
enddo

!      ------  Calculate SPH sum for h_sw:

h_sw = 0.0
v_sw = 0.0

current=> last
Do while (associated(current))
   if (current%Pint_type.eq.1) then
      i = current%pair_i 
      j = current%pair_j 
      h_sw (i)   = h_sw(i)   + mass(j)*current%w 
      h_sw (j)   = h_sw(j)   + mass(i)*current%w
      v_sw (:,i) = v_sw(:,i) + mass(j)*current%w * (vx(:,j)/rho(j)) 
      v_sw (:,j) = v_sw(:,j) + mass(i)*current%w * (vx(:,i)/rho(i))
   endif      
   current=>current%next
enddo

!      ------  Compute the normalized rho, rho=sum(rho)/sum(w)
               
do i = 1, npoin                     
   if (if_out_domain(i).eq.1) cycle   
   if (allocated(if_correct)) then
       if (if_correct(i).eq.1 ) then     ! for nodes close to vn0 boundaries                                         
           h_sw  (i) = h_sw  (i)/wi(i)   ! ** CHANGED MP 6 June 2017  
       endif
   elseif (nor_density) then            ! not always allocated. EXE works, debugger stops
       h_sw  (i) = h_sw  (i)/wi(i)       ! MP changed 06 Feb 18
   endif
   if (h_sw(i).LT.h_inf_SW) then         ! ** avoid n .ne.1 at water points
       v_sw(:,i) = 0.0
       h_sw(i)   = h_inf_SW
   else   
       v_sw(:,i) = v_sw(:,i)/wi(i)
   endif   
enddo

!      ------  Obtain h_DF = hs + hf  

h_DF = rho + h_sw

deallocate (wi)

END SUBROUTINE Get_s_at_w_Wir_SW
 
 
!  ** MP s_implicit 17.03.18   2021   *************************************+

!----------------------------------------------------------------------

      subroutine Second_Step_FS_SW_FEM 

!----------------------------------------------------------------------

 
implicit none

if (.NOT.allocated (gstif).AND.SPH_t_Integ_Alg.eq.7) then 
   allocate (gstif    (npoin,npoin)) ; gstif    = 0.0 
   allocate (vstar    (ndimn,npoin)) ; vstar    = 0.0  
   allocate (vn_plus1 (ndimn,npoin)) ; vn_plus1 = 0.0 
   allocate (grad_p   (ndimn,npoin)) ; grad_p   = 0.0 
   allocate (RHS_p    (npoin))       ; RHS_p    = 0.0
   allocate (iffix    (npoin))       ; iffix    = 0 
   allocate (fixed    (npoin))       ; fixed    = 0.0 
endif
ss_niac  = 0
ss_niac0 = 0   

if     (type_FS_problem.eq.5) then           ! -------------------------  TYPE 5 SHAO ------------------------
   if (.not.allocated(intmat_FS))  then      ! Only 1st time. Note that moving 2d nodes will need update
      Call Get_P_BCs_SW_fem                  !   read iffix and fixed
   endif
   Call AuxVars_FS_SW_fem  
   if (itimestep_sph.eq.1) rho0 = rho        !  
   Call Get_RHS_FS_SW_fem_5         !   Comp.(a) vstar <- vx (b) div_vstar, (c) p from intforces, (d) rhs_p
   Call Get_Laplacian_SW_fem 
   !rhs_p  =  rhs_p*10000
   !gstif2 = gstif2*10000 
   call solvegc2 (npoin, nelem_FS, nnode_FS, 1, npoin, 2, intmat_FS , &   ! p formulation from FS 2017   
                  iffix, fixed, gstif2, rhs_p, p, tolcg_MAIN) 
      !                                       ---                                        
   call Get_vn1_FS_SW_5

elseif (type_FS_problem.eq.4) then           ! -------------------------  TYPE 4  p  ------------------------
   if (.not.allocated(intmat_FS))  then      ! Only 1st time. Note that moving 2d nodes will need update
      Call Get_P_BCs_SW_fem                  !   read iffix and fixed
   endif
   Call AuxVars_FS_SW_fem  

   Call Get_RHS_FS_SW_fem         !   Comp.(a) vstar <- vx (b) div_vstar, (c) p from intforces, (d) rhs_p
   Call Get_Laplacian_SW_fem 
   rhs_p  =  rhs_p*10000
   gstif2 = gstif2*10000 
   call solvegc2 (npoin, nelem_FS, nnode_FS, 1, npoin, 2, intmat_FS , &   ! p formulation from FS 2017   
                  iffix, fixed, gstif2, rhs_p, p, tolcg_MAIN) 
      !                                       ---                                        
   call Get_vn1_FS_SW

elseif (type_FS_problem.eq.3) then
   if (.not.allocated(intmat_FS))  then      ! Only 1st time. Note that moving 2d nodes will need update
      Call Get_P_BCs_SW_fem                  !   read iffix and fixed
   endif
   Call AuxVars_FS_SW_fem 

   Call Get_RHS_FS_SW_fem         !   Comp.(a) vstar <- vx (b) div_vstar, (c) p from intforces, (d) rhs_p
   Call Get_Laplacian_SW_fem 
   rhs_p  =  rhs_p*10000
   gstif2 = gstif2*10000 
   call solvegc2 (npoin, nelem_FS, nnode_FS, 1, npoin, 2, intmat_FS , &   ! p formulation from FS 2017   
                  iffix, fixed, gstif2, rhs_p, rho, tolcg_MAIN) 
      !                                       -----              
   call Get_vn1_FS_SW
endif

if (type_FS_problem.eq.1) then
   Call Get_P_BCs_SW_fem                     !   read iffix and fixed
   Call AuxVars_FS_SW_fem                    !   Obtains nelem_FS, nnode_FS, intmat_FS (nnode, nelem)  
   Call Get_RHS_FS_SW_fem                    !   Comp.(a) vstar <- vx (b) div_vstar, (c) p from intforces, (d) rhs_p
   Call Get_Laplacian_SW_fem
   rhs_p  =  rhs_p*10000
   gstif2 = gstif2*10000
   call Get_vn1_FS_SW                             !   Gets grad pn+1 -> vn+1 calling AVR solver 
elseif (type_FS_problem.eq.2) then
   if (.not.allocated(intmat_FS))  then      ! Only 1st time. Note that moving 2d nodes will need update
      Call Get_P_BCs_SW_fem                  !   read iffix and fixed
      Call AuxVars_FS_SW_fem 
   endif 
   Call Get_RHS_FS_SW_fem                    !   Comp.(a) vstar <- vx (b) div_vstar, (c) p from intforces, (d) rhs_p
   Call Get_Laplacian_SW_fem
   rhs_p  =  rhs_p*10000
   gstif2 = gstif2*10000
   call Get_vn1_FS_SW                             !   Gets grad pn+1 -> vn+1 calling AVR solver 
endif

CONTAINS


!......................................................................

      subroutine Initial_alloc_FS_SW_FEM

!......................................................................

 
implicit none

END subroutine Initial_alloc_FS_SW_FEM

!......................................................................

      subroutine Get_RHS_FS_SW_fem_5

!......................................................................

 
implicit none
 
real   (irk) Le, he, Pe 
real   (irk) kii, kij, kji, kjj
real   (irk) mii, mij, mji, mjj 
real   (irk) div_vstar, fact, g, facti, factj
real   (irk) FiI, FiJ
integer(ink) i, j, ielem 
real    (irk), allocatable:: X_tmp   (:,:)  !   
real    (irk), allocatable:: Vx_tmp  (:,:)  !  
real    (irk), allocatable:: rho_tmp   (:)  !

 ! -------------------  SW only --

if (.not.allocated(X_tmp)) then
    allocate ( X_tmp  (ndimn, npoin) ); X_tmp  = 0.0
    allocate ( Vx_tmp (ndimn, npoin) ); Vx_tmp  = 0.0  
    allocate ( rho_tmp       (npoin) ); rho_tmp = 0.0 
endif   

rhs_p  = 0.0
g      = const (1)

if (type_FS_problem.eq.5) then

     X_tmp  (:,1:npoin) = x  (:,1:npoin)  ! we need x = x0 to get p
     Vx_tmp (:,1:npoin) = vx (:,1:npoin)
     x (:,1:npoin)      = x0 (:,1:npoin)
     vx (:,1:npoin)     = vx0 (:,1:npoin)
     call IntForces_SW (1)                 !  = FS code called with x, vx at the beginning of timestep  NO P comp  ** MP 22nd nov
     x0  (:,1:npoin)  = x      (:,1:npoin) ! we recover original x= and x      vx0 (:,1:npoin)  = vx (:,1:npoin)
     vx0 (:,1:npoin)  = vx     (:,1:npoin)
     x   (:,1:npoin)  = X_tmp  (:,1:npoin)
     vx  (:,1:npoin)  = Vx_tmp (:,1:npoin) 
        
     if (ndimn.eq.1) then    !   ------------------  NDIMN = 1 -----------------

        do ielem = 1, nelem_FS
           i  = intmat_FS(1,ielem)
           j =  intmat_FS(2,ielem)
           Le = (x(1,i)-x(1,j))*(x(1,i)-x(1,j))
           Le = Le**0.5
           mii =    Le/3
           mjj =    mii
           mij =    mii/2
           mji =    mii/2 
           FiI      = (rho0(i)-rho(i))/(dt_sph*dt_sph*rho(i))
           FiJ      = (rho0(j)-rho(j))/(dt_sph*dt_sph*rho(j))
           fact     =   (Le/6.) 
           rhs_p (i) =  rhs_p (i) + fact * (2*FiI +   FiJ) 
           rhs_p (j) =  rhs_p (j) + fact * (  FiI + 2*FiJ)
        enddo
        
        if (if_right_wall.NE.1) then     
            i  = intmat_FS(1,nelem_FS)    ! ***** FS activated and changed for dam break dry
            j  = intmat_FS(2,nelem_FS) 
            Le = geome_FS (3,nelem_FS) 
            he = rho(j)            
            if (itimestep_sph.eq.1) then 
               rhs_p (j) = rhs_p(j) + (1.0/ he) * (0 - p(j))/Le 
            else
             rhs_p (j) = rhs_p(j) + (1.0/ he) * grad_p(1,j)   !  MP chgd 5 feb 2021
            endif   
        endif
        if(if_left_wall.NE.1) then
           i  = intmat_FS(1,1)     
           j  = intmat_FS(2,1)
           Le = geome_FS (3,1) 
           he = rho(j)            
           if (itimestep_sph.eq.1) then 
               rhs_p (i) = rhs_p(i) - (1.0/ he) * (p(i) - 0.0)/Le 
           else
               rhs_p (i) = rhs_p(i) - (1.0/ he) * grad_p(1,i)   !  MP chgd 5 feb 2021
           endif  
        endif 

     endif    !   -------------end NDIMN = 1 -----------------
endif

END subroutine Get_RHS_FS_SW_fem_5

!......................................................................

      subroutine Get_vn1_FS_Sw_5  

!......................................................................

implicit none

integer(ink)  ielem, inode, i,j, nprer , namat, niter

real   (irk), allocatable:: grad_Pel (:,:), mmatL(:), Rhs(:)
integer(ink), allocatable:: bprer(:,:) 
real   (irk) Le, hi, g

nprer = 0
if (.not.allocated(mmatL))  then
    allocate (mmatL(npoin), grad_Pel (ndimn, nelem_FS), rhs(npoin) )
    allocate (bprer(2,nprer))
endif
mmatl = 0.0 ; grad_Pel = 0.0; Rhs = 0.0
g     = const(1) 

if (ndimn.eq.1) then    !   ------------------  NDIMN = 1 -----------------

   do ielem = 1, nelem_FS
      i  = intmat_FS(1,ielem)
      j =  intmat_FS(2,ielem)
      Le = geome_FS (3, ielem)    
      if (type_FS_problem.eq.5) then
          grad_Pel (1, ielem) = geome_Fs (1,ielem)*  p(i) + geome_FS (2, ielem)*  p(j)  !  sum dNi/dx*pi
      endif
      rhs (i)  = rhs (i) + grad_Pel(1,ielem)*Le/2
      rhs (j)  = rhs (j) + grad_Pel(1,ielem)*Le/2
      mmatL(i) = mmatL(i) + Le/2   ! lumped mass matrix
      mmatL(j) = mmatL(j) + Le/2  
   enddo 

   nprer = 0
   bprer = 0
   namat = ndimn
   niter = 3

   call solve_AVR  (npoin,      nelem_FS, nnode_FS, ndimn, &   !  ** MP 22nd nov   grad h <- grad P 
                    namat,      ngeom_FS,    nprer, niter, &
                    intmat_FS , geome_FS,    bprer, mmatl, &
                    rhs,        grad_p )    

   call get_grad_tools (p, 0, grad_p)    ! SPH computation...perhaps better than fem + avr
   i = i
 
   if (type_FS_problem.eq.5) then
     
      do i  = 1, npoin
         dvstar(1,i) = -(dt_sph/rho(i))*grad_p(1,i)
         vx(:,i)     = vx (:,i) + dvstar(1,i)
         x (:,i) = x0(:,i) + 0.5*dt_sph*(vx0(:,i) + vx(:,i))
      enddo
      
      call sum_density_Tools          !   obtains h at n+1 alternative, use pn+1

      if (if_left_wall.eq.1) then     !   when walls, vn+1 = 0
           i = intmat_FS (1,1)        !   MP 2021 29 Jan
           vx(1,i) = 0.0 
      endif
      if (if_right_wall.eq.1) then
          j = intmat_FS (2, nelem_FS)
          vx(1,j) = 0.0
      endif
      i = intmat_FS (1,1)
      j = intmat_FS (2, nelem_FS)
 !     vx(1,i) = 0.0
 !     vx(1,j) = 0.0

   endif 

endif

END subroutine  Get_vn1_FS_Sw_5  

!......................................................................

      subroutine Get_P_BCs_SW_fem 

!......................................................................

 implicit none

integer(ink) ipres,idest, npres
integer(ink),allocatable:: ifpre(:)
real   (irk),allocatable:: presc(:) 
real   (irk) xxx

1005 format(a60)

read (dat_file,1005) text                               
write(chk_file,1005) text
read (dat_file,*)    npres     
write(chk_file,*)    npres
if (npres.gt.0) then
   if (.not.allocated(ifpre))  then
       allocate (ifpre(npres),presc(npres))
   endif
   read (dat_file,1005) text                                
   write(chk_file,1005) text
   read (dat_file,*)    (ifpre(ipres), ipres=1,npres)                                
   write(chk_file,*)    (ifpre(ipres), ipres=1,npres) 
   read (dat_file,1005) text                                
   write(chk_file,1005) text
   read (dat_file,*)    (presc(ipres), ipres=1,npres)                                
   write(chk_file,*)    (presc(ipres), ipres=1,npres)    
   do ipres=1,npres
      idest = ifpre(ipres)
      iffix(idest) = 1
      fixed(idest) = presc(ipres)
   enddo
   deallocate (ifpre, presc)
endif
!                          --- to impose dp/dn = 0 and vn+1=0 BCs
read (dat_file,1005) text                               
write(chk_file,1005) text
read (dat_file,*)    if_left_wall, if_right_wall, X_left_wall, X_right_wall     
write(chk_file,*)    if_left_wall, if_right_wall, X_left_wall, X_right_wall

END subroutine  Get_P_BCs_SW_fem


!......................................................................

      subroutine AuxVars_FS_SW_fem

!......................................................................

 
implicit none

integer(ink) ielem, i, j, if_sort
integer(ink), allocatable:: indexx (:)    ! Chuan sorting
real   (irk), allocatable:: tempxx(:)  
real   (irk) Le
if (ndimn.eq.1) then

    if (.not.allocated (indexx)) then
       allocate (indexx(npoin), tempxx(npoin))
    endif

    if_sort = 0           !  1 if nodes disordered
    do i = 1, npoin-1
       if (x(1,i).gt.x(1,i+1)) then
          if_sort = 1
          exit
       endif
    enddo
    if (if_sort.eq.1) then
       write (*,*) ' *********  disordered FS **** '
       tempxx(1:npoin) = x(1,1:npoin)
       indexx = (/(i,i=1,npoin)/)     ! FS Chuan 
       call sort (tempxx, indexx, npoin)
       do ielem = 1,nelem_fs
          intmat_FS (1,ielem) = indexx (ielem)
          intmat_FS (2,ielem) = indexx (ielem+1)
       enddo
    endif

    nnode_FS = 2
    nelem_FS = npoin-1 
    if (.not.allocated (intmat_FS)) then
       allocate (intmat_FS(nnode_FS,nelem_FS))
    endif
    intmat_FS = 0
    do ielem = 1,nelem_FS
      intmat_FS (1, ielem) = ielem
      intmat_FS (2, ielem) = ielem+1
    enddo
    ngeom_FS = nnode_FS*ndimn + 1
    if (.not.allocated (geome_FS)) then
       allocate (geome_FS (ngeom_FS, nelem_FS))
    endif 
    geome_FS = 0.0
    do ielem = 1, nelem_FS
       i  = intmat_FS (1,ielem)
       j  = intmat_FS (2,ielem)
       Le = (x(1,i)-x(1,j))**2.
       Le = Le**0.5
       geome_FS(1,ielem) = -1./Le     ! This are dN1/dx =-1/L and dN2 = 1/Le
       geome_FS(2,ielem) = - geome_FS(1,ielem)
       geome_FS(3,ielem) = Le
    enddo
endif

END subroutine AuxVars_FS_SW_fem

!......................................................................

      subroutine Get_RHS_FS_SW_fem

!......................................................................

 
implicit none
 
real   (irk) Le, he, Pe , FiI, FiJ
real   (irk) kii, kij, kji, kjj
real   (irk) mii, mij, mji, mjj 
real   (irk) div_vstar, fact, g, facti, factj
integer(ink) i, j, ielem 
real    (irk), allocatable:: X_tmp   (:,:)  !   
real    (irk), allocatable:: Vx_tmp  (:,:)  !  
real    (irk), allocatable:: rho_tmp   (:)  !

 ! -------------------  SW only --

if (.not.allocated(X_tmp)) then
    allocate ( X_tmp  (ndimn, npoin) ); X_tmp  = 0.0
    allocate ( Vx_tmp (ndimn, npoin) ); Vx_tmp  = 0.0  
    allocate ( rho_tmp       (npoin) ); rho_tmp = 0.0 
endif   

rhs_p  = 0.0
g      = const (1)
X_tmp  (:,1:npoin) = x  (:,1:npoin)
Vx_tmp (:,1:npoin) = vx (:,1:npoin)
!rho_tmp  (1:npoin) = rho  (1:npoin)

x (:,1:npoin)      = x0 (:,1:npoin)
vx (:,1:npoin)     = vx0 (:,1:npoin)
!rho(:)             = rho0(:)

call IntForces_SW (1)                    !  = FS code called with x, vx at the beginning of timestep  NO P comp  ** MP 22nd nov
vstar (:,1:npoin) = vx_TMP (:,1:npoin)   !  changed, is not vx0 but vx, and now vx is v_n

if (ndimn.eq.1) then    !   ------------------  NDIMN = 1 -----------------

   do ielem = 1, nelem_FS
      i  = intmat_FS(1,ielem)
      j =  intmat_FS(2,ielem)
      Le = (x(1,i)-x(1,j))*(x(1,i)-x(1,j))
      Le = Le**0.5
      mii =    Le/3
      mjj =    mii
      mij =    mii/2
      mji =    mii/2 
 
      if (type_FS_problem.eq.1) then              !    stedy stateheat
         rhs_p (i) =  rhs_p (i) + Le/2  
         rhs_p (j) =  rhs_p (j) + Le/2 

      elseif (type_FS_problem.eq.2) then          !    transient heat
         rhs_p (i) =  rhs_p (i) + mii*p(i) + mij*p(j)  
         rhs_p (j) =  rhs_p (j) + mji*p(i) + mjj*p(j)

      elseif (type_FS_problem.eq.3) then         ! -------------------  SW only      

!         div_vstar = (vx(1,j) - vx(1,i))/Le
!         fact      =   (Le/4)*(1. - dt_sph*div_vstar) 
!        fact   =   (Le/2) * ( rho0(i) + rho0(j) )/ ( rho(i) + rho(j) )  ! averages in elemns
!        facti  =   (Le/2) *(rho(i)/rho0(i))                             ! direct form
!        factj  =   (Le/2) *(rho(i)/rho0(i))                             ! direct form 
!        rhs_p (i) =  rhs_p (i) + fact   !   ***** facti      we use averaged form 
!        rhs_p (j) =  rhs_p (j) + fact   !   ***** factj
!
!   ------------   Fi = int( Ni.Nj dA ) * hi/h0i
!
!        rhs_p (i) =  rhs_p (i) + mii*rho(i)/rho0(i) + mij*rho(j)/rho0(j)   ! deactivated MP 27th
!        rhs_p (j) =  rhs_p (j) + mji*rho(i)/rho0(i) + mjj*rho(j)/rho0(j) 

!   ------------   aternative, consistent with M in Laplacian

        fact      = (Le/2) *( rho(i) + rho(j) ) / ( rho0(i) + rho0(j) )    ! check ?
        rhs_p (i) =  rhs_p (i) + fact
        rhs_p (j) =  rhs_p (j) + fact
      
      elseif (type_FS_problem.eq.4) then        ! 21 Jan 2021  chg MP to include p formulation from   SW only      

        div_vstar = (vx(1,j) - vx(1,i))/Le
        fact      =   (Le/4)*(1. - 2.*dt_sph*div_vstar) 
     
         rhs_p (i) =  rhs_p (i) + fact  
         rhs_p (j) =  rhs_p (j) + fact

   !     FiI = (rho(i)/rho0(j))
   !     FiJ = (rho(i)/rho0(j))
   !     rhs_p(i) = rhs_p(i) - (Le/4.) + (Le/6.)*(2*FiI +   FiJ)
   !     rhs_p(j) = rhs_p(j) - (Le/4.) + (Le/6.)*(  FiI + 2*FiJ)

      endif 

   enddo
   !           ********************   changed for lake mp 30 0ct  dp/dn = 0
   !        ---------           14th Oct in fem taken out BCs 1 and npoin out of the loop in elements

   if (type_FS_problem.eq.3) then
       if (if_right_wall.NE.1) then     
           i  = intmat_FS(1,nelem_FS)     
           j  = intmat_FS(2,nelem_FS) 
           Le = geome_FS (3,nelem_FS)
           he = (rho(i) + rho(j))/2 
           rhs_p (j) = rhs_p(j) + dt_sph * dt_sph * g* (rho(j) - rho(i))/Le  
       endif
       if(if_left_wall.NE.1) then
           i  = intmat_FS(1,1)     
           j  = intmat_FS(2,1)
           Le = geome_FS (3,1) 
           rhs_p (i) = rhs_p(i) - dt_sph*dt_sph * g * (rho(j) - rho(i))/Le   
      endif 
   elseif (type_FS_problem.eq.4) then
       if (if_right_wall.NE.1) then     
           i  = intmat_FS(1,nelem_FS)    ! ***** FS activated and changed for dam break dry
           j  = intmat_FS(2,nelem_FS) 
           Le = geome_FS (3,nelem_FS) 
           he = (rho(i) + rho(j))/2 
           if (itimestep_sph.eq.1) then 
               rhs_p (j) = rhs_p(j) + (dt_sph * dt_sph / he) * (0 - p(j))/Le 
           else
               rhs_p (j) = rhs_p(j) + (dt_sph * dt_sph / he) * grad_p(1,j)   !  MP chgd 1 feb 2021
           endif   
       endif
       if(if_left_wall.NE.1) then
           i  = intmat_FS(1,1)    ! ***** FS activated and changed for dam break dry
           j  = intmat_FS(2,1)
           Le = geome_FS (3,1) 
           he = (rho(i) + rho(j))/2  
           rhs_p (j) = rhs_p(j) - (dt_sph * dt_sph / he) * (p(i) - 0)/Le  
      endif 
   endif

   x  (:,1:npoin) =  X_Tmp (:,1:npoin)   !  recover star state
   vx (:,1:npoin) = Vx_Tmp (:,1:npoin)
!   rho  (1:npoin) = rho_Tmp  (1:npoin)

endif    !   -------------end NDIMN = 1 -----------------

END subroutine Get_RHS_FS_SW_fem

!......................................................................

      subroutine Get_Laplacian_SW_fem 

!......................................................................


!      ------   Obtains laplacian matrix transient heat

implicit none
 
real   (irk) Le, he, factK, factM, g, Pe  
real   (irk) kii, kij, kji, kjj
real   (irk) mii, mij, mji, mjj 
integer(ink) i, j, ielem  
 
if (allocated (gstif2)) then
   deallocate (gstif2) 
   allocate ( gstif2 (nnode_FS,nnode_FS,nelem_FS))
else
   allocate ( gstif2 (nnode_FS,nnode_FS,nelem_FS))
endif
 
gstif2 = 0.0  
g      = const(1)       !     

if (ndimn.eq.1) then
   do ielem = 1, nelem_FS

      i  = intmat_FS(1,ielem)
      j =  intmat_FS(2,ielem)
      Le = (x(1,i)-x(1,j))*(x(1,i)-x(1,j))
      Le = Le**0.5    
      kij =  - 1/Le
      kji =  - 1/Le
      kii =    1/Le
      kjj =    1/Le
      mii =    Le/3
      mjj =    mii
      mij =    mii/2
      mji =    mii/2
 
      if (type_FS_problem.eq.1) then              !    stedy stateheat
          gstif2(1,1,ielem) = gstif2(1,1,ielem) - kii 
          gstif2(1,2,ielem) = gstif2(1,2,ielem) - kij 
          gstif2(2,1,ielem) = gstif2(2,1,ielem) - kji  
          gstif2(2,2,ielem) = gstif2(2,2,ielem) - kjj
      elseif (type_FS_problem.eq.2) then          !    transient heat
          gstif2(1,1,ielem) = gstif2(1,1,ielem) + dt_sph * kii + mii
          gstif2(1,2,ielem) = gstif2(1,2,ielem) + dt_sph * kij + mij 
          gstif2(2,1,ielem) = gstif2(2,1,ielem) + dt_sph * kji + mji 
          gstif2(2,2,ielem) = gstif2(2,2,ielem) + dt_sph * kjj + mjj
      elseif (type_FS_problem.eq.3) then          !    FS SW
!          Pe    = (p(i) + p(j)) / 2.
          he    = (rho(i) + rho(j) ) /2.
          factK =  (g*dt_sph**2 )
          factM =  1/he
          gstif2(1,1,ielem) = gstif2(1,1,ielem) + factK * kii  + factM * mii
          gstif2(1,2,ielem) = gstif2(1,2,ielem) + factK * kij  + factM * mij 
          gstif2(2,1,ielem) = gstif2(2,1,ielem) + factK * kji  + factM * mji 
          gstif2(2,2,ielem) = gstif2(2,2,ielem) + factK * kjj  + factM * mjj

      elseif (type_FS_problem.eq.4) then          !    FS SW  p based 21 Jan 2012 MP chgd
          he    = (rho(i) + rho(j) ) /2.
          Pe    = (p(i) + p(j)) / 2.
          factK =  (dt_sph**2 )/he 
          factM =  1./(2.*Pe)
          gstif2(1,1,ielem) = gstif2(1,1,ielem) + factK * kii  + factM * mii
          gstif2(1,2,ielem) = gstif2(1,2,ielem) + factK * kij  + factM * mij 
          gstif2(2,1,ielem) = gstif2(2,1,ielem) + factK * kji  + factM * mji 
          gstif2(2,2,ielem) = gstif2(2,2,ielem) + factK * kjj  + factM * mjj

      elseif (type_FS_problem.eq.5) then          !    FS SW  Shao
          he    = (rho(i) + rho(j) ) /2.
          factK =  1.00/he 
          gstif2(1,1,ielem) = gstif2(1,1,ielem) + factK * kii   
          gstif2(1,2,ielem) = gstif2(1,2,ielem) + factK * kij   
          gstif2(2,1,ielem) = gstif2(2,1,ielem) + factK * kji   
          gstif2(2,2,ielem) = gstif2(2,2,ielem) + factK * kjj   
      endif

   enddo
endif

END subroutine  Get_Laplacian_SW_fem

!......................................................................

      subroutine Get_vn1_FS_Sw 

!......................................................................

 
implicit none

integer(ink)  ielem, inode, i,j, nprer , namat, niter

real   (irk), allocatable:: grad_Pel (:,:), mmatL(:), Rhs(:)
integer(ink), allocatable:: bprer(:,:) 
real   (irk) Le, hi, g

nprer = 0
if (.not.allocated(mmatL))  then
    allocate (mmatL(npoin), grad_Pel (ndimn, nelem_FS), rhs(npoin) )
    allocate (bprer(2,nprer))
endif
mmatl = 0.0 ; grad_Pel = 0.0; Rhs = 0.0
g     = const(1) 

if (ndimn.eq.1) then    !   ------------------  NDIMN = 1 -----------------
   do ielem = 1, nelem_FS
      i  = intmat_FS(1,ielem)
      j =  intmat_FS(2,ielem)
      Le = geome_FS (3, ielem)   !  ** MP 22nd nov   grad h -> grad Pel
      if (type_FS_problem.eq.3) then
          grad_Pel (1, ielem) = geome_Fs (1,ielem)*rho(i) + geome_FS (2, ielem)*rho(j)  !  sum dNi/dx*hi
      elseif (type_FS_problem.eq.4) then
          grad_Pel (1, ielem) = geome_Fs (1,ielem)*  p(i) + geome_FS (2, ielem)*  p(j)  !  sum dNi/dx*pi
      endif
      rhs (i)  = rhs (i) + grad_Pel(1,ielem)*Le/2
      rhs (j)  = rhs (j) + grad_Pel(1,ielem)*Le/2
      mmatL(i) = mmatL(i) + Le/2   ! lumped mass matrix
      mmatL(j) = mmatL(j) + Le/2  
   enddo 

   nprer = 0
   bprer = 0
   namat = ndimn
   niter = 3

   call solve_AVR  (npoin,      nelem_FS, nnode_FS, ndimn, &   !  ** MP 22nd nov   grad h <- grad P 
                    namat,      ngeom_FS,    nprer, niter, &
                    intmat_FS , geome_FS,    bprer, mmatl, &
                    rhs,        grad_p )    

   call get_grad_tools (p, 0, grad_p)    ! SPH computation...perhaps better than fem + avr
   i = i

   if (type_FS_problem.eq.3) then
      do i  = 1, npoin
!                          hi = (2*p(i)/g)**0.5              !   ****  MP 7th nov 2017
!                          rho(i) = hi
         vx(:,i) = vstar(:,i) - (g * dt_sph )*grad_p (:,i)    !  ** MP 22nd nov   grad h -> grad Pel
         x (:,i) = x0   (:,i) + 0.5*dt_sph*(vx0(:,i) + vx(:,i))  
!                       vx(:,i) = vstar(:,i) - (dt_sph/hi)*grad_p(:,i)
      enddo
   
      if (if_left_wall.eq.1) then    !   when walls, vn+1 = 0
          i = intmat_FS (1,1)
          vx(1,i) = 0.0 
      endif
      if (if_right_wall.eq.1) then
          j = intmat_FS (2, nelem_FS)
          vx(1,j) = 0.0
      endif
  
   elseif (type_FS_problem.eq.4) then
     
      do i  = 1, npoin
         hi = (2*p(i)/g)**0.5              !   ****  MP 7th nov 2017
         rho(i) = hi
         vx(:,i) = vstar(:,i) - (dt_sph/hi)*grad_p(:,i)
         x (:,i) = x0(:,i) + 0.5*dt_sph*(vx0(:,i) + vx(:,i))
      enddo
      
      if (if_left_wall.eq.1) then    !   when walls, vn+1 = 0
           i = intmat_FS (1,1)        ! MP 2021 29 Jan
           vx(1,i) = 0.0 
      endif
      if (if_right_wall.eq.1) then
          j = intmat_FS (2, nelem_FS)
          vx(1,j) = 0.0
      endif
      i = intmat_FS (1,1)
      j = intmat_FS (2, nelem_FS)
 !     vx(1,i) = 0.0
 !     vx(1,j) = 0.0

   endif 

endif

END subroutine  Get_vn1_FS_Sw   

END subroutine Second_Step_FS_SW_FEM   ! ---------------------------------


!----------------------------------------------------------------------

      subroutine Second_Step_FS_SW_SFM 

!----------------------------------------------------------------------
 
implicit none

real   (irk), SAVE:: xminbg, xmaxbg, deltxbg,  Lbg
!real   (irk)  xminbg, xmaxbg, deltxbg,  Lbg

real    (irk), allocatable:: X_tmp   (:,:)  !  
real    (irk), allocatable:: Vx_tmp  (:,:)  !  
real    (irk), allocatable:: rho_tmp   (:)  !

if (.not.allocated(vn_plus1)) then
   allocate (vn_plus1 (ndimn,npoin)) ; vn_plus1 = 0.0 
   allocate (grad_p   (ndimn,npoin)) ; grad_p   = 0.0 
   allocate (vstar    (ndimn,npoin)) ; vstar    = 0.0 
   allocate (dvstar   (ndimn,npoin)) ; dvstar   = 0.0 
else
   vn_plus1 = 0.0; grad_p   = 0.0; vstar    = 0.0 ; dvstar   = 0.0
endif

if (type_FS_problem.eq.1) then   

    if (.not.allocated(X_tmp)) then
       allocate ( X_tmp  (ndimn, npoin) ); X_tmp  = 0.0
       allocate ( Vx_tmp (ndimn, npoin) ); Vx_tmp  = 0.0  
       allocate ( rho_tmp       (npoin) ); rho_tmp = 0.0 
    endif

!   Call AuxVars2_FS_SW_SFM        !  ***  to test, we are using the test sub 2 
   Call AuxVars_FS_SW_SFM      
   
   If (itimestep_SPH.eq.1) then   !  after computing npoin_FS...
!!      Call Init_P_BCs_SW_SFM 
   endif

!!   Call Get_RHS_FS_SW_SFM         !   Comp.(a) vstar <- vx (b) div_vstar, (c) p from intforces, (d) rhs_p
   
!!   Call Get_Laplacian_SW_SFM 
!!      rhs_p  =  rhs_p*10000
!!      gstif2 = gstif2*10000 
  
!!  call solvegc2 (npoin_FS, nelem_FS, nnode_FS, 1, npoin_FS, 2, intmat_FS , &   ! p formulation from FS 2017   
!!                 iffix, fixed, gstif2, rhs_p, p_FS, tolcg_MAIN) 
      !                                       ---                                        
!!  call Get_vn1_FS_SFM


endif

CONTAINS

!......................................................................

     subroutine AuxVars_FS_SW_SFM

!......................................................................
 
implicit none

real   (irk) xx, hsml_avg , dist, distxbg, distbg, distx2
integer(ink) ipoin, npaux, nelaux, ipoibgx, inx0, inx1, inx, dinx  
integer(ink) i,j, ie, ip, inode, ielem, ipoibg, ipoifs, ielbg, ielfs

!  --------------  Task 1  -- temporarily recover X* to get p  

X_tmp  (:,1:npoin ) = x  (:,1:npoin )
Vx_tmp (:,1:npoin ) = vx (:,1:npoin )
!rho_tmp  (1:npoin ) = rho  (1:npoin )

x (:,1:npoin)      = x0 (:,1:npoin)
vx (:,1:npoin)     = vx0 (:,1:npoin)
!rho(:)             = rho0(:)
call IntForces_SW (1)    !  in this way we get p(i) at time n using rho_n

! back
x0  (:,1:npoin) = x  (:,1:npoin ); x  (:,1:npoin ) = X_tmp  (:,1:npoin )
Vx0 (:,1:npoin) = Vx (:,1:npoin ); Vx (:,1:npoin ) = VX_tmp (:,1:npoin )

!      ------------  Task 2: set background grid  ----------------------------------------

xminbg  = 1.e10;  xmaxbg  = -1.e10; hsml_avg = 0
npaux  = 0
do ipoin = 1,npoin
   if (If_Out_Domain(ipoin).eq.1)  CYCLE
   npaux = npaux + 1
   if (x(1,ipoin).lt.xminbg) xminbg = x(1,ipoin)
   if (x(1,ipoin).gt.xmaxbg) xmaxbg = x(1,ipoin)
   hsml_avg = hsml_avg + hsml(ipoin)
enddo  
hsml_avg = 0.5*hsml_avg / npaux

Lbg     = (xmaxbg-xminbg)
xmaxbg = xmaxbg + 0.005*Lbg     !  expand BG mesh to avoid border nodes
xminbg = xminbg - 0.005*Lbg
Lbg    = xmaxbg -  xminbg

!    CHK TEST

do ip = 1, npoin
   vx(1,ip) = (x (1,ip) - xminbg)/Lbg
enddo

  
npoin_bg = (Lbg/hsml_avg) + 2 ; nelem_bg = npoin_bg - 1
deltxbg  = Lbg/(npoin_bg-1)

nnode_bg = 2
ngeom_bg =  nnode_bg * ndimn + 1
nvars_bg = 6       ! -----  n points linked from sph; 2 rho, 3 p, 4 u*, 5 grad P  6 sum(1/Rij) 

IF (.NOT.allocated (coord_BG)) then
   mpoin_bg = npoin_bg; melem_bg = nelem_bg
   allocate (coord_bg (ndimn   , mpoin_bg))  ; coord_bg     = 0.0
   allocate (intmat_bg(nnode_bg, mpoin_bg))  ; intmat_bg    = 0
   allocate (geome_bg (ngeom_bg, melem_bg))  ; geome_bg     = 0.0
   allocate (act_nodes_bg  (mpoin_bg ) )     ; act_nodes_bg = 0 
   allocate (act_elems_bg  ( melem_bg) )     ; act_elems_bg = 0 
   allocate ( vars_bg(nvars_bg, mpoin_bg))   ; vars_bg      = 0.0
ELSE
   if (npoin_bg.GT.mpoin_bg) then
      deallocate ( coord_bg )
      deallocate ( intmat_bg )
      deallocate ( geome_bg )
      deallocate ( act_nodes_bg ) 
      deallocate ( act_elems_bg )
      deallocate ( vars_bg )
      mpoin_bg = npoin_bg; melem_bg = nelem_bg
      allocate (coord_bg (ndimn   , mpoin_bg))  ; coord_bg     = 0.0
      allocate (intmat_bg(nnode_bg, mpoin_bg))  ; intmat_bg    = 0
      allocate (geome_bg (ngeom_bg, melem_bg))  ; geome_bg     = 0.0
      allocate (act_nodes_bg  (mpoin_bg ) )     ; act_nodes_bg = 0 
      allocate (act_elems_bg  ( melem_bg) )     ; act_elems_bg = 0 
      allocate ( vars_bg(nvars_bg, mpoin_bg))   ; vars_bg      = 0.0
   endif 
ENDIF

do ipoibgx = 1, npoin_bg
   coord_bg  (1, ipoibgx ) = xminbg + (ipoibgx -1)*deltxbg
enddo
do ie  = 1, nelem_bg
intmat_bg (1, ie ) = ie 
intmat_bg (2, ie ) = ie +1 
enddo

geome_bg = 0.0
do ie = 1, nelem_bg
   i  = intmat_bg (1,ie)
   j  = intmat_bg (2,ie)
   geome_bg(1,ie) = -1./deltxbg     ! This are dN1/dx =-1/L and dN2 = 1/Le
   geome_bg(2,ie) = - geome_bg(1,ie)
   geome_bg(3,ie) = deltxbg
enddo
 
!      ------------  Task 3: transfer vars to BG grids

vars_bg = 0.0
distbg   = deltxbg
DO ipoin = 1,npoin 
   xx = x(1,ipoin)    
   ipoibgx    = (xx-xminbg)/deltxbg + 1 
   dinx       = hsml(ipoin)/deltxbg 
   inx0       = max(1,ipoibgx-dinx)
   inx1       = min(npoin_bg ,ipoibgx + dinx + 1) 
   do inx  = inx0 , inx1
      ipoibgx = inx
      distx2  = ( x(1,ipoin) - coord_bg(1,ipoibgx) )
      distx2  = distx2*distx2
      dist    = max( 1.e-6, (distx2)**0.5 )
      dist    = dist/distbg                                              ! distance is normalized   
      vars_bg (1,ipoibgx) = vars_bg (1,ipoibgx) + 1.00                   !  
      vars_bg (2,ipoibgx) = vars_bg (2,ipoibgx) + rho     (ipoin)/dist   !   
      vars_bg (3,ipoibgx) = vars_bg (3,ipoibgx) + vx     (1,ipoin)/dist  ! vx is vstar after RK4 
      vars_bg (4,ipoibgx) = vars_bg (4,ipoibgx) + p       (ipoin)/dist   !  
      vars_bg (5,ipoibgx) = vars_bg (5,ipoibgx) + grad_p (1,ipoin)/dist  ! 
      vars_bg (6,ipoibgx) = vars_bg (6,ipoibgx) + 1.00/dist              !        
   enddo
ENDDO

DO ip  = 1, npoin_bg
   if (vars_bg(1,ip).ge.1.e-3) then
       vars_bg(2:5,ip) = vars_bg(2:5,ip)/vars_bg(6,ip)         ! Vstar stored in vars_bg(2,:)
   endif
ENDDO    

!  -----  Task 4.1 active nodes and elements in BG grid

act_nodes_bg = 0  !      ------  array with 0 (nde not active) or number of order  
npaux = 0         !              nr of active nodes 
do ip  = 1, npoin_bg
   if (Vars_bg(1,ip).GT.1.e-3) then
      npaux = npaux + 1
      act_nodes_bg (ip) =  npaux
   endif
enddo 

act_elems_bg = 0!      ------ array with 0 (elem not active) or number of order 
do ie = 1, nelem_bg
   do inode = 1, nnode_bg
      ip    = intmat_bg (inode, ie)
      if (act_nodes_bg(ip).eq.1) then
         nelaux = nelaux + 1
         act_elems_bg(ie) = 1
      endif
      if (act_elems_bg(ie).eq.1) EXIT
   enddo
enddo
nelaux = 0
do ie  = 1, nelem_bg
   if (act_elems_bg(ip).eq.1) then
       nelaux = nelaux + 1
       act_elems_bg(ie) = nelaux   
   endif
enddo 
                            
npoin_FS = npaux; nelem_FS = nelaux   ! we have the number of active nodes and elements in BG

!  -----  Task 4.2  Get FS grid info

Nnode_FS = nnode_bg
ngeom_FS = ngeom_bg
nvars_FS = nvars_bg            ! -----  n points linked from sph; 2 rho, 3 p, 4 u*, 5 grad P  6 sum(1/Rij) 

IF (.NOT.allocated (coord_FS)) then
   allocate (coord_FS (ndimn   , npoin_FS))   ; coord_FS     = 0.0
   allocate (intmat_FS(nnode_FS, npoin_FS))   ; intmat_FS    = 0
   allocate (geome_FS (ngeom_FS, nelem_FS))   ; geome_FS     = 0.0
   allocate (act_nodes_FS  (npoin_FS ) )      ; act_nodes_FS = 0 
   allocate (act_elems_FS  ( nelem_FS) )      ; act_nodes_FS = 0
   allocate ( vars_FS(nvars_FS, npoin_FS))    ; vars_FS      = 0.0
   
   allocate ( gstif2    (nnode_FS, nnode_FS, nelem_FS) )   ; gstif2    = 0.0
   allocate (iffix      (npoin_FS) )                   ; iffix     = 0
   allocate (fixed      (npoin_FS) )                   ; fixed     = 0.0
   allocate ( rhs_p     (npoin_FS) )                   ; rhs_p     = 0.0
   allocate ( rho_FS    (npoin_FS) )                   ; rho_FS    = 0.0
   allocate ( p_FS      (npoin_FS) )                   ; p_FS      = 0.0
   allocate ( vstar_FS  (ndimn, npoin_FS) )            ; vstar_FS  = 0.0
   allocate ( grad_p_FS (ndimn, npoin_FS) )            ; grad_p_FS = 0.0
   mpoin_FS = npoin_FS; melem_FS = nelem_FS

ELSE

   if (npoin_FS.GT.mpoin_FS) then
      deallocate ( coord_FS, intmat_FS, geome_FS, act_nodes_FS, act_elems_FS, vars_FS  )
      deallocate ( rhs_p,  rho_FS, p_FS, vstar_FS, grad_p_FS, gstif2, iffix, fixed  )
      allocate (coord_FS (ndimn   , npoin_FS))   ; coord_FS     = 0.0
      allocate (intmat_FS(nnode_FS, npoin_FS))   ; intmat_FS    = 0
      allocate (geome_FS (ngeom_FS, nelem_FS))   ; geome_FS     = 0.0
      allocate (act_nodes_FS  (npoin_FS ) )      ; act_nodes_FS = 0 
      allocate (act_elems_FS  ( nelem_FS) )      ; act_nodes_FS = 00
      allocate ( vars_FS(nvars_FS, npoin_FS))    ; vars_FS      = 0.0

      allocate (iffix         (npoin_FS)  )                   ; iffix     = 0
      allocate (fixed         (npoin_FS)  )                   ; fixed     = 0.0
      allocate ( rhs_p     (npoin_FS) )                       ; rhs_p     = 0.0
      allocate ( rho_FS    (npoin_FS) )                       ; rho_FS    = 0.0
      allocate ( p_FS      (npoin_FS) )                       ; p_FS      = 0.0
      allocate ( vstar_FS  (ndimn, npoin_FS) )                ; vstar_FS  = 0.0
      allocate ( grad_p_FS (ndimn, npoin_FS) )                ; grad_p_FS = 0.0
      allocate ( gstif2    (nnode_FS, nnode_FS, nelem_FS) )   ; gstif2    = 0.0
      mpoin_FS = npoin_FS; melem_FS = nelem_FS
   endif 

ENDIF

!  -----  Task 4.3  ------  obtain cross ref arrays from FS to BG  act_nodes_FS and act_elems_fs

do ipoibg = 1, npoin_bg
   if ( act_nodes_bg(ipoibg).GT.0) then
      ipoifs = act_nodes_bg(ipoibg)
      act_nodes_FS (ipoifs) = ipoibg
   endif
enddo  

do ielbg = 1, nelem_bg
   if ( act_elems_bg(ielbg).GT.0) then
      ielfs = act_elems_bg(ielbg)
      act_elems_FS (ielfs) = ielbg
   endif
enddo  
 
!  -----  Task 4.4    transfer info from BG to FS  

do ipoifs = 1, npoin_FS
   ipoibg = act_nodes_FS (ipoifs)
   coord_FS(:, ipoifs) = coord_bg (:, ipoibg) 
   Vars_FS (:, ipoifs) = Vars_bg  (:, ipoibg)  
enddo
do ielfs = 1, nelem_FS
   ielbg = act_elems_FS (ielfs)
   intmat_FS (:, ielfs)  = intmat_bg (:, ielbg) 
   geome_FS  (:, ielfs)  = geome_bg  (:, ielbg)  
enddo

END subroutine AuxVars_FS_SW_SFM


!......................................................................

     subroutine AuxVars2_FS_SW_SFM    !   **** test with fixed grid 1D

!......................................................................
 
implicit none

real   (irk) xx, hsml_avg , dist, distxbg, distbg, distx2
integer(ink) ipoin, npaux, ipoibgx, inx0, inx1, inx, dinx, ipoibg 
integer(ink) i,j, ie, ip

!  ----------------------   get aux vars 

X_tmp  (:,1:npoin ) = x  (:,1:npoin )
Vx_tmp (:,1:npoin ) = vx (:,1:npoin )
!rho_tmp  (1:npoin ) = rho  (1:npoin )

x (:,1:npoin)      = x0 (:,1:npoin)
vx (:,1:npoin)     = vx0 (:,1:npoin)
!rho(:)             = rho0(:)
call IntForces_SW (1)    !  in this way we get p(i) at time n using rho_n

x0  (:,1:npoin)  = x      (:,1:npoin) ! we recover original x= and x      vx0 (:,1:npoin)  = vx (:,1:npoin)
vx0 (:,1:npoin)  = vx     (:,1:npoin)
x   (:,1:npoin)  = X_tmp  (:,1:npoin)
vx  (:,1:npoin)  = Vx_tmp (:,1:npoin) 

!      ------------  Task 1: set backgroud grid  ----------------------------------------


if (.NOT.allocated (coord_BG)) then

   xminbg  = 1.e10;  xmaxbg  = -1.e10; hsml_avg = 0
   npaux  = 0
   do ipoin = 1,npoin
      if (If_Out_Domain(ipoin).eq.1)  CYCLE
      npaux = npaux + 1
      if (x(1,ipoin).lt.xminbg) xminbg = x(1,ipoin)
      if (x(1,ipoin).gt.xmaxbg) xmaxbg = x(1,ipoin)
      hsml_avg = hsml_avg + hsml(ipoin)
   enddo  
   hsml_avg = 0.5*hsml_avg / npaux
   ! xminbg = xminbg - hsml_avg/2. ; xmaxbg = xmaxbg + hsml_avg/2.
   
   Lbg      = (xmaxbg-xminbg)
   npoin_bg = npoigx; nelem_bg = npoin_bg - 1
   deltxbg  = Lbg/(npoin_bg-1) 
    
!  npoin_bg = (Lbg/hsml_avg) + 2 ; nelem_bg = npoin_bg - 1
!  deltxbg  = Lbg/(npoin_bg-1)  
!   npoin_bg = npoigx; nelem_bg = npoin_bg - 1             !   use topo grid, which is fixed
!   Lbg      = xmaxg - xming
!   deltxbg  = deltxg
   
   nnode_bg = 2
   ngeom_bg =  nnode_bg * ndimn + 1
   nvars_bg = 6       ! -----  n points linked from sph; 2 rho, 3 p, 4 u*, 5 grad P  6 sum(1/Rij) 

   allocate (coord_bg (ndimn   , npoin_bg))  ; coord_bg     = 0.0
   allocate (intmat_bg(nnode_bg, npoin_bg))  ; intmat_bg    = 0
   allocate (geome_bg (ngeom_bg, nelem_bg))  ; geome_bg     = 0.0
   allocate (act_nodes_bg  (npoin_bg ) )     ; act_nodes_bg = 0 
   allocate (act_elems_bg  ( nelem_bg) )     ; act_elems_bg = 0 
   allocate ( vars_bg(nvars_bg, npoin_bg))   ; vars_bg      = 0.0
   mpoin_bg = npoin_bg; melem_bg = nelem_bg
    
   do ipoibgx = 1, npoin_bg
      coord_bg  (1, ipoibgx ) = xminbg + (ipoibgx -1)*deltxbg
   enddo
   do ie  = 1, nelem_bg
   intmat_bg (1, ie ) = ie 
   intmat_bg (2, ie ) = ie +1 
   enddo
   
   geome_bg = 0.0
   do ie = 1, nelem_bg
      i  = intmat_bg (1,ie)
      j  = intmat_bg (2,ie)
      geome_bg(1,ie) = -1./deltxbg     ! This are dN1/dx =-1/L and dN2 = 1/Le
      geome_bg(2,ie) = - geome_bg(1,ie)
      geome_bg(3,ie) = deltxbg
   enddo

endif


!  -----  Task 2 FS concentrated grid... provisional transfer 

npoin_Fs = npoin_bg 
nelem_FS = nelem_bg
Nnode_FS = nnode_bg
ngeom_FS = ngeom_bg
nvars_FS = nvars_bg            ! -----  n points linked from sph; 2 rho, 3 p, 4 u*, 5 grad P  6 sum(1/Rij) 

IF (.NOT.allocated (coord_FS)) then
   allocate (coord_FS (ndimn   , npoin_FS))   ; coord_FS     = 0.0
   allocate (intmat_FS(nnode_FS, npoin_FS))   ; intmat_FS    = 0
   allocate (geome_FS (ngeom_FS, nelem_FS))   ; geome_FS     = 0.0
   allocate (act_nodes_FS  (npoin_FS ) )      ; act_nodes_FS = 0 
   allocate (act_elems_FS  ( nelem_FS) )      ; act_nodes_FS = 0
   allocate ( vars_FS(nvars_FS, npoin_FS))    ; vars_FS      = 0.0
   
   allocate (iffix         (npoin_FS)  )                   ; iffix     = 0
   allocate (fixed         (npoin_FS)  )                   ; fixed     = 0.0
   allocate ( rhs_p     (npoin_FS) )                       ; rhs_p     = 0.0
   allocate ( rho_FS    (npoin_FS) )                       ; rho_FS    = 0.0
   allocate ( p_FS      (npoin_FS) )                       ; p_FS      = 0.0
   allocate ( vstar_FS  (ndimn, npoin_FS) )                ; vstar_FS  = 0.0
   allocate ( grad_p_FS (ndimn, npoin_FS) )                ; grad_p_FS = 0.0
   allocate ( gstif2    (nnode_FS, nnode_FS, nelem_FS) )   ; gstif2    = 0.0
   mpoin_FS = npoin_FS; melem_FS = nelem_FS
ENDIF

coord_FS     = coord_bg
intmat_FS    = intmat_bg
geome_FS     = geome_bg
act_nodes_FS = act_nodes_bg
act_elems_FS = act_elems_bg
 
!      ------------  Task 3: transfer vars to FS grid and BG grids

vars_bg = 0.0
distbg   = Lbg/(npoin_bg-1) 
DO ipoin = 1,npoin 
   xx = x(1,ipoin)    
   ipoibgx    = (xx-xminbg)/deltxbg + 1 
   dinx       = hsml(ipoin)/deltxbg 
   inx0       = max(1,ipoibgx-dinx)
   inx1       = min(npoin_bg ,ipoibgx + dinx + 1) 
   do inx  = inx0 , inx1
      ipoibgx = inx
      distx2  = ( x(1,ipoin) - coord_bg(1,ipoibgx) )
      distx2  = distx2*distx2
      dist    = max( 1.e-6, (distx2)**0.5 )
      dist    = dist/distbg                                              ! distance is normalized   
      vars_bg (1,ipoibgx) = vars_bg (1,ipoibgx) + 1.00                   !  
      vars_bg (2,ipoibgx) = vars_bg (2,ipoibgx) + rho     (ipoin)/dist   !   
      vars_bg (3,ipoibgx) = vars_bg (3,ipoibgx) + vx     (1,ipoin)/dist  ! vx is vstar after RK4 
      vars_bg (4,ipoibgx) = vars_bg (4,ipoibgx) + p       (ipoin)/dist   !  
      vars_bg (5,ipoibgx) = vars_bg (5,ipoibgx) + grad_p (1,ipoin)/dist  ! 
      vars_bg (6,ipoibgx) = vars_bg (6,ipoibgx) + 1.00/dist              !        
   enddo
ENDDO

DO ip  = 1, npoin_bg
   if (vars_bg(1,ip).ge.1.e-3) then
       vars_bg(2:5,ip) = vars_bg(2:5,ip)/vars_bg(6,ip)         ! Vstar stored in vars_bg(2,:)
   endif
ENDDO    

vars_FS = vars_bg

END subroutine AuxVars2_FS_SW_SFM


!......................................................................

      subroutine Init_P_BCs_SW_SFM 

!......................................................................

implicit none

integer(ink) ipres,idest, npres
integer(ink),allocatable:: ifpre(:)
real   (irk),allocatable:: presc(:) 
real   (irk) xxx

1005 format(a60)

read (dat_file,1005) text                               
write(chk_file,1005) text
read (dat_file,*)    npres     
write(chk_file,*)    npres
!                          --- to impose dp/dn = 0 and vn+1=0 BCs
read (dat_file,1005) text                               
write(chk_file,1005) text
read (dat_file,*)    if_left_wall, if_right_wall, X_left_wall, X_right_wall     
write(chk_file,*)    if_left_wall, if_right_wall, X_left_wall, X_right_wall

END subroutine  Init_P_BCs_SW_SFM 


!......................................................................

      subroutine Get_RHS_FS_SW_SFM

!......................................................................

 
implicit none
 
real   (irk) Le, he, Pe , FiI, FiJ
real   (irk) div_vstar, fact, facti, factj
integer(ink) i, j, ielem 

 ! -------------------  SW only --

rhs_p  = 0.0

if (ndimn.eq.1) then    !   ------------------  NDIMN = 1 -----------------

   do ielem = 1, nelem_FS
      i  = intmat_FS(1,ielem)
      j =  intmat_FS(2,ielem)
      Le = (coord_FS(1,i)-coord_FS(1,j))*(coord_FS(1,i)-coord_FS(1,j))
      Le = Le**0.5

      if (type_FS_problem.eq.1) then        ! 21 Jan 2021  chg MP to include p formulation from   SW only      
         div_vstar = ( vars_FS (3,i) - vars_FS (3,j))/Le
         fact      =   (Le/4)*(1. - 2.*dt_sph*div_vstar)  
         rhs_p (i) =  rhs_p (i) + fact  
         rhs_p (j) =  rhs_p (j) + fact 
      endif 

   enddo
   
   !        ---------           

   if (type_FS_problem.eq.1) then
       if (if_right_wall.NE.1) then     
           i  = intmat_FS(1,nelem_FS)    ! ***** FS activated and changed for dam break dry
           j  = intmat_FS(2,nelem_FS) 
           Le = geome_FS (3,nelem_FS) 
           he = (vars_FS (2,i)+vars_FS (2,j))/2.                                   ! (rho(i) + rho(j))/2 
           if (itimestep_sph.eq.1) then 
               rhs_p (j) = rhs_p(j) + (dt_sph * dt_sph / he) * (0. - p_FS(j))/Le 
           else
               rhs_p (j) = rhs_p(j) + (dt_sph * dt_sph / he) * vars_FS (5,j)       !  grad_p(1,j)  (i)----(j)| 
           endif   
       endif
       if(if_left_wall.NE.1) then
           i  = intmat_FS(1,1)    ! ***** FS activated and changed for dam break dry
           j  = intmat_FS(2,1)
           Le = geome_FS (3,1) 
           he = (vars_FS (2,i)+vars_FS (2,j))/2.                                  !  (rho(i) + rho(j))/2        
           if (itimestep_sph.eq.1) then 
               rhs_p (i) = rhs_p(i) + (dt_sph * dt_sph / he) * (p_FS(i) - 0.)/Le 
           else
               rhs_p (i) = rhs_p(i) + (dt_sph * dt_sph / he) * (- vars_FS (5,i) ) !  -grad_p(1,i)  |(i)----(j) 
           endif  
       endif 

   endif

endif    !   -------------end NDIMN = 1 -----------------

END subroutine Get_RHS_FS_SW_SFM

!......................................................................

      subroutine Get_Laplacian_SW_SFM 

!......................................................................


!      ------   Obtains laplacian matrix transient heat

implicit none
 
real   (irk) Le, he, factK, factM, g, Pe  
real   (irk) kii, kij, kji, kjj
real   (irk) mii, mij, mji, mjj 
integer(ink) i, j, ielem  
 
!if (allocated (gstif2)) then    ! ****  already done in auxvars. Dimensioned to mélem. GE. nelem_FS
!   deallocate (gstif2) 
!   allocate   (gstif2 (nnode_FS,nnode_FS,nelem_FS))
!else
!   allocate ( gstif2 (nnode_FS,nnode_FS,nelem_FS))
!endif
 
gstif2 = 0.0  
g      = const(1)       !     

if (ndimn.eq.1) then
   do ielem = 1, nelem_FS

      i  = intmat_FS(1,ielem)
      j =  intmat_FS(2,ielem)
      Le = (coord_FS(1,i)-coord_FS(1,j))*(coord_FS(1,i)-coord_FS(1,j))  !   (x(1,i)-x(1,j))*(x(1,i)-x(1,j))
      Le = Le**0.5    
      kij =  - 1/Le
      kji =  - 1/Le
      kii =    1/Le
      kjj =    1/Le
      mii =    Le/3
      mjj =    mii
      mij =    mii/2
      mji =    mii/2
 
      if (type_FS_problem.eq.1) then          !    FS SW  p based 21 Jan 2012 MP chgd
          he    = ( vars_FS (2,i) + vars_FS (2,j) )/2.      ! (rho(i) + rho(j) ) /2.
          Pe    = ( vars_FS (4,i) + vars_FS (4,j) )/2.      ! (p(i) + p(j)) / 2.
          factK =  (dt_sph**2 )/he 
          factM =  1./(2.*Pe)
          gstif2(1,1,ielem) = gstif2(1,1,ielem) + factK * kii  + factM * mii
          gstif2(1,2,ielem) = gstif2(1,2,ielem) + factK * kij  + factM * mij 
          gstif2(2,1,ielem) = gstif2(2,1,ielem) + factK * kji  + factM * mji 
          gstif2(2,2,ielem) = gstif2(2,2,ielem) + factK * kjj  + factM * mjj
      endif

   enddo
endif

END subroutine  Get_Laplacian_SW_SFM

!......................................................................

      subroutine Get_vn1_FS_SFM 

!......................................................................

 
implicit none

integer(ink)  ipoin , iel_fs, i, j
 
real   (irk) g, xx1, Le, N1, N2

g     = const(1) 

if (type_FS_problem.eq.1) then
   if (ndimn.eq.1) then    !   ------------------  NDIMN = 1 -----------------     
      Le  = (coord_FS (1,2)-coord_FS (1,1))
      do ipoin = 1, npoin
         if (If_Out_Domain(ipoin).eq.1)  CYCLE 
         i   = (( (x(1,ipoin) - coord_FS (1,1)) / (coord_FS (1,2)-coord_FS (1,1)) ) +1 ) 
         i   = min(npoin_FS - 1, i)         ! max (i, 1); i = min (i, npoin_FS-1)
         i   = max (1, i)
         j   = i + 1                        ! min (i +1 , npoin_FS) 
         iel_fs = i
         N1     = ( coord_FS (1,j)- x(1,ipoin) ) / (coord_FS (1,j) - coord_FS (1,i)  )
         N2     = ( x(1,ipoin) - coord_FS (1,i) ) / (coord_FS (1,j) - coord_FS (1,i)  )
         grad_p (1, ipoin) = p_fs(I)*geome_Fs (1,iel_fs) + p_fs(J)*geome_Fs (2,iel_fs)
         p     (ipoin)     = N1* p_fs(I)   + N2 * p_fs(J)  
         rho   (ipoin)     = (2.*p(ipoin)/g)**0.5
         vstar (1,ipoin)   =  vx (1,ipoin)
         vx    (1,ipoin)   =  vstar (1,ipoin) - (dt_sph/g) * grad_p (1, ipoin)
         x     (1,ipoin)   =  x0    (1,ipoin) + (dt_sph/2.)* ( vx0(1,ipoin) + vx (1,ipoin) )
      enddo
   endif

   if (if_left_wall.eq.1) then    !   when walls, vn+1 = 0
       i = intmat_FS (1,1)        ! MP 2021 29 Jan
       vx(1,i) = 0.0 
   endif
   if (if_right_wall.eq.1) then
      j = intmat_FS (2, nelem_FS)
      vx(1,j) = 0.0
   endif 

   Call Vn0_BCs_SW ( dt_sph)

endif 

END subroutine  Get_vn1_FS_SFM   




END subroutine Second_Step_FS_SW_SFM   ! ---------------------------------



!-------------------------------------------------------------------

       Subroutine Update_Const_SW (icode, Cindex, Cval, Cval0)
       
!-------------------------------------------------------------------

!      ------   when running a set of cases, updtates the varying const or recovers it
!               MP CHGD 20 march 2020
implicit none

integer(ink) icode, Cindex, icase, jcase, j_Cindex
real   (irk) vj_Val0
real   (irk) Cval, Cval0
real   (irk) dens, denss, densw, poros 
       
IF (icode.eq.1) THEN
   Cval0        = const(Cindex)
   const(Cindex) = Cval
ELSEIF (icode.eq.2) THEN
   icase = OC_mat (nproblems_index,1)
   DO jcase = 1, ncases_mat
      if (icase.eq.jcase) then
         if (Cindex.eq.2.and.const(17).gt.1.e-4 ) then      
             dens  = Cval                                !  input Cval here is the current value of dens                           
             denss = const (17)                          !  from density, if possible get n
             densw = const (18)                          !  and store it to be printed
             dens  = Cval
             poros = (denss-dens)/(denss-densw)
             Cval  = poros                               ! Cval was dens and now is porosity
         endif                                           !  
         const(Cindex) = Cval0                           ! we store back in const the initial value
      elseif (icase.NE.jcase) then
         j_Cindex = Ccases_mat (jcase,1)
         vj_Val0  = Ccases_mat (jcase,2)          
         if (j_Cindex.eq.2.and.const(17).gt.1.e-4 ) then      
             dens  = const  (2)                          ! Ccases_mat (jcase,2) ! *******vj_Val0  !  input Cval here is the current value of dens                           
             denss = const (17)                          !  from density, if possible get n
             densw = const (18)                          !  and store it to be printed
             poros = (denss-dens)/(denss-densw)
             OC_mat (nproblems_index, jcase +3)= poros   ! we store (output to csv) porosity
         endif                                           ! j_Cindex = 2 density var -> porosity 
      endif                                              ! cases icase=jcase or <>
   enddo                                                 ! jcase 
ENDIF

End Subroutine Update_Const_SW


!-------------------------------------------------------------------

       Subroutine Update_Const_SW_BAK (icode, Cindex, Cval, Cval0)
       
!-------------------------------------------------------------------

!      ------   when running a set of cases, updtates the varying const or recovers it
!               MP CHGD 20 march 2020
implicit none

integer(ink) icode, Cindex
real   (irk) Cval, Cval0
real   (irk) dens, denss, densw, poros 
       
if (icode.eq.1) then
   Cval0        = const(Cindex)
   const(Cindex) = Cval
elseif (icode.eq.2) then
   if (Cindex.eq.2.and.const(17).gt.1.e-4 ) then      
       dens  = Cval                                !  input Cval here is the current value of dens                           
       denss = const (17)                          !  from density, if possible get n
       densw = const (18)                          !  and store it to be printed
       dens  = Cval
       poros = (denss-dens)/(denss-densw)
       Cval  = poros                               ! Cval was dens and now is porosity
   endif                                           !  
   const(Cindex) = Cval0                           ! we store back in const the initial value
endif

End Subroutine Update_Const_SW_BAK

!-------------------------------------------------------------------

       Subroutine Stop_case_SW
       
!-------------------------------------------------------------------

!      ------   when running a set of cases, updtates the varying const or recovers it
!               MP CHGD 25 march 2020
implicit none

integer(ink) ipoin, idimn, npts_out, i0 , icw , is_target
real   (irk) v_tmp, v_tmp_s, xmax_tmp, xmax_tmp_s
       
ic_stop_cases = 0        ! comes from driver time vars, when 1 means reached stop condition
npts_out      = 0        ! number of points out of domain
v_tmp         = 0.0      ! average velocity of active nodes
v_tmp_s       = 0.0 
xmax_tmp      = 0.0      ! 1D problems, max X reached by the leading soil node
xmax_tmp_s    = 0.0

do ipoin = 1, npoin 
   if (if_out_domain(ipoin).eq.1) then
      npts_out = npts_out + 1
      CYCLE
   endif
   if (ipoin.le.npoin_s) then
      v_tmp_s = v_tmp + abs (vx(1,ipoin))
      xmax_tmp_s = max ( x(1,ipoin), xmax_tmp_s)
   endif
   v_tmp    = v_tmp + abs (vx(1,ipoin))
   xmax_tmp = max ( x(1,ipoin), xmax_tmp)
enddo

v_tmp = v_tmp/npoin
if (npoin_s.gt.0) then 
   v_tmp_s = v_tmp_s/npoin_s
endif

if( ic_cases_mat.GE.1.AND.itimestep.GT.1) then           !  CHGD MP 6th June 2022
   if (v_tmp_s .LT.vstop_cases.AND.itimestep.GT.10)  ic_stop_cases = 1
   if (npts_out.GT.npoin/5)              ic_stop_cases = 1
   if (time_sph.GT.timestop_cases)  ic_stop_cases = 1
   if (dt      .LT.dtimestop_cases) ic_stop_cases = 1
   if (ndimn.eq.1) OC_mat(nproblems_index, 4 + ncases_mat ) = xmax_tmp_s  
else
   if (time_sph.GE.time_end)        ic_stop_cases = 1
endif

if (ic_stop_cases.EQ.1) then
   i0 = i0
endif

! and now, time to store in OC_mat, positions 3 + ncases_mat +2 +...

if (ic_stop_cases.eq.1.AND.ic_cases_win.GT.0) then 
   do icw = 1, ncase_win
      is_target = info_cw(icw) % is_target                 !  CHGD MP 20 July 2022
      ! i0  = 3 + ncases_mat + 2 + 6*(icw-1)               !  CHGD MP 20 July 2022
      i0  =   3 + ncases_mat + 2 + (icw-1) * nwin_vars     !  all windows can be target
      OC_MAT (nproblems_index, i0 + 1 ) = T0_cw     (icw) 
      OC_MAT (nproblems_index, i0 + 2 ) = Tf_cw     (icw) 
      OC_MAT (nproblems_index, i0 + 3 ) = hpeak_cw  (icw) 
      OC_MAT (nproblems_index, i0 + 4 ) = vnpeak_cw (icw) 
      OC_MAT (nproblems_index, i0 + 5 ) = qpeak_cw  (icw) 
      OC_MAT (nproblems_index, i0 + 6 ) = QTot_cw   (icw)
      if (allocated(dstminpeak_cw).and.is_target.eq.1) then      ! CHGD MP 20 July 2022
         OC_MAT (nproblems_index, i0 + 7 ) = dstminpeak_cw(icw)  ! distance activated when target
      endif
   enddo       
endif

End Subroutine Stop_case_SW


!-------------------------------------------------------------------

       Subroutine Close_files_SW
       
!-------------------------------------------------------------------

!      ------    

implicit none

if (ic_stop_cases.EQ.1) then 
   close (dat_file)
   close (chk_file)
   close (res_file1)
   close (res_file2)
   close (res_file3)
   close (res_file4)   !  deactivate for debugging  ********************  DEBUG
   close (res_file5)
   close (res_file6)
   close (res_file7)  ! Saeidmt 12Nov2023
endif

End Subroutine Close_files_SW

!-------------------------------------------------------------------

       Subroutine Out_Plot_Trigger_SW 
       
!-------------------------------------------------------------------

!      ------  Output for TRIGGER   MP March 2022 April 2022
    
implicit none

integer (ink) islope, i_stats, iproblem, ipoig, izone, index
integer (ink), allocatable:: Slope_in_range(:)           ! MP june 2022
                                                         ! to avoid plotting pots with slopes out of ran
real    (irk) PoF, FoS_mean, FoS_2, FoS_Var, Fos_sigma, FoS_CoeffVar 
real    (irk) Zgrad, dSlope, Zx, Zy
real    (irk) xx, yy, zp                                 ! Saeidmt 30March2022   Adding .CSV for GiS
1007 format (*(G0.7,:,","))                              ! Saeidmt 30March2022   Adding .CSV for GiSge 

!   -----  Task 1 obtain statistics for each problem if ic_MCarlo:NE.6 FOSM

if (.NOT.allocated(slope_in_range)) allocate (slope_in_range(npoig))
slope_in_range = 1

IF (ic_MCarlo.NE.6) THEN  ! avoid if ic_MCarlo:NE.6 FOSM
    Do islope = 1, N_slopes
    Do izone  = 1, n_GZones
       PoF = 0.0; FoS_mean = 0.0; FoS_2 = 0.0; FoS_Var = 0.0; Fos_sigma = 0.0; FoS_CoeffVar = 0.0; 
       Do iproblem = 1, n_MCarlo_problems                  ! each zone nMC problems 
          index = (izone-1)*n_MCarlo_problems + iproblem                   
          if (FoS_Slope (index, islope) .LT.1 ) poF = PoF + 1     !  Saeid remark 20th May 2022 
          FoS_mean = FoS_mean + FoS_Slope (index,islope)
          FoS_2    = FoS_2    + FoS_Slope (index,islope)*FoS_Slope (index,islope)
       enddo
       FoS_mean     = FoS_mean/ n_MCarlo_problems         !   mean
       FoS_2        = FoS_2  /  n_MCarlo_problems          ! E[x**2]
       FoS_Var      = FoS_2 - FoS_mean*FoS_mean   ! sigma
       Fos_sigma    = (FoS_Var)**0.5              ! coeff of variability
       if (FoS_mean.GT.1.e-6) then 
          FoS_CoeffVar = Fos_sigma / FoS_mean   
       else
          FoS_CoeffVar = 999.
       endif 
       index = (izone-1) * N_stats  !!!   *****
       FoS_stats (islope, index + 1) = PoF / n_MCarlo_problems
       FoS_stats (islope, index + 2) = FoS_mean  
       FoS_stats (islope, index + 3) = Fos_sigma
       FoS_stats (islope, index + 4) = FoS_CoeffVar
       FoS_stats (islope, index + 5) = (FoS_mean-1.)/FoS_sigma
    ENDDO
    ENDDO
ENDIF   

!   -----  Task 2 obtain stats at DEM points

Do ipoig = 1, npoig
   izone = topol (16, ipoig)
   Zx = topol (2,ipoig)
   Zy = topol (3,ipoig)
   Zgrad = (Zx * Zx + Zy*Zy)
   if (Zgrad.GT.0.001) then
      Zgrad = Zgrad ** 0.5
   else
      Zgrad = 0.0
   endif
   if (Zgrad.LT.slope(1).OR.Zgrad.GT.slope(N_slopes)) then
      slope_in_range(ipoig) = 0   
   endif
   dslope =  (slope(N_slopes) - slope(1))/ (N_slopes-1)
   islope =  (Zgrad - slope(1))/dslope + 1
   if ( islope.LT.1 ) islope = 1
   if ( islope.GT.N_slopes ) islope = N_slopes
   index = (izone-1) * N_stats  !!!   *****
   FoS_DEM (ipoig, 1) = FoS_stats (islope,index+1)
   FoS_DEM (ipoig, 2) = FoS_stats (islope,index+2)
   FoS_DEM (ipoig, 3) = FoS_stats (islope,index+3)
   FoS_DEM (ipoig, 4) = FoS_stats (islope,index+4)
   FOS_DEM (ipoig, 5) = FoS_stats (islope,index+5)
enddo

!      ------  Task 3  print in gid.res

write(gid_res,*) 'Result "PoF" "Prob of Failure" ',time,'   Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (slope_in_range(ipoig).eq.0) then
      cycle
   endif
   write(gid_res,*) ipoig,  FoS_DEM (ipoig,1)    
enddo
write(gid_res,*) ' End Values ' 

write(gid_res,*) 'Result "FoS mean" "Average Fos" ',time,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (slope_in_range(ipoig).eq.0) then
      cycle
   endif
   write(gid_res,*) ipoig, FoS_DEM (ipoig,2)    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "FoS sigma" "sigma" '     ,time,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (slope_in_range(ipoig).eq.0) then
      cycle
   endif
   write(gid_res,*) ipoig, FoS_DEM (ipoig,3)    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "FoS variab" "coeff variab" ',time,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (slope_in_range(ipoig).eq.0) then
      cycle
   endif
   write(gid_res,*) ipoig, FoS_DEM (ipoig,4)    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "Reliab index" "reliab" ',time,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (slope_in_range(ipoig).eq.0) then
      cycle
   endif
   write(gid_res,*) ipoig, FoS_DEM (ipoig,5)    
enddo
write(gid_res,*) ' End Values '

do ipoig = 1, npoig                                 ! Saeidmt 30March2022   Adding .CSV for GiS
   xx    = coorg(1,ipoig)
   yy    = coorg(2,ipoig)
   zp    = topol(1,ipoig)
   write(gis_csv1,1007) ipoig, xx, yy, zp, FoS_DEM (ipoig,1), FoS_DEM (ipoig,2), FoS_DEM (ipoig,3), FoS_DEM (ipoig,4), FoS_DEM (ipoig,5)
enddo
do islope = 1, N_slopes                                 ! Saeidmt 30March2022   Adding .CSV for GiS
   write(gis_csv2,1007) slope(islope), FoS_stats (islope, index + 1), FoS_stats (islope, index + 2), FoS_stats (islope, index + 3), FoS_stats (islope, index + 4), FoS_stats (islope, index + 5), FoS_stats (islope, index + 6)
enddo

End Subroutine Out_Plot_Trigger_SW  


!-------------------------------------------------------------------

       Subroutine Get_Prob_Propag_SW_BAK 
       
!-------------------------------------------------------------------

!      ------  omputes FOSM for propagation   MP June 2022  
    
implicit none

integer (ink)  i0, icase, itarget, idest, ntarget
real    (irk)  Fi_plus, Fi_minus, Sigma_i, DFidXi, Sigma_Pr2 , dXi


i0 = 3 + ncases_Mat+2 + 6*(ncase_win-1)
ntarget = 7

DO itarget = 1, ntarget
   Fi_plus = 0.0;  Fi_minus = 0.0; Sigma_i = 0.0; DFidXi = 0.0; Sigma_Pr2 = 0.0 
   DO icase = 1, ncases_mat
      idest    = 2*icase
      Fi_plus  = OC_Mat ( idest,   i0 + itarget)
      Fi_minus = OC_Mat ( idest+1, i0 + itarget)
      Sigma_i  = FOSMms (3,icase)
      if (Sigma_i.GT.1.e-4) then
         dXi    = Sigma_i /10.
         DFidXi = (Fi_plus-Fi_minus)/dXi
      else
         DFidXi = 0.0
      endif
      Sigma_Pr2 = Sigma_Pr2 + ( DFidXi*Sigma_i) **2.  
   ENDDO
   Sigma_Pr2 = Sigma_Pr2**0.5
   FOSM_Pr2 (1,itarget) = OC_Mat (1, i0 + itarget) 
   FOSM_Pr2 (2,itarget) = Sigma_Pr2
ENDDO

End Subroutine Get_Prob_Propag_SW_BAK  

!-------------------------------------------------------------------

       Subroutine Get_Prob_Propag_SW 
       
!-------------------------------------------------------------------

!      ------  omputes FOSM for propagation   MP June 2022  
    
implicit none

integer (ink)  i0, icase, itarget, idest, ntarget, iproblem, iwin
real    (irk)  Fi_plus, Fi_minus, Sigma_i, DFidXi, Sigma_Pr2 , dXi
real    (irk)  mean,  secondM  

! i0 = 3 + ncases_Mat+2 + 6*(ncase_win-1)  ! CHGD MP 20 July 2022
ntarget = nwin_vars

If (ic_cases_mat.EQ.20) THEN   !    FOSM for propagation
  
   if (.NOT.allocated(FOSM_Pr2))   then                !  MP 20 July 2022
      allocate ( FOSM_Pr2 (2, ncase_win*nwin_vars) )   !  to use more than 1 eindow
      FOSM_Pr2 = 0.0
   endif
   
   DO iwin = 1, ncase_win      ! CH MP 20 July 2022
      i0 = 3 + ncases_Mat+2 + (iwin-1)*nwin_vars     ! CHGD MP 20 July 2022
      DO itarget = 1, ntarget              ! itarget is the variable considered in the window iwin
         Fi_plus = 0.0;  Fi_minus = 0.0; Sigma_i = 0.0; DFidXi = 0.0; Sigma_Pr2 = 0.0 
         DO icase = 1, ncases_mat
            idest    = 2*icase
            Fi_plus  = OC_Mat ( idest,   i0 + itarget)
            Fi_minus = OC_Mat ( idest+1, i0 + itarget)
            Sigma_i  = FOSMms (3,icase)
            if (Sigma_i.GT.1.e-4) then
               dXi    = Sigma_i /10.
               DFidXi = (Fi_plus-Fi_minus)/(2.*dXi)
            else
               DFidXi = 0.0
            endif
            Sigma_Pr2 = Sigma_Pr2 + ( DFidXi*Sigma_i) **2.  
         ENDDO
         Sigma_Pr2 = Sigma_Pr2**0.5
         FOSM_Pr2 (1,itarget) = OC_Mat (1, i0 + itarget) 
         FOSM_Pr2 (2,itarget) = Sigma_Pr2
      ENDDO
   ENDDO

ELSE   !   this is for MCarlo
 
   if (.NOT.allocated(MCarlo_Pr2)) then                !  MP 20 July 2022                
      allocate ( MCarlo_Pr2 (2, ncase_win*nwin_vars) ) !  to use more than 1 eindow
      MCarlo_Pr2 = 0.0
   endif
   
   DO iwin = 1, ncase_win      ! CH MP 20 July 2022
      i0 = 3 + ncases_Mat+2 + (iwin-1)*nwin_vars     ! CHGD MP 20 July 2022
      DO itarget = 1, ntarget
         mean    = 0.0
         secondM = 0.0
         do iproblem = 1, nproblems
            mean    = mean + OC_mat (iproblem,i0 + itarget)
            secondM = secondM + (OC_mat (iproblem,i0 + itarget))**2
         enddo
         mean    = mean / nproblems
         secondM = secondM / nproblems
         Sigma_Pr2 = secondM - mean*mean
         MCarlo_Pr2 (1,itarget) = mean
         MCarlo_Pr2 (2,itarget) = Sigma_Pr2**0.5
      ENDDO
   ENDDO
ENDIF 

End Subroutine Get_Prob_Propag_SW  

!-------------------------------------------------------------------

       Subroutine Get_pixel_hVars_SW 
       
!-------------------------------------------------------------------

!      ------  omputes FOSM for propagation   MP June 2022  
    
implicit none 

integer (ink) ipoin, ipoigx, ipoigy, ipoig, ieleg, n1, n2, n3, n4
integer (ink) Lnode (4), inode, iproblem, n_pr_pts
real    (irk) OCM (npixel_HVars)
real    (irk) xx, yy, hh, vv, qq

iproblem = nproblems_index      ! simpler to use

if (.NOT.allocated (FOSM_pixel) ) then
   allocate ( FOSM_pixel(npoig) )
   FOSM_pixel(:)%if_alloc                = 0  ! state initial not allocated. Made 1 when allocating 
   FOSM_pixel(:)%if_all_FOSM             = 0  ! vor FOSM_Pr alloc when not allocated
endif                           

DO ipoin = 1, npoin
   
   Call Check_Out_Domain_SW 
   if ( if_Out_Domain(ipoin).eq.1) CYCLE

   xx = X (1,ipoin)  
   if (ndimn.eq.2) then
       yy = X (2,ipoin)
   else
      yy = 0.0
   endif
   
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
   
   ieleg     = (ipoigy-1)*(npoigx-1) + ipoigx
   ipoig     = (ipoigy-1)*npoigx + ipoigx
   lnode (1) = ipoig
   lnode (2) = ipoig + 1
   lnode (3) = ipoig + 1 + npoigx
   lnode (4) = lnode (3) -1 

   DO inode = 1, 4
                                                   ! trick to allocate %OCMM wif(allocated)/o asking 
      ipoig = lnode(inode)                         ! which does not work... introduce the field if_alloc
      if (FOSM_pixel(ipoig)  % if_alloc.EQ.0) then ! Initialized to zero when allocating FOSM_pixel
         allocate (FOSM_pixel(ipoig)  % OCMM(nproblems, npixel_hVars) )
         allocate (FOSM_pixel(ipoig)  % nr_pts_pix(nproblems) )
         FOSM_pixel(ipoig)  % OCMM       = 0
         FOSM_pixel(ipoig)  % nr_pts_pix = 0
         FOSM_pixel(ipoig)  % if_alloc   = 1         ! now we change the estate to 1 (allocated)
      endif

      OCM = FOSM_pixel(ipoig)  % OCMM (iproblem,:)
      FOSM_pixel(ipoig)% nr_pts_pix(iproblem) = FOSM_pixel(ipoig) % nr_pts_pix(iproblem)  +  1 
      hh    = rho (ipoin)
      if (allocated (h_DF) ) then
         hh = h_Df (ipoin)
      endif
      vv = vx (1,ipoin)* vx (1,ipoin)
      if (ndimn.eq.2) then
         vv = vv + vx (2,ipoin)* vx (2,ipoin)
      endif
      vv = vv**0.5
      qq = vv*hh

      OCM (1) = OCM (1) + hh
      OCM (3) = OCM (3) + vv
      OCM (5) = OCM (5) + qq
      if (hh.GT.OCM (2)) OCM (2) = hh
      if (vv.GT.OCM (4)) OCM (4) = vv
      if (qq.GT.OCM (6)) OCM (6) = qq

      FOSM_pixel(ipoig) % OCMM (iproblem,:) = OCM(:)

   enddo
ENDDO

END Subroutine Get_pixel_hVars_SW

!-------------------------------------------------------------------

       Subroutine Get_pixel_FOSM_SW 
       
!-------------------------------------------------------------------

!      ------  omputes FOSM for propagation   MP June 2022  
    
implicit none 

integer (ink) ipoig, nptsx , ipixVar, icase, idest, n_pr_pts, iproblem
integer (ink) if_Void
real    (irk) Fi_plus,Fi_minus,Sigma_i , DFidXi, Sigma_Pr2 , DXi 
real    (irk) OCM (nproblems, npixel_HVars), FOSM_prr (2,npixel_HVars)

DO ipoig = 1, npoig

   if_Void = 0
   if (FOSM_pixel(ipoig) % if_alloc.EQ.1) then   ! wecompute it for cells  not void

      OCM   = FOSM_pixel(ipoig) % OCMM

      Do iproblem = 1, nproblems  ! normalize vars 2,4 and 6 (h,v,q avg)
         nptsx = FOSM_pixel(ipoig) % nr_pts_pix(iproblem)
         if ( nptsx.eq.0)  then 
             nptsx   = 1    ! just in case.The if_allocated should be enough
             if_Void = 1
         endif
         do ipixVar = 1, npixel_HVars, 2
            OCM (iproblem,ipixVar) = OCM (iproblem,ipixVar) /nptsx
         enddo
      enddo
      
      if (if_Void.Eq.1) then
         CYCLE
      endif

      if (FOSM_pixel(ipoig) % if_all_FOSM.eq.0) then
         allocate   (FOSM_pixel(ipoig) % FOSM_pr (2,npixel_HVars))
         FOSM_pixel(ipoig) % FOSM_pr = 0.0
         FOSM_pixel(ipoig) % if_all_FOSM = 1
      endif

      deallocate (FOSM_pixel(ipoig) % OCMM)  
       
      FOSM_Prr = 0.0
         
      DO ipixVar = 1, npixel_HVars          
         Fi_plus = 0.0;  Fi_minus = 0.0; Sigma_i = 0.0; DFidXi = 0.0; Sigma_Pr2 = 0.0 
         DO icase = 1, ncases_mat
            idest    = 2*icase
            Fi_plus  = OCM ( idest,   ipixVar)
            Fi_minus = OCM ( idest+1, ipixVar)
            Sigma_i  = FOSMms (3,icase)
            if (Sigma_i.GT.1.e-4) then
               dXi    = Sigma_i /10.
               DFidXi = (Fi_plus-Fi_minus)/(2.*dXi)
            else
               DFidXi = 0.0
            endif
            Sigma_Pr2 = Sigma_Pr2 + ( DFidXi*Sigma_i) **2.  
         ENDDO
         Sigma_Pr2 = Sigma_Pr2**0.5
         FOSM_Prr (1,ipixVar) = OCM (1, ipixVar) !! perhaps OC_MAT ??? ***
         FOSM_Prr (2,ipixVar) = Sigma_Pr2
      ENDDO
      FOSM_pixel(ipoig) % FOSM_pr (:,:) = FOSM_Prr (:,:) 
      FOSM_pixel(ipoig) % if_alloc = 0
   ENDIF
    
ENDDO

END Subroutine Get_pixel_FOSM_SW



!-------------------------------------------------------------------

       Subroutine Plot_pixel_RES_SW
       
!-------------------------------------------------------------------

!      ------  Output FOSM for propagation   MP June 2022  
    
implicit none 

integer (ink) ipoig, nptsx , ipixVar, icase, idest
real    (irk) time_dummy     


time_dummy = 0.0

write(gid_res,*) 'Result "FoS av hh" "hh" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (FOSM_pixel(ipoig) % if_all_FOSM.EQ.1) then
      write(gid_res,*) ipoig, FOSM_pixel(ipoig) % FOSM_pr (1, 1)  
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "FoS av vv" "vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (FOSM_pixel(ipoig) % if_all_FOSM.EQ.1) then
      write(gid_res,*) ipoig, FOSM_pixel(ipoig) % FOSM_pr (1, 3)  
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "FoS av qq" "qq" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (FOSM_pixel(ipoig) % if_all_FOSM.EQ.1) then
      write(gid_res,*) ipoig, FOSM_pixel(ipoig) % FOSM_pr (1, 5)  
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sigma FoS av hh" "sigma hh" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (FOSM_pixel(ipoig) % if_all_FOSM.EQ.1) then
      write(gid_res,*) ipoig, FOSM_pixel(ipoig) % FOSM_pr (2, 1)  
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sigma FoS av vv" "sigma vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (FOSM_pixel(ipoig) % if_all_FOSM.EQ.1) then
      write(gid_res,*) ipoig, FOSM_pixel(ipoig) % FOSM_pr (2, 3)  
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sigma FoS av qq" "sigma qq" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (FOSM_pixel(ipoig) % if_all_FOSM.EQ.1) then
      write(gid_res,*) ipoig, FOSM_pixel(ipoig) % FOSM_pr (2, 5)  
   endif    
enddo
write(gid_res,*) ' End Values '

End subroutine Plot_pixel_RES_SW 


!-------------------------------------------------------------------

       Subroutine Get_MC_pixel_Vars_SW 
       
!-------------------------------------------------------------------

!      ------  Computes FOSM for propagation   MP June 2022  
    
implicit none 

integer (ink) ipoin, ipoigx, ipoigy, ipoig, ieleg, n1, n2, n3, n4
integer (ink) Lnode (4), inode, iproblem, n_pr_pts
real    (irk) OCM (npixel_HVars)
real    (irk) xx, yy, hh, vv, qq, dens, hvuln, Pvuln  !  MP 22 sept 2022 
real    (irk) Intens, Suscept,  Vuln                  !  MP 22 sept 2022 

if (.NOT.allocated (MC_pixel) ) then
   allocate ( MC_pixel(npoig) )
   MC_pixel(:)%if_OCMC_alloc     = 0  ! state initial not allocated. Made 1 when allocating 
   MC_pixel(:)%nr_pts_pix        = 0 
   MC_pixel(:)%if_MC_Pr_alloc    = 0  ! state initial not allocated. Made 1 when allocating 
   MC_pixel(:)%nr_problem_pix    = 0   
endif                           

DO ipoin = 1, npoin
   
   Call Check_Out_Domain_SW 
   if ( if_Out_Domain(ipoin).eq.1) CYCLE

   xx = X (1,ipoin)  
   if (ndimn.eq.2) then
       yy = X (2,ipoin)
   else
      yy = 0.0
   endif
   
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
   
   ieleg     = (ipoigy-1)*(npoigx-1) + ipoigx
   ipoig     = (ipoigy-1)*npoigx + ipoigx
   lnode (1) = ipoig
   lnode (2) = ipoig + 1
   lnode (3) = ipoig + 1 + npoigx
   lnode (4) = lnode (3) -1 

   DO inode = 1, 4
                                                   ! trick to allocate %OCMM wif(allocated)/o asking 
      ipoig = lnode(inode)                         ! which does not work... introduce the field if_alloc
      if (MC_pixel(ipoig)%if_OCMC_alloc.EQ.0) then ! Initialized to zero when allocating FOSM_pixel
         allocate (MC_pixel(ipoig)  % OCMC    (npixel_hVars) )
         MC_pixel(ipoig)  % OCMC              = 0.
         MC_pixel(ipoig)  % nr_pts_pix        = 0
         MC_pixel(ipoig)  % if_OCMC_alloc     = 1   ! now we change the estate to 1 (allocated)
      endif
      
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.0) then                !  aug 8 MP
         allocate (MC_pixel(ipoig)  % MC_Pr (4,npixel_hVars) )       ! aug 22 MP           
         MC_pixel(ipoig)  % MC_Pr             = 0.
         MC_pixel(ipoig)  % nr_problem_pix    = 0
         MC_pixel(ipoig)  % if_MC_Pr_alloc    = 1   ! now we change the estate to 1 (allocated)
      endif


      OCM =  MC_pixel(ipoig)  % OCMC (:)
      MC_pixel(ipoig)% nr_pts_pix =  MC_pixel(ipoig) % nr_pts_pix +  1 
      hh    = rho (ipoin)
      if (allocated (h_DF) ) then
         hh = h_Df (ipoin)
      endif
      vv = vx (1,ipoin)* vx (1,ipoin)
      if (ndimn.eq.2) then
         vv = vv + vx (2,ipoin)* vx (2,ipoin)
      endif
      vv = vv**0.5
      qq = vv*hh
!     ------  If topol(16,ipoig).NE.1  it means 
      OCM (1) = OCM (1) + hh
      OCM (3) = OCM (3) + vv
!      OCM (5) = OCM (5) + qq
      if (hh.GT.OCM (2)) OCM (2) = hh
      if (vv.GT.OCM (4)) OCM (4) = vv
!      if (qq.GT.OCM (6)) OCM (6) = qq
      if (topol(16, ipoig).GT.1.e-3)  then
         dens  = const(2)
         hvuln = hh 
         Pvuln = 0.5*dens*vv*vv/1000.
         call Get_Assets_Vulnerability (ipoig, hvuln, Pvuln, Intens, Suscept,  Vuln )   !!  Vuln modulus
         OCM (5) = OCM (5) + Intens
         OCM (7) = OCM (7) + Suscept
         OCM (9) = OCM (9) + Vuln
         if (Intens.GT.OCM  (6)) OCM  (6) = Intens
         if (Intens.GT.OCM  (8)) OCM  (8) = Suscept
         if (Intens.GT.OCM (10)) OCM (10) = Vuln 
      endif

      MC_pixel(ipoig) % OCMC (:) = OCM(:)

   enddo
ENDDO


End Subroutine Get_MC_pixel_Vars_SW 


!-------------------------------------------------------------------

       Subroutine Norm_MC_pixel_Vars_SW (offset)
       
!-------------------------------------------------------------------

!      ------  Normalizes the values at every pixel after timestep loop 
    
implicit none 

integer (ink), intent(IN), OPTIONAL:: offset       ! MPAug 2022 to plot every problem
logical if_offset
real    (irk) mean, sigma2, sigmatmp
integer (ink) ipoig, nptsx , ipixVar, nptxel


if_offset  = PRESENT (offset)

DO ipoig = 1, npoig

!  ...............   Task 1: normalize OCMC

   IF (MC_pixel(ipoig)%if_OCMC_alloc.EQ.1) THEN
      
      nptsx = MC_pixel(ipoig)  % nr_pts_pix 
      if ( nptsx.eq.0)  nptsx = 1    ! just in case.The if_allocated should be enough
      do ipixVar = 1, npixel_HVars, 2
          MC_pixel(ipoig) % OCMC(ipixVar) = MC_pixel(ipoig) % OCMC(ipixVar) /nptsx
      enddo

!  ...............   Task 2: compute MC_pr  Row 1 accum h, row 2 h**2
   
      do ipixVar = 1, npixel_HVars 
         MC_pixel(ipoig) % MC_Pr(1,ipixVar) = MC_pixel(ipoig) % MC_Pr(1,ipixVar) &
                                            + MC_pixel(ipoig) % OCMC(ipixVar)
         MC_pixel(ipoig) % MC_Pr(2,ipixVar) = MC_pixel(ipoig) % MC_Pr(2,ipixVar) &
                                            + (MC_pixel(ipoig) % OCMC(ipixVar))**2
      enddo

!  ...............   Task 3: update nr_problem_pix = +1

      MC_pixel(ipoig)%nr_problem_pix = MC_pixel(ipoig)%nr_problem_pix + 1 

!  ...............   Task 4 set back to zero OCMC and % nr_pts_pix

      MC_pixel(ipoig)  % nr_pts_pix = 0
      MC_pixel(ipoig) % OCMC(:)     = 0.0
      MC_pixel(ipoig)%if_OCMC_alloc = 0
      deallocate (MC_pixel(ipoig) % OCMC)

!  ...............   Task 5 statistics if offset eq 1

      if (if_offset) then  
         nptxel =  MC_pixel(ipoig)  % nr_problem_pix  
         if (nptxel.EQ.0) nptxel = 1
         DO ipixVar = 1, npixel_HVars 
            mean              = MC_pixel(ipoig)  % MC_Pr (1,ipixVar)/nptxel
            sigma2            = MC_pixel(ipoig)  % MC_Pr (2,ipixVar)/nptxel
            MC_pixel(ipoig)  % MC_Pr (3,ipixVar) = mean
            sigmatmp = sigma2 - mean*mean
         !   if (abs(sigmatmp).lt.1.e-6) then    ! MP11 11 2022 to avoin¡d nans
            if (sigmatmp.lt.1.e-6) then
               sigmatmp = 0.0
            endif
            MC_pixel(ipoig)  % MC_Pr (4,ipixVar) = (sigmatmp)**0.5 
         enddo
      endif
   
   ENDIF

ENDDO

End Subroutine Norm_MC_pixel_Vars_SW 


!-------------------------------------------------------------------

       Subroutine Stats_MC_pixel_Vars_SW 
       
!-------------------------------------------------------------------

!      ------  Computes FOSM for propagation   MP June 2022  
    
implicit none 

integer (ink) ipoig, ipixVar, nptxel
real    (irk) mean, sigma2, sigmatmp

DO ipoig = 1, npoig

   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then   ! we compute it for cells  not void
       nptxel =  MC_pixel(ipoig)  % nr_problem_pix  
       if (nptxel.EQ.0) nptxel = 1
       DO ipixVar = 1, npixel_HVars 
          mean              = MC_pixel(ipoig)  % MC_Pr (1,ipixVar)/nptxel
          sigma2            = MC_pixel(ipoig)  % MC_Pr (2,ipixVar)/nptxel
          MC_pixel(ipoig)  % MC_Pr (1,ipixVar) = mean
          sigmatmp = sigma2 - mean*mean
          if (abs(sigmatmp).lt.1.e-6) then
             sigmatmp = 0.0
          endif
          MC_pixel(ipoig)  % MC_Pr (2,ipixVar) = (sigmatmp )**0.5   
       enddo
   endif

ENDDO

End Subroutine Stats_MC_pixel_Vars_SW 


!-------------------------------------------------------------------

       Subroutine Plot_MC_pixel_Vars_SW_old (offset)
       
!-------------------------------------------------------------------

!      ------  Computes FOSM for propagation   MP June 2022  
    

!      ------  Output FOSM for propagation   MP June 2022  
    
implicit none 

integer (ink), intent(IN), OPTIONAL:: offset       ! MPAug 2022 to plot every problem
logical if_offset

integer (ink) ipoig, nptsx , ipixVar, icase, idest 
real    (irk) time_dummy     
integer (ink) time_dum

character (72) text

1000 format (A72)

time_dummy = nproblems_index    !       0.0
time_dum   = nproblems_index 
idest      = 0
if_offset  = PRESENT (offset)

rewind (gid_res)
read (gid_res,1000) text
if (if_offset) then
   idest = 2  
!   rewind (gid_res)
!   read (gid_res,1000) text 
endif


!  -----------------------------  avg h, v , intens, suscept, Vulnerability 

write(gid_res,*) 'Result "MC av hh" "01hh" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,1)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "MC av vv" "02vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,3)   
   endif    
enddo
write(gid_res,*) ' End Values '

if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "Intens" "03intens" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,5)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Avsuscept" "04Suscept" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,7)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "av Vuln" "05Vuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,9)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif

!  ----------------------------- peak h, v and q 

write(gid_res,*) 'Result "MC peak hh" "11hhpk" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,2)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "MC peak vv" "12vvpeak" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,4)   
   endif    
enddo
write(gid_res,*) ' End Values '


if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "pkk Intens" "13PKintns" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,6)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "pkk suscept" "14PkSuscept" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,8)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "PKk Vuln" "15pkVuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,10)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif


!  -----------------------------  sigma of avg h, v , intens, suscept, Vulnerability 

write(gid_res,*) 'Result "Sav hh" "06Sgmhh" ',time_dum ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,1)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sigma av vv" "07sigma vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,3)   
   endif    
enddo
write(gid_res,*) ' End Values '

if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "sgm av Intens" "08sgmintens" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,5)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Sgm avSuscp" "09sgmScpt" ',time_dum ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,7)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Sgm Vuln" "10sgmVuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,9)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif

!  -----------------------------  sigma of peak h, v , intens, suscept, Vulnerability 

write(gid_res,*) 'Result "sgm peak hh" "16Sg peak hh" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,2)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sgm peak vv" "17Sg peak vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,4)   
   endif    
enddo
write(gid_res,*) ' End Values '

if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "SpkInt" "18s pkint" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,6)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "S pk sc" "19Sg PkSc" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,8)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Sgpk Vuln" "20Sgpk Vuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,10)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif

End Subroutine Plot_MC_pixel_Vars_SW_old 


!-------------------------------------------------------------------

       Subroutine Plot_MC_pixel_Vars_SW  (offset)
       
!-------------------------------------------------------------------

!      ------  Computes MC for propagation   MP June 2022  
    

!      ------  Output Mc for propagation   MP June 2022  
    
implicit none 

integer (ink), intent(IN), OPTIONAL:: offset       ! MPAug 2022 to plot every problem
logical if_offset

integer (ink) ipoig, nptsx , ipixVar, icase, idest, ivar , idimn      !  MP2023 11 09
real    (irk) time_dummy 
real    (irk)  temp1    
integer (ink) time_dum
integer (ink) ip_Asset,  ip_activ  
integer (ink) NP_A        ! debugg tool                                        !  MP 02 10 2022

real  (irk) cost_per_m2, pix_area, asset_cost  ! added new feature MP 2023 11 09

real    (irk), allocatable,SAVE:: Avg_Asset_Vals (:,:)  ! avd vals h,v, I, S & Vuln @ assets (10,n_Assets)
integer( ink), allocatable,SAVE:: NP_Asset_Vals  (:)    ! Nr of pts in assets(n_Assets)

character (72) text

1000 format (A72)

time_dummy = nproblems_index    !       0.0
time_dum   = nproblems_index 
idest      = 0
if_offset  = PRESENT (offset)

pix_area   = deltxg * deltyg    !  2023 11 09 for cost estimation

if (ic_Vuln.EQ.1) then                       !  Allocate and set to zero  MP 06 10 2022
  if (.NOT.allocated(Avg_Asset_Vals)) then
      allocate ( Avg_Asset_Vals (10,n_Assets)) ; Avg_Asset_Vals = 0.0
      allocate ( NP_Asset_Vals     (n_Assets)) ; NP_Asset_Vals  = 0 
   endif
endif

rewind (gid_res)
read (gid_res,1000) text
if (if_offset) then
   idest = 2  
!   rewind (gid_res)
!   read (gid_res,1000) text 
endif

rewind (QGIS_asset)   !  MP 2023 11 09
rewind (QGIS_avg)
rewind (QGIS_peaks)

!  ==============  we will plot first averages at all Assets 1..N:Assets if ic_vuln=1

IF (ic_Vuln.eq.1) THEN     !   add values to Avg_assets, count number of vaues 
   
   Avg_Asset_Vals = 0.0; NP_Asset_Vals = 0   !   Bug detected MP August 2023

   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
         NP_Asset_Vals  (ip_Asset) =   NP_Asset_Vals  (ip_Asset) +1                
         do ivar = 1, 10,2
            Avg_Asset_Vals (ivar,ip_Asset) =  Avg_Asset_Vals (ivar,ip_Asset) &
                                            +  MC_pixel(ipoig)  % MC_Pr(1+idest,ivar)               
         enddo               
         do ivar = 1, 10,2
            Avg_Asset_Vals (ivar+1,ip_Asset) =  Avg_Asset_Vals (ivar+1,ip_Asset) &
                                            +  MC_pixel(ipoig)  % MC_Pr(2+idest,ivar)               
         enddo   
      endif
   enddo
    
   Do ip_Asset = 1, N_Assets    !   Normalize to obtain averages 
      if (Asset_Info (ip_Asset)% np_in_Asset.GT.0) then
         Avg_Asset_Vals (:,ip_Asset) = Avg_Asset_Vals (:,ip_Asset) & 
                                      /Asset_Info (ip_Asset)% np_in_Asset
      else
          Avg_Asset_Vals (:,ip_Asset) = 0.0  
      endif
  !     if (NP_Asset_Vals  (ip_Asset).GT.0) then   CHANGED MP 11 11 2022
  !       Avg_Asset_Vals (:,ip_Asset) = Avg_Asset_Vals (:,ip_Asset)/NP_Asset_Vals  (ip_Asset)
  !    else
  !        Avg_Asset_Vals (:,ip_Asset) = -999
  !     endif
   enddo

   !   -----------------   write blocks

   ! (1) average of heights @ ASSETS
   
   write(gid_res,*) 'Result "MChh AsT" "001hA" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
         write(gid_res,*) ipoig, Avg_Asset_Vals (1,ip_Asset)   
      endif    
   enddo
   write(gid_res,*) ' End Values ' 

   ! (2) average of velocities @ ASSETS 
   
   write(gid_res,*) 'Result "MC v AsT" "002vA" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
          write(gid_res,*) ipoig, Avg_Asset_Vals (3,ip_Asset)   
      endif    
   enddo
   write(gid_res,*) ' End Values ' 

   ! (3) average of Intensities @ ASSETS 
   
   write(gid_res,*) 'Result "MC I AsT" "003IA" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
          write(gid_res,*) ipoig, Avg_Asset_Vals (5,ip_Asset)   
      endif    
   enddo
   write(gid_res,*) ' End Values '  

   ! (4) average of Susceptibilities @ ASSETS 
   
   write(gid_res,*) 'Result "MC S AsT" "004SA" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
          write(gid_res,*) ipoig, Avg_Asset_Vals (7,ip_Asset)   
      endif    
   enddo
   write(gid_res,*) ' End Values '  

   ! (5) average of Vulnerbilities @ ASSETS 
   
   write(gid_res,*) 'Result "MC VulAsT" "005VnSA" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
          write(gid_res,*) ipoig, Avg_Asset_Vals (9,ip_Asset)   
      endif    
   enddo
   write(gid_res,*) ' End Values ' 

   ! (6) Sigma of Vulnerbilities  @ ASSETS
   
   write(gid_res,*) 'Result "SgVulAsT" "010VnSA" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
          write(gid_res,*) ipoig, Avg_Asset_Vals (10,ip_Asset)   
      endif    
   enddo
   write(gid_res,*) ' End Values ' 

!   Do ip_Asset = 1, N_Assets    !   UN-Normalize ********  NOT USED LATER     MP 2023 11 09 deactivated
!      Avg_Asset_Vals (:,ip_Asset) = Avg_Asset_Vals (:,ip_Asset)*NP_Asset_Vals  (ip_Asset)
!   enddo

ENDIF


!  -----------------------------  avg h, v , intens, suscept, Vulnerability 

write(gid_res,*) 'Result "MC av hh" "01hh" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,1)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "MC av vv" "02vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,3)   
   endif    
enddo
write(gid_res,*) ' End Values '

if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "Intens" "03intens" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,5)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Avsuscept" "04Suscept" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,7)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "av Vuln" "05Vuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,9)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif

!  ----------------------------- peak h, v and q 

write(gid_res,*) 'Result "MC peak hh" "11hhpk" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,2)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "MC peak vv" "12vvpeak" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,4)   
   endif    
enddo
write(gid_res,*) ' End Values '


if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "pkk Intens" "13PKintns" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,6)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "pkk suscept" "14PkSuscept" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,8)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "PKk Vuln" "15pkVuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(1+idest,10)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   ! (1) print assets to QGIS_asset file  MP 2023 11 09 
   
   9991 format (2x, i6, 8(2x,F10.4),2x, G10.4)
   write(QGIS_asset,*) ' ipoin   X     Y    havg   vavg  Intens Suscept Vuln  Sigma Vuln cost@asset USDs '
   Do ipoig = 1, npoig
      ip_Asset = topol (16, ipoig)
      ip_activ = MC_pixel(ipoig) % if_MC_Pr_alloc
      if (ip_asset.NE.0.AND. ip_activ.NE.0) then
         Cost_per_m2 = Asset_Info (ip_Asset)%Cost_per_m2 (1)
         asset_cost = Cost_per_m2 * pix_area *  Avg_Asset_Vals (9,ip_Asset)* Asset_Info (ip_Asset)% np_in_Asset 
         write(QGIS_asset,9991) ipoig, (coorg(1,ipoig) ), (coorg(2,ipoig) ), Avg_Asset_Vals (1,ip_Asset) , &
                                 Avg_Asset_Vals (3,ip_Asset), Avg_Asset_Vals (5,ip_Asset), Avg_Asset_Vals (7,ip_Asset), &
                                 Avg_Asset_Vals (9,ip_Asset) , Avg_Asset_Vals (1,ip_Asset), asset_cost  
      endif    
   enddo 

   ! (2) print averages of h,v,I,susc, Vulnand sigma Vuln  to QGIS_avg file  MP 2023 11 09 
   
   write(QGIS_avg,*) ' ipoin   X     Y    havg   vavg  Intens Suscept Vuln  Sigma Vuln '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(QGIS_avg,9991) ipoig, (coorg(1,ipoig) ), (coorg(2,ipoig) ), MC_pixel(ipoig)%MC_Pr(1+idest,1)  , &
                                 MC_pixel(ipoig)%MC_Pr(1+idest,3)  , MC_pixel(ipoig)%MC_Pr(1+idest,5)  ,  &
                                 MC_pixel(ipoig)%MC_Pr(1+idest,7)  , MC_pixel(ipoig)%MC_Pr(1+idest,9)  ,  &
                                 MC_pixel(ipoig)%MC_Pr(2+idest,9)    
      endif    
   enddo 

endif


!  -----------------------------  sigma of avg h, v , intens, suscept, Vulnerability 

write(gid_res,*) 'Result "Sav hh" "06Sgmhh" ',time_dum ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,1)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sigma av vv" "07sigma vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,3)   
   endif    
enddo
write(gid_res,*) ' End Values '

if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "sgm av Intens" "08sgmintens" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,5)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Sgm avSuscp" "09sgmScpt" ',time_dum ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,7)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Sgm Vuln" "10sgmVuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,9)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif

!  -----------------------------  sigma of peak h, v , intens, suscept, Vulnerability 

write(gid_res,*) 'Result "sgm peak hh" "16Sg peak hh" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,2)   
   endif    
enddo
write(gid_res,*) ' End Values '

write(gid_res,*) 'Result "sgm peak vv" "17Sg peak vv" ',time_dummy ,'  Scalar OnNodes "where" '
write(gid_res,*) ' Values '
Do ipoig = 1, npoig
   if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
      write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,4)   
   endif    
enddo
write(gid_res,*) ' End Values '

if (ic_vuln.eq.1) then

   write(gid_res,*) 'Result "SpkInt" "18s pkint" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,6)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "S pk sc"  "19Sg PkSc" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,8)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

   write(gid_res,*) 'Result "Sgpk Vuln" "20Sgpk Vuln" ',time_dummy ,'  Scalar OnNodes "where" '
   write(gid_res,*) ' Values '
   Do ipoig = 1, npoig
      if (MC_pixel(ipoig) % if_MC_Pr_alloc.EQ.1) then
         write(gid_res,*) ipoig, MC_pixel(ipoig)  % MC_Pr(2+idest,10)   
      endif    
   enddo
   write(gid_res,*) ' End Values '

endif
!  ************************************************************

End Subroutine Plot_MC_pixel_Vars_SW 


END MODULE SPH_SW_2019