module SPH_FEM_variable_types_2019

implicit none

integer,parameter::k_min=selected_int_kind(1),                               &
                   k_short=selected_int_kind(4),                             &
                   k_int=kind(0),                                            &
                   k_long=selected_int_kind(range(0)+1),                     &
                   k_real=kind(0.0),                                         &
                   k_double=kind(0d0),                                       &
                   k_quad=selected_real_kind(precision(0d0)+1)
                   
!  integer,parameter::ink=k_int,irk=k_double
integer,parameter::ink=k_int,irk=k_real

real    (irk),parameter   :: zero  = 0.0  !  definition to avoid type incompatibility
integer (ink),parameter   :: zeroI = 0    !  definition to avoid type incompatibility
integer (ink),parameter   :: oneI  = 1    !  definition to avoid type incompatibility
integer (ink),parameter   :: twoI  = 2    !  definition to avoid type incompatibility

real    (irk),parameter   :: very_small = 1.e-35 	!  definition of a very small quantity
real    (irk), parameter  :: very_large = 1.e35  	!  idem. very large                   


end module SPH_FEM_variable_types_2019 
