!--------------------------------------------------------


                 MODULE GFL_Global_Vars_2019 
                 

!--------------------------------------------------------

! comes from geoflow 2005
 
                       
USE SPH_FEM_variable_types_2019

implicit none

real   (irk), allocatable:: coord    (:,:)      ! -- coords
!real   (irk), allocatable:: topol    (:,:)      ! -- topo data: Z, Zx Zy Zxx...(nodes)
real   (irk), allocatable:: topel    (:,:)      ! -- topo data: Z, Zx Zy Zxx...(elems)
real   (irk), allocatable:: geome    (:,:)      ! --
real   (irk), allocatable:: mmatl    (:)        ! --
real   (irk), allocatable:: geomb    (:,:)      ! --
real   (irk), allocatable:: rtanr    (:,:)      ! --
real   (irk), allocatable:: rtang    (:,:)      ! --
real   (irk), allocatable:: unkno    (:,:)      ! -- Variable de campo
real   (irk), allocatable:: unkne    (:,:)      ! -- Variable de campo
real   (irk), allocatable:: unpre    (:,:)      ! --
real   (irk), allocatable:: delun    (:,:)      ! --
real   (irk), allocatable:: rhs0     (:,:)      ! --
real   (irk), allocatable:: rhs1     (:,:)      ! --
real   (irk), allocatable:: rhs2     (:,:)      ! --
real   (irk), allocatable:: gelem    (:,:)      ! --
real   (irk), allocatable:: gpoin    (:,:)      ! --
real   (irk), allocatable:: fluxp    (:,:)      ! --
real   (irk), allocatable:: fluyp    (:,:)      ! --
real   (irk), allocatable:: patma    (:,:)      ! --
real   (irk), allocatable:: patmel   (:,:)      ! --
real   (irk), allocatable:: constF   (:)        ! --
real   (irk), allocatable:: twind    (:)        ! --
real   (irk), allocatable:: datm     (:)        ! --
real   (irk), allocatable:: wave     (:)        ! --
real   (irk), allocatable:: deltel   (:)        ! -- optimal dt at elems
real   (irk), allocatable:: denom    (:)        ! --

integer (ink), allocatable:: intmat  (:,:)      ! -- connectivity (elems)
integer (ink), allocatable:: intmab  (:,:)      ! -- connect. b.sides
integer (ink), allocatable:: btanr   (:)        ! --
integer (ink), allocatable:: ltang   (:)        ! --
integer (ink), allocatable:: bconl   (:)        ! -- list of bound. nodes
integer (ink), allocatable:: bprer   (:,:)      ! --
integer (ink), allocatable:: babhz   (:)        ! --
integer (ink), allocatable:: bload   (:)        ! --
integer (ink), allocatable:: babso   (:)        ! -- list of B abs. nodes
integer (ink), allocatable:: lprpoi  (:)        ! -- list output points
integer (ink), allocatable:: lweir   (:)        ! -- list of weir nodes
integer (ink), allocatable:: iactnd  (:)        ! -- index active nodes
integer (ink), allocatable:: iactel  (:)        ! -- index active elements
integer (ink), allocatable:: bindex  (:)        ! -- index for boundary nodes

! ---------  general -----

integer (ink) nsour      ! --  Algorithm for sources: 0 nothing 2 RK split
                         !     we transfer to it GFL_t_Integ_Alg
		          
integer (ink) ndimnF     ! --
integer (ink) namat      ! --
integer (ink) nelem      ! --
integer (ink) npoinF     ! --        
integer (ink) nnode      ! --
integer (ink) ngaus      ! --
integer (ink) ntotg      ! --
		         ! --
! integer (ink) ntopo      ! -- defined in topo module
integer (ink) nconstF    ! --
integer (ink) ngeom      ! --
integer (ink) ngeob      ! --
integer (ink) nelemb     ! --
integer (ink) neleb      ! --
integer (ink) nboun      ! --
		          
integer (ink) niter      ! --
		          
integer (ink) nunpre     ! --
integer (ink) nbour      ! --
integer (ink) nprer      ! --          
integer (ink) nload	     ! -- defined in topo module
integer (ink) nabso      ! --
integer (ink) nabshz     ! --
integer (ink) ntypbc     ! --

integer (ink) csmooth	 ! --  
		        
integer (ink) nprpoi     ! --  Number of points at which we want output

real   (irk) time_cpu0, time_cpuf		  !    variable para el tiempo de cálculo

integer(ink) datg_file, chkg_file
integer(ink) gidg_msh,  gidg_res, gid3g_msh, gid3g_res 	 ! files for input/output
integer(ink) resg_file1,resg_file2, resg_file3           ! files for input/output

real(irk) pi

character(60) text                    ! general purpose char string


!  ---------   SW only -----

integer (ink) icpwpF         ! --  1 if we are using pore press, 0 otherwise
			     
integer (ink) ntype          ! --     
integer (ink) icoutp	     ! --
			      
integer (ink) nweir	                        ! -- number of weir BCs points
real    (irk) cweir(2)                          ! -- 2 cts: v = c1 * href^c2
real   (irk), allocatable:: zweir    (:)        ! -- z of weir pts

integer (ink) npatma	     ! --
integer (ink) ntwind	     ! --
integer (ink) ndatm    	     ! --
integer (ink) nwave 	     ! --
			      
integer (ink) nacti   	     ! --
			      
real    (irk) comp	         ! --
real    (irk) swl   
			     
			     
			     
END MODULE GFL_Global_Vars_2019 
			     
			     
			     
			     
			     
			     
			     