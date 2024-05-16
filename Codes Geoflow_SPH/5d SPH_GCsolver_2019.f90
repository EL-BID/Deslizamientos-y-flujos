MODULE GCsph_solvers_2019


USE SPH_FEM_variable_types_2019

implicit none
PRIVATE

Public :: solvegc, solvegc2, solve_AVR

CONTAINS


! --------------------------------------------------------- 

       subroutine solvegc2 ( npoin , nelem, nnode,ndofn,   &
                             ntotv , nevab,                &
                             intmat, iffix, fixed,         &
                             gstif , rcg  ,  ficg, tolCG   )  

! --------------------------------------------------------- 

!      ------  Solves the pressure laplacian by PJCG

implicit none

integer (ink),INTENT   (IN):: npoin !  nodes in the fem mesh
integer (ink),INTENT   (IN):: nelem !  elements 
integer (ink),INTENT   (IN):: nnode !  nodes per element 
integer (ink),INTENT   (IN):: ndofn !  dofs per node 
integer (ink),INTENT   (IN):: ntotv !  total nr of dofs = npoin*ndofn
integer (ink),INTENT   (IN):: nevab !  elem vars = nnode*ndofn
real    (irk),INTENT   (IN):: tolCG !  tolerance for CG iterations convergence

integer (ink),INTENT   (IN):: intmat (:,:) ! intmat(nnode,nelem)  conect.

real    (irk),INTENT   (IN):: fixed    (:) ! fixed (ntotv) value of pwp
real    (irk),INTENT   (IN):: gstif(:,:,:) ! gstif (nevab,nevab,nelemT)
                                           !  glb.K
real    (irk),INTENT  (OUT):: ficg     (:)  ! solution found by CG

integer (ink),INTENT(INOUT):: iffix    (:) ! iffix(ntotv) = 1 if presc.
                                           ! code changes it to -999 exception
real    (irk),INTENT(INOUT):: rcg      (:) ! rcg  (ntotv) RHS
                                           ! rhs at the beginning, code uses it
real    (irk),allocatable:: r1cg (:)   ! r1cg (ntotv) used in CG
real    (irk),allocatable:: scg  (:)   ! scg  (ntotv) used in CG
real    (irk),allocatable:: pcg  (:)   ! pcg  (ntotv) used in CG
real    (irk),allocatable:: q1cg (:)   ! q1cg (ntotv) used in CG
real    (irk),allocatable:: apcg (:)   ! apcg (ntotv) used in CG

integer (ink) miter     !  max number of iters. allowed
integer (ink) iiter     !  auxiliar counter of iterations

real    (irk) Zero      !  Aux zero, useful for calls if double precision
real    (irk) rnorm0, rnorm1, rnorre
real    (irk) rnormX
integer (ink) itotv

Zero = 0.0

if (.not.allocated (r1cg)) then   
    allocate (r1cg (ntotv)); r1cg = zero
    allocate  (scg (ntotv));  scg = zero  
    allocate  (pcg (ntotv));  pcg = zero  
    allocate (q1cg (ntotv)); q1cg = zero
    allocate (apcg (ntotv)); apcg = zero
endif

miter = 10*ntotv
zero  = 0.0

!      ------  Initializes

ficg   = zero
q1cg   = zero

!      ------  Obtains preconditioning matrix  q1cg(i,i)=q1cg(i)

call  getqmat                

!      ------  Initializes

call initcg  (rnorm0)

!      ------  avoid solving K.dp = 0 if fixed is 0.0

rnormX    = zero
do itotv  = 1, ntotv
   rnormX = rnormX + fixed(itotv)*fixed(itotv)
enddo

rnormX = rnormX + rnorm0
rnormX = sqrt(rnormX)
if (rnormX.lt.tolCG) then

   ficg = zero
   
else   

!      ------  Main CG iteration loop
 
   do iiter = 1,miter
 
      call jcg (rnorm0,rnorm1)

      rnorre = abs(rnorm1/rnorm0)

      !  write(*,*) ' iter = ',iiter,' norm.rel.= ',rnorre, &
      !             ' res= ',rnorm1

      if((rnorm1/rnorm0).le.tolCG)  EXIT

   enddo

where (iffix.eq.-999) iffix = 0

endif

CONTAINS

!---------------------------------------------------------- 

      subroutine getqmat  

!----------------------------------------------------------- 


!   ------- Get  q1 = diag(k)-1 matriz


implicit none

integer (ink) ielem, inode, idofn, ievab, lnode,nposn, itotv

real    (irk) acc 

do  ielem = 1,nelem 
    do  inode = 1,nnode 
        do  idofn = 1,ndofn
            ievab = (inode-1)*ndofn + idofn
            acc   = gstif(ievab,ievab,ielem)
            lnode = intmat (inode,ielem)
            nposn = (lnode-1)*ndofn + idofn
            q1cg(nposn) = q1cg(nposn) + acc
        enddo
    enddo
 enddo

do itotv = 1,ntotv
   if(abs(q1cg(itotv)).le.(very_small)) then	! ++++  MP
       write(*,*) ' var=0 at dofn ',itotv
       if (iffix(itotv).eq.0) iffix(itotv)=-999
   else
       q1cg(itotv) = 1./q1cg(itotv)
   endif
 enddo

END SUBROUTINE getqmat 

!------------------------------------------------- 

      subroutine initcg  (rnorm0)

!-------------------------------------------------- 


!      ------  Initialize

implicit none

integer (ink) itotv
real    (irk) rnorm0 

!      ------  initialize all vectors that matters

r1cg = zero
scg  = zero
pcg  = zero
apcg = zero

!      ------  Inits fi & pcg

do  itotv = 1,ntotv
    ficg(itotv) = fixed(itotv)
    pcg (itotv) = ficg (itotv)
enddo

call kbyp ( 0 )

!      ------  Inits rcg,scg,pcg and get norm of r0

rnorm0   = zero
do itotv = 1,ntotv
    rcg(itotv) = rcg(itotv) - apcg(itotv)
    scg(itotv) = q1cg(itotv) * rcg(itotv)
    pcg(itotv) = scg(itotv)
    rnorm0 = rnorm0 + rcg(itotv)*rcg(itotv)
enddo

!      ------  zero prescribed dof components

do itotv=1,ntotv
   if(iffix(itotv).ne.0) then
      pcg  (itotv) = zero
      rcg  (itotv) = zero
      scg  (itotv) = zero
   endif
enddo


END SUBROUTINE initcg 


!-------------------------------------------------- 
      subroutine jcg (rnorm0,rnorm1)

!--------------------------------------------------- 
 
!      ------  GC solver

implicit none

integer (ink) itotv 
real    (irk) rnorm0, rnorm1, prod1, prod2, prod3
real    (irk) alfa  ,   beta 

!      ------  Get k * p

call       kbyp ( 1 )

!      ------  get alfa = (r * s/p *ap)

prod1 =  zero
prod2 =  zero

do  itotv = 1,ntotv
    if (iffix(itotv).eq.0) then
       prod1 = prod1 + rcg(itotv)*scg(itotv)
       prod2 = prod2 + pcg(itotv)*apcg(itotv)
    endif
enddo
alfa = prod1/prod2

!      ------  Computes  fi = fi + alfa*p

prod3  = zero
rnorm1 = zero
do  itotv = 1,ntotv
    ficg (itotv) = ficg (itotv) + alfa * pcg (itotv)
    r1cg (itotv) =  rcg (itotv) - alfa * apcg(itotv)
    rnorm1     = rnorm1 + r1cg (itotv) * r1cg(itotv)
    scg(itotv) =          q1cg (itotv) * r1cg(itotv)
    prod3      = prod3 +  r1cg (itotv) * scg (itotv)
enddo

!      ------  Checks convergence

if((rnorm1/rnorm0).gt.tolcg) then
    beta = prod3/prod1                     ! --- beta = (r1*s1/r0*s0)
    do  itotv = 1,ntotv                    ! --- Nuevo p & r
        pcg (itotv) = scg (itotv) + beta*pcg(itotv)
        rcg (itotv) = r1cg(itotv)
    enddo
endif 


END SUBROUTINE jcg


!------------------------------------------------------ 

      subroutine  kbyp ( iflag )

!------------------------------------------------------- 
 
!      ------     ap = k * p  

implicit none

integer (ink) iflag, itotv, ielem  
integer (ink) inode, idofn, ievab, lnode, nposn
integer (ink) jnode, jdofn, jevab, knode, mposn 
real    (irk) sum


do itotv =1,ntotv
   apcg(itotv) = zero
enddo

do ielem = 1,nelem
   do inode = 1,nnode
      do idofn = 1,ndofn
         ievab = (inode-1)*ndofn + idofn
         lnode = intmat(inode,ielem)
         nposn = (lnode-1)*ndofn + idofn
         do jnode = 1,nnode
            do jdofn = 1,ndofn
               jevab = (jnode-1)*ndofn + jdofn
               knode = intmat(jnode,ielem)
               mposn = (knode-1)*ndofn + jdofn
               sum   = gstif(ievab,jevab,ielem)*pcg(mposn)
               if(iflag.eq.1) then
                  if(iffix(nposn)+iffix(mposn).ne.0) sum= zero
               endif
               apcg(nposn) = apcg(nposn) + sum
            enddo
         enddo 
      enddo
   enddo
 enddo

END SUBROUTINE kbyp  


END SUBROUTINE  solvegc2



! --------------------------------------------------------- 

       subroutine solvegc (  npoin , ntotv ,       &
                             iffix, fixed, gstif,  &
                             rcg  ,  ficg, tolCG   )  

! --------------------------------------------------------- 


!      ------  Solves the pressure laplacian by PJCG

implicit none

integer (ink),INTENT   (IN):: npoin !  nodes in the fem mesh
integer (ink),INTENT   (IN):: ntotv !  total nr of dofs = npoin*ndofn
real    (irk),INTENT   (IN):: tolCG !  tolerance for CG iterations convergence


real    (irk),INTENT   (IN):: fixed    (:) ! fixed (ntotv) value of pwp
real    (irk),INTENT   (IN):: gstif(:,:)   ! gstif (npoin1,npoin1)
                                           !  glb.K
real    (irk),INTENT  (OUT):: ficg     (:)  ! solution found by CG

integer (ink),INTENT(INOUT):: iffix    (:) ! iffix(ntotv) = 1 if presc.
                                           ! code changes it to -999 exception
real    (irk),INTENT(INOUT):: rcg      (:) ! rcg  (ntotv) RHS
                                           ! rhs at the beginning, code uses it

real    (irk),allocatable:: r1cg (:)   ! r1cg (ntotv) used in CG
real    (irk),allocatable:: scg  (:)   ! scg  (ntotv) used in CG
real    (irk),allocatable:: pcg  (:)   ! pcg  (ntotv) used in CG
real    (irk),allocatable:: q1cg (:)   ! q1cg (ntotv) used in CG
real    (irk),allocatable:: apcg (:)   ! apcg (ntotv) used in CG

integer (ink) miter     !  max number of iters. allowed
integer (ink) iiter     !  auxiliar counter of iterations


real    (irk) Zero      !  Aux zero, useful for calls if double precision
real    (irk) rnorm0, rnorm1, rnorre
real    (irk) rnormX
integer (ink) itotv

Zero = 0.0

if (.not.allocated (r1cg)) then   
    allocate (r1cg (ntotv)); r1cg = zero
    allocate  (scg (ntotv));  scg = zero  
    allocate  (pcg (ntotv));  pcg = zero  
    allocate (q1cg (ntotv)); q1cg = zero
    allocate (apcg (ntotv)); apcg = zero
endif

miter = 10*ntotv
zero  = 0.0

!      ------  Initializes

ficg   = zero
q1cg   = zero

!      ------  Obtains preconditioning matrix  q1cg(i,i)=q1cg(i)

call  getqmat                

!      ------  Initializes

call initcg  (rnorm0)

!      ------  avoid solving K.dp = 0 if fixed is 0.0

rnormX    = zero
do itotv  = 1, ntotv
   rnormX = rnormX + fixed(itotv)*fixed(itotv)
enddo

rnormX = rnormX + rnorm0
rnormX = sqrt(rnormX)
if (rnormX.lt.tolCG) then

   ficg = zero
   
else   

!      ------  Main CG iteration loop
 
   do iiter = 1,miter
 
      call jcg (rnorm0,rnorm1)

      rnorre = abs(rnorm1/rnorm0)

      !  write(*,*) ' iter = ',iiter,' norm.rel.= ',rnorre, &
      !             ' res= ',rnorm1

      if((rnorm1/rnorm0).le.tolCG)  EXIT

   enddo

where (iffix.eq.-999) iffix = 0

endif


CONTAINS

!---------------------------------------------------------- 

      subroutine getqmat  

!----------------------------------------------------------- 


!   ------- Get  q1 = diag(K)-1 matriz


implicit none

integer (ink) ipoin, itotv

Do ipoin =1,npoin
    q1cg(ipoin) = gstif(ipoin,ipoin)
Enddo


do itotv = 1,ntotv
   if(abs(q1cg(itotv)).le.(1.e-9)) then	! ++++  MP
       write(*,*) ' var=0 at dofn ',itotv
       if (iffix(itotv).eq.0) iffix(itotv)=-999
   else
       q1cg(itotv) = 1./q1cg(itotv)
   endif
 enddo

END SUBROUTINE getqmat


!------------------------------------------------- 

      subroutine initcg  (rnorm0)

!-------------------------------------------------- 


!      ------  Initialize

implicit none

integer (ink) itotv
real    (irk) rnorm0 

!      ------  initialize all vectors that matters

r1cg = zero
scg  = zero
pcg  = zero
apcg = zero

!      ------  Inits fi & pcg

do  itotv = 1,ntotv
    ficg(itotv) = fixed(itotv)
    pcg (itotv) = ficg (itotv)
enddo

call kbyp ( 0 )

!      ------  Inits rcg,scg,pcg and get norm of r0

rnorm0   = zero
do itotv = 1,ntotv
    rcg(itotv) = rcg(itotv) - apcg(itotv)
    scg(itotv) = q1cg(itotv) * rcg(itotv)
    pcg(itotv) = scg(itotv)
    rnorm0 = rnorm0 + rcg(itotv)*rcg(itotv)
enddo

!      ------  zero prescribed dof components

do itotv=1,ntotv
   if(iffix(itotv).ne.0) then
      pcg  (itotv) = zero
      rcg  (itotv) = zero
      scg  (itotv) = zero
   endif
enddo


END SUBROUTINE initcg 

!-------------------------------------------------- 
      subroutine jcg (rnorm0,rnorm1)

!--------------------------------------------------- 
 
!      ------  GC solver

implicit none

integer (ink) itotv 
real    (irk) rnorm0, rnorm1, prod1, prod2, prod3
real    (irk) alfa  ,   beta 

!      ------  Get k * p

call       kbyp ( 1 )

!      ------  get alfa = (r * s/p *ap)

prod1 =  zero
prod2 =  zero

do  itotv = 1,ntotv
    if (iffix(itotv).eq.0) then
       prod1 = prod1 + rcg(itotv)*scg(itotv)
       prod2 = prod2 + pcg(itotv)*apcg(itotv)
    endif
enddo
alfa = prod1/prod2

!      ------  Computes  fi = fi + alfa*p

prod3  = zero
rnorm1 = zero
do  itotv = 1,ntotv
    ficg (itotv) = ficg (itotv) + alfa * pcg (itotv)
    r1cg (itotv) =  rcg (itotv) - alfa * apcg(itotv)
    rnorm1     = rnorm1 + r1cg (itotv) * r1cg(itotv)
    scg(itotv) =          q1cg (itotv) * r1cg(itotv)
    prod3      = prod3 +  r1cg (itotv) * scg (itotv)
enddo

!      ------  Checks convergence

if((rnorm1/rnorm0).gt.tolcg) then
    beta = prod3/prod1                     ! --- beta = (r1*s1/r0*s0)
    do  itotv = 1,ntotv                    ! --- Nuevo p & r
        pcg (itotv) = scg (itotv) + beta*pcg(itotv)
        rcg (itotv) = r1cg(itotv)
    enddo
endif 


END SUBROUTINE jcg

!------------------------------------------------------ 

      subroutine  kbyp ( iflag )

!------------------------------------------------------- 
 
!      ------     ap = K * p  

implicit none

integer (ink) iflag, itotv  
integer (ink) ipoin
integer (ink) jpoin
real    (irk) sum


do itotv =1,ntotv
   apcg(itotv) = zero
enddo

Do ipoin = 1, npoin
    Do jpoin = 1, npoin
        sum   = gstif(ipoin,jpoin)*pcg(jpoin)
        if(iflag.eq.1) then
            if(iffix(ipoin)+iffix(jpoin).ne.0) sum= zero
        endif
        apcg(ipoin) = apcg(ipoin) + sum
     Enddo
Enddo

END SUBROUTINE kbyp  



END SUBROUTINE solvegc 


! -------------------------------------------------------------------

       subroutine solve_AVR     (npoin, nelem, nnode, ndimn, &
                                 namat, ngeom, nprer, niter, &
                                 intmat,geome, bprer, mmatl, &
                                 rhs  , Fi_pres )

! -------------------------------------------------------------------


!      npoin                   number of nodes in mesh
!      nelem                   number of elements
!      nnode                   nodes per element
!      ndimn                   dimension
!      namat                   nr of dofn in the unknowns array
!      ngeom                   dimension of geome, store Ni,j and area
!      niter                   AVR number of iterations
!      nprer                   number of prescribed DOFn
!      intmat (nnode,nelem)    Connectivity matrix
!      geome  (ngeom,nelem)    N1x..Nnx,N1y..Nny,...area (2d is 2Area)
!      bprer      (2,nprer)    node and dofn prescribed for var iprer
!      mmatl        (npoin)    Lumped mass matrix
!
!      rhs    (namat,npoin)    Input : Right hand side of the eqns system
!      Fi_pres(namat,npoin)    In/Out: prescribed values and solution

!      Fi     (namat,npoin)    Aux   : Solution of the system  (local)
!      M_by_Fi(namat,npoin)    auxiliar array M*Fi             (local)
!      mmat   (nnode,nnode)    Lumped mass matrix              (local) 
!
!      We solve  M.Fi = rhs  with prescribed dofn given in bprer
!         where  M  is the mass matrix M = int(Ni Nj)
!
!         algorithm is:
!                    Fi_n+1 = Fi_n + inv(Ml)* ( rhs - M.Fi_n)
!         stored as:
!                    Fi = Fi + inv(Ml)* ( rhs - M_by_Fi)
!
!      In the case of TG algorithm what we solve is
!
!                    dFi = dFi + inv(Ml)* ( rhs - M_by_dFi) 
!
!    ** Output is (i) Fi and Fi_pres, where we merge fi and fi_pres **
 
implicit none

integer(ink) npoin, nelem, nnode, ndimn,namat, ngeom, niter, nprer
integer(ink) bprer (2,nprer)
integer(ink) iiter, iamat, ielem, inode, ipoin,jnode, jpoin   
integer(ink) iunkn, ipres

integer(ink) intmat(nnode,nelem)

real   (irk) geome (ngeom,nelem), mmatl (npoin)
real   (irk) rhs (namat,npoin) , Fi_Pres (namat,npoin)

real   (irk) Fi(namat,npoin), M_by_Fi (namat,npoin)  		! INTERNALS
real   (irk) mmat(nnode,nnode)

real   (irk) cm, coe

!      ------  Get Mass matrix

if (ndimn.eq.1.or.ndimn.eq.2) then  !  *** check Mmat for tetrahedra
   mmat = 1.0
   do inode = 1,nnode
      mmat(inode,inode) = 2.0
   enddo
endif

Fi =  0.0

DO iiter = 1,niter

   IF (iiter.eq.1) THEN
       do ipoin = 1,npoin
          cm    = mmatl(ipoin)
          do iamat = 1,namat
             Fi(iamat,ipoin) = rhs(iamat,ipoin)/cm
          enddo
       enddo
   ELSE
       M_by_Fi =  0.0
       
       do ielem = 1,nelem 
          if (nnode.eq.3) then 
             coe = geome (7,ielem)/24.   !  2D triangles
          elseif (nnode .eq.2) then
             coe = geome (3,ielem)/6.    !  1D linear
          elseif (nnode .eq.4) then
             write(*,*) ' solve...check 3D  '
             pause
             stop
          endif
          do inode = 1,nnode 
             ipoin = intmat (inode,ielem)
             do jnode = 1,nnode 
                jpoin = intmat (jnode,ielem)
                cm    = mmat(inode,jnode)
                do iamat = 1,namat 
                   M_by_Fi(iamat,ipoin) = M_by_Fi (iamat,ipoin)      & 
                                        + cm * coe * Fi (iamat,jpoin)
                enddo
             enddo
          enddo
       enddo       
       
       do ipoin = 1,npoin 
          cm = mmatl (ipoin)
          do iamat = 1,namat 
             Fi (iamat,ipoin) = Fi (iamat,ipoin) & 
                                + (rhs (iamat,ipoin)-M_by_Fi(iamat,ipoin))/cm
          enddo
       enddo

   ENDIF

!      ------  Apply BC's  prescribed values

   IF (nprer .ne.0) THEN
      do ipres = 1,nprer 
         ipoin = bprer (1,ipres)
         iunkn = bprer (2,ipres) 
         Fi (iunkn,ipoin) = Fi_pres (iunkn,ipoin)
      enddo
   ENDIF
   
 

ENDDO                   !      ------  ends loop in iterations

!        transcribe Fi onto Fi_pres

Fi_pres  = Fi 


END SUBROUTINE solve_AVR 



END MODULE GCsph_solvers_2019
