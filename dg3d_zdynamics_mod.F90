#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg3d_zdynamics_mod
!=======================================================================================================!
! Author: R. D. Nair (rnair@ucar.edu) 10/28/2014 , 03/07/2015 
!=======================================================================================================!
  use kinds,           only: real_kind, longdouble_kind 
  use dimensions_mod,  only: ne, np, nlev
  use element_mod,     only: element_t
  use control_mod, only : test_case
  use derivative_mod,  only: derivative_t
  use quadrature_mod,  only: quadrature_t, gausslobatto, jacobi, gauss, legendre
  use physical_constants, only : rearth , rrearth, g, dd_pi
  use dg_userdata_mod, only : dg_userdata_t
  use dg_core_mod, only : sphere2contra, contra2sphere 
!=======================================================================================================!
   implicit none
!------------------------------------------------------------------------------------------------
!        For vertical z-coordinates 
!------------------------------------------------------------------------------------------------
   Integer, public, parameter :: vnp = 4   !! = np as a default,  may be changed
   Integer, public, parameter :: vnel = nlev/vnp   ! Note: nlev should be divisible by vnp 

   real (kind=real_kind), public, parameter :: z_top = 12000.0D0   ! default, dcmip200 
  !real (kind=real_kind), public, parameter :: z_top = 10000.0D0   ! for nhgw dcmip-31 tests 
  !real (kind=real_kind), public, parameter :: z_top = 30000.0D0   ! for mtn dcmip-21,22 
  !real (kind=real_kind), public, parameter :: z_top = 44000.0D0   ! for JW_bcl dcmip-4x     
  !real (kind=real_kind), public, parameter :: z_top = 44307.0D0   ! for 3d-Rossby Haurwitz  
  !real (kind=real_kind), public, parameter :: z_top = 30000.0D0   ! for 3d-Rossby Haurwitz  

   real (kind=real_kind), public, parameter :: z_bot = 0.0D0
   real (kind=real_kind), public :: z_wel  ! width
   real (kind=real_kind), public :: vder(vnp,vnp),  vgp(vnp), vgw(vnp)
   real (kind=real_kind), public :: z_grd(nlev)     ! 1d vertical grid
   real (kind=real_kind), public :: z_gsg(vnp,vnel) ! Gauss  2d grid (GL or GLL) 
   real (kind=real_kind), public :: ime_gl(vnp,2)   ! interpolation matrix GL
   real (kind=real_kind), public :: hed_gl(vnp,2)   ! boundary test-fn for GL
   real (kind=real_kind), public :: rad_der(vnp,2)  ! differentiated Radau edge values 
   real (kind=real_kind), public :: z_wts(vnp,nlev) ! vertical interpolation weights   
   real (kind=real_kind), public :: zsponge_coef(nlev)! vertical sponge coeffeicients    
   real (kind=real_kind), public :: zsponge_level   ! z-level from which sponge layer starts   
   integer, public :: vindx(nlev)                   ! vertical uniform output/grid index

  !Thermodynamic constants 
   real (kind=real_kind), public, parameter  :: Cp_by_Cv = 1004.0D0/717.0D0   ! gama = cp/cv       
   real (kind=real_kind), public, parameter  :: CC_0     = 27.5629410929719D0 !rd**gama * p_0**(-rd/cv) 

   Character(len=15), public :: zgrid_type 
   Character(len=15), public :: test_type 
   Integer, private, parameter :: ibot=1, itop=2 
!=======================================================================================================! 
  public  :: vertical_GLL_zgrid, vertical_FD_zgrid ,vertical_GL_zgrid, edge_GL_1d
  public  :: Zsponge_function_top
  public  :: print_vertical_zinfo
  public  :: Vertical_Euler_Forcing
  public  :: vertical_solve_update         
  public  :: vertical_advect_fdsrc , vertical_advect_dgsrc
  public  :: vertical_z_integral, vertical_unif_interpolate 
  public  :: Compute_Pressure, Compute_SoundSpeed
  public  :: Vertically_ExplicitRK3_EulerSolver            

  private :: copyto_1d_elem, copyfrom_1d_elem, copyfrom_1dsys_elem
  private :: DG1d_GLL_rhs, dg1d_rk3_solver, fd1d_rhs_src, fd1d_rk3_solver
  private :: DG1d_compute_rhs, Apply_TopSponge_Layer
  private :: Vertical_DGEuler_Solver,Vertical_LXF_Flux, Compute_Vgrad_term, TopBot_noflux_BC
  private :: Recover_Primary_NHvars, Compute_EulerSys_rhs
  private :: Decompose2_Column_elements, Gather_TimeIndep_DatColumn, Copyfrom_Column2_elem2D
!=======================================================================================================!
     Contains


!=======================================================================================================!
 Subroutine Vertically_ExplicitRK3_EulerSolver(ie,userdata,dtv,neqn)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie,neqn
    type (dg_userdata_t) , intent(inout) :: userdata
    type (element_t)     , pointer :: elem(:)
    real (kind=real_kind), intent(in)    :: dtv      
    real (kind=real_kind),dimension(vnp,vnel)  :: rho,w,u1,u2,pt
    real (kind=real_kind),dimension(vnp,vnel)  :: prp,rhp,prb,rhb,rtb 
    real (kind=real_kind),dimension(vnp,vnel)  :: sg13,sg23, speed 
    real (kind=real_kind),dimension(np,np,2,nlev):: uv
    real (kind=real_kind),dimension(np,np,2,nlev) :: contuv
    real (kind=real_kind), dimension(2,2,np,np) :: A_inv, Amat 
    real (kind=real_kind), dimension(vnp,vnel,neqn) :: vrhs, sp_vec0 
    real (kind=real_kind), dimension(vnp,vnel,neqn) :: rk0_vec, rk1_vec, rhs
    real (kind=real_kind) :: s1, trm, sg,sgv, ginv(2,2)  
    integer :: i,j, ke, k, eq 

       elem  => userdata%elem

 !! Extracting primary variables from storage 

    A_inv(:,:,:,:) = elem(ie)%ainv(:,:,:,:)
     Amat(:,:,:,:) = elem(ie)%amat(:,:,:,:)

 ! Converting (u,v) into contravariant velocity for the vertical Euer solver 
    do k =1, nlev 
     contuv(:,:,:,k) = sphere2contra(elem(ie)%state%uv(:,:,:,k),A_inv)
    enddo 

 !! For each element,  vertical column-wise ODE solver 

   do j =1, np
   do i =1, np
   
     ginv(:,:) = elem(ie)%g_ijinv(:,:,i,j)

  !Transferring data to 1D vertical elements 
     Call  Decompose2_Column_elements(ie,i,j,userdata,contuv,rho,w,u1,u2,pt,rhb,prb,rtb)
     Call  Gather_TimeIndep_DatColumn(ie,i,j,neqn,userdata,sg,sgv,sg13,sg23,sp_vec0)


 !! Constructing "U" in the ODE dU/dt = RHS,  U=[ro',u1ro,u2ro,wro,thp']*sg 
     Call Construct_SateVec_NHvars(neqn,sg,u1,u2,rho,w,pt,rhb,prb,rtb,rhp,prp,speed,rk0_vec)

  !!RK_stage-1
     Call Vertical_DGEuler_solver(neqn,sg,sgv,ginv,sg13,sg23,rk0_vec,speed, &
                                              rhp,u1,u2,w,rho,pt,prp,rhs)
    do eq=1,neqn 
      rk1_vec(:,:,eq) = rk0_vec(:,:,eq) + 1.0*dtv * rhs(:,:,eq)
    enddo

     Call Recover_Primary_NHvars(neqn,sg,prb,rhb,rtb,rk1_vec,u1,u2,rho,w,pt,rhp,prp,speed)

  !!RK_stage-2
     Call Vertical_DGEuler_solver(neqn,sg,sgv,ginv,sg13,sg23,rk1_vec,speed, &
                                              rhp,u1,u2,w,rho,pt,prp,rhs)
     do eq=1,neqn 
      rk1_vec(:,:,eq) = (3.0D0*rk0_vec(:,:,eq) + rk1_vec(:,:,eq) + dtv * rhs(:,:,eq)) *0.25D0
     enddo

     Call Recover_Primary_NHvars(neqn,sg,prb,rhb,rtb,rk1_vec,u1,u2,rho,w,pt,rhp,prp,speed)

  !!RK_stage-3
     Call Vertical_DGEuler_solver(neqn,sg,sgv,ginv,sg13,sg23,rk1_vec,speed, &
                                              rhp,u1,u2,w,rho,pt,prp,vrhs)
     do eq=1,neqn 
      rk1_vec(:,:,eq) = (rk0_vec(:,:,eq) + 2.0D0*(rk1_vec(:,:,eq) + dtv * rhs(:,:,eq)) ) / 3.0D0
     enddo

   Call Recover_Primary_NHvars(neqn,sg,prb,rhb,rtb,rk1_vec,u1,u2,rho,w,pt,rhp,prp,speed)

!! Updating the changes by the vertical Euler solver on NH primary variables 

   Call Copyfrom_Column2_elem2D(ie,i,j,userdata,contuv,rho,w,u1,u2,pt)

  end do 
  end do 

 ! Converting back to spherical (u,v) components
    do k =1, nlev 
     elem(ie)%state%uv(:,:,:,k) = contra2sphere(contuv(:,:,:,k),Amat)
    enddo 

 End Subroutine Vertically_ExplicitRK3_EulerSolver

!=======================================================================================================!
 Subroutine Copyfrom_Column2_elem2D(ie,i,j,userdata,contuv,rho,w,u1,u2,pt)
!=======================================================================================================!
    Implicit None
    integer, intent(in)  :: ie,i,j 
    type (dg_userdata_t) , intent(inout) :: userdata
    real (kind=real_kind), dimension(vnp,vnel), intent(in)  :: rho,w,u1,u2,pt
    real (kind=real_kind), dimension(np,np,2,nlev), intent(inout) :: contuv
    type (element_t), pointer :: elem(:)

    integer :: k,el,lev 

       elem  => userdata%elem

      do el=1,vnel
      do  k=1,vnp
       lev = k + (el-1) * vnp    !index change 
         elem(ie)%state%rho(i,j,lev) =  rho(k,el)
          contuv(i,j,1,lev) = u1(k,el)
          contuv(i,j,2,lev) = u2(k,el)
          elem(ie)%state%w3d(i,j,lev) = w(k,el) 
          elem(ie)%state%pt(i,j,lev) =  pt(k,el) 
      end do
      end do

 End Subroutine Copyfrom_Column2_elem2D     

!=======================================================================================================!
 Subroutine Decompose2_column_elements(ie,i,j,userdata,contuv,rho,w,u1,u2,pt,rhb,prb,rtb)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie,i,j 
    type (dg_userdata_t) , intent(in) :: userdata
    real (kind=real_kind),dimension(np,np,2,nlev), intent(in) :: contuv
    real (kind=real_kind),dimension(vnp,vnel),intent(out)  :: rho,w,u1,u2,pt, prb,rtb,rhb
    type (element_t)     , pointer :: elem(:)

    integer :: k,el,lev 

       elem  => userdata%elem

      do el=1,vnel
      do  k=1,vnp
       lev = k + (el-1) * vnp    !index change 
         rho(k,el) = elem(ie)%state%rho(i,j,lev)
          u1(k,el) = contuv(i,j,1,lev)
          u2(k,el) = contuv(i,j,2,lev)
           w(k,el) = elem(ie)%state%w3d(i,j,lev)
          pt(k,el) = elem(ie)%state%pt(i,j,lev)
         rhb(k,el) = elem(ie)%state%rhb(i,j,lev)
         prb(k,el) = elem(ie)%state%prb(i,j,lev)
         rtb(k,el) = elem(ie)%state%rtb(i,j,lev)
      end do
      end do

 End Subroutine Decompose2_column_elements

!=======================================================================================================!
 Subroutine Gather_TimeIndep_DatColumn(ie,i,j,neqn,userdata,sg,sgv,sg13,sg23,sponge_vec)            
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie,i,j,neqn 
    type (dg_userdata_t) , intent(in) :: userdata
    real (kind=real_kind),dimension(vnp,vnel),intent(out)  :: sg13,sg23
    real (kind=real_kind),dimension(vnp,vnel,neqn),intent(out)  :: sponge_vec 
    real (kind=real_kind),intent(out)  :: sg,sgv      
    type (element_t)     , pointer :: elem(:)

    integer :: k,el,lev 

       elem  => userdata%elem

        sg   = elem(ie)%sg(i,j)
        sgv  = elem(ie)%sgv(i,j)

      do el=1,vnel
      do  k=1,vnp
       lev = k + (el-1) * vnp    !index change 
         sponge_vec(k,el,:) = elem(ie)%state%vec0(i,j,:,lev)
           sg13(k,el) = elem(ie)%sgv13(i,j,lev)
           sg23(k,el) = elem(ie)%sgv23(i,j,lev)
      end do
      end do

 End Subroutine Gather_TimeIndep_DatColumn
!!=======================================================================================================!
! Subroutine Vertically_ExplicitRK3_EulerSolver(ie,userdata,dtv,neqn)
!!=======================================================================================================!
!    Implicit None
!    integer, intent(in)               :: ie,neqn
!    type (dg_userdata_t) , intent(inout) :: userdata
!    type (element_t)     , pointer :: elem(:)
!    real (kind=real_kind), intent(in)    :: dtv      
!    real (kind=real_kind),dimension(np,np,nlev)  :: rho3d,w3d,the,prs
!    real (kind=real_kind),dimension(np,np,2,nlev):: uv
!    real (kind=real_kind),dimension(np,np,2,nlev) :: contuv
!    real (kind=real_kind), dimension(2,2,np,np) :: A_inv, Amat 
!    real (kind=real_kind), dimension(vnp,vnel,neqn) :: vrhs
!    real (kind=real_kind), dimension(np,np,neqn,nlev) :: rk0_vec, rk1_vec, rhs
!    real (kind=real_kind), dimension(np,np) :: sg
!    real (kind=real_kind) :: s1, trm
!    integer :: i,j, ke, k, eq 
!
!       elem  => userdata%elem
!
! !! Extracting primary variables from storage 
!
!        do k=1,nlev
!          do j=1,np
!          do i=1,np
!             rho3d(i,j,k) = elem(ie)%state%rho(i,j,k) 
!              uv(i,j,1,k) = elem(ie)%state%uv(i,j,1,k) 
!              uv(i,j,2,k) = elem(ie)%state%uv(i,j,2,k) 
!               w3d(i,j,k) = elem(ie)%state%w3d(i,j,k) 
!               the(i,j,k) = elem(ie)%state%pt(i,j,k) 
!          end do 
!          end do 
!        end do 
!
!    A_inv(:,:,:,:) = elem(ie)%ainv(:,:,:,:)
!     Amat(:,:,:,:) = elem(ie)%amat(:,:,:,:)
!       sg(:,:)   = elem(ie)%sg(:,:)
!
  ! Converting to contravariant velocity for the vertical Euer solver 
!    do k =1, nlev 
!     contuv(:,:,:,k) = sphere2contra(uv(:,:,:,k),A_inv)
!    enddo 
!
! !! Constructing "U" in the ODE dU/dt = RHS,  U=[ro',u1ro,u2ro,wro,thp']*sg 
!   Call Construct_SateVec_NHvars(ie,userdata,neqn,contuv,rho3d,w3d,the,prs,rk0_vec)

!  !!RK_stage-1
!   Call Compute_EulerSys_rhs(ie,userdata,neqn,rk0_vec,contuv,rho3d,w3d,the,prs,rhs)
!
!   do k=1,nlev
!    do eq=1,neqn 
!      rk1_vec(:,:,eq,k) = rk0_vec(:,:,eq,k) + 1.0*dtv * rhs(:,:,eq,k)
!    enddo
!   enddo

!  !!RK_stage-2
!   Call Recover_Primary_NHvars(ie,userdata,neqn,rk1_vec,contuv,rho3d,w3d,the,prs)  
!   Call Compute_EulerSys_rhs(ie,userdata,neqn,rk1_vec,contuv,rho3d,w3d,the,prs,rhs)
!
!   do eq=1,neqn 
!    do k=1,nlev
!      rk1_vec(:,:,eq,k) = (3.0D0*rk0_vec(:,:,eq,k) + rk1_vec(:,:,eq,k) + dtv * rhs(:,:,eq,k)) *0.25D0
!    enddo
!   enddo
!
!  !!RK_stage-3
!   Call Recover_Primary_NHvars(ie,userdata,neqn,rk1_vec,contuv,rho3d,w3d,the,prs)  
!   Call Compute_EulerSys_rhs(ie,userdata,neqn,rk1_vec,contuv,rho3d,w3d,the,prs,rhs)
!
!   do eq=1,neqn 
!    do k=1,nlev
!      rk1_vec(:,:,eq,k) = (rk0_vec(:,:,eq,k) + 2.0D0*(rk1_vec(:,:,eq,k) + dtv * rhs(:,:,eq,k)) ) / 3.0D0
!    enddo
!   enddo

!   Call Recover_Primary_NHvars(ie,userdata,neqn,rk1_vec,contuv,rho3d,w3d,the,prs)  
!
!!! Updating the changes by the vertical Euler solver on NH primary variables 
!
! ! Converting back to spherical coordinates 
!    do k =1, nlev 
!     uv(:,:,:,k) = contra2sphere(contuv(:,:,:,k),Amat)
!    enddo 
!        do k=1,nlev
!          do j=1,np
!          do i=1,np
!              elem(ie)%state%rho(i,j,k) = rho3d(i,j,k)
!              elem(ie)%state%uv(i,j,1,k) = uv(i,j,1,k)
!              elem(ie)%state%uv(i,j,2,k) = uv(i,j,2,k)
!              elem(ie)%state%w3d(i,j,k) = w3d(i,j,k)
!              elem(ie)%state%pt(i,j,k) = the(i,j,k)
!          end do 
!          end do 
!        end do 
!
! End Subroutine Vertically_ExplicitRK3_EulerSolver

!=======================================================================================================!
  Subroutine Construct_SateVec_NHvars(neqn,sg,u1,u2,rho,w,pt,rhb,prb,rtb,rhp,prp,speed,vstate_vec)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: neqn
    real (kind=real_kind), intent(in)  :: sg
    real (kind=real_kind), dimension(vnp,vnel), intent(in)  :: rho,u1,u2,w,pt,rhb,rtb,prb  
    real (kind=real_kind), dimension(vnp,vnel), intent(out)  :: rhp, prp, speed 
    real (kind=real_kind), dimension(vnp,vnel,neqn), intent(out)  :: vstate_vec

    real (kind=real_kind), dimension(vnp,vnel) :: rhog
    real (kind=real_kind) :: th,ro, pr 
    integer :: i,j,ke 


!!  State vector from primitive varaibles for each vertial(1d) element 

       rhog(:,:) = rho(:,:) *sg
       vstate_vec(:,:,1) = (rho(:,:) - rhb(:,:)) *sg
       vstate_vec(:,:,2) = u1(:,:) * rhog(:,:)
       vstate_vec(:,:,3) = u2(:,:) * rhog(:,:)
       vstate_vec(:,:,4) =  w(:,:) * rhog(:,:)
       vstate_vec(:,:,5) = (rho(:,:) * pt(:,:) - rtb(:,:)) * sg

    do ke=1,vnel 
      do j=1,vnp
        th = pt(j,ke)
        ro = rho(j,ke)
        pr =  Compute_Pressure(ro,th)
        prp(j,ke) = pr - prb(j,ke) 
        rhp(j,ke) = ro - rhb(j,ke) 
        speed(j,ke) = Compute_SoundSpeed(pr,ro)
      enddo
    enddo

  End Subroutine Construct_SateVec_NHvars

!=======================================================================================================!
  Subroutine Recover_Primary_NHvars(neqn,sg,prb,rhb,rtb,rk_dat,u1,u2,rho,w,pt,rhp,prp,speed)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: neqn
    real (kind=real_kind), dimension(vnp,vnel,neqn), intent(in)  :: rk_dat     
    real (kind=real_kind), dimension(vnp,vnel), intent(in)  :: rhb,rtb,prb  
    real (kind=real_kind), dimension(vnp,vnel), intent(out)  :: rho,u1,u2,w,pt
    real (kind=real_kind), dimension(vnp,vnel), intent(out)  :: rhp, prp, speed 

    real (kind=real_kind), dimension(vnp,vnel) :: rhog
    real (kind=real_kind) :: sg,th,ro, pr 
    integer :: i,j,ke


!! Primary varaibles from the state vector at each rk_stage, for  column element 

        rho(:,:) = rk_dat(:,:,1)/sg  + rhb(:,:)
       rhog(:,:) = rho(:,:) * sg
         u1(:,:) = rk_dat(:,:,2) / rhog(:,:)
         u2(:,:) = rk_dat(:,:,3) / rhog(:,:)
          w(:,:) = rk_dat(:,:,4) / rhog(:,:)
         pt(:,:) =(rk_dat(:,:,5) / sg + rtb(:,:) )/rho(:,:)

 !! Computed Pressure, speed of sound as  needed for rhs evaluation 

    do ke=1,vnel 
      do j=1,vnp
        th = pt(j,ke)
        ro = rho(j,ke)
        pr =  Compute_Pressure(ro,th)
        prp(j,ke) = pr - prb(j,ke) 
        rhp(j,ke) = ro - rhb(j,ke) 
        speed(j,ke) = Compute_SoundSpeed(pr,ro)
      enddo
    enddo

  End Subroutine Recover_Primary_NHvars

!=======================================================================================================!
!  Subroutine Recover_Primary_NHvars(ie,userdata,neqn,rk_dat,contuv,rho3d,w3d,the,prs)
!=======================================================================================================!
!    Implicit None
!    integer, intent(in)               :: ie,neqn
!    type (dg_userdata_t) , intent(in) :: userdata
!    type (element_t)     , pointer :: elem(:)
!    real (kind=real_kind), dimension(np,np,neqn,nlev), intent(in)  :: rk_dat      
!    real (kind=real_kind), dimension(np,np,nlev), intent(out)  :: rho3d,w3d,the, prs
!    real (kind=real_kind), dimension(np,np,2,nlev), intent(out)  :: contuv
!
!    real (kind=real_kind), dimension(np,np) :: sg, rhog
!    real (kind=real_kind) :: ro, pt  
!    integer :: i,j,k
!
!       elem  => userdata%elem
!
!       sg(:,:)   = elem(ie)%sg(:,:)
!
!!! Primary varaibles from the state vector at each rk_stage, for an element 
!
!    do k =1, nlev
!       rho3d(:,:,k) = rk_dat(:,:,1,k)/sg(:,:)  + elem(ie)%state%rhb(:,:,k)
!       rhog(:,:) = rho3d(:,:,k) * sg(:,:) 
!       contuv(:,:,1,k) = rk_dat(:,:,2,k) / rhog(:,:)
!       contuv(:,:,2,k) = rk_dat(:,:,3,k) / rhog(:,:)
!            w3d(:,:,k) = rk_dat(:,:,4,k) / rhog(:,:)
!            the(:,:,k) =(rk_dat(:,:,5,k) / sg(:,:) + elem(ie)%state%rtb(:,:,k) )/rho3d(:,:,k) 
!    end do
!
! !! Computed Pressure needed for rhs evaluation 
!
!  do k=1,nlev
!    do j=1,np
!    do i=1,np
!      pt = the(i,j,k)
!      ro = rho3d(i,j,k)
!      prs(i,j,k) =  Compute_Pressure(ro,pt)
!    enddo
!    enddo
!  enddo
!
!  End Subroutine Recover_Primary_NHvars

!=======================================================================================================!
!  Subroutine Construct_SateVec_NHvars(ie,userdata,neqn,contuv,rho3d,w3d,the,prs,vstate_vec)
!=======================================================================================================!
!    Implicit None
!    integer, intent(in)               :: ie,neqn 
!    type (dg_userdata_t) , intent(in) :: userdata
!    type (element_t)     , pointer :: elem(:)
!    real (kind=real_kind), dimension(np,np,nlev), intent(in)  :: rho3d,w3d,the
!    real (kind=real_kind), dimension(np,np,2,nlev), intent(in):: contuv
!    real (kind=real_kind), dimension(np,np,nlev), intent(out)  :: prs      
!    real (kind=real_kind), dimension(np,np,neqn,nlev), intent(out)  :: vstate_vec 
!
!    real (kind=real_kind), dimension(2,2,np,np) :: a_inv
!    real (kind=real_kind), dimension(np,np) :: sg, rhog 
!    real (kind=real_kind) :: pt,ro
!    integer :: i,j,k 
!
!       elem  => userdata%elem
!       sg(:,:) = elem(ie)%sg(:,:)
!

!!  State vector from primitive varaibles for each element 

!    do k =1, nlev 
!       rhog(:,:) = rho3d(:,:,k) *sg(:,:) 
!       vstate_vec(:,:,1,k) = (rho3d(:,:,k) - elem(ie)%state%rhb(:,:,k)) *sg(:,:) 
!       vstate_vec(:,:,2,k) = contuv(:,:,1,k) * rhog(:,:) 
!       vstate_vec(:,:,3,k) = contuv(:,:,2,k) * rhog(:,:) 
!       vstate_vec(:,:,4,k) = w3d(:,:,k) * rhog(:,:)
!       vstate_vec(:,:,5,k) = (rho3d(:,:,k) * the(:,:,k) - elem(ie)%state%rtb(:,:,k)) * sg(:,:) 
!    end do 
!
!    do k=1,nlev
!    do j=1,np
!    do i=1,np
!      pt = the(i,j,k)
!      ro = rho3d(i,j,k)
!      prs(i,j,k) =  Compute_Pressure(ro,pt)
!    enddo
!    enddo
!    enddo
!

!  End Subroutine Construct_SateVec_NHvars

!=======================================================================================================!
  Subroutine Compute_EulerSys_rhs(ie,userdata,neqn,svec,contuv,rho3d,w3d,the,prs,zrhs_sys)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie,neqn 
    type (dg_userdata_t) , intent(in) :: userdata
    type (element_t)     , pointer :: elem(:)
    real (kind=real_kind),dimension(np,np,neqn,nlev), intent(in)  :: svec
    real (kind=real_kind),dimension(np,np,nlev), intent(in)  :: rho3d,w3d,the, prs
    real (kind=real_kind),dimension(np,np,2,nlev), intent(in):: contuv
    real (kind=real_kind), intent(out) :: zrhs_sys(np,np,neqn,nlev)

    real (kind=real_kind), dimension(2,2,np,np) ::  a_inv, g_inv 
    real (kind=real_kind), dimension(2,2) ::   gij      
    real (kind=real_kind), dimension(vnp,vnel,neqn) :: vrhs, vec  
    real (kind=real_kind), dimension(vnp,vnel)   :: u1,u2,w,rho,pt,prp,rhp,rtp,wt 
    real (kind=real_kind), dimension(vnp,vnel)   :: sgv13, sgv23, speed , pr 
    real (kind=real_kind), dimension(nlev)   :: d1d          
    real (kind=real_kind), dimension(neqn,nlev)   :: svec_old, svec_new , spng_src 
    real (kind=real_kind)  :: sg, sgh, sgv, ginv(2,2) 
    integer :: i,j,k,ke , eq 

    ! vgw & vder are vertical Gauss weights and Derivative matrix 
    ! initially activated  in dg3d_test_mod.F90 
       
    elem  => userdata%elem

  !Copy to vertical 1d element coordinate location 

    do j=1,np
    do i=1,np

        sg   = elem(ie)%sg(i,j)
        sgv  = elem(ie)%sgv(i,j) 
        ginv(:,:) = elem(ie)%g_ijinv(:,:,i,j)

        svec_old(:,:) =  elem(ie)%state%vec0(i,j,:,:)
        svec_new(:,:) = svec(i,j,:,:) 

        d1d(:) = elem(ie)%sgv13(i,j,:) 
        sgv13(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) = elem(ie)%sgv23(i,j,:) 
        sgv23(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) =  prs(i,j,:) - elem(ie)%state%prb(i,j,:)
        prp(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) =  rho3d(i,j,:) - elem(ie)%state%rhb(i,j,:)
        rhp(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) =  rho3d(i,j,:) * the(i,j,:) - elem(ie)%state%rtb(i,j,:)
        rtp(:,:) =  copyto_1d_elem(d1d(:))

             u1(:,:) = copyto_1d_elem(contuv(i,j,1,:))
             u2(:,:) = copyto_1d_elem(contuv(i,j,2,:))
              w(:,:) = copyto_1d_elem(w3d(i,j,:))
            rho(:,:) = copyto_1d_elem(rho3d(i,j,:))
             pt(:,:) = copyto_1d_elem(the(i,j,:))
             pr(:,:) = copyto_1d_elem(prs(i,j,:))
   
        !State vector in 1d element structure 
         do eq = 1, neqn 
             vec(:,:,eq) = copyto_1d_elem(svec_new(eq,:))
         enddo 

         do ke =1, vnel 
         do k =1, vnp 
           speed(k,ke) = Compute_SoundSpeed(pr(k,ke),rho(k,ke))
         enddo 
         enddo 

     ! spng_src(:,:) =  Apply_TopSponge_Layer(neqn,svec_old,svec_new)


     ! 1D vertical DG solver and update 

      if ((zgrid_type == "GL").or.(zgrid_type == "GLL")) Then 
          Call Vertical_DGEuler_solver(neqn,sg,sgv,Ginv,sgv13,sgv23,vec,speed,rhp,u1,u2,w,rho,pt,prp,vrhs) 
      endif 
     ! Copy back to vertical 1d array location (nlev) for each eqn 
           zrhs_sys(i,j,:,:) = copyfrom_1dsys_elem(neqn,vrhs(:,:,:))

     ! Adding sponge effect as source term 
     !    zrhs_sys(i,j,:,:) = zrhs_sys(i,j,:,:) + spng_src(:,:)
    end do
    end do

  End Subroutine Compute_EulerSys_rhs

!=======================================================================================================!
  subroutine Vertical_Euler_Forcing(ie,userdata,neqn,uv,rho3d,w3d,the,prs,zrhs_sys)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie,neqn 
    type (dg_userdata_t) , intent(in) :: userdata
    type (element_t)     , pointer :: elem(:)
    real (kind=real_kind),dimension(np,np,nlev), intent(in)  :: rho3d,w3d,the,prs     
    real (kind=real_kind),dimension(np,np,2,nlev), intent(in):: uv
    real (kind=real_kind),dimension(np,np,2,nlev) :: contuv
    real (kind=real_kind),dimension(np,np,nlev) :: wtl       
    real (kind=real_kind), intent(out) :: zrhs_sys(np,np,neqn,nlev)

    real (kind=real_kind), dimension(2,2,np,np) ::  a_inv, g_inv 
    real (kind=real_kind), dimension(vnp,vnel,neqn) :: vrhs , svec 
    real (kind=real_kind), dimension(vnp,vnel)   :: u1,u2,w,rho,pt,prp,rhp,rtp,wt 
    real (kind=real_kind), dimension(vnp,vnel)   :: sgv13, sgv23, speed , pr 
    real (kind=real_kind), dimension(nlev)   :: d1d          
    real (kind=real_kind), dimension(neqn,nlev)   :: svec_old, svec_new , spng_src 
    real (kind=real_kind)  :: sg, sgh, sgv, ginv(2,2)
    integer :: i,j,k,ke,eq   

    ! Note: This for using the vertical effect for the 3D system   as a forcing, 
    ! and only one ODE solver required for this case (2D + 1D discretization) 
       
    elem  => userdata%elem

    a_inv(:,:,:,:) = elem(ie)%ainv(:,:,:,:)

    do k =1, nlev 
     contuv(:,:,:,k) = sphere2contra(uv(:,:,:,k),a_inv)
    enddo 

!   Initial state 

  !Copy to vertical 1d element coordinate location 
    do j=1,np
    do i=1,np

    ginv(:,:) = elem(ie)%g_ijinv(:,:,i,j)

        sg   = elem(ie)%sg(i,j)
        sgv  = elem(ie)%sgv(i,j) 
      svec_old(:,:) =  elem(ie)%state%vec0(i,j,:,:)

        d1d(:) = elem(ie)%sgv13(i,j,:) 
        sgv13(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) = elem(ie)%sgv23(i,j,:) 
        sgv23(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) =  prs(i,j,:) - elem(ie)%state%prb(i,j,:)
        prp(:,:) =  copyto_1d_elem(d1d(:))

        d1d(:) =  rho3d(i,j,:) - elem(ie)%state%rhb(i,j,:)
        rhp(:,:) =  copyto_1d_elem(d1d(:))
          svec_new(1,:) = d1d(:) *sg 
           d1d(:)  = rho3d(i,j,:)  *sg 
          svec_new(2,:) = contuv(i,j,1,:) * d1d(:) 
          svec_new(3,:) = contuv(i,j,2,:) * d1d(:) 
          svec_new(4,:) = w3d(i,j,:) * d1d(:) 

        d1d(:) =  rho3d(i,j,:) * the(i,j,:) - elem(ie)%state%rtb(i,j,:)
        rtp(:,:) =  copyto_1d_elem(d1d(:))
          svec_new(5,:) = d1d(:)  * sg 

             u1(:,:) = copyto_1d_elem(contuv(i,j,1,:))
             u2(:,:) = copyto_1d_elem(contuv(i,j,2,:))
              w(:,:) = copyto_1d_elem(w3d(i,j,:))
            rho(:,:) = copyto_1d_elem(rho3d(i,j,:))
             pt(:,:) = copyto_1d_elem(the(i,j,:))
             pr(:,:) = copyto_1d_elem(prs(i,j,:))
   
        ! state vector in 1d element structure 
         do eq = 1, neqn 
             svec(:,:,eq) = copyto_1d_elem(svec_new(eq,:))
         enddo

         do ke =1, vnel 
         do k =1, vnp 
           speed(k,ke) = Compute_SoundSpeed(pr(k,ke),rho(k,ke))
         enddo 
         enddo 
       spng_src(:,:) =  0.0D0  

     ! spng_src(:,:) =  Apply_TopSponge_Layer(neqn,svec_old,svec_new)

     ! 1D vertical DG solver and update 

      if ((zgrid_type == "GL").or.(zgrid_type == "GLL")) Then 
          Call vertical_DGEuler_solver(neqn,sg,sgv,ginv,sgv13,sgv23,svec,speed,rhp,u1,u2,w,rho,pt,prp,vrhs) 
      endif 
     ! Copy back to vertical 1d array location (nlev) for each eqn 
          zrhs_sys(i,j,:,:) = copyfrom_1dsys_elem(neqn,vrhs(:,:,:))

     ! Adding sponge effect as source term 
          zrhs_sys(i,j,:,:) = zrhs_sys(i,j,:,:) + spng_src(:,:)
    end do
    end do

  End subroutine Vertical_Euler_Forcing
!=======================================================================================================!
   Function Apply_TopSponge_Layer(neqn,vec_old,vec_new) result(spng_src)
!=======================================================================================================!
    Implicit None
    integer, intent(in)  :: neqn
    real (kind=real_kind), dimension(neqn,nlev), intent(in) ::  vec_old, vec_new 
    real (kind=real_kind) :: spng_src(neqn,nlev)
    real (kind=real_kind) :: tow           
    integer :: i,k,eq,ke

!     tow = 1.0D0 / 25.0D0 
      tow = 0.20 
    spng_src(:,:) = 0.0D0 

    do k=1,nlev 

     if (z_grd(k) > zsponge_level) then 
       do eq=1,neqn -1
         spng_src(eq,k) = -tow * zsponge_coef(k) * (vec_new(eq,k) - vec_old(eq,k)) 
       enddo
     endif  

    enddo

   End Function Apply_TopSponge_Layer
!=======================================================================================================!
  function Zsponge_function_top(z_up) result(spg_fn)     
!=======================================================================================================!
    Implicit None
    real (kind=real_kind), intent(in) :: z_up           
    real (kind=real_kind)             :: spg_fn(nlev)
    real (kind=real_kind)  :: zk
    integer :: k 

! dcmip specification 
    do  k=1,nlev  
      if (z_grd(k) > z_up) then 
       zk = ((z_grd(k) - z_up)/(z_top - z_up) )
       spg_fn(k) = (sin(DD_PI*0.5D0 * zk))**4 
      else 
       spg_fn(k) = 0.0D0 
      endif 
    end do 
 end function Zsponge_function_top

!=======================================================================================================!
   Subroutine  Vertical_DGEuler_Solver(neqn,sg,sgv,ginv,sgv13,sgv23,vec,speed,rhp,u1,u2,w,rho,pt,prp,vrhs)
!=======================================================================================================!
    Implicit None
    integer, intent(in) :: neqn       
    real (kind=real_kind), intent(in) :: sg, sgv, ginv(2,2)  
    real (kind=real_kind), dimension(vnp,vnel), intent(in) :: sgv13,sgv23, speed 
    real (kind=real_kind), dimension(vnp,vnel,neqn), intent(in) :: vec
    real (kind=real_kind), dimension(vnp,vnel), intent(inout) :: u1,u2,w
    real (kind=real_kind), dimension(vnp,vnel), intent(inout) :: rho,pt,prp,rhp
    real (kind=real_kind), dimension(vnp,vnel,neqn), intent(out) ::  vrhs 

    real (kind=real_kind), dimension(vnp,vnel)::  wt
    real (kind=real_kind), dimension(2,vnel) :: sgv13_int, sgv23_int, speed_int,fj_max 
    real (kind=real_kind), dimension(2,vnel) :: u1_int,u2_int,w_int,pt_int,rho_int,wt_int,prp_int
    real (kind=real_kind), dimension(2,vnel) :: u1_ext,u2_ext,w_ext,pt_ext,rho_ext,wt_ext,prp_ext

    real (kind=real_kind), dimension(vnp,vnel,neqn) :: flz,src, grad_term, num_flux
    real (kind=real_kind), dimension(2,0:vnel+1,neqn) :: vec_int, vec_ext, flz_int, flz_ext  
    real (kind=real_kind), dimension(2,0:vnel+1) ::  fj_int
    real (kind=real_kind), dimension(vnp) :: vmm        

    real (kind=real_kind) ::  sgh,sgrho,sgrho2(2), sgrhowt,sghprp,  sgrhowt2(2), sghpp(2)  
    real (kind=real_kind) ::  g13,g23 , g13_int(2), g23_int(2) , am
    integer :: i,j,k , ke, eq 

 !vertical inverted mass-matrix 
   do k=1,vnp
      vmm(k) =  2.0D0/ (z_wel * vgw(k) )
   enddo 

     sgh = sg / sgv  

     am = 1.0D0

 !The state vector form primary varibales 
     do ke=1, vnel
     do  k=1, vnp  
      !sgrho = sg * rho(k,ke) 
      !vec(k,ke,1) = sg * rhp(k,ke) 
      !vec(k,ke,2) = sgrho *u1(k,ke)
      !vec(k,ke,3) = sgrho *u2(k,ke)
      !vec(k,ke,4) = sgrho * w(k,ke)
      !vec(k,ke,5) = sg * rtp(k,ke)  
        wt(k,ke) =  (w(k,ke) + am*sgv13(k,ke) * u1(k,ke) +  am*sgv23(k,ke) * u2(k,ke)) /sgv 
     enddo 
     enddo 

    ! wt(1,1)  =  0.0D0 
 !vertical Euler fluxes 
     do ke=1, vnel
     do  k=1, vnp  
           sgrhowt = sg * rho(k,ke) * wt(k,ke) 
           sghprp  = sgh * prp(k,ke)
       flz(k,ke,1) = sgrhowt 
       flz(k,ke,2) = sgrhowt *u1(k,ke) + sghprp*(sgv13(k,ke)*ginv(1,1) + sgv23(k,ke)*ginv(1,2))
       flz(k,ke,3) = sgrhowt *u2(k,ke) + sghprp*(sgv13(k,ke)*ginv(2,1) + sgv23(k,ke)*ginv(2,2)) 
       flz(k,ke,4) = sgrhowt *w(k,ke)  + sghprp 
       flz(k,ke,5) = sgrhowt *pt(k,ke)

     enddo 
     enddo 

!     maybe sourced with horizontal dynamics 
      src(:,:,:) = 0.0D0 
      src(:,:,4) =  -rhp(:,:) * g *sg 
  !!  src(1,1,4) =  -rhp(1,1) * g *sg *0.2 

   !Extrapolate edge GL values 
     do ke=1, vnel
      sgv13_int(:,ke) = edge_GL_1d(sgv13(:,ke))
      sgv23_int(:,ke) = edge_GL_1d(sgv23(:,ke))
      speed_int(:,ke) = edge_GL_1d(speed(:,ke)) 

      rho_int(:,ke) = edge_GL_1d(rho(:,ke))
       u1_int(:,ke) = edge_GL_1d(u1(:,ke))
       u2_int(:,ke) = edge_GL_1d(u2(:,ke))
        w_int(:,ke) = edge_GL_1d(w(:,ke))
       pt_int(:,ke) = edge_GL_1d(pt(:,ke))
      prp_int(:,ke) = edge_GL_1d(prp(:,ke))
     end do


     do ke=1, vnel
      wt_int(:,ke) = (w_int(:,ke) + am*sgv13_int(:,ke)*u1_int(:,ke) + am*sgv23_int(:,ke)*u2_int(:,ke)) /sgv 
     end do

      wt_int(ibot,1)  =  0.0D0 
      wt_int(itop,vnel)  = 0.0D0 

    ! Construct edge state vector 

    !u1_int(ibot,1)  =  0.0D0 
    !u2_int(ibot,1)  =  0.0D0 
    !!w_int(ibot,1)  =  0.0D0 

    !sgv13_int(ibot,1)  = sgv13(ibot,1) 
    !sgv23_int(ibot,1)  = sgv23(ibot,1) 

    !if (sgv13_int(ibot,1) /= 0.0D0 )  u1_int(ibot,1)  = -w_int(ibot,1) / sgv13_int(ibot,1) 
    !if (sgv23_int(ibot,1) /= 0.0D0 )  u2_int(ibot,1)  = -w_int(ibot,1) / sgv23_int(ibot,1) 
      

    do eq = 1, neqn
     do ke=1, vnel
        vec_int(:,ke,eq) = edge_GL_1d(vec(:,ke,eq))
     end do 
    end do 

   ! do ke=1, vnel
   !  sgrho2(:) = sg * rho_int(:,ke) 
   !  vec_int(:,ke,1) = edge_GL_1d(vec(:,ke,1))
   !  vec_int(:,ke,2) = sgrho2(:) *u1_int(:,ke)
   !  vec_int(:,ke,3) = sgrho2(:) *u2_int(:,ke)
   !  vec_int(:,ke,4) = sgrho2(:) *w_int(:,ke)
   !  vec_int(:,ke,5) = edge_GL_1d(vec(:,ke,5))
   ! enddo 

    do ke=1,vnel
     fj_int(:,ke) = speed_int(:,ke) /sgv  + abs(wt_int(:,ke)) 
    !fj_int(:,ke) = speed_int(:,ke)   + abs(wt_int(:,ke)) 
    end do 

    !fj_int(ibot,1) = 0.0
     fj_int(itop,0) = fj_int(ibot,1) 
     fj_int(ibot,vnel+1) = fj_int(itop,vnel) 

    do ke=1,vnel
     fj_max(ibot,ke) = max(fj_int(ibot,ke), fj_int(itop,ke-1)) 
     fj_max(itop,ke) = max(fj_int(itop,ke), fj_int(ibot,ke+1)) 
    end do 

    !For the lowest boundary (wt_int ==> 0)  
    ! w_int(ibot,1) =  -(u1_int(ibot,1)*sgv13_int(ibot,1) + u2_int(ibot,1)*sgv23_int(ibot,1))
                            
    ! Construct edge fluxes 
     do ke=1, vnel
           sgrhowt2(:) = sg * rho_int(:,ke) * wt_int(:,ke)
              sghpp(:) = sgh *prp_int(:,ke)
       flz_int(:,ke,1) = sgrhowt2(:) 
       flz_int(:,ke,2) = sgrhowt2(:) * u1_int(:,ke) + sghpp(:)*(ginv(1,1)*sgv13_int(:,ke) + &
                                                                ginv(1,2)*sgv23_int(:,ke)) 
       flz_int(:,ke,3) = sgrhowt2(:) * u2_int(:,ke) + sghpp(:)*(ginv(2,1)*sgv13_int(:,ke) + &
                                                                ginv(2,2)*sgv23_int(:,ke)) 
       flz_int(:,ke,4) = sgrhowt2(:) * w_int(:,ke)  + sghpp(:)
       flz_int(:,ke,5) = sgrhowt2(:) * pt_int(:,ke)

     enddo

    ! vec_int(ibot,1,4) =  sg*w_int(ibot,1)*rho_int(ibot,1)
    ! flz_int(ibot,1,4) =  sgh* prp_int(ibot,1) !!- vec_int(ibot,1,1)*g *z_wel 

    Call  TopBot_noflux_BC(sgv13_int,sgv23_int,vec_int,flz_int) 

    do ke=1,vnel
      vec_ext(ibot,ke,:) = vec_int(itop,ke-1,:)
      vec_ext(itop,ke,:) = vec_int(ibot,ke+1,:)

      flz_ext(ibot,ke,:) = flz_int(itop,ke-1,:)
      flz_ext(itop,ke,:) = flz_int(ibot,ke+1,:)
    enddo


     num_flux(:,:,:) =  Vertical_LxF_Flux(fj_max,vec_int,vec_ext,flz_int,flz_ext)

     grad_term(:,:,:) = Compute_Vgrad_term(flz)  

    ! Vertical DG RHS 

     do eq=1,neqn    
      do ke=1,vnel
       do k=1,vnp
        vrhs(k,ke,eq) = (grad_term(k,ke,eq)  - num_flux(k,ke,eq)) * vmm(k) + src(k,ke,eq)  
       !vrhs(k,ke,eq) = 0.0  
       enddo
      enddo
     enddo

   End Subroutine  Vertical_DGEuler_Solver

!=======================================================================================================!
   Subroutine  TopBot_noflux_BC(sg13_int,sg23_int,vec_int,flz_int)
!=======================================================================================================!
    Implicit None 
    real (kind=real_kind), dimension(2,vnel), intent(in)  :: sg13_int, sg23_int       
    real (kind=real_kind), dimension(2,0:vnel+1,5), intent(inout)  :: vec_int,  flz_int
    real (kind=real_kind) :: num_flux(np,vnel,5)
    real (kind=real_kind) :: vec2, vec3, vec_up,vec_dn , lxf_ibot, lxf_itop , gu1,gu2 
    integer :: k,eq,ke


        flz_int(:,0,:) = 0.0D0
        flz_int(itop,0,:) =  flz_int(ibot,1,:)  
!       flz_int(itop,0,4) =  flz_int(ibot,1,4) 


        vec_int(:,0,:) = 0.0D0
        vec_int(itop,0,:) = vec_int(ibot,1,:)
        vec_int(itop,0,4) = -vec_int(ibot,1,4) 

        gu1 =  sg13_int(ibot,1)
        gu2 =  sg23_int(ibot,1)

!   If ((gu1 /= 0.0D0).or.(gu2 /= 0.0D0)) then 
!     vec2 = vec_int(ibot,1,2) 
!     vec3 = vec_int(ibot,1,3) 
!     vec_int(ibot,1,4) = -(vec2*gu1 + vec3*gu2) !+ vec_int(ibot,1,1)*g 
!     vec_int(itop,0,4) =  vec_int(ibot,1,4)  
!   endif 

        flz_int(:,vnel+1,:) = 0.0D0
        flz_int(ibot,vnel+1,4) =  flz_int(itop,vnel,4)

        vec_int(:,vnel+1,:) = 0.0D0
        vec_int(ibot,vnel+1,:) = vec_int(itop,vnel,:)
        vec_int(ibot,vnel+1,4) = -vec_int(itop,vnel,4)

   End Subroutine  TopBot_noflux_BC

!=======================================================================================================!
   Function Vertical_LXF_Flux(fj_max,vec_int,vec_ext,flz_int,flz_ext) result(num_flux) 
!=======================================================================================================!
    Implicit None 
    real (kind=real_kind), dimension(2,0:vnel+1,5) :: vec_int, vec_ext, flz_int, flz_ext
    real (kind=real_kind), intent(in) :: fj_max(2,vnel)
    real (kind=real_kind) :: num_flux(np,vnel,5)
    real (kind=real_kind) :: flz_up,flz_dn, vec_up,vec_dn , lxf_ibot, lxf_itop 
    integer :: k,eq,ke

      do eq=1, 5          
       do ke=1, vnel
           flz_up = flz_int(ibot,ke,eq)
           flz_dn = flz_ext(ibot,ke,eq)
           vec_up = vec_int(ibot,ke,eq)
           vec_dn = vec_ext(ibot,ke,eq)
           lxf_ibot = 0.5D0 * ((flz_up + flz_dn) - fj_max(ibot,ke) *(vec_up - vec_dn))

           flz_up = flz_ext(itop,ke,eq)
           flz_dn = flz_int(itop,ke,eq)
           vec_up = vec_ext(itop,ke,eq)
           vec_dn = vec_int(itop,ke,eq)
           lxf_itop = 0.5D0 * ((flz_up + flz_dn) - fj_max(itop,ke) *(vec_up - vec_dn))

         do k = 1, vnp 
           num_flux(k,ke,eq) =  hed_gl(k,2) * lxf_itop - hed_gl(k,1) * lxf_ibot 
         end do 
       end do
      end do

   End Function Vertical_LXF_Flux

!=======================================================================================================!
   Function Compute_Vgrad_term(flz) result(zgrad)
!=======================================================================================================!
    Implicit None
    real (kind=real_kind), dimension(vnp,vnel,5), intent(in) ::  flz
    real (kind=real_kind) :: zgrad(vnp,vnel,5)
    real (kind=real_kind) :: sm,flux(vnp)
    integer :: i,k,eq,ke

    do eq=1,5 
    do ke=1,vnel

     flux(:) = flz(:,ke,eq)

    !Gradient weak (DG) 
      do k=1,vnp
        sm = 0.0D0
       do i=1,vnp
         sm = sm + flux(i) *vder(i,k)  * vgw(i)
       enddo
        zgrad(k,ke,eq) = sm
      enddo

    enddo
    enddo

   End Function Compute_Vgrad_term

!=======================================================================================================!
  subroutine vertical_advect_dgsrc(ie,userdata,w3d,rho3d,srcdg)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie
    type (dg_userdata_t) , intent(in) :: userdata
    real (kind=real_kind), intent(in)    :: w3d(np,np,nlev)
    real (kind=real_kind), intent(in) :: rho3d(np,np,nlev)
    real (kind=real_kind), intent(out) :: srcdg(np,np,nlev)

    real (kind=real_kind), dimension(vnp,vnel)   :: w, rho, src 
    real (kind=real_kind), dimension(nlev)   :: d1d          
    integer :: i,j,k

    ! vgw & vder are vertical Gauss weights and Derivative matrix 
    ! generated in dg3d_test_mod.F90 

       
    do j=1,np
    do i=1,np

     ! Copy to vertical 1d element coordinate location 

              w(:,:) =  copyto_1d_elem(w3d(i,j,:))
            rho(:,:) =  copyto_1d_elem(rho3d(i,j,:))

     ! 1D vertical DG solver and update 
     !if (zgrid_type == "GL") then

      if ((zgrid_type == "GL").or.(zgrid_type == "GLL")) Then 
           src(:,:) = DG1d_compute_rhs(w,rho) 
         ! src(:,:) = DG1d_GLL_rhs(vgw,vder,w,rho)   !testing 
      endif 


     ! Copy back to vertical 1d array location (nlev) 
           srcdg(i,j,:) = copyfrom_1d_elem(src)

    end do
    end do

  end subroutine vertical_advect_dgsrc

!=======================================================================================================!
   subroutine vertical_solve_update(ie,userdata,dtime,w3d,rho3d)
!=======================================================================================================!
!  Used for the time-Split version , For advection problem
    Implicit None 
    integer, intent(in)               :: ie
    type (dg_userdata_t) , intent(in) :: userdata
    real (kind=real_kind), intent(in)    :: dtime
    real (kind=real_kind), intent(in)    :: w3d(np,np,nlev)
    real (kind=real_kind), intent(inout) :: rho3d(np,np,nlev)

    real (kind=real_kind), dimension(vnp,vnel)   :: w, rho          
    real (kind=real_kind), dimension(nlev)   :: d1d, vrho, vw  
    integer :: i,j,k 

   !!Note: Grid choice FD/GL shoulbe be made at dg3d_tests_mod 

    if (zgrid_type == "FD") Then 
      do j=1,np
      do i=1,np

        vw(:) = w3d(i,j,:) 
        vrho(:) = rho3d(i,j,:) 

        Call FD1d_rk3_solver(dtime,vw,vrho) 

        rho3d(i,j,:) = vrho(:)  
      end do
      end do
    endif 

    if ((zgrid_type == "GL").or.(zgrid_type == "GLL")) Then 
      do j=1,np
      do i=1,np
     ! Copy to vertical 1d element coordinate location 
              w(:,:) =  copyto_1d_elem(w3d(i,j,:))
            rho(:,:) =  copyto_1d_elem(rho3d(i,j,:))

     ! 1D vertical DG solver and update 
          Call dg1d_rk3_solver(dtime,w,rho) 
         
     ! Copy back to vertical 1d array location (nlev) 
          rho3d(i,j,:) = copyfrom_1d_elem(rho) 

      end do
      end do
    endif 

  end subroutine vertical_solve_update 

!=======================================================================================================!
   function fd1d_rhs_src(win,rho) result(src) 
!=======================================================================================================!
    Implicit None 

    real (kind=real_kind), intent(in) :: win(nlev)
    real (kind=real_kind), intent(in) :: rho(nlev)

    real (kind=real_kind), dimension(vnp,vnel) :: w,ro 
    real (kind=real_kind), dimension(nlev) :: src 
    real (kind=real_kind), dimension(0:nlev+1) :: flz
    real (kind=real_kind) :: s1,s2,alpha, sm, dz, dz2, dz12 ,wk
    integer :: i,j,k
     
  !! An FD approximation of vertical 1d advection. 
  !! The grid line (z_grd) is assumed to be incresing from bottom to top.

       flz = 0.0D0 
       src = 0.0D0 

     do  k = 1, nlev 
       flz(k) = rho(k) *win(k) 
     enddo 

   ! boundary conditions 
       flz(1) = 0.0D0 
       flz(nlev) = 0.0D0 
        
       dz2 = 2.0D0 *  z_wel 
       dz12 = 12.0D0 *  z_wel 

   do  k = 1, nlev 
     if (k <= 2) then    !bottom edge 
      src(k) = (flz(k+3) - 6.0D0*flz(k+2) + 18.0D0*flz(k+1) -10.0D0*flz(k) - 3.0D0* flz(k-1)) / dz12  !4th-order

     elseif ((k > 2) .and. (k < nlev-2)) then
      src(k) = (flz(k-2) + 8.0D0*(-flz(k-1) + flz(k+1)) - flz(k+2)) / dz12  !4th-order

     elseif (k >= nlev-2) then     !top edge 
      src(k) = (3.0D0*flz(k+1) + 10.0D0*flz(k) - 18.0D0*flz(k-1) + 6.0D0*flz(k-2) - flz(k-3)) / dz12  !4th-order

     else     ! default/dumb 
      src(k) = (flz(k+1) - flz(k-1)) / dz2      !2nd-order
     endif 
    enddo 

      src(1) = 0.0D0 
      src(nlev) = 0.0D0 

    !! dU/dt = - src = RHS 
     src(:) = - src(:)   ! sign flip is necessary to keep this term at the ODE-RHS


   End function fd1d_rhs_src

!=======================================================================================================!
   function DG1d_compute_rhs(win,rhoin) result(rhs) 
!=======================================================================================================!
    Implicit None 

   !Integer, parameter :: ibot=1, itop=2
    real (kind=real_kind), dimension(vnp,vnel), intent(in) :: win
    real (kind=real_kind), dimension(vnp,vnel), intent(inout) ::  rhoin

    real (kind=real_kind), dimension(vnp,vnel) :: rhs,w , rho, flz
    real (kind=real_kind), dimension(2,0:vnel+1) ::  rhoe_ext,we_ext,rhoe_int,we_int 
    real (kind=real_kind), dimension(2,0:vnel+1) ::  flz_ext, flz_int 
    real (kind=real_kind), dimension(2,vnel) :: lxf, jumpf
    real (kind=real_kind), dimension(vnp) :: flux, grad, cflx
    real (kind=real_kind) :: s1,s2,alpha, sm, dl , dz2 
    real (kind=real_kind) :: rho_up,rho_dn,f_up,f_dn 
    integer :: i,j,k , ke 

  !! Note:  ibot ==> left edge & itop ==> right edge of 1d vertical  element in a FV sense, 
  !! when counting from the element_index = 1 (bottom) to vnel (top) 

  rhoin(1,1) = 0.0D0 
  rhoin(vnp,vnel) = 0.0D0 

    do ke=1,vnel 
    do k=1,vnp  
     w(k,ke) = win(k,ke)
     rho(k,ke) = rhoin(k,ke)
    flz(k,ke) = w(k,ke)*rho(k,ke) 
    enddo 
    enddo 

     we_int(:,:) = 0.0D0 
     we_ext(:,:) = we_int(:,:)     
    rhoe_int(:,:) = 0.0D0 

     dz2 = 2.0D0 / z_wel 

  w(1,1) = 0.0D0 
  w(vnp,vnel) = 0.0D0 
!  rho(1,1) = 0.0D0 
!  rho(vnp,vnel) = 0.0D0 

   if (zgrid_type =="GL") then 
    do ke=1,vnel 
     rhoe_int(:,ke) = edge_gl_1d(rhoin(:,ke))
     we_int(:,ke) = edge_gl_1d(win(:,ke))
    enddo 

   elseif (zgrid_type == "GLL") then
    do ke=1,vnel
     we_int(ibot,ke) = w(1,ke)
     we_int(itop,ke) = w(vnp,ke)
     rhoe_int(ibot,ke) = rho(1,ke)
     rhoe_int(itop,ke) = rho(vnp,ke)
   enddo
  endif 

    !no-flux boundary condition 
    we_int(ibot,1) = 0.0D0 
    we_int(itop,vnel) = 0.0D0 
!   rhoe_int(ibot,1) = 0.0D0 
!   rhoe_int(itop,vnel) = 0.0D0 

     we_int(itop,0) = - we_int(ibot,1) 
     we_int(ibot,vnel+1) = - we_int(itop,vnel) 
    rhoe_int(:,0) = 0.0D0 
    rhoe_int(:,vnel+1) = 0.0D0 

    do ke=1,vnel 
     flz_int(itop,ke) = we_int(itop,ke)*rhoe_int(itop,ke) 
     flz_int(ibot,ke) = we_int(ibot,ke)*rhoe_int(ibot,ke) 
    enddo 

  ! flz_int(itop,vnel) = 0.0D0       
  ! flz_int(ibot,1) = 0.0D0       

     flz_int(itop,0) = - flz_int(ibot,1) 
     flz_int(ibot,vnel+1) = - flz_int(itop,vnel) 

    do ke=1,vnel 
       we_ext(ibot,ke) =  we_int(itop,ke-1) 
       we_ext(itop,ke) =  we_int(ibot,ke+1) 
     rhoe_ext(ibot,ke) = rhoe_int(itop,ke-1) 
     rhoe_ext(itop,ke) = rhoe_int(ibot,ke+1) 
      flz_ext(ibot,ke) = flz_int(itop,ke-1) 
      flz_ext(itop,ke) = flz_int(ibot,ke+1) 
    enddo 

    !Lax-Fred fluxes at each edge 

  do ke=1,vnel 
     alpha = max( abs(we_int(ibot,ke)),abs(we_ext(ibot,ke)) )
      f_up = flz_int(ibot,ke) 
      f_dn = flz_ext(ibot,ke) 
      rho_up = rhoe_int(ibot,ke) 
      rho_dn = rhoe_ext(ibot,ke) 
     lxf(ibot,ke) = 0.5D0 * ((f_up + f_dn) - alpha*(rho_up - rho_dn)) 

     jumpf(ibot,ke) = lxf(ibot,ke) - f_up 

     alpha = max( abs(we_int(itop,ke)),abs(we_ext(itop,ke)) )
      f_dn = flz_int(itop,ke) 
      f_up = flz_ext(itop,ke) 
      rho_up = rhoe_ext(itop,ke) 
      rho_dn = rhoe_int(itop,ke) 
     lxf(itop,ke) = 0.5D0 * ((f_up + f_dn) - alpha*(rho_up - rho_dn)) 

     jumpf(itop,ke) = lxf(itop,ke) - f_dn 
  enddo 
 
   ! DG discretization 

  do ke=1,vnel 
     flux(:) = flz(:,ke) 

    do k=1,vnp   
      sm = 0.0D0 
    !Gradient weak (DG) 
      do i=1,vnp 
       sm = sm + flux(i) *vder(i,k)  * vgw(i) 
      enddo 
    !Gradient strong (SD) 
   !    do i=1,vnp 
   !      sm = sm + flux(i) * vder(k,i) 
   !    enddo 
      grad(k) = sm 
    enddo 

   !Flux correction for SD 
   ! do k=1,vnp   
   !   cflx(k) = rad_der(k,1)*jumpf(itop,ke) + rad_der(k,2)*jumpf(ibot,ke) 
   ! enddo 

   if (zgrid_type == "GL") then 
    !rhs evaluation of the ODE 
    do k=1,vnp   
     rhs(k,ke) =  (grad(k) - (hed_gl(k,2)*lxf(itop,ke) - hed_gl(k,1)*lxf(ibot,ke)) ) * dz2 /vgw(k)  !DG 
   ! rhs(k,ke) = -(grad(k) + cflx(k)) * dz2    ! For SD 
    enddo 

   else if (zgrid_type == "GLL") then 
    do k=1,vnp
      s1 = 0.0D0
      s2 = 0.0D0
     if (k == 1) s1 = 1.0D0
     if (k == vnp) s2 = 1.0D0
     rhs(k,ke) = (grad(k) - (s2* lxf(itop,ke) - s1*lxf(ibot,ke)) ) * dz2/vgw(k)  
    enddo
   endif 

  enddo 

 rhs(1,1) = 0.0D0 
 rhs(vnp,vnel) = 0.0D0 

 End function DG1d_compute_rhs 

!=======================================================================================================!
   subroutine vertical_advect_fdsrc(ie,userdata,w3d,rho3d,srcfd)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie
    type (dg_userdata_t) , intent(in) :: userdata
    real (kind=real_kind), intent(in) :: w3d(np,np,nlev)
    real (kind=real_kind), intent(in) :: rho3d(np,np,nlev)
    real (kind=real_kind), intent(out) :: srcfd(np,np,nlev)

    real (kind=real_kind), dimension(nlev)   :: w, rho, src
    integer :: i,j,k

    ! vgw & vder are vertical Gauss weights and Derivative matrix 
    ! generated in dg3d_test_mod.F90 

    do j=1,np
    do i=1,np

     ! Copy to vertical 1d coordinate location 
             w(:) =  w3d(i,j,:)
           rho(:) =  rho3d(i,j,:)

     ! 1D vertical DG solver and update 
           src(:) =  fd1d_rhs_src(w,rho)

     ! Copy back to vertical 1d array location (nlev) 
           srcfd(i,j,:) = src(:)

    end do
    end do

  end subroutine vertical_advect_fdsrc

!=======================================================================================================!
   function DG1d_GLL_rhs(gw,der,win,rho) result(rhs) 
!=======================================================================================================!
    Implicit None 

   !Integer, parameter :: ibot=1,itop=2
    real (kind=real_kind), intent(in) :: der(vnp,vnp), gw(vnp) 
    real (kind=real_kind), intent(in) :: win(vnp,vnel)
    real (kind=real_kind), intent(in) :: rho(vnp,vnel)

    real (kind=real_kind), dimension(vnp,vnel) :: rhs,w 
    real (kind=real_kind), dimension(2,0:vnel+1) ::  rhoe_ext,we_ext,rhoe_int,we_int 
    real (kind=real_kind), dimension(2,0:vnel+1) ::  flz_ext, flz_int 
    real (kind=real_kind), dimension(2,vnel) :: lxf
    real (kind=real_kind) :: flux(vnp), grad(vnp)
    real (kind=real_kind) :: s1,s2,alpha, sm, dl 
    real (kind=real_kind) :: rho_up,rho_dn,f_up,f_dn 
    integer :: i,j,k , ke 

     w(:,:) = win(:,:) 
     we_int(:,:) = 0.0D0 
     we_ext(:,:) = we_int(:,:)     

    w(1,1) = 0.0D0 
    w(vnp,vnel) = 0.0D0 

    do ke=1,vnel 
     we_int(itop,ke) = w(1,ke) 
     we_int(ibot,ke) = w(vnp,ke) 
     rhoe_int(itop,ke) = rho(1,ke) 
     rhoe_int(ibot,ke) = rho(vnp,ke) 
    enddo 

    !no-flux boundary condition 
     we_int(itop,1) = 0.0D0 
     we_int(ibot,vnel) = 0.0D0 


    do ke=1,vnel 
     flz_int(itop,ke) = we_int(itop,ke)*rhoe_int(itop,ke) 
     flz_int(ibot,ke) = we_int(ibot,ke)*rhoe_int(ibot,ke) 
    enddo 

    rhoe_int(:,0) = 0.0D0 
    rhoe_int(:,vnel+1) = 0.0D0 
    flz_int(:,0) = 0.0D0 
    flz_int(:,vnel+1) = 0.0D0 
    

    do ke=1,vnel
       we_ext(ibot,ke) =  we_int(itop,ke-1)
       we_ext(itop,ke) =  we_int(ibot,ke+1)
     rhoe_ext(ibot,ke) = rhoe_int(itop,ke-1)
     rhoe_ext(itop,ke) = rhoe_int(ibot,ke+1)
      flz_ext(ibot,ke) = flz_int(itop,ke-1)
      flz_ext(itop,ke) = flz_int(ibot,ke+1)
    enddo


    !Lax-Fred fluxes at each edge 
  do ke=1,vnel 

     alpha = max(we_int(ibot,ke),we_ext(ibot,ke))
      f_up = flz_int(ibot,ke)
      f_dn = flz_ext(ibot,ke)
      rho_up = rhoe_int(ibot,ke)
      rho_dn = rhoe_ext(ibot,ke)
      lxf(ibot,ke) = 0.5D0 * ((f_up + f_dn) - alpha*(rho_up - rho_dn))

      alpha = max(we_int(itop,ke),we_ext(itop,ke))
      f_up = flz_ext(itop,ke)
      f_dn = flz_int(itop,ke)
      rho_up = rhoe_ext(itop,ke)
      rho_dn = rhoe_int(itop,ke)
      lxf(itop,ke) = 0.5D0 * ((f_up + f_dn) - alpha*(rho_up - rho_dn))


  enddo 
 
   ! DG discretization 

  do ke=1,vnel 
     !Gradient 
      do k=1,vnp   
        flux(k) = w(k,ke)*rho(k,ke) 
      enddo 
    do k=1,np   
      sm = 0.0D0 
       do i=1,vnp 
        sm = sm + flux(i) *der(i,k) * gw(i) 
       enddo 
      grad(k) = sm 
    enddo 

     !rhs for dg1d 
    do k=1,np   
      s1 = 0.0D0 
      s2 = 0.0D0 
      if (k == 1) s1 = 1.0D0 
      if (k == np) s2 = 1.0D0 

     rhs(k,ke) = (grad(k) - (s2* lxf(itop,ke) - s1*lxf(ibot,ke)) ) * 2.0D0/(gw(k) * z_wel) 
    enddo 

  enddo 
 end function DG1d_GLL_rhs

!=======================================================================================================!
 subroutine dg1d_rk3_solver(dtime,w,rho)
!=======================================================================================================!
    Implicit None 

    real (kind=real_kind), intent(in)    :: dtime
    real (kind=real_kind), intent(in)    :: w(vnp,vnel)
    real (kind=real_kind), intent(inout) :: rho(vnp,vnel)

    real (kind=real_kind), dimension(vnp,vnel) :: rk0, rk1, rhs 
    real (kind=real_kind) :: s1, trm           
    integer :: ke, k 

      rk0(:,:) = rho(:,:)

   rhs(:,:) = 0.0D0 

   rhs(:,:) = DG1d_compute_rhs(w,rk0) 
  

  !RK_stage-1
   do ke=1,vnel
    do k=1,vnp
      rk1(k,ke) = rk0(k,ke) + dtime * rhs(k,ke)
    enddo
   enddo

  rhs(:,:) = DG1d_compute_rhs(w,rk1) 

  !RK_stage-2
   do ke=1,vnel
    do k=1,vnp
      rk1(k,ke) = (3.0D0*rk0(k,ke) + rk1(k,ke) + dtime * rhs(k,ke)) *0.25D0 
    enddo
   enddo

   rhs(:,:) = DG1d_compute_rhs(w,rk1) 

  !RK_stage-3
   do ke=1,vnel
    do k=1,vnp
      rho(k,ke) = (rk0(k,ke) + 2.0D0*(rk1(k,ke) + dtime * rhs(k,ke)) ) / 3.0D0 
    enddo
   enddo
 end subroutine dg1d_rk3_solver

!=======================================================================================================!
 subroutine fd1d_rk3_solver(dtime,w,rho)
!=======================================================================================================!
    Implicit None

    real (kind=real_kind), intent(in)    :: dtime
    real (kind=real_kind), intent(in)    :: w(nlev)
    real (kind=real_kind), intent(inout) :: rho(nlev)

    real (kind=real_kind), dimension(nlev) :: rk0, rk1, rhs
    real (kind=real_kind) :: s1, trm
    integer :: ke

      rk0(:) = rho(:)

   rhs(:) = 0.0D0

   rhs(:) = fd1d_rhs_src(w,rho)


  !RK_stage-1
   do ke=1,nlev
      rk1(ke) = rk0(ke) + dtime * rhs(ke)
   enddo

   rhs(:) = fd1d_rhs_src(w,rk1)

  !RK_stage-2
   do ke=1,nlev
      rk1(ke) = (3.0D0*rk0(ke) + rk1(ke) + dtime * rhs(ke)) *0.25D0
   enddo

   rhs(:) = fd1d_rhs_src(w,rk1)

  !RK_stage-3
   do ke=1,nlev
      rho(ke) = (rk0(ke) + 2.0D0*(rk1(ke) + dtime * rhs(ke)) ) / 3.0D0
   enddo

 end subroutine fd1d_rk3_solver

!=======================================================================================================!
     subroutine print_vertical_zinfo(ie,userdata)
!=======================================================================================================!
    Implicit None
    integer, intent(in)               :: ie
    type (dg_userdata_t) , intent(in) :: userdata

    integer :: i,j,k

   if (userdata%hybrid%par%masterproc .and. (ie == userdata%nets)) then

         !  if (zgrid_type == "GLL") then     
         !    print*, '   --- Vertical GLL grid type--'
         !     elseif (zgrid_type == "FD") then 
         !    print*, '   --- Vertical FD  grid type--'
         !   endif

     if ((zgrid_type == "GLL").or.(zgrid_type == "GL"))  then 
       print*,'-------- Vertical (z) Vital Statistics (RDN) ------------'
       print*,'      '
       print*,'    Vertical Grid Type:  ', zgrid_type 
       print*,'    Test case, class  :  ', test_case
       print*,'    Test name / type  :  ', test_type      

       print*,'    Number of vertical 1d-elements and Gauss points within '
       print*, 'vnel= ', vnel,  ' vnp = ', vnp
       print*,'    Gauss points and weights '
         do j = 1, vnp
          write(*,*) j, vgp(j), vgw(j)
         enddo
       print*,'    Derivative matrix (vgp x vgp) '
         do j = 1, vnp
          write(*,*) vder(:,j)
         enddo

         if (zgrid_type == "GL") then 
          print*,'    Interpolation  matrix (vnp x 2) '
          write(*,*) ime_gl(:,1)
          write(*,*) ime_gl(:,2)

          print*,'    Testfunction at the boundary    '
          write(*,*) hed_gl(:,1)
          write(*,*) hed_gl(:,2)
       print*,'      '
          print*,'    Radau edge derivatives (Lt/Rt)  '
          write(*,*) rad_der(:,1)
          write(*,*) rad_der(:,2)
         endif 

       print *,' -------------------------------'
       print *,' Z-grid parameters: z_top,         z_bot,         z_wel'
       write(*,*)  z_top, z_bot, z_wel

       print *,' ---------------------------------------------'
       print *,'              Vertical Z(Gauss) grid          '
       print *,' ---------------------------------------------'
          do i = 1, vnel
          do j = 1, vnp
              k = (i-1)*vnp + j
           print*, k, z_grd(k)
          enddo
          enddo
       print *,' ---------------------------------------------'

     else

       print *,'       Vertical Z-(Finite Difference) Grid    '
       print *,' Z-grid parameters: z_top,         z_bot,         z_wel'
       write(*,*)  z_top, z_bot, z_wel
       print *,' ---------------------------------------------'
          do i = 1, nlev
            print*, i, z_grd(i)
          enddo
       print *,' ---------------------------------------------'
      endif 
   endif

 end subroutine print_vertical_zinfo

!=======================================================================================================!
     subroutine vertical_FD_zgrid(userdata)
!=======================================================================================================!
    Implicit None
    type (dg_userdata_t) , intent(in) :: userdata
    real (kind=real_kind) ::  dz, z_edg(nlev) 
    integer :: i,j,k

    ! zgrid_type = "FD"   ! to be selected at dg3d_tests_mod 

     ! A simple equi-spaced FD grid  (debug purpose) 
     ! Monotonically increasing from bottom (k=1) to top (k=nlev) 

           z_wel = (z_top - z_bot)/dble(nlev-1)

     !Monotonically increasing z-grid 
        do  k = 1, nlev
           z_edg(k) = z_bot +  z_wel * dble(k-1)
        enddo


        ! Reverse order  
        do k = 1, nlev
         ! z_grd(nlev + 1-k) = z_edg(k)     ! index reversed 
           z_grd(k) = z_edg(k)
        enddo


     End subroutine vertical_FD_zgrid
!!------------------------------------------------------
    Function edge_GL_1d(f1) result(gle)
     Implicit None
     real (kind=real_kind), intent(in), Dimension(vnp) :: f1
     real (kind=real_kind) ::  v1, v2, gle(2)
       Integer :: k
 
 
         v1 = 0.0D0
         v2 = 0.0D0

         do k=1, vnp
          v1  = v1 + ime_gl(k,1)*f1(k)
          v2  = v2 + ime_gl(k,2)*f1(k)
         enddo
 
        !! [gl(-1), gl(+1)] 
         gle(1) = v1   !left/bottom edge 
         gle(2) = v2   !right/top  edge
 
     End Function edge_GL_1d

!=======================================================================================================!
     subroutine vertical_GL_zgrid(userdata)
!=======================================================================================================!
    Implicit None
    type (dg_userdata_t) , intent(in) :: userdata
    type (quadrature_t) :: gp,gs 
    real (kind=real_kind) ::  dz, z_edg(vnel+1) , z_rvs(nlev), ugrid(nlev) 
    real (kind=real_kind) ::  der(vnp,vnp)
    real (kind=longdouble_kind), dimension(0:vnp)  ::  leg,ldr,jac, djac,ldr1 
    real (kind=real_kind), dimension(vnp) ::  wts, bfn
    real (kind=longdouble_kind) :: x, val, zero, alpha,beta , one , tval(2)
    real (kind=real_kind) :: zm, dif,pt,xi,s1, deno, tol,edg(2) 
    integer :: i,j,k, kj, deg 

!    zgrid_type = "GL"   ! to be selected at dg3d_tests_mod 

     zero = 0.0_longdouble_kind
      one = 1.0_longdouble_kind

      tol = 1.0D-12 

 ! Legendre poly from Jacobi-poly by setting 
     alpha = zero 
     beta = zero 

       deg = vnp-1 !highest degree of the polynomial space

! Gauss points and weights in the  vertical dircetion

           gs = gauss(vnp)
       vgp(:) = gs%points(:)    !vertical GL points 
       vgw(:) = gs%weights(:)   !vertical GL weights 

       do k = 0, vnp - 1
         x  = vgp(vnp - k) 
         Call  jacobi(vnp,x,alpha,beta,jac,djac)  !"vnp"-degree Legendre poly
         leg(vnp-k) =  jac(vnp) !legendre poly @ GL pts
         ldr(vnp-k) = djac(vnp) !legendre deri of dgree "vnp" poly  @ GL pts
         ldr1(vnp-k) = djac(deg) !legendre deri of  degree "deg" poly  
       enddo 

         s1 = 0.5D0 * ((-1.0D0)**vnp) 

        do k = 1, vnp
         rad_der(k,1) = 0.5D0 *(ldr(k) + ldr1(k))  !left / ibot 
         rad_der(k,2) =    s1 *(ldr(k) - ldr1(k))  !right / itop 
        enddo 

      !  do k = 0, vnp - 1
      !     x  = vgp(vnp - k) 
      !    Call  jacobi(deg+1,x,alpha,beta,jac,djac)  !"vnp"-degree Legendre poly
      !    leg(vnp-k) =  jac(vnp) !legendre poly @ GL pts
      !    ldr(vnp-k) = djac(vnp) !legendre deri @ GL pts
      !  enddo 

 ! GL derivative  ematrix 

       do i = 1, vnp
         do j = 1, vnp
           if (i == j) then
              vder(i,j) =  vgp(i) / (1.0D0 - vgp(i)*vgp(i) )
             else
              vder(i,j) =  ldr(i)/ (ldr(j) *(vgp(i) - vgp(j)) )
          endif
        enddo
       enddo

  ! Edge target value for GL grid interpolation 

      tval(1) = -one 
      tval(2) =  one 

      do k = 1, 2
            x = tval(k)  !target grid 
           leg(:) =   legendre(x,vnp)    !Legendre polynomial value at "x" 
          do j = 1, vnp
            deno = (x - vgp(j)) * ldr(j)
            ime_gl(j,k) = leg(vnp) / deno
          enddo
      enddo

!  Test function at the boundary [-1,1] 

       do k = 1, vnp
         bfn(:) = 0.0D0 
         bfn(k) = 1.0D0 
         edg(:) = edge_GL_1d(bfn) 
         hed_gl(k,1) = edg(1) 
         hed_gl(k,2) = edg(2) 
       enddo


! Vertical GL-type grid generation  (increasing upward)  

          z_wel = (z_top - z_bot)/dble(vnel)
        do  k = 1, vnel+1
          z_edg(k) = z_bot +  z_wel * dble(k-1)
        enddo

! Non-repeated (GL) physical points at the element edges 

          do k = 1, vnel
           do j = 1, vnp
              kj = (k-1)*vnp + j
                z_grd(kj) = (z_edg(k) + z_edg(k+1) +  z_wel *vgp(j) )/2.0D0
                z_rvs(kj) = z_grd(kj)     
           enddo
          enddo

         z_grd(:) = z_rvs(:) !out file 

         do k = 1, vnel
           do j = 1, vnp
              kj = (k-1)*vnp + j
              z_gsg(j,k) = z_rvs(kj)     
           enddo
          enddo

! Pre-compute weights for vertical uniform data (post-processing) interpolation 
! onto a uniform resolution "ugrid" with "nlev" points 

         dz = (z_top - z_bot)/dble(nlev-1)
        do  k = 1, nlev
           ugrid(k) = z_bot +  dz * dble(k-1)
        enddo

        do  k = 1, nlev
            pt = ugrid(k)
         vindx(k) = int(pt/z_wel) + 1   ! element index of target point
        enddo
         vindx(1) = 1
         vindx(nlev) = vnel

  ! GL interpolation with basis functions 
  ! Create interpolation weights associated with GL grid 

       do k=1, nlev
        kj = vindx(k)
        zm = (z_edg(kj) + z_edg(kj+1))/2.0D0
        xi = 2.0D0* (ugrid(k) - zm)/z_wel    !quadrature grid [-1,1] 

         wts(:) = 0.0D0
        do j = 1, vnp
         dif =  xi - vgp(j)
         if (abs(dif) < tol) then
         wts(j) = 1.0D0
            exit    !from inner loop
         else
          x = xi + zero  !target point  (not longdouble!) 
          leg(:) =   legendre(x,vnp)    !Legendre polynomial value at "x" 
          deno = dif * ldr(j)
          wts(j) = leg(vnp) / deno
         endif
        enddo
         z_wts(:,k) = wts(:) 
       enddo

     End subroutine vertical_GL_zgrid

!=======================================================================================================!
     subroutine vertical_GLL_zgrid(userdata)
!=======================================================================================================!
 ! Monotonically increasing z-grid (from k=1 (bottom) to top k=nlev) 
    Implicit None
    type (dg_userdata_t) , intent(in) :: userdata
    type (quadrature_t) ::  gll 
    real (kind=real_kind) ::  dz, z_edg(vnel+1) , z_rvs(nlev)
    real (kind=real_kind) ::  der(vnp,vnp)
    real (kind=longdouble_kind), dimension(vnp,vnp)  ::  leg,dvv
    real (kind=longdouble_kind) :: zero,one,four
    integer :: i,j,k, kj

!    zgrid_type = "GLL"   ! to be selected at dg3d_tests_mod 

     zero = 0.0_longdouble_kind
      one = 1.0_longdouble_kind
     four = 4.0_longdouble_kind

! Verical Gauss points and weights 
           gll = gausslobatto(vnp)
       vgp(:) = gll%points(:)
       vgw(:) = gll%weights(:)

!  if (vnp == np) then       ! Homme default  
! Vertical derivative matrix 
!       der(:,:) = userdata%deriv%Dvv(:,:)
!      vder(:,:) = Transpose(der(:,:))
!  else

! specified "vnp" number of  GLL points, weights, derivative matrix 

      do i=1,vnp
       leg(:,i) = legendre(gll%points(i),vnp-1)
     end do

    dvv(:,:) = zero
    do j=1,vnp
       do i=1,j-1
          dvv(j,i) = (one/(vgp(i)-vgp(j)))*leg(vnp,i)/leg(vnp,j)
       end do
       dvv(j,j) = zero
       do i=j+1,vnp
          dvv(j,i) = (one/(vgp(i)-vgp(j)))*leg(vnp,i)/leg(vnp,j)
       end do
    end do

    dvv(vnp,vnp) = + vnp*(vnp-1)/four
    dvv(1,1)     = - vnp*(vnp-1)/four

        der(:,:) = dvv(:,:)
       vder(:,:) = Transpose(der(:,:))
      !vder(:,:) = dvv(:,:)
!  endif 

     ! debug data 
     !gllp(1) = -1.00000000000000D0
     !gllp(2) = -0.447213595499958D0
     !gllp(3) =  0.447213595499958D0
     !gllp(4) =  1.00000000000000D0
     !vgw(1) = 0.166666666666667D0
     !vgw(2) = 0.833333333333333D0
     !vgw(3) = 0.833333333333333D0
     !vgw(4) = 0.166666666666667D0
     !vgp(:) = gllp(:) 

     z_wel = (z_top - z_bot)/dble(vnel)
           do  k = 1, vnel+1
            z_edg(k) = z_bot +  z_wel * dble(k-1)
           enddo

        ! repeated physical points at the element edges 
          do k = 1, vnel
           do j = 1, vnp
              kj = (k-1)*vnp + j
                z_grd(kj) = (z_edg(k) + z_edg(k+1) +  z_wel *vgp(j) )/2.0D0
                z_rvs(kj) = z_grd(kj) 
           enddo
          enddo

         z_grd(:) = z_rvs(:) !out file 

         do k = 1, vnel
           do j = 1, vnp
              kj = (k-1)*vnp + j
                   z_gsg(j,k) = z_rvs(kj)     !Noreverse 
           enddo
          enddo

End subroutine vertical_GLL_zgrid

!=======================================================================================================!
  Function  Compute_SoundSpeed(prs,rho) result(speed)
!=======================================================================================================!
    Implicit None
    real (kind=real_kind) :: rho, speed, prs

     speed =  sqrt(Cp_by_Cv * prs / rho)

   End Function Compute_SoundSpeed
!=======================================================================================================!
  Function  Compute_Pressure(rho,the) result(prs)
!=======================================================================================================!
    Implicit None
    real (kind=real_kind) :: rho, the, prs

     prs = CC_0 * (rho * the)**Cp_by_Cv 

   End Function Compute_Pressure

!=======================================================================================================!
  Function vertical_unif_interpolate(ff) result(fout)
!=======================================================================================================!
    Implicit None
    real (kind=real_kind), intent(in) :: ff(np,np,nlev)
    real (kind=real_kind) :: fout(np,np,nlev)
    real (kind=real_kind) :: fe(vnp,vnel)
    real (kind=real_kind) :: f1d(nlev)
    real (kind=real_kind) :: s1
    integer :: i,j,k,el,l

 ! Vertical GL interpolation onto a uniform grid using 
 ! the pre-computed 1d weights "z_wts" 

      if (zgrid_type =="GL") then

         do j=1,np    
         do i=1,np   

        ! element-wise decomposition 
         fe(:,:) = copyto_1d_elem(ff(i,j,:)) 

       do k=1, nlev
         el = vindx(k)
         s1 = 0.0D0 
         do l=1, vnp 
          s1 = s1 + fe(l,el) * z_wts(l,k) 
         enddo 
       fout(i,j,k) = s1 
       enddo 
    end do 
    end do 

    else   
      fout = ff     !No interpolation 
    endif 

 End Function vertical_unif_interpolate

!=======================================================================================================!
  Function vertical_z_integral(f1d) result(v_int)
    Implicit None
    real (kind=real_kind), intent(in) :: f1d(nlev)
    real (kind=real_kind)             :: fe(vnp,vnel)
    real (kind=real_kind) :: v_int, s1, dz 
    integer :: k,el,l

      dz = z_wel * 0.5D0 

  if (zgrid_type =="FD") then 
       
   !Cell-centered averaging 
        s1 = 0.0D0 
      do  k = 1, nlev-1
        s1 = s1 + (f1d(k) + f1d(k+1))*0.5D0 
      enddo 
        v_int = s1*z_wel  !uniform vertical grid with dz = z_wel 

   ! To normalize the integral 
      v_int = v_int / (z_top - z_bot) 

  else     !high-order methods 

  ! element-wise decomposition 
    do el=1,vnel
    do  k=1,vnp
     l = k + (el-1) * vnp
     fe(k,el) = f1d(l)
    end do
    end do

      v_int = 0.0D0 

   !Vertical integration 
    do el=1,vnel
      s1 = 0.0D0 
    do  k=1,vnp
      s1 = s1 + fe(k,el) * vgw(k) 
    end do
    ! v_int = v_int + s1*dz       
      v_int = v_int + s1* 0.5D0     !normalized 
    end do

   endif 

 end function vertical_z_integral 

!=======================================================================================================!
  function copyto_1d_elem(f1d) result(fe)       
    Implicit None
    real (kind=real_kind), intent(in) :: f1d(nlev)
    real (kind=real_kind)             :: fe(vnp,vnel)
    integer :: k,el,l 

    do el=1,vnel
    do  k=1,vnp   
     l = k + (el-1) * vnp  
     fe(k,el) = f1d(l) 
    end do 
    end do 
 end function copyto_1d_elem
!=======================================================================================================!
  function copyfrom_1d_elem(fe) result(f1d)
    Implicit None
    real (kind=real_kind),intent(in) :: fe(vnp,vnel)
    real (kind=real_kind)             :: f1d(nlev)
    integer :: k,el,l

    do el=1,vnel
    do  k=1,vnp
     l = k + (el-1) * vnp
     f1d(l) = fe(k,el)
    end do
    end do

 end function copyfrom_1d_elem

!=======================================================================================================!
  function copyfrom_1dsys_elem(neqn,fe) result(ff_column)
    Implicit None
    integer,intent(in) :: neqn            
    real (kind=real_kind),intent(in) :: fe(vnp,vnel,neqn)
    real (kind=real_kind)             :: ff_column(neqn,nlev)
    integer :: k,el,l, eq 

    do eq =1,neqn 

      do el=1,vnel
      do  k=1,vnp
         l = k + (el-1) * vnp
         ff_column(eq,l) = fe(k,el,eq)
      end do
      end do

    end do

 end function copyfrom_1dsys_elem
!=======================================================================================================!
End module dg3d_zdynamics_mod
