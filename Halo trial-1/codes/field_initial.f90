module initial_field

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
  use alpha_profile
  use velocity_profile
!
  implicit none
!
  double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq
  double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br
  double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
  double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi



! 
contains
    subroutine field_initialization
        call construct_alpha_profile
        call construct_velocity_profile
        ! integer :: 

        ! do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys ! CHANGED!
        !   B_R(i) = 0.0001*(1.0-x(i)**2.)*exp(-x(i)**2.)
        !   B_phi(i) = 0.0
        !   B_eq(i) = exp(-radius/R - x(i)**2./2.)
        !   alpha_Br(i) = B_r(i)*alpha_cap/(1+(B_r(i)**2+B_phi(i)**2)/B_eq(i)**2)
        !   alpha_Bphi(i) = B_phi(i)*alpha_cap/(1+(B_r(i)**2+B_phi(i)**2)/B_eq(i)**2)
        !   Uz_Br(i) = B_r(i)*U_z_cap
        !   Uz_Bphi(i) = B_phi(i)*U_z_cap
        ! end do

        B_r = 0.0001*(1.0-(x**2./h))*exp((-x**2./h))                       !CHANGE: the old implementation is given above
        B_phi = 0.0                                                        !CHANGE: put a /h wherever there is x !CONCERN: figure this out
        B_eq = exp(-radius/R - x**2./(2.*h))
        alpha_Br = B_r*alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        alpha_Bphi = B_phi*alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        Uz_Br = B_r*U_z_cap
        Uz_Bphi = B_phi*U_z_cap


    end subroutine field_initialization
!
  end module initial_field