module initial_field

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
!
  implicit none
!
  double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi, B_eq
  double precision, dimension(nx) :: alpha_Br, alpha_Bphi, Uz_Br, Uz_Bphi,d_alpha_Br
  double precision, dimension(nx) :: d2_alpha_Br,d_alpha_Bphi,d2_alpha_Bphi,d_Uz_Br
  double precision, dimension(nx) :: d2_Uz_Br,d_Uz_Bphi,d2_Uz_Bphi, old_Br, old_Bphi
  double precision, parameter :: alp=dt/dx**2  !alpha for finite difference scheme
  double precision, parameter :: beta= dt/dx !beta for fin diff 



! 
contains
    subroutine field_initialization
      
        B_r = 0.0001*(1.0-x**2.)*exp(-x**2.)
        B_phi = 0.0
        B_eq = exp(-radius/R - x**2./2.)
        alpha_Br = B_r*alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        alpha_Bphi = B_phi*alpha_cap/(1+(B_r**2+B_phi**2)/B_eq**2)
        Uz_Br = B_r*U_z_cap
        Uz_Bphi = B_phi*U_z_cap


    end subroutine field_initialization
!
  end module initial_field