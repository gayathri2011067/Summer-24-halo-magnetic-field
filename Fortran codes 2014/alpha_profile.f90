module alpha_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use omega_profile
  use spatial_derivatives
!
  implicit none
!
 
  double precision, dimension(nx) :: alpha_cap
  double precision, dimension(nx) :: alpha_fz,d_alpha,d2_alpha

 

! 
contains
    subroutine construct_alpha_profile
        character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
        alpha_cap = x 
        alpha_fz = omega_fzr*(l**2/h)*(x/h)*exp(1/2 - (x/h)**2)*sqrt(2.)
        call impose_boundary_conditions(alpha_fz, ghost_zone_type2)
        call spatial_derivative(alpha_fz,6,d_alpha,d2_alpha)


    end subroutine construct_alpha_profile
!
  end module alpha_profile