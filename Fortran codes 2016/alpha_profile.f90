module alpha_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use omega_profile
!
  implicit none
!
 
  double precision, dimension(nx) :: alpha_cap
  double precision, dimension(nx) :: alpha_fz

 

! 
contains
    subroutine construct_alpha_profile
        alpha_cap = x 
        alpha_fz = omega_fzr*(l**2/h)*(x/h)*exp(1/2 - (x/h)**2)*sqrt(2.)

    end subroutine construct_alpha_profile
!
  end module alpha_profile