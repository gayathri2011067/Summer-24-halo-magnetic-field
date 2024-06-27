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
 
  double precision, dimension(nx) :: alpha_k
  double precision, dimension(nx) :: alpha_cap,d_alpha_cap

 

! 
contains
    subroutine construct_alpha_profile
        character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
        alpha_cap = sin(pi*x)
        alpha_k = alpha_0*alpha_cap
        d_alpha_cap = cos(pi*x)
  

    end subroutine construct_alpha_profile
!
  end module alpha_profile