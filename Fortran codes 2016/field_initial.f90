module initial_field

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use eta_profile
!
  implicit none
!
  double precision, dimension(nx) :: B_r, B_phi, dBr, d2Br, dBphi, d2Bphi
  


! 
contains
    subroutine field_initialization
      
        B_r = sin(pi*x)
        B_phi = 0.0

    end subroutine field_initialization
!
  end module initial_field