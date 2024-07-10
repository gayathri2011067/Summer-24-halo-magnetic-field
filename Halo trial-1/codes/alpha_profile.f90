module alpha_profile

  use parameters
  use time_grid
  use physical_grid
  use make_a_grid
  use omega_profile
  ! use spatial_derivatives
!
  implicit none
!
 
  double precision, dimension(nx) :: alpha_k
  double precision, dimension(nx) :: alpha_cap,d_alpha_cap,d2alpha_cap, alpha_cap2

 

! 
contains
    subroutine construct_alpha_profile
        character(len=30) :: ghost_zone_type2 = 'relative anti-symmetric'
        ! integer :: i                                                          

        ! do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys !  CHANGED!
        !     alpha_cap(i) = sin(pi*x(i))
        !     alpha_k(i) = alpha_0*alpha_cap(i)
        !     d_alpha_cap(i) = cos(pi*x(i))
        ! end do
        ! alpha_cap = sin(pi*x)
        ! alpha_k = alpha_0*alpha_cap
        ! d_alpha_cap = cos(pi*x)
        !CHANGE: the above was first halo implementation, now using just the alpha equation.
            alpha_cap = sin((pi*x)/h)          !CONCERN: h need to be the disc height, !TRY: h = 1.0
            alpha_k = alpha_0*alpha_cap
            d_alpha_cap = cos((pi*x)/h)


      
  

    end subroutine construct_alpha_profile
!
  end module alpha_profile