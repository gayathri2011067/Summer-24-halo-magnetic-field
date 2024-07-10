module spatial_derivatives

    use time_grid
    use physical_grid
    use make_a_grid

    implicit none
 
contains

    subroutine impose_boundary_conditions(fun, ghost_zone_type)

    character(len=30), intent(in) :: ghost_zone_type
    double precision, dimension(nx), intent(inout) :: fun
    integer :: i
    
    if (nxghost/=0) then
      select case (ghost_zone_type)
      case ('symmetric')
        do i = 1, nxghost
          fun(i)=  fun(2*(nxghost+1)-i)  !Symmetric     about z=-h    Neumann   BC on Br, Bp: dBrdz=dBpdz=0 at z=-h
          fun(nx+1-i)=  fun(nx+1-2*(nxghost+1)+i)
        end do
      case ('anti-symmetric')
        do i = 1, nxghost
            fun(i)= -fun(2*(nxghost+1)-i)  !Antisymmetric about z=-h    Dirichlet BC on Br, Bp
            fun(nx+1-i)= -fun(nx+1-2*(nxghost+1)+i)
        end do
      case ('relative anti-symmetric')
        do i = 1, nxghost
            fun(i) = 2.*fun(nxghost+1) -fun(2*(nxghost+1)-i)  !Relative antisymmetric   (Specify alpha_m in ghost zones)
            fun(nx+1-i) = 2.*fun(nx+1-(nxghost+1)) - fun(nx+1-2*(nxghost+1)+i)
        end do
      case default
        print *, 'Invalid boundary condition'
        return
      end select
    end if

    end subroutine impose_boundary_conditions

        
  subroutine spatial_derivative(fun, fd_order, first_derivative, second_derivative)
    
    integer, intent(in) :: fd_order
    double precision, dimension(nx), intent(in) :: fun
    double precision, dimension(nx), intent(inout) :: first_derivative, second_derivative
    integer :: i

  
   

    ! ! CHANGED !
    ! Compute finite difference derivatives
    ! select case (fd_order)
    ! case (2)
    !   do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys
    !     first_derivative(i) = (fun(i+1) - fun(i-1)) / (2.0 * dx)
    !     second_derivative(i) = (fun(i-1) - 2.0 * fun(i) + fun(i+1)) / (dx ** 2)
    !   end do
    ! case (4)
    !   do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys
    !     first_derivative(i) = (fun(i-2) - 8.0 * fun(i-1) + 8.0 * fun(i+1) - fun(i+2)) / (12.0 * dx)
    !     second_derivative(i) = (-fun(i-2) + 16.0 * fun(i-1) - 30.0 * fun(i) &
    !     + 16.0 * fun(i+1) - fun(i+2)) / (12.0 * (dx ** 2))
    !   end do
    ! case (6)
    !   do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys
    !     first_derivative(i) = (-fun(i-3) + 9.0 * fun(i-2) - 45.0 * fun(i-1) + 45.0 * fun(i+1) &
    !     - 9.0 * fun(i+2) + fun(i+3)) / (60.0 * dx)
    !     second_derivative(i) = (2.0 * fun(i-3) - 27.0 * fun(i-2) + 270.0 * fun(i-1) &
    !     - 490.0 * fun(i) + 270.0 * fun(i+1) - 27.0 * fun(i+2) + 2.0 * fun(i+3)) / (180.0 * (dx ** 2))
    !   end do
    ! ! case (8)
    ! !   do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys
    ! !     first_derivative(i) = (3.0 * fun(i-4) - 32.0 * fun(i-3) + 168.0 * fun(i-2) &
    ! !     - 672.0 * fun(i-1) + 672.0 * fun(i+1) - 168.0 * fun(i+2) &
    ! !     + 32.0 * fun(i+3) - 3.0 * fun(i+4)) / (840.0 * dx)
    ! !     second_derivative(i) = (-9.0 * fun(i-4) + 128.0 * fun(i-3) - 1008.0 * fun(i-2) &
    ! !     + 8064.0 * fun(i-1) - 14350.0 * fun(i) + 8064.0 * fun(i+1) &
    ! !     - 1008.0 * fun(i+2) + 128.0 * fun(i+3) - 9.0 * fun(i+4)) / (5040.0 * (dx ** 2))
    ! !   end do
    ! ! case (10)
    ! !   do i=1+nxghost+nxvacuum, 1 + nxghost + nxvacuum + nxphys
    ! !     first_derivative(i) = (-2.0 * fun(i-5) + 25.0 * fun(i-4) - 150.0 * fun(i-3) &
    ! !     + 600.0 * fun(i-2) - 2100.0 * fun(i-1) + 2100.0 * fun(i+1) - 600.0 * fun(i+2) &
    ! !     + 150.0 * fun(i+3) - 25.0 * fun(i+4) + 2.0 * fun(i+5)) / (2520.0 * dx)
    ! !     second_derivative(i) = (8.0 * fun(i-5) - 125.0 * fun(i-4) + 1000.0 * fun(i-3) &
    ! !     - 6000.0 * fun(i-2) + 42000.0 * fun(i-1) - 73766.0 * fun(i) + 42000.0 * fun(i+1) &
    ! !     - 6000.0 * fun(i+2) + 1000.0 * fun(i+3) - 125.0 * fun(i+4) + 8.0 * fun(i+5)) / (25200.0 * (dx ** 2))
    ! !   end do !NOTE: The 10th  and 8th order FD scheme is not implemented yet, but can be added after changing the number of ghost zones.


    ! case default
    !   print *, 'Invalid order of finite difference'
    !   return
    ! end select


    select case (fd_order)    !CHANGE: old code above
    case (2)
      do i = 1+nxghost, nxphys+nxghost
        first_derivative(i) = (fun(i+1) - fun(i-1)) / (2.0 * dx)
        second_derivative(i) = (fun(i-1) - 2.0 * fun(i) + fun(i+1)) / (dx ** 2)
      end do
    case (4)
      do i = 1+nxghost, nxphys+nxghost
        first_derivative(i) = (fun(i-2) - 8.0 * fun(i-1) + 8.0 * fun(i+1) - fun(i+2)) / (12.0 * dx)
        second_derivative(i) = (-fun(i-2) + 16.0 * fun(i-1) - 30.0 * fun(i) &
        + 16.0 * fun(i+1) - fun(i+2)) / (12.0 * (dx ** 2))
      end do
    case (6)
      do i = 1+nxghost, nxphys+nxghost
        first_derivative(i) = (-fun(i-3) + 9.0 * fun(i-2) - 45.0 * fun(i-1) + 45.0 * fun(i+1) &
        - 9.0 * fun(i+2) + fun(i+3)) / (60.0 * dx)
        second_derivative(i) = (2.0 * fun(i-3) - 27.0 * fun(i-2) + 270.0 * fun(i-1) &
        - 490.0 * fun(i) + 270.0 * fun(i+1) - 27.0 * fun(i+2) + 2.0 * fun(i+3)) / (180.0 * (dx ** 2))
      end do
    ! case (8)
    !   do i = 1+nxghost, nxphys+nxghost
    !     first_derivative(i) = (3.0 * fun(i-4) - 32.0 * fun(i-3) + 168.0 * fun(i-2) &
    !     - 672.0 * fun(i-1) + 672.0 * fun(i+1) - 168.0 * fun(i+2) &
    !     + 32.0 * fun(i+3) - 3.0 * fun(i+4)) / (840.0 * dx)
    !     second_derivative(i) = (-9.0 * fun(i-4) + 128.0 * fun(i-3) - 1008.0 * fun(i-2) &
    !     + 8064.0 * fun(i-1) - 14350.0 * fun(i) + 8064.0 * fun(i+1) &
    !     - 1008.0 * fun(i+2) + 128.0 * fun(i+3) - 9.0 * fun(i+4)) / (5040.0 * (dx ** 2))
    !   end do
    ! case (10)
    !   do i = 1+nxghost, nxphys+nxghost
    !     first_derivative(i) = (-2.0 * fun(i-5) + 25.0 * fun(i-4) - 150.0 * fun(i-3) &
    !     + 600.0 * fun(i-2) - 2100.0 * fun(i-1) + 2100.0 * fun(i+1) - 600.0 * fun(i+2) &
    !     + 150.0 * fun(i+3) - 25.0 * fun(i+4) + 2.0 * fun(i+5)) / (2520.0 * dx)
    !     second_derivative(i) = (8.0 * fun(i-5) - 125.0 * fun(i-4) + 1000.0 * fun(i-3) &
    !     - 6000.0 * fun(i-2) + 42000.0 * fun(i-1) - 73766.0 * fun(i) + 42000.0 * fun(i+1) &
    !     - 6000.0 * fun(i+2) + 1000.0 * fun(i+3) - 125.0 * fun(i+4) + 8.0 * fun(i+5)) / (25200.0 * (dx ** 2))
    !   end do !NOTE: The 10th  and 8th order FD scheme is not implemented yet, but can be added after changing the number of ghost zones.


    case default
      print *, 'Invalid order of finite difference'
      return
    end select

  end subroutine spatial_derivative
      
end module spatial_derivatives
