program run_all

    use alpha_profile
    use velocity_profile
    use initial_field
    use parameters
    use time_grid
    use physical_grid
    use make_a_grid
    use eta_profile
    use omega_profile
    use equations
    use spatial_derivatives
    use timestepping
  
    implicit none

    integer :: i, j
    character(len=20) :: filename, xfile, omegafile, alphafile, Br_ini_file,B_phi_ini_file,&
    
    B_r_final_file,B_phi_final_file, time_file, first_der,second_der

    ! Call subroutine to construct the physical grid
    call construct_grid

    ! Call subroutine to construct eta profile
    call construct_eta_profile
    call construct_omega_profile
    call construct_alpha_profile
    call field_initialization


    ! Define the output file name
    filename = 'eta_fz_values.txt'
    xfile= 'z_values.txt'
    omegafile= 'omega_values.txt'
    alphafile= 'alpha_values.txt'
    Br_ini_file='Br_ini.txt'
    B_phi_ini_file='B_phi_ini.txt'
    B_r_final_file='Br_final.txt'
    B_phi_final_file='B_phi_final.txt'
    time_file='time.txt'
    first_der='first_derivative.txt'
    second_der='second_derivative.txt'




    ! Open the file for writing
    open(unit=10, file=filename)
    open(unit=17, file=xfile)
    open(unit=18, file=omegafile)
    open(unit=19, file=alphafile)
    open(unit=20, file=Br_ini_file)
    open(unit=21, file=B_phi_ini_file)
    open(unit=22, file=B_r_final_file)
    open(unit=23, file=B_phi_final_file)
    open(unit=24, file=time_file)
    open(unit=25, file=first_der)
    open(unit=26, file=second_der)
    ! Write the values to the file
    do i = 1, nx
        write(10, '(F12.8)') eta_fz(i)
        write(17, '(F12.8)') x(i)
        write(18, '(F12.8)') omega_fzr(i)
        write(19, '(F12.8)') alpha_fz(i)
        write(20, '(F12.8)') B_r(i)
        write(21, '(F12.8)') B_phi(i)

    end do

    ! Close the file
    close(10)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)

    print *, 'eta_fz values have been saved to ', filename
    print *, 'z values have been saved to ', xfile
    print *, 'omega values have been saved to ', omegafile
    print *, 'alpha values have been saved to ', alphafile
    print *, 'Br_ini values have been saved to ', Br_ini_file
    print *, 'B_phi_ini values have been saved to ', B_phi_ini_file
    

  call field_initialization
  call RK4
    write(25, '(F12.8)') dBr
    write(26, '(F12.8)') d2Br

    close(22)
    close(23)
    close(24)
    close(25)
    close(26)
  print*, 'n2=', n2
  print*, 'dt=', dt
  print*, 'Nt=', Nt
  print*, 'n1=', n1
  print*, 'total_t=', total_t
  print*, 't=', t
  print*, 'first=', first
end program run_all
