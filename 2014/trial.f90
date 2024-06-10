program main
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none

    ! Define constants
    real(wp), parameter :: kpc_to_cm = 3.085677581e21_wp
    real(wp), parameter :: km_to_cm = 1.0e5_wp
    real(wp), parameter :: s_to_Gyr = 3.16887646e-17_wp
    real(wp), parameter :: pi = 3.14159265358979323846_wp

    ! Define parameters
    integer, parameter :: T = 17
    integer, parameter :: Nt = 7000 * T
    integer, parameter :: Nz = 101

    real(wp) :: radius_dim, r_d_dim, h_d_dim, l_dim, you, eta_dim, h_dim, td_dim
    real(wp) :: B_0_dim, omega_0_dim, r_omega_dim, omega_dim, alpha_0_dim, G_dim
    real(wp) :: U_0_dim, k_dim, R_dim, z_i_dim, z_f_dim, z_i, z_f, dr, dt
    real(wp), allocatable :: r_arr(:), B_eq(:), U_z(:), alpha_tilda(:)
    real(wp) :: eta, omega_0, r_omega, h_d, r_d, l, U_0, k, B_0, R, radius, h, td, alpha_0, G
    integer :: i, r0_index

    ! Assign values to parameters
    radius_dim = 4.0_wp * kpc_to_cm
    r_d_dim = 10.0_wp * kpc_to_cm
    h_d_dim = 0.35_wp * kpc_to_cm
    l_dim = 0.1_wp * kpc_to_cm
    you = 10.0_wp * km_to_cm / s_to_Gyr
    eta_dim = (1.0_wp / 3.0_wp) * l_dim * you
    h_dim = 0.381_wp * kpc_to_cm
    td_dim = (h_dim**2) / (eta_dim / (kpc_to_cm**2) * s_to_Gyr)
    B_0_dim = 8.2e-6_wp
    omega_0_dim = 127.0_wp * km_to_cm / kpc_to_cm
    r_omega_dim = 2.0_wp * kpc_to_cm

    write(*,*) "eta_dim = ", eta_dim / (kpc_to_cm**2 / s_to_Gyr)

    omega_dim = omega_0_dim * (1.0_wp + (radius_dim / r_omega_dim)**2)**(-0.5_wp)
    alpha_0_dim = (l_dim**2) * omega_dim / h_dim
    G_dim = -45.6_wp * km_to_cm / kpc_to_cm

    U_0_dim = 1.0_wp * km_to_cm / s_to_Gyr
    k_dim = 0.1_wp * km_to_cm * kpc_to_cm / s_to_Gyr

    R_dim = 20.0_wp * kpc_to_cm
    z_i_dim = -h_dim
    z_f_dim = h_dim
    z_i = z_i_dim / h_dim
    z_f = z_f_dim / h_dim

    allocate(r_arr(Nz), B_eq(Nz), U_z(Nz), alpha_tilda(Nz))
    r_arr = [(z_i + (z_f - z_i) * real(i-1, wp) / (Nz-1), i = 1, Nz)]

    eta = eta_dim / eta_dim
    omega_0 = omega_0_dim * (h_dim / eta_dim)
    r_omega = r_omega_dim / z_f_dim
    h_d = h_d_dim / h_dim
    r_d = r_d_dim / h_dim
    l = l_dim / h_dim

    U_0 = (U_0_dim * eta_dim / h_dim**2)
    k = (k_dim * h_dim**2 / eta_dim)
    B_0 = B_0_dim / B_0_dim
    R = R_dim / h_dim
    radius = radius_dim / h_dim
    h = h_dim / h_dim
    td = td_dim * eta_dim / h_dim**2
    B_eq = B_0 * exp(-(radius/R) - r_arr**2 / (2.0_wp * h**2))

    U_z = r_arr
    alpha_0 = alpha_0_dim * h_dim / eta_dim
    G = G_dim * h_dim**2 / eta_dim

    alpha_tilda = sin(pi * r_arr)

    ! Write data to a file
    open(unit=10, file='output_data.txt', status='replace')
    write(10, '(A)') "r B_r B_phi alpha_m"
    do i = 1, Nz
        write(10, '(F10.6, 3X, F10.6, 3X, F10.6, 3X, F10.6)') r_arr(i), B_eq(i), 0.0_wp, 0.0_wp
    end do
    close(10)

    deallocate(r_arr, B_eq, U_z, alpha_tilda)

end program main

