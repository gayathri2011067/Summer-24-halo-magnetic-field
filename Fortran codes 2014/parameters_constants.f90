module parameters

  implicit none

!  CONSTANT VALUES

  double precision, parameter :: pi= 3.14156295358
  double precision, parameter :: g_mp=   1.67262158d-24
  double precision, parameter :: cm_kpc= 3.08567758d21
  double precision, parameter :: km_kpc= 3.08567758d16
  double precision, parameter :: cm_km=  1.d5
  double precision, parameter :: s_Gyr=  1.d9*365.25d0*24*3600
  double precision, parameter :: G_muG= 1.d6

!***************************************************************************************************************

! WITH DIMENSION

double precision, parameter :: radius_dim = 4.d0 ! kpc
double precision, parameter :: r_d_dim = 10.d0 ! kpc
double precision, parameter :: h_d_dim = 0.35d0 ! kpc
double precision, parameter :: eta_dim = 0.1d26*s_Gyr/(cm_kpc**2.) ! cm2/s --> kpc2/Gyr
double precision, parameter :: h_dim = h_d_dim !kpc !CHANGE: if disc flaring
!double precision, parameter :: h_dim = -eqxn- !NOTE: DISC FLARING
double precision, parameter :: t_d_dim = h_dim**2./eta_dim ! Gyr
double precision, parameter :: omega_0_dim = 127.*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
double precision, parameter :: r__omega_dim = 2. ! kpc
double precision, parameter :: l_dim = 0.1 ! kpc
double precision, parameter :: omega_dim = omega_0_dim*(1.+(radius_dim/r__omega_dim)**2.)**(-0.5) ! 1/Gyr
! double precision, parameter :: G_dim= -omega_dim ! 1/Gyr !NOTE: from paper
double precision, parameter :: G_dim = -45.6*s_Gyr/km_kpc   ! km/s.kpc --> 1/Gyr
double precision, parameter :: alpha_0_dim = (l_dim**2.)*omega_dim/(h_dim) ! kpc/Gyr
double precision, parameter :: U_0_dim = 1.*s_Gyr/km_kpc ! km/s --> kpc/Gyr
double precision, parameter :: k_dim = 0.1*s_Gyr/km_kpc !km.kpc/s --> kpc**2/Gyr
double precision, parameter :: R_dim = 20.!kpc
double precision, parameter :: z_i_dim = -h_dim!kpc
double precision, parameter :: z_f_dim = +h_dim !kpc

!***************************************************************************************************************

! DIMENSIONLESS

double precision, parameter :: radius = radius_dim/h_dim
double precision, parameter :: r_d = r_d_dim/h_dim
double precision, parameter :: h_d = h_d_dim/h_dim
double precision, parameter :: eta = eta_dim/(eta_dim) !NOTE: use profile for thick disc
double precision, parameter :: h = h_dim/h_dim
double precision, parameter :: t_d = t_d_dim/(h_dim**2./eta_dim)
double precision, parameter :: omega_0 = omega_0_dim*h_dim/(eta_dim) 
double precision, parameter :: r__omega = r__omega_dim/h_dim
double precision, parameter :: l = l_dim/h_dim
double precision, parameter :: omega = omega_dim*h_dim/(eta_dim)
double precision, parameter :: G = G_dim*(h_dim**2/(eta_dim))
double precision, parameter :: alpha_0 = alpha_0_dim*h_dim/eta_dim
double precision, parameter :: U_0 = U_0_dim*h_dim/eta_dim !NOT_SURE!CHECK: old code might be wrong
double precision, parameter :: k = k_dim/(eta_dim) !NOT_SURE
double precision, parameter :: R = R_dim/h_dim
double precision, parameter :: z_i = z_i_dim/h_dim
double precision, parameter :: z_f = z_f_dim/h_dim

!*******************************************************************************************************************************

double precision, parameter :: R_alpha = alpha_0
double precision, parameter :: R_omega = G
double precision, parameter :: R_k = k
double precision, parameter :: R_U = U_0

!********************************************************************************************************************************














end module parameters