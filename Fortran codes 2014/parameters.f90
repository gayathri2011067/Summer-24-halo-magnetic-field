module parameters
implicit none    
  

  !GRID PARAMETERS
  integer,parameter :: Nz = 51               !spatial resolution
  integer,parameter :: Nt = 100000           !time resolution


  !_________________________________________________________________________


  !  CONSTANTS FOR UNIT CONVERSIONS
  double precision, parameter :: pi = 3.14156295358
  double precision, parameter :: g_mp =   1.67262158d-24
  double precision, parameter :: cm_kpc = 3.08567758d21
  double precision, parameter :: km_kpc = 3.08567758d16
  double precision, parameter :: cm_km =  1.d5
  double precision, parameter :: s_Gyr =  1.d9*365.25d0*24*3600
  double precision, parameter :: G_muG = 1.d6


  !__________________________________________________________________________



  !  DIMENSIONLESS UNITS
  double precision, parameter :: eta = 1.
  double precision, parameter :: h =  1.
  double precision, parameter :: B0 =  1.
  double precision, parameter :: td =  1.

  !  DISC PARAMETERS
  double precision, parameter :: radius = 10. 
  double precision, parameter :: U0 = 2.38

  !  DISC PARAMETERS
  double precision, parameter :: R0 = 8.5
  double precision, parameter :: z0 = 0.1
  double precision, parameter :: rho0 = 0.020
  double precision, parameter :: h0 = 0.1



  !_____________________________________________________________________________

  !  B_eq
  double precision, parameter :: B_eq = 3.25d-6
  

  !______________________________________________________________________________

  !  ALPHA PROFILE
  double precision, parameter :: alpha = 0.1





end module parameters