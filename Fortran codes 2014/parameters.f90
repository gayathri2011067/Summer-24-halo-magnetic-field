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


!  SWITCHES

  logical :: Damp=       .false.  !Set to 0 for FOSA, 1 for minimal tau approximation
  logical :: Alg_quench= .false.  !Works with dyn_quench=F; Set to 1 for algebraic quenching (alpha= alpha_k/(1+Emag/Beq^2))
  logical :: Dyn_quench= .true.  !Works with alg_quench=F; Set to 1 for dynamical quenching (d_alpha_m/dt eqn incl in sim)
  logical :: Alp_sin=    .true.  !Set to 1 to get a sinusoidal alpha profile; Set to 0 for a linear alpha profile.
  logical :: Alp_squared=.true.  !Set to 1 to include alpha^2 effect; set to 0 to use alpha-omega approximation equations
  logical :: Shear=      .true.  !Set to 1 to include Omega effect in dynamo, 0 for alpha^2 dynamo

!***************************************************************************************************************


!  DIMENSIONLESS UNITS

  double precision, parameter :: etat=1.
  double precision, parameter :: h0=   1.
  double precision, parameter :: B0=  1.
  double precision, parameter :: Bseed=1.d-3  !Amplitude of seed magnetic field, as a fraction of Beq

!***************************************************************************************************************


!  DISC PARAMETERS

  double precision, parameter :: r0_kpc=10.  !Fiducial radius
  double precision, parameter :: V0_kms=250.  !Circular rotation speed at r=r0_kpc in km/s
  double precision, parameter :: r_kpc=4.  !Radius at which simulation is performed

!***************************************************************************************************************


!  TURBULENCE

  double precision, parameter :: l_kpc= 0.1  !Used in Krause's law; Size of largest turbulent eddies, in parsecs
  double precision, parameter :: h0_kpc= 0.5  !The designated unit of length, corresp to a half-disc thickness at fiducial radius r=r0
  double precision, parameter :: h0_km=h0_kpc*km_kpc  !The designated unit of length, corresp to a typical half-disk thickness, in km
  double precision, parameter :: v_turb_kms= 10.  !Typical turbulent velocity; together with h0_kpc it determines the unit of time for the simulation
  double precision, parameter :: etat_cm2s= l_kpc*cm_kpc*v_turb_kms*cm_km/3  !Turbulent diffusivity in units of cm^2/s
  double precision, parameter :: td0_Gyr= h0_kpc**2/etat_cm2s/s_Gyr*cm_kpc*cm_kpc  !Typical vertical turbulent diffusion timescale in units of Gyr
  double precision, parameter :: td0_s= td0_Gyr*s_Gyr  !Typical vertical turbulent diffusion timescale in units of seconds
  double precision, parameter :: nH0_cm3=0.1  !Number density in cm^{-3} at r=0 (only used to get physical value for B0)
  double precision, parameter :: rho0_gcm3=nH0_cm3*g_mp  !Density in g/cm^3 (only used to get physical value for B0)
  double precision, parameter :: B0_muG= sqrt(4*pi*rho0_gcm3)*v_turb_kms*cm_km*G_muG
  double precision, parameter :: req_kpc= 20.  !Scale length of equiparition field in kpc

!***************************************************************************************************************


! U, SCALE HEIGHT, ALPHA AND OMEGA

!------------------------------------------------------------------------------------------------------------------
!  VELOCITY U_z
  double precision, parameter :: U0_kms=0.  !Vertical mean velocity in km/s
!------------------------------------------------------------------------------------------------------------------
!   SCALE HEIGHT
  double precision, parameter :: r_D_kpc= 10.  !Relevant only if Flaring=1; characteristic radius of the hyperbolically varying scale height in kpc
  double precision, parameter :: h_D_kpc= h0_kpc/(1.+(r0_kpc/r_D_kpc)**2)**(1./2)  !Relevant only if Flaring=1; amplitude of the scale height in kpc
  double precision, parameter :: h_kpc= h_D_kpc*(1.+(r_kpc/r_D_kpc)**2)**(1./2)  !Disk scale height at r_kpc where simulation is performed
  double precision, parameter :: td_Gyr= h_kpc**2/etat_cm2s/s_Gyr*cm_kpc*cm_kpc  !Typical vertical turbulent diffusion timescale in units of Gyr
  double precision, parameter :: td_s= td_Gyr*s_Gyr  !Typical vertical turbulent diffusion timescale in units of seconds
!------------------------------------------------------------------------------------------------------------------
!  ROTATION CURVE
  double precision, parameter :: r_om_kpc= 2.  !Relevant only if Om_Brandt =1; r_om is the characteristic radius of the Brandt profile in kpc
  double precision, parameter :: om0_kmskpc= (1.+(r0_kpc/r_om_kpc)**2)**(1./2)*V0_kms/r0_kpc  !om0 of Brandt profile in physical units
  double precision, parameter :: om_kmskpc= om0_kmskpc/(1. +(r_kpc/r_om_kpc)**2)**(1./2)
  double precision, parameter :: G_kmskpc=-om0_kmskpc*(r_kpc/r_om_kpc)**2/(1. +(r_kpc/r_om_kpc)**2)**(3./2)
!------------------------------------------------------------------------------------------------------------------
!  ALPHA EFFECT
  double precision, parameter :: C_alp= 1.
  double precision, parameter :: alp_k_mean_kms= C_alp*l_kpc**2/h_kpc*om_kmskpc  !maximum alpha_k in km/s
!------------------------------------------------------------------------------------------------------------------
!  DIMENSIONLESS PARAMETERS
  double precision, parameter :: l= l_kpc/h0_kpc*h0
  double precision, parameter :: h= h_kpc/h0_kpc*h0
  double precision, parameter :: tau=(l/h0)**2/3  !tau in units of td=h^2/etat
  double precision, parameter :: tautilde=1.  !Ratio of tau (eddy turnover time) to correlation time of the turbulence
  double precision, parameter :: U0=U0_kms/h0_km*td0_s  !Vertical mean velocity in normalized units h0/td0; 
  double precision, parameter :: G=G_kmskpc/km_kpc*td0_s
!  double precision, parameter :: R_omega= -18.75*25./26  !R_omega=Gh^2/etat
  double precision, parameter :: R_omega= G*h**2/etat  !R_omega=Gh^2/etat
  double precision :: tau_mta=0.
!  -18.75 corresponds to 25 km/s/kpc used in SSSB06. Value in Chamandy et al 2012 (paper I) is -18.75*25./26. 
!  double precision, parameter :: R_alpha= 1.5  
  double precision, parameter :: alp_k_mean= alp_k_mean_kms/h0_km*td0_s
  double precision, parameter :: R_alpha= alp_k_mean*h/etat
!  R_alpha=alp0*h/etat; Ampl of alpha effect, in normalized units h/td; 0.75 corresponds to 0.5 km/s used in SSSB06
  double precision, parameter :: Dyn=R_omega*R_alpha
  integer :: n_alp=  1  !Relevant only if Module Alp_sin eq 1. Sets the rate of variation of alpha in z
  double precision :: factau= 1.0 !1.  !Multiplication factor: tau_mta=factau*tau
  double precision :: Rm_inv= 0.  !1.e-5	;Inverse magnetic Reynolds number
!------------------------------------------------------------------------------------------------------------------
!  DYNAMICAL QUENCHING
  !  New Vishniac flux
  double precision, parameter :: fac_NV=0.  !fudge factor (f in Sharanya's notes) in the new Vishniac flux <b^2>=|fac|*Beq^2; 
  !  Set to 0. for no New Vishniac flux. May be +/-.
  !  Fickian diffusive flux
  double precision, parameter :: kappa=0.3  !Fickian diffusion coefficient
!------------------------------------------------------------------------------------------------------------------



end module parameters