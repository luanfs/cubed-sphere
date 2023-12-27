module constants
!======================================================
!
! This module contains all the constants needed for the other modules
!
! Based on module 'constants' from iModel
! https://github.com/pedrospeixoto/iModel
!======================================================
implicit none
save 
public

!---------------------------------------------------
!Kind attributions
!---------------------------------------------------

integer, parameter :: i2  = selected_real_kind(2,20)
integer, parameter :: i4  = selected_real_kind(6,20)
integer, parameter :: i8  = selected_real_kind(14,40)

integer, parameter :: r4  = selected_real_kind(6,37)  
integer, parameter :: r8  = selected_real_kind(12,100)

!---------------------------------------------------
! General Parameters
!---------------------------------------------------

!Pi 
real(kind=8), parameter :: pi   = dacos(-1.d0)
real(kind=8), parameter :: pi2  = 2.d0*pi
real(kind=8), parameter :: pio2 = pi/2.d0
real(kind=8), parameter :: piby2 = pi*0.5d0
real(kind=8), parameter :: pio4 = pi*0.25d0

!Degrees to radians coversion (multiply to obtain conversion)
real(kind=8), parameter :: deg2rad = pi / 180.d0

!Radians to Degrees coversion (multiply to obtain conversion)
real(kind=8), parameter :: rad2deg = 1.d0/deg2rad

!Very small real number (aprox 1.10e-7)
!real(kind=8), parameter :: eps = epsilon(1.)
real(kind=8), parameter :: eps = epsilon(1.)

!Very very small real number (aprox 1.10e-16)
real(kind=8), parameter :: eps2 = epsilon(pi)

!---------------------------------------------------
! Physical Parameters
!---------------------------------------------------

! Earth mean radius (meters)
real(kind=8), parameter :: erad     = 6371200.d0
real(kind=8), parameter :: rearth   = 6371200.d0
real(kind=8), parameter :: eradi    = 1.d0/erad
real(kind=8), parameter :: unitspharea    = 4.d0*pi

! Half length of cube edge
real(kind=8), parameter :: acube = erad/dsqrt(3.d0)

!Gravitational accerlaration of the Earth (m/s^2)
real(kind=8), parameter  :: grav    = 9.8066499999999994d0
real(kind=8), parameter  :: gravity = grav
real(kind=8), parameter  :: gravi   = 1.d0/grav
real(kind=8), parameter  :: gravo2   = grav*0.5d0

! Angular velocity of the Earth (rot/s)
real (kind=8), parameter :: omega   = 7.2921d0*10.d0**(-5)
real (kind=8), parameter :: rotatn   = 7.2921d0*10.d0**(-5)

!Days to seconds
real (kind=8), parameter :: day2sec = 86400.d0
real (kind=8), parameter :: sec2day = 1.d0/86400.d0

! Dry air gas constant [J/(kg K)]
real(kind=8), parameter :: rdry     = 287.   

! Dry air spec heat at const P [J/(kg K)]
real(kind=8), parameter :: cp       = 1004.  

! Dry air spec heat at const vol [J/(kg K)]
real(kind=8), parameter :: cv       = 717.              

! Water vapor gas constant [J/(kg K)]
real(kind=8), parameter :: rvap     = 461.               

! Reference pressure [Pa]
real(kind=8), parameter :: p00      = 1.e5            

! 0 Celsius temperature [K]
real(kind=8), parameter :: t00      = 273.15         

!---------------------------------------------------
! Directories
!--------------------------------------------------
!Output data directory
character (len=60), parameter::  datadir= "data/"
!Grid data directory
character (len=60), parameter::  griddir= "grid/"
!Directory for ploting thing using GMT
character (len=60), parameter::  gmtdir = "gmt/"
!Sources directory
character (len=60), parameter::  srcdir = "src/"
!Parameter files  directory
character (len=60), parameter::  pardir = "par/"
!Reference solutions directory - data will be read from here
character (len=60), parameter::  refdir = "ref/"

!---------------------------------------------------
! Global constant variables
!---------------------------------------------------
!Flag for verbose output
logical :: showonscreen=.false.

!Simulation to be done
integer (i4) :: simulcase

!Number of faces on a cube
integer (i4) :: nbfaces = 6

!Interior grid indexes
integer (i4) :: i0, iend
integer (i4) :: j0, jend

! Number of ghost cells in a direction
integer(i4) :: nghost
integer(i4) :: n0, nend
integer(i4) :: hs

!Lat/lon grid size (for ploting)
integer (i4) :: n_lon = 1440
integer (i4) :: n_lat = 720

end module constants
