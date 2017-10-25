#Veeshan Narinesingh, Computational Methods in Physics Fall '17
#Professor Ari Maller

#Final Project
#Numerical Solution to the Barotropic Vorticity Equation

#Note:This code has been adapted from Peter Lynch's MatLab code which can be
#accessed here https://maths.ucd.ie/met/msc/fezzik/MatLab/matlab.html


#    d              g                       d
#   -- (Del^2-F)w + - J(w,Del^2(w)) + beta* -- w = 0.
#   dt              f                      dx
#
#     With a mean zonal flow a term
#     Ubar*d((Del^2w)/dx is added.


from numpy import linspace,pi,sin,cos,sqrt,meshgrid,size

RAFILT=1        # Robert-Asselin time filter.
RA_COEFF=0.0001 # Robert-Asselin filter coefficient.
FORWARD=0       # Forward time-step once per day.
FCLIP=0         # Fourier clipping of small scales.
ARAKAWA=1       # Energy/Enstrophy conserving Jacobian.
HOMOCLINIC=0    # Singular solution on Homoclinic orbit.
Ubar=0          # Mean Zonal flow (zero by default).
Hbar=10^4       # Scale Height (m).


#Initial Conditions Specifications for a pseudo-real 500 mb flow

DAYLEN=1            # Forecast length in days.
NX = 61
NY = 21   # Spatial resolution
DELTA_t = 1/12      # Timestep in hours
Ubar = 50          # Mean zonal wind (m/s).
Hbar=5500          # Mean Height (m) for 500mb surface.

# Part 1. Set constants and domain.

daylen=DAYLEN             #  Total time of integration (in days).
tlen = daylen*24*60*60    #  Change to seconds.
Delta_t = DELTA_t         #  Time-step (in hours).
Delta_t = Delta_t*60*60   #  Change to seconds.
nt = tlen/Delta_t         #  Number of time-steps.
t = linspace(0,nt,nt+1)*Delta_t        #  time variable.
time = t/(24*60*60)       #  time in days (for plots).
nspd = (24*60*60)/Delta_t #  time steps per day.
nx = NX
ny = NY      # Number of points in each direction
nxny = nx*ny



print('Grid size, nx=',nx,' ny=',ny)
print('Timesteps per day',nspd)

# Calculate the Coriolis parameter and beta parameter
Rearth = (4*10**7)/(2*pi);  # Radius of the Earth (meters). ###!!!!!!!!!!!!!NOT EXACT TRY CHANGING TO EXACT
Omega = 2*pi/(24*60*60) # Angular velocity of the Earth
phi0=45*(pi/180)        # Latitude the calculation is centered on
fcor0 = 2*Omega*sin(phi0) #Coiolis parameter


# Calculate the Rossby Radius of Deformation.
grav = pi**2 # gravitational acceleration in ms^-2 !!!!!!!!!!!!!!!!!!!!!not exact try changing to exact
L_R = sqrt(grav*Hbar)/fcor0  # Rossby Radius
F = 1/(L_R**2)                 # Factor in B.V. Equation

# Specify the domain size (adjusted below).
xlen = Rearth             # East-West Length of the Domain.
ylen = Rearth/3           # North-South Length of the Domain.

Delta_x = xlen/(nx)        #  Grid length in x-direction
Delta_y = ylen/(ny)        #  Grid length in y-direction
D_ratio = Delta_y/Delta_x #  Grid length ratio

# Define the grid to have redundant rows east and north.


x = linspace(0,nx,nx+1)*Delta_x
y = linspace(0,ny,ny+1)*Delta_y

(XMESH, YMESH) = meshgrid(x,y);
XX = XMESH; YY=YMESH;
