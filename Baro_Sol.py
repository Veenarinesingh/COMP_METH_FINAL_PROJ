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


from numpy import linspace,pi,sin,cos,sqrt,meshgrid,size,shape,zeros,transpose,fix,mean,array
from random import seed,random
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from numpy.fft import fft2
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
NX = 61             # Set spatial resolution
NY = 21
DELTA_t = 1/12      # Timestep in hours
Ubar = 50           # Mean zonal wind (m/s).
Hbar=5500           # Mean Height (m) for 500mb surface.

# Part 1. Set constants and domain.

daylen=DAYLEN             #  Total time of integration (in days).
tlen = daylen*24*60*60    #  Change to seconds.
Delta_t = DELTA_t         #  Time-step (in hours).
Delta_t = Delta_t*60*60   #  Change to seconds.
nt = tlen/Delta_t         #  Number of time-steps.
t = linspace(0,nt,nt+1)*Delta_t        #  time variable.
time = t/(24*60*60)       #  time in days (for plots).
nspd = (24*60*60)/Delta_t #  time steps per day.

#Set grid
nx = NX                   # Number of points in each direction
ny = NY
nxny = nx*ny



print('Grid size, nx=',nx,' ny=',ny)
print('Timesteps per day',nspd)

# Calculate the Coriolis parameter and beta parameter
Rearth = (4*10**7)/(2*pi);  # Radius of the Earth (meters). ###!!!!!!!!!!!!!NOT EXACT TRY CHANGING TO EXACT
Omega = 2*pi/(24*60*60) # Angular velocity of the Earth
phi0=45*(pi/180)        # Latitude the calculation is centered on
fcor0 = 2*Omega*sin(phi0) #Coiolis parameter
beta0 = 2*Omega*cos(phi0)/Rearth;


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
XX = transpose(XMESH); YY= transpose(YMESH);

# Section 2. Define the Initial Fields.
#  w is the dependent variable. w_0 is the Initial field.
#  Pedlosky gives equation for streamfunction. For more
#  ergonomic scales we use the geopotential height.
#
#  Note that w does NOT include the part due to the
#  mean zonal flow. THis must be added if required.
#  w is periodic in both directions. Z is not.


#set seed for same initial conditions every time
seed(a=2)

#randomize seed for different initial conditions
seed()


Z_0 = zeros((nx+1,ny+1));
 # Set the incoming values.
Nwavex = 3
Nwavey = 3
Nwaves = (2*Nwavex+1)*(2*Nwavey+1)


for kwave in range(-Nwavex,Nwavex+1):
    for lwave in range(-Nwavey,Nwavey+1):
         Amplitude = (2000.0*2*(random()-0.5)) / Nwaves
         phase= 2*pi*(random()-0.5)
         term = cos(2*pi*(kwave*(XX/xlen)+lwave*(YY/ylen)+phase))
         Z_0 = Z_0 + Amplitude*term

#  Add a constant to give typical 500mb values.
Zplus_0 = Z_0 + Hbar

#Add in the zonal mean flow.
Ztotal_0 = Zplus_0 - (fcor0/grav)*Ubar*YY##################################################ASK BOOTH WHAT THIS IS

XM, YM = XX/(10**6), YY/(10**6)



plt.figure()
CS = plt.contourf(XM, YM, Z_0)
plt.colorbar()
plt.title('500 mbar Geopotential Perturbation (m)')


# Plot the field including the mean flow


vecwmin = Ztotal_0.min();
vecwmax = Ztotal_0.max();
vecwmean = (vecwmax+vecwmin)/2;
vecwspan = (vecwmax-vecwmin);
vecw = linspace(vecwmean-vecwspan,vecwmean+vecwspan,21);


plt.figure()
plt.contourf(XM, YM, Ztotal_0,vecw)
plt.colorbar()
plt.title('500 mbar Geopotential Height (m)')


#generate initial streamfunction



w_0=Z_0#*(grav/fcor0)  #  w_0 is perturbation height (excl. Hbar and Ubar-terms).

#  Add the mean zonal flow ( -Ubar*y ).
wtotal_0 = w_0 + Hbar - (fcor0/grav)*Ubar*YY;

XM = XX*10**(-6)
YM=YY*10**(-6);   #   Rescaled axes.


# Save specific Fourier components for plotting and analysis.
R=zeros((nx,ny))
a1=[]
a2=[]
a3=[]
a1star=[]
a2star=[]
a3star=[]

[XXin,YYin] = meshgrid(x[1:nx],y[1:ny]);
XXin=transpose(XXin)
YYin=transpose(YYin)

R=w_0[0:nx,0:ny]
W_hat = fft2(R)
W_hat_0 = W_hat



a1.append(2*(W_hat[nx-1,ny-1]) /nxny)
a2.append(2*(W_hat[nx-3,1]) /nxny)
a3.append(2*(W_hat[4,0])/nxny)

a1star.append(2*(W_hat[1,1]) / nxny)
a2star.append(2*(W_hat[3,ny-1]) / nxny)
a3star.append(2*(W_hat[nx-4,1]) / nxny)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter3D(XM, YM, w_0, cmap='Greens')
surf = ax.plot_surface(XM, YM, w_0, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.title('Initial Stream Function');


#% Section 3. Integrate the BPV Equation in time

#%%%   Time-stepping is by leapfrog method.
#%%%   First step is forward.
#%%%   Define Q = (Del^2 - F)w. The BPVE is
#%%%   (d/dt)Q + J(w,Del^2(w)) + beta*(d/dx)w
#%%%              + Ubar*(d/dx)Del^2(w) = 0.
#%%%   We approximate the time derivative by
#%%%     ( Q(n+1)-Q(n-1))/(2*Delta_t)
#%%%   and the remaining terms by centered differences:
#%%%      R(n) = - ( J(w,Del^2(w)) + beta*(d/dx)w
#%%%                + Ubar*(d/dx)Del^2(w))
#%%%   Then the value of Q at the new time (n+1)*Delta_t is:
#%%%      Q(n+1) =  Q(n-1) + 2*Delta_t * R(n)
#%%%
#   When we have Q(n+1), we have to solve a Helmholtz
#   equation to get w(n+1). Then the cycle is repeated.

# Define working arrays to have correct size.

dwdx = zeros((nx+1,ny+1))
dwdy = zeros((nx+1,ny+1))
gradsq = zeros((nx+1,ny+1))
d2wdx2 = zeros((nx+1,ny+1))
d2wdy2 = zeros((nx+1,ny+1))
laplac = zeros((nx+1,ny+1))
dlapdx = zeros((nx+1,ny+1))
dlapdy = zeros((nx+1,ny+1))
wdldx = zeros((nx+1,ny+1))
wdldy = zeros((nx+1,ny+1))
dwdxl = zeros((nx+1,ny+1))
dwdyl = zeros((nx+1,ny+1))
dwdldydx = zeros((nx+1,ny+1))
dwdldxdy = zeros((nx+1,ny+1))
ddwdxldy = zeros((nx+1,ny+1))
ddwdyldx = zeros((nx+1,ny+1))
Jac1 = zeros((nx+1,ny+1))
Jac2 = zeros((nx+1,ny+1))
Jac3 = zeros((nx+1,ny+1))
Jarakawa = zeros((nx+1,ny+1))



#%%%%% Start of main time-stepping loop %%%%%%%%%%

w = w_0;

Energy=[]
Enstrophy=[]
NLsize=[]
CFL_nonlin=[]
R=zeros((nx,ny))

for n in range(1,4):

#Take derivatives using finite-difference method, setting periodic boundary
#conditions in x and y

#x derivative

    dwdx[1:nx,0:ny+2]= (w[2:nx+1,0:ny+1]-w[0:nx-1,0:ny+1])/(2*Delta_x);
    dwdx[0,0:ny+1] = (w[1,0:ny+1]-w[nx-1,0:ny+1])/(2*Delta_x);
    dwdx[nx,0:ny+1] = dwdx[0,0:ny+1];

# y-derivative of w
    dwdy[0:nx+1,1:ny] = (w[0:nx+1,2:ny+1]-w[0:nx+1,0:ny-1])/(2*Delta_y);

    dwdy[0:nx+1,0] = (w[0:nx+1,1]-w[0:nx+1,ny-1])/(2*Delta_y);
    dwdy[0:nx+1,ny] = dwdy[0:nx+1,0];


# Square of the gradient of w
    gradsq = dwdx**2+dwdy**2;

# Second x-derivative of w
    d2wdx2[1:nx,0:ny+1] = (w[2:nx+1,0:ny+1]+w[0:nx-1,0:ny+1]-2*w[1:nx,0:ny+1])/(Delta_x**2)
    d2wdx2[0,0:ny+1] = (w[1,0:ny+1]+w[nx-1,0:ny+1]-2*w[0,0:ny+1])/(Delta_x**2)
    d2wdx2[nx,0:ny+1] = d2wdx2[0,0:ny+1]

# Second y-derivative of w
    d2wdy2[0:nx+1,1:ny] = (w[0:nx+1,2:ny+1]+w[0:nx+1,0:ny-1]-2*w[0:nx+1,1:ny])/(Delta_y**2)
    d2wdy2[0:nx+1,0] = (w[0:nx+1,1]+w[0:nx+1,ny-1]-2*w[0:nx+1,0])/(Delta_y**2)
    d2wdy2[0:nx+1,ny] = d2wdy2[0:nx+1,0]

    laplac = d2wdx2+d2wdy2;

# x-derivative of laplacian
    dlapdx[1:nx,0:ny+1] = (laplac[2:nx+1,0:ny+1]-laplac[0:nx-1,0:ny+1])/(2*Delta_x)
    dlapdx[0,0:ny+1] = (laplac[1,0:ny+1]-laplac[nx-1,0:ny+1])/(2*Delta_x)
    dlapdx[nx,0:ny+1] = dlapdx[0,0:ny+1]


#y-derivative of laplacian
    dlapdy[0:nx+1,1:ny] = (laplac[0:nx+1,2:ny+1]-laplac[0:nx+1,0:ny-1])/(2*Delta_y)
    dlapdy[0:nx+1,0] = (laplac[0:nx+1,1]-laplac[0:nx+1,ny-1])/(2*Delta_y)
    dlapdy[0:nx+1,ny] = dlapdy[0:nx+1,0]

    Jacobi = dwdx*dlapdy - dwdy*dlapdx

#Compute the Arakawa Jacobian.

    Jac1 = Jacobi;

    wdldx = w*dlapdx
    wdldy = w*dlapdy

    dwdldydx[1:nx,0:ny+1] = (wdldy[2:nx+1,0:ny+1]-wdldy[0:nx-1,0:ny+1])/(2*Delta_x);
    dwdldydx[0,0:ny+1] = (wdldy[1,0:ny+1]-wdldy[nx-1,0:ny+1])/(2*Delta_x);
    dwdldydx[nx,0:ny+1] = dwdldydx[0,0:ny+1];


    dwdldxdy[0:nx+1,1:ny] = (wdldx[0:nx+1,2:ny+1]-wdldx[0:nx+1,0:ny-1])/(2*Delta_y)
    dwdldxdy[0:nx+1,0] = (wdldx[0:nx+1,1]-wdldx[0:nx+1,ny-1])/(2*Delta_y)
    dwdldxdy[0:nx+1,ny] = dwdldxdy[0:nx+1,0]

    Jac2 = dwdldydx - dwdldxdy

    dwdxl = dwdx*laplac
    dwdyl = dwdy*laplac

    ddwdxldy[0:nx+1,1:ny] = (dwdxl[0:nx+1,2:ny+1]-dwdxl[0:nx+1,0:ny-1])/(2*Delta_y)
    ddwdxldy[0:nx+1,0] = (dwdxl[0:nx+1,1]-dwdxl[0:nx+1,ny-1])/(2*Delta_y)
    ddwdxldy[0:nx+1,ny] = ddwdxldy[0:nx+1,0]

    ddwdyldx[1:nx,0:ny+1] = (dwdyl[2:nx+1,0:ny+1]-dwdyl[0:nx-1,0:ny+1])/(2*Delta_x);
    ddwdyldx[0,0:ny+1] = (dwdyl[1,0:ny+1]-dwdyl[nx-1,0:ny+1])/(2*Delta_x)
    ddwdyldx[nx,0:ny+1] = ddwdyldx[0,0:ny+1]

    Jac3 = ddwdxldy - ddwdyldx

    Jarakawa = (1/3)*(Jac1+Jac2+Jac3)

#%%% Use the energy and enstrophy preserving Jacobian.

    Jacobi = Jarakawa;

#  Compute the function to be stepped forward.
    Q_n = laplac - F*w

    #  First time through the loop:
    if n==1:
        Dt = Delta_t/2;
        Q_nm1 = Q_n;

        rmeshvec=linspace(0,nx-1,nx)
        smeshvec=linspace(0,ny-1,ny)
        [rmesh,smesh] = meshgrid(rmeshvec,smeshvec)

        rr = transpose(rmesh)
        ss =transpose(smesh)
        C_rs = 2*(cos(2*pi*rr/nx)-1)/Delta_x**2+2*(cos(2*pi*ss/ny)-1)/Delta_y**2-F

    if FORWARD==1 and fix(n/nspd)*nspd==n:
        print('Forward timestep once per day \n')
        Dt = Delta_t/2
        Q_nm1 = Q_n


    #Calculate the energy and enstrophy integrals, the conserved quantities
    Rgsq=zeros((nx,ny))
    Rgsq[0:nx,0:ny] = gradsq[0:nx,0:ny]
    Rwsq=zeros((nx,ny))
    Rwsq[0:nx,0:ny] = w[0:nx,0:ny]**2
    Energy.append(0.5 * mean(mean(Rgsq+F*Rwsq)))
    Rgsq[0:nx,0:ny] = laplac[0:nx,0:ny]
    Rwsq[0:nx,0:ny] = w[0:nx,0:ny]
    Enstrophy.append(0.5 * mean(mean((Rgsq-F*Rwsq)**2)))

    #  Estimate the size of the nonlinear terms

    NonL = mean(mean(abs(Jacobi)))
    Beta = mean(mean(abs(beta0*dwdx)))
    NLsize.append(NonL/Beta)

    umax = abs(dwdy).max()
    vmax = abs(dwdx).max()
    maxx=array((umax,vmax))
    VMAX = maxx.max()
    CFL_nonlin.append(abs(VMAX*Delta_t/Delta_x))

    Q_np1 = Q_nm1 - (2*Dt)*((grav/fcor0)*Jacobi + beta0*dwdx + Ubar*dlapdx)

    if RAFILT==1:
        RA_coeff = RA_COEFF
        Q_n = Q_n + RA_coeff*(Q_nm1+Q_np1-2*Q_n)

#   Section 3.3: Solve the Helmholtz Equation (Del^2-F)w = R.

#  Compute the fft of the right hand side
#  (strip off additional row and column).
    R[0:nx,0:ny] = Q_np1[0:nx,0:ny]
    R_hat = fft2(R)

#Compute the transform of the solution
    W_hat = R_hat/C_rs

#  Compute the Wave Power of the components

    a1.append(2*(W_hat[nx-1,ny-1]) / nxny)
    a2.append(2*(W_hat[nx-3,1]) / nxny)
    a3.append(2*(W_hat[4,0]) / nxny)
    a1star.append(2*(W_hat[1,1]) / nxny)
    a2star.append(2*(W_hat[3,ny-1]) / nxny)
    a3star.append(2*(W_hat[nx-4,1]) / nxny)

    
#  Fourier filtering 
   if(FCLIP)
     nfilt = 10;
     nfilt=min([nfilt,(nx+1)/2-1,(ny+1)/2-1]);
     mask(1:nx,1:ny) = ones(nx,ny);
     nx1 = 2+nfilt; nx2 = nx-nfilt;
     ny1 = 2+nfilt; ny2 = ny-nfilt; 
     mask(nx1:nx2,1:ny)=zeros(nx2-nx1+1,ny);
     mask(1:nx,ny1:ny2)=zeros(nx,ny2-ny1+1);
     W_hat = W_hat.*mask;
   end


plt.show()
