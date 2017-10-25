
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RAFILT=1;        % Robert-Asselin time filter.
RA_COEFF=0.0001; % Robert-Asselin filter coefficient.
FORWARD=0;       % Forward time-step once per day.
FCLIP=0;         % Fourier clipping of small scales.
ARAKAWA=1;       % Energy/Enstrophy conserving Jacobian.
HOMOCLINIC=0;    % Singular solution on Homoclinic orbit.
Ubar=0;          % Mean Zonal flow (zero by default).
Hbar=10^4;       % Scale Height (m).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  DEFINE PARAMETERS WHICH CHANGE OFTEN  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify type of Initial Conditions.
  DAYLEN=1;          %  Forecast length in days.
  NX = 61;  NY = 21;  %  Spatial resolution
  DELTA_t = 1/12;      %  Timestep in hours.

     ICtype =  0;



  Ubar = 50;          % Mean zonal wind (m/s).
  Hbar=5500;          % Mean Height (m) for 500mb surface.


%%%% Section 1. Define constants and domain.   

daylen=DAYLEN;             %  Total time of integration (in days).
tlen = daylen*24*60*60;    %  Change to seconds.
Delta_t = DELTA_t;         %  Time-step (in hours).
Delta_t = Delta_t*60*60;   %  Change to seconds.
nt = tlen/Delta_t;         %  Number of time-steps.
t = (0:nt)*Delta_t;        %  time variable.
time = t/(24*60*60);       %  time in days (for plots).
nspd = (24*60*60)/Delta_t;  %  time steps per day.
nx = NX;  ,  ny = NY;      % Number of points in each direction
nxny = nx*ny;

fprintf('Grid size, nx=%i  ny=%i \n',nx,ny)
fprintf('Timesteps per day, nspd=%i \n',nspd)

%%%  Calculate the Coriolis parameter and beta parameter
Rearth = (4*10^7)/(2*pi);   %  Radius of the Earth (metres).
Omega = 2*pi/(24*60*60);
phi0=45*(pi/180);
fcor0 = 2*Omega*sin(phi0);
beta0 = 2*Omega*cos(phi0)/Rearth;

%%%  Calculate the Rossby Radius of Deformation.
grav = pi^2;     %    m s^(-2)
L_R = sqrt(grav*Hbar)/fcor0;  % Rossby Radius
F = 1/L_R^2;                   % Factor in BPV equation

%%%  Specify the domain size (adjusted below).
xlen = Rearth;             % East-West Length of the Domain.  
ylen = Rearth/3;           % North-South Length of the Domain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta_x = xlen/(nx);        %  Grid length in x-direction
Delta_y = ylen/(ny);        %  Grid length in y-direction
D_ratio = Delta_y/Delta_x; %  Grid length ratio

% Define the grid to have redundant rows east and north.
x = (0:nx)*Delta_x;
y = (0:ny)* Delta_y;

[XMESH, YMESH ] = meshgrid(x,y);
XX = XMESH'; YY=YMESH';
 
