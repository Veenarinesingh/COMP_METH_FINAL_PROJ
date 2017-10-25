
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                     %%%
%%%  Solve the Barotropic Potential Vorticity Equation  %%%
%%%                                                     %%%
%%%    d              g                       d         %%%
%%%   -- (Del^2-F)w + - J(w,Del^2(w)) + beta* -- w = 0. %%%
%%%   dt              f                      dx         %%%
%%%                                                     %%%
%%%     With a mean zonal flow a term                   %%%
%%%     Ubar*d((Del^2w)/dx is added.                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                     %%%
%%%   The Barotropic PVE is approximated by centered    %%%
%%%   finite differences. The domain has nx*ny points.  %%%
%%%                                                     %%%
%%%   The boundary conditions are periodic:             %%%
%%%     w(nx+1,:) = w(1,:) and w(:,ny+1) = w(:,1)       %%%
%%%                                                     %%%
%%%   Arrays are of size (nx+1)*(ny+1) so that the      %%%
%%%   easternmost column and northernmost row are       %%%
%%%   redundant. This is for computational convenience. %%%
%%%                                                     %%%
%%%   For the Fourier Transforms, the redundant row     %%%
%%%   and column must be removed (this is expected      %%%
%%%   by the functions fft2 and ifft2). The redundant   %%%
%%%   information is filled in after the inverse        %%%
%%%   transformation.                                   %%%
%%%                                                     %%%
%%%   The time scheme is the leapfrog method.           %%%
%%%   (the first step is a forward step).               %%%
%%%                                                     %%%
%%%   The quantity to be stepped forward is             %%%
%%%           Q = Del^2(w)-F*w.                         %%%
%%%   Thus, after each time-step, it is necessary       %%%
%%%   to solve a Helmholtz equation to get w.           %%%
%%%   This is done by the Fourier Transform method.     %%%
%%%   (See Numerical Recipes, Chapter 19 (Sec. 19.4).   %%%
%%%                                                     %%%
%%%  Note: Stability for linear equation is determined  %%%
%%%  by Rossby wave phase speeds. The permissible time  %%%
%%%  step is quite long. For larger amplitudes, when    %%%
%%%  the nonlinear terms become bigger, so does the     %%%
%%%  wind speed, and it is the size of the gradient of  %%%
%%%  the stream function which determines the timestep. %%%
%%%                                                     %%%
%%%  Fourier clipping of smaller scales is possible.    %%%
%%%  With the Arakawa energy and enstrophy conserving   %%%
%%%  Jacobian, this is not required. Also, a forward    %%%
%%%  step once per day may be used, but is not required.%%%
%%%  The Robert-Asselin time filter is also available.  %%%
%%%                                                     %%%
%%%   A variety of initial conditions may be specified  %%%
%%%   determined by the parameter ICtype. Additional    %%%
%%%   types of IC may easily be added.                  %%%
%%%                                                     %%%
%%%   1 March, 2005: Modification made to allow for     %%%
%%%                  a mean zonal flow Ubar.            %%%
%%%                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                     %%%
%%% Author: Peter Lynch, Met Eireann, Glasnevin, Dublin %%%
%%% Email:  Peter.Lynch@met.ie                          %%%
%%%                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  diary DIARY %  SAVE OUTPUT IN DIARY IF REQUIRED.

clear    %    Clear the memory.
clf      %    clear the display.

% define the display colours.
whitebg([1 1 0.50]);  %  Beige background.
colormap('jet');    

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

     ICtype =  0
     ;   %%
     %ICtype =  2;   %%
     
     %ICtype =  7;   %%
     %ICtype =  3;   %%

%%%%%%%%%%%%%%%%% ICtype==0: Pseudo-real 500 mb flow.
if(ICtype==0)
  Ubar = 50;          % Mean zonal wind (m/s).
  Hbar=5500;          % Mean Height (m) for 500mb surface.
end
%%%%%%%%%%%%%%%%% Pure wave 1.
if(ICtype==1)
  DAYLEN=100;  % time for about two cycles.
  NX = 61;  NY = 21;  %  Spatial resolution
%%%%%%%%%%%%%%%%% Pure wave 2.
  AMP3 = 0.00; AMP2 = 0.00; AMP1 = 1.00;
  AMPSCALE = 1.0*10^5;
  PSI1 = 0; PSI2 = 0; PSI3 = 0;
  DELTA_t = 12;
end
if(ICtype==2)
  DAYLEN=100;  % time for about 3.5 cycles.
  NX = 61;  NY = 21;  %  Spatial resolution
  AMP3 = 0.00; AMP2 = 1.00; AMP1 = 0.00;
  AMPSCALE = 1.0*10^5;
  PSI1 = 0; PSI2 = 0; PSI3 = 0;
  DELTA_t = 12;
end
%%%%%%%%%%%%%%%%% Pure wave 3.
if(ICtype==3)
  DAYLEN=100;  % time for about 5.5 cycles.
  NX = 61;  NY = 21;  %  Spatial resolution
  AMP3 = 1.00; AMP2 = 0.00; AMP1 = 0.00;
  AMPSCALE = 1.0*10^5;
  PSI1 = 0; PSI2 = 0; PSI3 = 0;
  DELTA_t = 12;
end
%%%%%%%%%%%%%%%%% Large amplitude triad: 
%%%  Most energy leaves the triad.
if(ICtype==4)
  DAYLEN=100;    % 
  NX = 61;  NY = 21;  %  Spatial resolution
  AMP3 = 1.00; AMP2 = 0.01; AMP1 = 0.00;
  AMPSCALE = 10.0*10^5;
  PSI1 = 0; PSI2 = 0; PSI3 = 0;
  DELTA_t = 1;
end
%%%%%%%%%%%%%%%%% Homoclinic orbit. 
%%%%%%%%% Doesn't work. Hard to simulate unstable mode.
if(ICtype==5)
  HOMOCLINIC = 1;
  DAYLEN = 400;    
  NX = 61;  NY = 21;  %  Spatial resolution
  AMP3 =0; AMP2 = 1; AMP1 = 1;
  AMPSCALE = 0.25;
  PSI1 = 0; PSI2 = 0; PSI3 = pi;
  DELTA_t = 12;
end
%%%%%%%%%%%%%%%%% Precessing triad. <<<<<<< Put Correct values.
% See program ST_Link for definition of ICs.
if(ICtype==6)
  DAYLEN=6;    % time for three pulses with scale 0.5
  NX = 61;  NY = 21;  %  Spatial resolution
  AMP3 = 1.00; AMP2 = 0.01; AMP1 = 0.00;
  AMPSCALE = 1.0*10^6;
  PSI1 = 0; PSI2 = 0; PSI3 = 0;
  DELTA_t = 3;
end
%%%%%%%%%%%%%%%%% PRECESSION: Starfish pattern. 
% See program ST_Link for definition of ICs.
if(ICtype==7)
  DAYLEN=4800;
     DAYLEN=800;  % Small value for testing (one pulsation).
  NX = 61;  NY = 21;  %  Spatial resolution
  AMPS = [ 0.26847915 0.15402274 1.15552765 ]*1.0e+011;
  PHAS = [ 0.00000    0.00000    1.57079633 ];
    AMP1 = AMPS(1); AMP2 = AMPS(2); AMP3 = AMPS(3);
    PSI1 = PHAS(1); PSI2 = PHAS(2); PSI3 = PHAS(3);
  AMPSCALE = 3.4*10^(-7);
  DELTA_t = 12;
end
%%%%%%%%%%%%%%%%% 2D Oscillations. 
% See program ST_Link for definition of ICs.
if(ICtype==8)
  DAYLEN = 2400;
  NX = 46;  NY = 16;  %  Spatial resolution
  AMPS = [0.21094790   0.21178127   1.15552765]*1.0e+011;
  PHAS = [ 0.00000      0.00000     1.57079633];
    AMP1 = AMPS(1); AMP2 = AMPS(2); AMP3 = AMPS(3);
    PSI1 = PHAS(1); PSI2 = PHAS(2); PSI3 = PHAS(3);
  AMPSCALE = 3.4*10^(-7);
  DELTA_t = 12;
end
%%%%%%%%%%%%%%%%% Large amplitude triads: Predictability. 
%%%  Most energy leaves the triad.
%%%  Save fields and plot for regular and flipped phases.
%%%  (save w_0 wq4 on first run and load on second).
if(ICtype==9)
  DAYLEN=4;
  NX = 61;  NY = 21;  %  Spatial resolution
  AMP3 = 1.00; AMP2 = 0.01; AMP1 = 0.01;
  AMPSCALE = 60*10^(5);
  firstrun=1; PSI1 =  0; PSI2 =  0; PSI3 = 0;  %  Regular phases.
  % firstrun=0; PSI1 = pi; PSI2 = pi; PSI3 = 0;  %  Flipped phases.
  DELTA_t = 1;
end

%%%%%%%%%%%%%%%%% Period doubling. 
% Vary gscale. 1.0: Period one. 2.5: Period 2. 5.0: Chaos.
%%%   Do: "save PD100 time NRG1 NRG2 NRG3 a1 a2 a3", etc.
%%%   to save the relevant values in Pdouble for plotting.
if(ICtype==10)
   DAYLEN = 10000; 
   %  DAYLEN = 2000;    %  Small value for testing.
  NX = 31;  NY = 11;  %  Spatial resolution
  AMPS = [ 0.26847915 0.15402274 1.15552765 ]*1.0e+011;
  PHAS = [ 0.00000    0.00000    1.57079633 ];
    AMP1 = AMPS(1); AMP2 = AMPS(2); AMP3 = AMPS(3);
    PSI1 = PHAS(1); PSI2 = PHAS(2); PSI3 = PHAS(3);
  AMPSCALE = 3.5*10^(-7);
  AMPSCALE = 3.4*10^(-7);
  DELTA_t = 12;
  RA_COEFF=0.0002;  % Set RA_filter coefficient.
  % define the scaling factor gscale and apply it.
      gscale = 2.5;
      AMPSCALE = gscale*AMPSCALE;
      DAYLEN = round(DAYLEN/gscale);
      RA_COEFF = gscale*RA_COEFF;
end

%%%%%%%%%%%%%%%%% Period doubling. TESTING in August, 2004.
% Vary gscale. 1.0: Period one. 2.5: Period 2. 5.0: Chaos.
%%%   Do: "save PD100 time NRG1 NRG2 NRG3 a1 a2 a3", etc.
%%%   to save the relevant values in Pdouble for plotting.
if(ICtype==11)
   DAYLEN = 10000; 
      %       DAYLEN = 2000;    %  New value for testing.
  NX = 31;  NY = 11;  %  Spatial resolution
  AMPS = [ 0.26847915 0.15402274 1.15552765 ]*1.0e+011;
  PHAS = [ 0.00000    0.00000    1.57079633 ];
    AMP1 = AMPS(1); AMP2 = AMPS(2); AMP3 = AMPS(3);
    PSI1 = PHAS(1); PSI2 = PHAS(2); PSI3 = PHAS(3);
  AMPSCALE = 3.5*10^(-7);
  AMPSCALE = 3.4*10^(-7);
  DELTA_t = 12;
  RA_COEFF=0.0002;  % Set RA_filter coefficient.
       %   DAYLEN = 1000;
          PHAS = [ 0.00    0.00   pi/2 ];
          PSI1 = PHAS(1); PSI2 = PHAS(2); PSI3 = PHAS(3);
  % define the scaling factor gscale and apply it.
  gscale = 2.5;   %
     %  gscale = 1.0;
      AMPSCALE = gscale*AMPSCALE;
      DAYLEN = round(DAYLEN/gscale);
      RA_COEFF = gscale*RA_COEFF;
        
   end
% gscale = 2.65: Periot two; peaks tending towards equality
% gscale = 2.70: Periot two; peaks tending towards equality
% gscale = 2.80: Period two; peaks tending towards equality
% gscale = 2.90: Period two: switch over (also seen above).
% gscale = 3.00: BSBSBSSBSBMMBS-SBSB
% gscale = 3.10: BSBSBBSBSBSS
% gscale = 3.20: BS BS SB SS BS SS ...
% gscale = 3.30: Big and small peaks; Difficult to interpret.

% gscale = 2.50: period two.
% gscale = 2.60: period two.
% gscale = 3.50: Hint of three peaks each biger than the previous one.
% gscale = 3.60: Difficult to interpret.
% gscale - 4.00: BSS BSS BBS BS...
% gscale = 4.50: Difficult to interpret: many peaks similar size.
% gscale = 5.00: Looks fairly chaotic.
% gscale = 10.00: Clearly chaotic.


%%%%%%%%%%%%%%%%%
if(ICtype>11)
  fprintf('ICtype > 10 not permitted. \n'); return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  END DEFINE PARAMETERS WHICH CHANGE OFTEN  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
x = (0:nx)* Delta_x;
y = (0:ny)* Delta_y;

[XMESH, YMESH ] = meshgrid(x,y);
XX = XMESH'; YY=YMESH';
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Section 2. Define the Initial Fields.
%  w is the dependent variable. w_0 is the Initial field.
%  Pedlosky gives equation for streamfunction. For more
%  ergonomic scales we use the geopotential height.    
%
%  Note that w does NOT include the part due to the
%  mean zonal flow. THis must be added if required.
%  w is periodic in both directions. Z is not.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ICtype==0)    %% Pseudo-real 500mb geopotential field.

   % Define the field size.
   Z_0 = zeros(nx+1,ny+1);
   
   %    Set the seed of the random number generator
   %    to its initial value (same results each run).
            rand('state',0)
   %    Set the seed of the random number generator
   %    to yield different results each run.
        rand('state',sum(100*clock))
   
   % Set the incoming values.
   Nwavex = 3;
   Nwavey = 3;
   Nwaves = (2*Nwavex+1)*(2*Nwavey+1);
   for kwave=-Nwavex:Nwavex
   for lwave=-Nwavey:Nwavey
      Amplitude = (2000.0*2*(rand(1)-0.5)) / Nwaves;
      phase= 2*pi*(rand(1)-0.5);
      term = cos(2*pi*(kwave*(XX/xlen)+lwave*(YY/ylen)+phase));
      Z_0 = Z_0 + Amplitude*term;
   end
   end
   
   %%%%%%%%%%%%%%%%
   %%%   Add a constant to give typical 500mb values.
   Zplus_0 = Z_0 + Hbar; 

   % Add in the zonal mean flow.
   Ztotal_0 = Zplus_0 - (fcor0/grav)*Ubar*YY;
   
   % Plot the field
   
   figure
   
   XM = XX/10^6; YM = YY/10^6;
   contourf(XM,YM,Z_0)
   
   title('PERTURBATION GEOPOTENTIAL HEIGHT')
   colorbar
   
  
   
   fprintf('Press RETURN to continue \n');  pause
   
   figure
   
else  %%%%%%%%%%%%   More exotic Initial Conditions.

%%%   First, define the wave parameters:
%  Third wave has zonal wavenumber four. 
%  Second wave has zonal wavenumber three. 
%  First wave has zonal wavenumber one. 

   k3 = (2*pi/xlen)*4;   
   k2  = -0.75*k3;
   k1  = -0.25*k3;

   Lenx(3) = 2*pi/k3;
   Lenx(2) = 2*pi/k2;
   Lenx(1) = 2*pi/k1;

%  Solve the equations for the Resonant Triad.
   Aq = k3;
   Bq = k3*(k1^2+k2^2+2*F)+(k1+k2)*(k3^2+F);
   Cq = k3*(k2^2+F)*(k1^2+F)+k2*(k1^2+F)*(k3^2+F)+k1*(k3^2+F)*(k2^2+F);
   
   lsq = ( -Bq+sqrt(Bq^2-4*Aq*Cq) ) / (2*Aq);
   l = sqrt(lsq);
   l3 = 0; l2 = l; l1 = -l;
   
   Leny(1) = 2*pi/l1;
   Leny(2) = 2*pi/l2;
   Leny(3) = 2*pi/(l3+eps);  %  l3=0. Avoid zero division.

%  Calculate the interaction coefficients
   K1SQ = k1^2+l1^2; 
   K2SQ = k2^2+l2^2; 
   K3SQ = k3^2+l3^2;

   B1 = 0.5*(K2SQ-K3SQ)*(k2*l3-k3*l2); mu1 = sqrt(abs(B1/(K1SQ+F)));
   B2 = 0.5*(K3SQ-K1SQ)*(k3*l1-k1*l3); mu2 = sqrt(abs(B2/(K2SQ+F)));
   B3 = 0.5*(K1SQ-K2SQ)*(k1*l2-k2*l1); mu3 = sqrt(abs(B3/(K3SQ+F)));
   SumB = (B1+B2+B3)/(abs(B1)+abs(B2)+abs(B3));

   kwave  = [ k1      k2      k3      ];
   lwave  = [ l1      l2      l3      ];
   KSQ    = [ K1SQ    K2SQ    K3SQ    ];
   MU     = [ mu1     mu2     mu3     ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Adjust domain size for periodicity in y, i.e., to
%  ensure an exact number of waves in domain.
   ywvs = abs(ylen/Leny(1));
   yadj = ywvs/round(ywvs);
   ylen = ylen/yadj;
   aspect = ylen/xlen;     % Aspect ratio of domain

   Xwaves = abs( xlen ./ Lenx' );
   Ywaves = abs( ylen ./ Leny' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %  Get the Frequencies and phase-speeds.
   sigma1 = - beta0*k1/(K1SQ+F); cphase1 = sigma1/k1;
   sigma2 = - beta0*k2/(K2SQ+F); cphase2 = sigma2/k2;
   sigma3 = - beta0*k3/(K3SQ+F); cphase3 = sigma3/k3;
   
   tau1 = Lenx(1)/cphase1; 
   tau2 = Lenx(2)/cphase2; 
   tau3 = Lenx(3)/cphase3; 
   period1 = abs(tau1/(24*60*60)); %  period in days.
   period2 = abs(tau2/(24*60*60)); %  period in days.
   period3 = abs(tau3/(24*60*60)); %  period in days.
      
   sigma  = [ sigma1  sigma2  sigma3  ];
   cphase = [ cphase1 cphase2 cphase3 ];
   period = [ period1 period2 period3 ];

% Discrete form of the frequencies, phase-speeds and periods.

   C_kl_1 = -2*((cos(k1*Delta_x)-1)/Delta_x^2+(cos(l1*Delta_y)-1)/Delta_y^2);
   C_kl_2 = -2*((cos(k2*Delta_x)-1)/Delta_x^2+(cos(l2*Delta_y)-1)/Delta_y^2);
   C_kl_3 = -2*((cos(k3*Delta_x)-1)/Delta_x^2+(cos(l3*Delta_y)-1)/Delta_y^2);
   sins1 = -(beta0*k1/(C_kl_1+F))*(sin(k1*Delta_x)/(k1*Delta_x));
   sins2 = -(beta0*k2/(C_kl_2+F))*(sin(k2*Delta_x)/(k2*Delta_x));
   sins3 = -(beta0*k3/(C_kl_3+F))*(sin(k3*Delta_x)/(k3*Delta_x));
   Sigma1 = asin(sins1*Delta_t)/Delta_t;
   Sigma2 = asin(sins2*Delta_t)/Delta_t;
   Sigma3 = asin(sins3*Delta_t)/Delta_t;

   Cphase1 = Sigma1/k1; Tau1 = Lenx(1)/Cphase1; 
   Cphase2 = Sigma2/k2; Tau2 = Lenx(2)/Cphase2; 
   Cphase3 = Sigma3/k3; Tau3 = Lenx(3)/Cphase3; 
   Period1 = abs(Tau1/(24*60*60)); %  period in days.
   Period2 = abs(Tau2/(24*60*60)); %  period in days.
   Period3 = abs(Tau3/(24*60*60)); %  period in days.

   Sigma  = [ Sigma1  Sigma2  Sigma3  ];
   Cphase = [ Cphase1 Cphase2 Cphase3 ];
   Period = [ Period1 Period2 Period3 ];

% Replace Discrete form of the frequencies.

   sigma1c= sigma1; sigma2c= sigma2; sigma3c= sigma3;
   sigma1 = Sigma1; sigma2 = Sigma2; sigma3 = Sigma3;

   sigmac = [ sigma1c sigma2c sigma3c ];
   Sigma  = [ Sigma1  Sigma2  Sigma3  ]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the resonance conditions apply.
   sum_k = sum(kwave);
   sum_l = sum(lwave);
   sum_sigma = sum(sigma);

   cmax=max(abs(cphase));
   CFL_rossby = abs(cmax*Delta_t/Delta_x);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Now set up the components and combine them.
   
%%%  w_0 = Amp * cos(k*x+l*y-sigma*t+psi)

% First, Scale coefficients from streamfunction to height:
   AMP = [ AMP1, AMP2, AMP3 ]; %  Before conversion.
   Amp1 = (fcor0/grav)*AMP1;
   Amp2 = (fcor0/grav)*AMP2; 
   Amp3 = (fcor0/grav)*AMP3; 

%%%   Now scale the component coefficients:
   Amp = [ Amp1, Amp2, Amp3 ]; %  Before re-scaling.
   gamma = AMPSCALE;
   Amp1 = gamma*Amp1; Amp2 = gamma*Amp2; Amp3 = gamma*Amp3;
   Amp = [ Amp1, Amp2, Amp3 ];  %  After re-scaling.

%    Special procedures for the Homoclinic orbit only.
     if(HOMOCLINIC) 
       Amp1 = mu1*Amp1;
       Amp2 = mu2*Amp2;
       Amp3 = mu3*Amp3;
       Askale =  gamma / (max(abs([Amp1,Amp2,Amp3]))); 
       Amp1 = Askale*Amp1; Amp2 = Askale*Amp2; Amp3 = Askale*Amp3;
       Amp = [ Amp1, Amp2, Amp3 ]  % Homoclinic orbit.
       eta = sqrt( 0.5*(Amp1^2+Amp2^2+2*Amp3^2) ); eta
     end

   psi1 = PSI1; psi2 = PSI2; psi3 = PSI3;
   Psi    = [ psi1    psi2    psi3    ];
 
%  Construct the initial field from the components.
   w_01 = cos(k1*XX+l1*YY+psi1);  %  mesh(XX,YY,w_01); %pause
   w_02 = cos(k2*XX+l2*YY+psi2);  %  mesh(XX,YY,w_02); %pause
   w_03 = cos(k3*XX+l3*YY+psi3);  %  mesh(XX,YY,w_03); %pause
   w_0 = Amp1*w_01+Amp2*w_02+Amp3*w_03;

   Z_0 = w_0;  Ztotal_0 = w_0;
   XM = XX/10^6; YM = YY/10^6;

end   %    End of construction of initial streamfunction.

%%%%%%%%%%%
% Plot the field including the mean flow
vecwmin = min(min(Ztotal_0));
vecwmax = max(max(Ztotal_0));
vecwmean = (vecwmax+vecwmin)/2;
vecwspan = (vecwmax-vecwmin);
vecw = linspace(vecwmean-vecwspan,vecwmean+vecwspan,21);
contourf(XM,YM,Ztotal_0,vecw)
title('TOTAL GEOPOTENTIAL HEIGHT')
colorbar
fprintf('Press RETURN to continue \n');  pause
figure

w_0 = Z_0;  %  w_0 is perturbation height (excl. Hbar and Ubar-terms).

%  Add the mean zonal flow ( -Ubar*y ).
wtotal_0 = w_0 + Hbar - (fcor0/grav)*Ubar*YY;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   XM = XX*10^(-6); YM=YY*10^(-6);   %   Rescaled axes.
%     mesh(XM,YM,w_0); % pause

   plotnn = 50;
   if(plotnn==50 & ICtype>0)
     zmax=1; zmin=-zmax;
     subplot(2,2,1); surf(XM,YM,w_01); title('WAVE 1');
       axis([0,7,0,3,zmin,zmax]); shading interp;
     subplot(2,2,2); surf(XM,YM,w_02); title('WAVE 2');
       axis([0,7,0,3,zmin,zmax]); shading interp;
     subplot(2,2,3); surf(XM,YM,w_03); title('WAVE 3');
       axis([0,7,0,3,zmin,zmax]); shading interp
       xlabel('X'); ylabel('Y'); 
     w0max = max(max(w_0));
     subplot(2,2,4); surf(XM,YM,w_0/w0max); title('INITIAL FIELD'); 
       axis([0,7,0,3,zmin,zmax]); shading interp;
       xlabel('X'); ylabel('Y'); 
     fprintf('Press RETURN to continue \n');  pause
     %print -depsc trinit.eps;
   end

   plotnn = 60;
   if(plotnn==60 & ICtype>0)
     subplot(2,2,1); pcolor(XM,YM,w_01); title('WAVE 1');
       axis('off'); shading interp;
     subplot(2,2,2); pcolor(XM,YM,w_02); title('WAVE 2');
       axis('off'); shading interp;
     subplot(2,2,3); pcolor(XM,YM,w_03); title('WAVE 3');
       axis('off'); shading interp
       xlabel('X'); ylabel('Y'); 
     w0max = max(max(w_0));
     subplot(2,2,4); pcolor(XM,YM,w_0/w0max); 
       title('INITIAL FIELD'); 
       axis('off'); shading interp;
       xlabel('X'); ylabel('Y'); 
     fprintf('Press RETURN to continue \n');  pause
     %print -depsc tpinit.eps;
   end

% Save specific Fourier components for plotting and analysis.
   [XXin,YYin] = meshgrid(x(1:nx),y(1:ny));
   XXin = XXin'; YYin=YYin';
   R = w_0(1:nx,1:ny);
   W_hat = (fft2(R));
   W_hat_0 = W_hat;  

   a1(1) = 2*(W_hat(1+nx-1,1+ny-1)) / nxny; 
   a2(1) = 2*(W_hat(1+nx-3,1+1   )) / nxny;
   a3(1) = 2*(W_hat(1+4   ,1+0   )) / nxny;

   a1star(1) = 2*(W_hat(1+1   ,1+1   )) / nxny; 
   a2star(1) = 2*(W_hat(1+3   ,1+ny-1)) / nxny;
   a3star(1) = 2*(W_hat(1+nx-4,1+0   )) / nxny;

%  Plot the Initial Conditions
subplot(1,1,1);
mesh(XM,YM,w_0);
xlabel('x'); ylabel('y'); zlabel('w');
title('Initial Stream Function'); 
drawnow
fprintf('Press RETURN to continue \n');

pause


figure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Section 3. Integrate the BPV Equation in time

%%%   Time-stepping is by leapfrog method.
%%%   First step is forward.
%%%   Define Q = (Del^2 - F)w. The BPVE is
%%%   (d/dt)Q + J(w,Del^2(w)) + beta*(d/dx)w 
%%%              + Ubar*(d/dx)Del^2(w) = 0.
%%%   We approximate the time derivative by
%%%     ( Q(n+1)-Q(n-1))/(2*Delta_t)    
%%%   and the remaining terms by centered differences:
%%%      R(n) = - ( J(w,Del^2(w)) + beta*(d/dx)w 
%%%                + Ubar*(d/dx)Del^2(w)) 
%%%   Then the value of Q at the new time (n+1)*Delta_t is:
%%%      Q(n+1) =  Q(n-1) + 2*Delta_t * R(n)   
%%%
%%%   When we have Q(n+1), we have to solve a Helmholtz
%%%   equation to get w(n+1). Then the cycle is repeated.

% Define working arrays to have correct size.

dwdx = zeros(nx+1,ny+1);
dwdy = zeros(nx+1,ny+1);
gradsq = zeros(nx+1,ny+1);
d2wdx2 = zeros(nx+1,ny+1);
d2wdy2 = zeros(nx+1,ny+1);
laplac = zeros(nx+1,ny+1);
dlapdx = zeros(nx+1,ny+1);
dlapdy = zeros(nx+1,ny+1);
wdldx = zeros(nx+1,ny+1);
wdldy = zeros(nx+1,ny+1);
dwdxl = zeros(nx+1,ny+1);
dwdyl = zeros(nx+1,ny+1);
dwdldydx = zeros(nx+1,ny+1);
dwdldxdy = zeros(nx+1,ny+1);
ddwdxldy = zeros(nx+1,ny+1);
ddwdyldx = zeros(nx+1,ny+1);
Jac1 = zeros(nx+1,ny+1);
Jac2 = zeros(nx+1,ny+1);
Jac3 = zeros(nx+1,ny+1);
Jarakawa = zeros(nx+1,ny+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Start of main time-stepping loop %%%%%%%%%%

w = w_0;

fprintf('  Total Number of Steps: %i \n', nt)

for n=1:nt
   
   if(fix(n/10)*10==n | n==nt) 
      fprintf('  Starting step number: %i \n', n)
   end
      
%%%    Section 3.1: Compute the derivatives, Laplacian and Jacobian.

% x-derivative of w
dwdx(2:nx,1:ny+1) = (w(3:nx+1,1:ny+1)-w(1:nx-1,1:ny+1))/(2*Delta_x);
dwdx(1,1:ny+1) = (w(2,1:ny+1)-w(nx,1:ny+1))/(2*Delta_x);
dwdx(nx+1,1:ny+1) = dwdx(1,1:ny+1);

% y-derivative of w
dwdy(1:nx+1,2:ny) = (w(1:nx+1,3:ny+1)-w(1:nx+1,1:ny-1))/(2*Delta_y);
dwdy(1:nx+1,1) = (w(1:nx+1,2)-w(1:nx+1,ny))/(2*Delta_y);
dwdy(1:nx+1,ny+1) = dwdy(1:nx+1,1);

% Square of the gradient of w
gradsq = dwdx.^2+dwdy.^2;

% Second x-derivative of w
d2wdx2(2:nx,1:ny+1) = (w(3:nx+1,1:ny+1)+w(1:nx-1,1:ny+1)-2*w(2:nx,1:ny+1))/(Delta_x^2);
d2wdx2(1,1:ny+1) = (w(2,1:ny+1)+w(nx,1:ny+1)-2*w(1,1:ny+1))/(Delta_x^2);
d2wdx2(nx+1,1:ny+1) = d2wdx2(1,1:ny+1);

% Second y-derivative of w
d2wdy2(1:nx+1,2:ny) = (w(1:nx+1,3:ny+1)+w(1:nx+1,1:ny-1)-2*w(1:nx+1,2:ny))/(Delta_y^2);
d2wdy2(1:nx+1,1) = (w(1:nx+1,2)+w(1:nx+1,ny)-2*w(1:nx+1,1))/(Delta_y^2);
d2wdy2(1:nx+1,ny+1) = d2wdy2(1:nx+1,1);

laplac = d2wdx2+d2wdy2;

% x-derivative of laplacian
dlapdx(2:nx,1:ny+1) = (laplac(3:nx+1,1:ny+1)-laplac(1:nx-1,1:ny+1))/(2*Delta_x);
dlapdx(1,1:ny+1) = (laplac(2,1:ny+1)-laplac(nx,1:ny+1))/(2*Delta_x);
dlapdx(nx+1,1:ny+1) = dlapdx(1,1:ny+1);

% y-derivative of laplacian
dlapdy(1:nx+1,2:ny) = (laplac(1:nx+1,3:ny+1)-laplac(1:nx+1,1:ny-1))/(2*Delta_y);
dlapdy(1:nx+1,1) = (laplac(1:nx+1,2)-laplac(1:nx+1,ny))/(2*Delta_y);
dlapdy(1:nx+1,ny+1) = dlapdy(1:nx+1,1);

Jacobi = dwdx.*dlapdy - dwdy.*dlapdx;

%%%%%    Compute the Arakawa Jacobian.

Jac1 = Jacobi;

wdldx = w.*dlapdx;
wdldy = w.*dlapdy;

dwdldydx(2:nx,1:ny+1) = (wdldy(3:nx+1,1:ny+1)-wdldy(1:nx-1,1:ny+1))/(2*Delta_x);
dwdldydx(1,1:ny+1) = (wdldy(2,1:ny+1)-wdldy(nx,1:ny+1))/(2*Delta_x);
dwdldydx(nx+1,1:ny+1) = dwdldydx(1,1:ny+1);

dwdldxdy(1:nx+1,2:ny) = (wdldx(1:nx+1,3:ny+1)-wdldx(1:nx+1,1:ny-1))/(2*Delta_y);
dwdldxdy(1:nx+1,1) = (wdldx(1:nx+1,2)-wdldx(1:nx+1,ny))/(2*Delta_y);
dwdldxdy(1:nx+1,ny+1) = dwdldxdy(1:nx+1,1);

Jac2 = dwdldydx - dwdldxdy;

dwdxl = dwdx.*laplac;
dwdyl = dwdy.*laplac;

ddwdxldy(1:nx+1,2:ny) = (dwdxl(1:nx+1,3:ny+1)-dwdxl(1:nx+1,1:ny-1))/(2*Delta_y);
ddwdxldy(1:nx+1,1) = (dwdxl(1:nx+1,2)-dwdxl(1:nx+1,ny))/(2*Delta_y);
ddwdxldy(1:nx+1,ny+1) = ddwdxldy(1:nx+1,1);

ddwdyldx(2:nx,1:ny+1) = (dwdyl(3:nx+1,1:ny+1)-dwdyl(1:nx-1,1:ny+1))/(2*Delta_x);
ddwdyldx(1,1:ny+1) = (dwdyl(2,1:ny+1)-dwdyl(nx,1:ny+1))/(2*Delta_x);
ddwdyldx(nx+1,1:ny+1) = ddwdyldx(1,1:ny+1);

Jac3 = ddwdxldy - ddwdyldx;

Jarakawa = (1/3)*(Jac1+Jac2+Jac3);
 
%%% Use the energy and enstrophy preserving Jacobian.
if(ARAKAWA) Jacobi = Jarakawa; end  

%%%%%%%%%%%%%
%  Compute the function to be stepped forward.
   Q_n = laplac - F*w;
   
%  First time through the loop:
   if(n==1)
     Dt = Delta_t/2;
     Q_nm1 = Q_n;

     %% [rmesh,smesh] = meshdom(0:nx-1,0:ny-1);
     %% rr = rmesh'; ss =flipud(smesh)'; 
     [rmesh,smesh] = meshgrid(0:nx-1,0:ny-1);
     rr = rmesh'; ss =smesh'; 
     C_rs = 2*(cos(2*pi*rr/nx)-1)/Delta_x^2+2*(cos(2*pi*ss/ny)-1)/Delta_y^2 - F;
   end

%  Forward step once per day.
   if(FORWARD & fix(n/nspd)*nspd==n)
     fprintf('Forward timestep once per day \n');
     Dt = Delta_t/2;
     Q_nm1 = Q_n;
   end

%  Calculate the energy and enstrophy integrals
   Rgsq(1:nx,1:ny) = gradsq(1:nx,1:ny);
   Rwsq(1:nx,1:ny) = w(1:nx,1:ny).^2;
   E(n) = 0.5 * mean(mean(Rgsq+F*Rwsq));
   Rgsq(1:nx,1:ny) = laplac(1:nx,1:ny);
   Rwsq(1:nx,1:ny) = w(1:nx,1:ny);
   S(n) = 0.5 * mean(mean((Rgsq-F*Rwsq).^2));

%  Estimate the size of the nonlinear terms
   NonL = mean(mean(abs(Jacobi)));
   Beta = mean(mean(abs(beta0*dwdx)));
   NLsize(n) = NonL/Beta;

   umax = max(max(abs(dwdy)));
   vmax = max(max(abs(dwdx)));
   VMAX = max(umax,vmax);
   CFL_nonlin(n) = abs(VMAX*Delta_t/Delta_x);

%%%    Section 3.2: Step forward one time-step (leapfrog scheme).

   Q_np1 = Q_nm1 - (2*Dt)*((grav/fcor0)*Jacobi + beta0*dwdx + Ubar*dlapdx);

%  Apply the Robert-Asselin filter if required
   if ( RAFILT==1 )
      RA_coeff = RA_COEFF; 
      Q_n = Q_n + RA_coeff*(Q_nm1+Q_np1-2*Q_n);
   end

%%%    Section 3.3: Solve the Helmholtz Equation (Del^2-F)w = R.

%  Compute the fft of the right hand side
%  (strip off additional row and column).
   R(1:nx,1:ny) = Q_np1(1:nx,1:ny);
   R_hat = fft2(R);

%  Compute the transform of the solution
   W_hat = R_hat ./ C_rs ;
 
%  Compute the Wave Power of the components

     a1(n+1) = 2*(W_hat(1+nx-1,1+ny-1)) / nxny; 
     a2(n+1) = 2*(W_hat(1+nx-3,1+1   )) / nxny;
     a3(n+1) = 2*(W_hat(1+4   ,1+0   )) / nxny;

     a1star(n+1) = 2*(W_hat(1+1,   1+1)) / nxny; 
     a2star(n+1) = 2*(W_hat(1+3,1+ny-1)) / nxny;
     a3star(n+1) = 2*(W_hat(1+nx-4,1+0)) / nxny;
 
%  Fourier filtering 
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

%  Compute the inverse transform to get the solution at (n+1)*Delta_t.
   w_new = real(ifft2(W_hat)); % We assume w is real.
   w(1:nx,1:ny) = w_new;
   w(nx+1,1:ny) = w(1,1:ny);     % Fill in additional column at east.
   w(1:nx+1,ny+1)=w(1:nx+1,1);   % Fill in additional row at north.

%  Add the term for the zonal mean flow.
   wtotal = w + Hbar - (fcor0/grav)*Ubar*YY;

% Plot the solution at (about) ten stages over the integration.
        % hold on
        % nplot=10; 
        % np=fix(nt/nplot) ; if(np==0) np=1; end
        % if(fix(n/np)*np==n) mesh(XX,YY,w); pause; end
  % np=1; if(fix(n/np)*np==n) mesh(XX,YY,w); pause; end

%  Save particular values at each time-step. 
   w_center(n) = w(fix(nx/2),fix(ny/2));

%  Save an east-west mid cross-section each time-step. 
%%%      w_section(1:nx+1,n) = w(1:nx+1,fix(ny/2));
 
%  Shuffle the fields at the end of each time-step
   Dt = Delta_t;
   Q_nm1 = Q_n;

%  Save the fields at quarterpoints of the integration
   if(n==1*nt/4) wq1=w; end
   if(n==2*nt/4) wq2=w; end
   if(n==3*nt/4) wq3=w; end
   if(n==4*nt/4) wq4=w; end

     % subplot(2,1,1); contourf(XM,YM,wtotal,vecw); 
     % subplot(2,1,2); pcolor(XM,YM,wtotal);
     % axis('off'); shading interp;
     % pcolor(XM,YM,wtotal);% shading interp;
       if(ICtype==0)
         contourf(XM,YM,wtotal,vecw); drawnow;
         title('Intermediate Stream Function'); 
         colorbar
     %   fprintf('Press RETURN to continue \n'); 
     %   pause(0.1)
       else
         contourf(XM,YM,w,vecw); drawnow;
         
         
         
         title('Intermediate Stream Function'); 
         colorbar
       end
end

%%%%%%%% End of the time-stepping loop  %%%%%%%%%

if(ICtype==0) 
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the final solution
mesh(XM,YM,w);
title('Final Stream Function');
fprintf('Press RETURN to continue \n');  pause

% Fill in particular values at initial time. 
w_center(nt+1:-1:2) = w_center(nt:-1:1);
w_center(1) = w_0(fix(nx/2),fix(ny/2));

hold off

plot(t',w_center);
title('Central Value');
fprintf('Press RETURN to continue \n');  pause

%%%   Plot cross-sections at various times.

%  plot(w_section(:,1:fix(nt/10):nt))
%  title('x-sections');
%  pause

%  Plot the Fourier coefficients at start and end.
  W_hat_N = W_hat;   
  mesh(rr,ss,abs(W_hat_0));
  title('FFT(0)          '); hold on; 
  fprintf('Press RETURN to continue \n'); pause
  mesh(rr,ss,abs(W_hat_N)); 
  title('          FFT(N)'); hold off; 
  fprintf('Press RETURN to continue \n'); pause

% Plot the Amplitudes of the components 
  plotnn=800;
  if(plotnn==80)
    plot(time,abs(a1),':m',time,abs(a2),'--r',time,abs(a3),'-b')
    hgt = max(abs([a1 a2 a3]));
    axis([0,1.1*DAYLEN,0,1.1*hgt])
    title('Component Amplitudes')
    %%  legend('a1','a2','a3',-1)
    if(ICtype==7)
      text(4900,0.41,'a3'); text(4900,0.10,'a1'); text(4900,0.05,'a2')  
      %print -depsc tramps.eps; 
    end
    if(ICtype==10) 
      text(4050,0.80,'a3'); text(4050,0.62,'a1'); text(4050,0.3,'a2')
      %print -depsc trdubl.eps; 
    end
    pause; %% legend off
  end

% Plot the phases of the components
  plotnn=850;
  if(plotnn==85)
    % phi_1 = phase(a1)+(sigma1)*t;
    % phi_2 = phase(a2)+(sigma2)*t;
    % phi_3 = phase(a3)+(sigma3)*t;
    ang_1 = angle(a1)+(sigma1)*t;
    ang_2 = angle(a2)+(sigma2)*t;
    ang_3 = angle(a3)+(sigma3)*t;
    phi_1 = unwrap(ang_1); phi_1 = phi_1-phi_1(1);
    phi_2 = unwrap(ang_2); phi_2 = phi_2-phi_2(1);
    phi_3 = unwrap(ang_3); phi_3 = phi_3-phi_3(1);
    plot(time,phi_1,time,phi_2,time,phi_3);
    title('Component Phases'); pause
  end

%  Plot Energy and Enstrophy  Integrals.
   E(nt+1) = E(nt);  %  Fill in the last value.
   S(nt+1) = S(nt);  %  Fill in the last value.
   NLsize(n+1) = NLsize(nt); %  last value.
   CFL_nonlin(n+1) = CFL_nonlin(nt); %  last value.
   plotnn=88;
   if(plotnn==88)
     plot(time,E); title('Total Energy');
     fprintf('Press RETURN to continue \n'); pause
     plot(time,S); title('Total Enstrophy'); 
     fprintf('Press RETURN to continue \n'); pause
%    plot(time,NLsize); title('Nonlin Size');    pause
%    plot(time,CFL_nonlin); title('CFL_nonlin'); pause
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moreplots = 0;
if(moreplots==0) return; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Calculate and plot Energy Components.
   NRG1 = 0.25*(K1SQ+F)*abs(a1).^2;
   NRG2 = 0.25*(K2SQ+F)*abs(a2).^2;
   NRG3 = 0.25*(K3SQ+F)*abs(a3).^2;
   NRG  = NRG1+NRG2+NRG3;
   EnCor = NRG(1)/E(1);
   E = EnCor*E;
   plotnn=89;
   if(plotnn==89)
     plot(time,NRG1,time,NRG2,time,NRG3,time,NRG,'c',time,E,'m')
     title('ENERGY')
     if(plotnn==891) print -depsc trnrgs.eps; end
     pause
   end

%  Calculate and plot Enstrophy Components.
   ENS1 = (K1SQ+F)*NRG1;
   ENS2 = (K2SQ+F)*NRG2;
   ENS3 = (K3SQ+F)*NRG3;
   ENS = ENS1+ENS2+ENS3;
   StCor = ENS(1)/S(1);
   S = StCor*S;
   plotnn=900;
   if(plotnn==90)
     plot(time,ENS1,time,ENS2,time,ENS3,time,ENS,'c',time,S,'m')
     title('ENSTROPHY')
     pause
   end

% Plot the phase difference and primary energy.
  plotnn=854;
  if(plotnn==854)
    NRG12 = (NRG1+NRG2); NRG12 = NRG12/max(NRG12);
    plot(time,(180/pi)*(phi_1-phi_2)/2,'r',time,180*NRG12,'b');
    grid on;
    axis([0,time(nt+1)*(5000/4800),0,200])
    title(' Azimuth;  Energy in Waves 1 and 2')
    if(plotnn==854) print -depsc trfaze.eps; end
    pause
  end

%%%  Plot the solution at selected times.
   plotnn = 3011;
   if(plotnn==301)
     figure('Position',[300,150,300,600]);  
     whitebg([ 1 1 0.5]);
     zmax=1.1*max(max(w_0)); zmin=-zmax;
     subplot(5,1,1); surf(XM,YM,w_0); 
       title('t=0'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,1,2); surf(XM,YM,wq1); 
       title('t=T/4'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,1,3); surf(XM,YM,wq2); 
       title('t=T/2'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,1,4); surf(XM,YM,wq3); 
       title('t=3T/4'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,1,5); surf(XM,YM,wq4);
       axis([0,7,0,3,zmin,zmax])
       title('t=T'); shading interp;
       xlabel('X'); ylabel('Y'); 
     pause
     print -depsc trsurf.eps;
   end

   plotnn = 3021;
   if(plotnn==302)
     figure('Position',[400,175,300,600]);  
     whitebg([ 1 1 0.5]);
     colormap('jet');
     subplot(5,1,1); pcolor(XM,YM,w_0); axis('off');
       ylabel(''); title('t=0'); shading interp;
     subplot(5,1,2); pcolor(XM,YM,wq1);  axis('off');
       ylabel(''); title('t=T/4'); shading interp;
     subplot(5,1,3); pcolor(XM,YM,wq2);  axis('off');
       ylabel(''); title('t=T/2'); shading interp;
     subplot(5,1,4); pcolor(XM,YM,wq3);  axis('off');
       ylabel(''); title('t=3T/4'); shading interp;
     subplot(5,1,5); pcolor(XM,YM,wq4);
       xlabel('X'); ylabel('Y'); title('t=T'); shading interp;
     pause
     print -depsc tcolor.eps;
   end     

%%%  New plots: solution at five times, in two forms.
   plotnn=202; 
   if(plotnn==200)
     figure('Position',[300,150,550,600]);  
     whitebg([ 1 1 0.5]); colormap('jet');
     zmax=1.1*max(max(w_0)); zmin=-zmax;
     subplot(5,2,1); surf(XM,YM,w_0); 
       title('t=0'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,2,3); surf(XM,YM,wq1); 
       title('t=T/4'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,2,5); surf(XM,YM,wq2); 
       title('t=T/2'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,2,7); surf(XM,YM,wq3); 
       title('t=3T/4'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(5,2,9); surf(XM,YM,wq4);
       axis([0,7,0,3,zmin,zmax])
       title('t=T'); shading interp;
       xlabel('X'); ylabel('Y'); %%  zlabel('Z'); 
  
     subplot(5,2,2); pcolor(XM,YM,w_0); axis('off');
       ylabel(''); title('t=0'); shading interp;
     subplot(5,2,4); pcolor(XM,YM,wq1);  axis('off');
       ylabel(''); title('t=T/4'); shading interp;
     subplot(5,2,6); pcolor(XM,YM,wq2);  axis('off');
       ylabel(''); title('t=T/2'); shading interp;
     subplot(5,2,8); pcolor(XM,YM,wq3);  axis('off');
       ylabel(''); title('t=3T/4'); shading interp;
     subplot(5,2,10); pcolor(XM,YM,wq4);
       xlabel('X'); ylabel('Y'); title('t=T'); shading interp;
     pause
     print -depsc tryeah.eps;
   end

%%%  New plots: solution at three times, in two forms.
   plotnn=2501; 
   if(plotnn==250)
     figure('Position',[300,150,550,600]);  
     whitebg([ 1 1 0.5]); colormap('jet');
     zmax=1.1*max(max(w_0)); zmin=-zmax;
     subplot(3,2,1); surf(XM,YM,w_0); 
       title('t=0'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(3,2,3); surf(XM,YM,wq1); 
       title('t=T/4'); shading interp;
       axis([0,7,0,3,zmin,zmax]); axis('off');
     subplot(3,2,5); surf(XM,YM,wq4);
       axis([0,7,0,3,zmin,zmax])
       title('t=T'); shading interp;
       xlabel('X'); ylabel('Y'); %%  zlabel('Z'); 
  
     subplot(3,2,2); pcolor(XM,YM,w_0); axis('off');
       ylabel(''); title('t=0'); shading interp;
     subplot(3,2,4); pcolor(XM,YM,wq1);  axis('off');
       ylabel(''); title('t=T/4'); shading interp;
     subplot(3,2,6); pcolor(XM,YM,wq4);
       xlabel('X'); ylabel('Y'); title('t=T'); shading interp;
     pause
     if(ICtype==7) print -depsc trsixp.eps; end
   end

   if(ICtype==9)
     if(firstrun)
       w_0reg=w_0; wq4reg=wq4;
       save OBVERS w_0reg wq4reg
       w_0flp=w_0; wq4flp=wq4; % adjust in second run.
     else
       load OBVERS % w_0reg wq4reg
       w_0flp=w_0; wq4flp=wq4; 
     end
     figure('Position',[400,175,600,400]);  
     whitebg([1 1 0.5]); colormap('jet'); 
     subplot(2,2,1); pcolor(XM,YM,w_0reg); axis('off');
       ylabel(''); title('t=0'); shading interp;
       text(3.1,1.0,'+')
     subplot(2,2,3); pcolor(XM,YM,wq4reg);
       title('t=T'); shading interp
       text(3.1,1.0,'+')
     subplot(2,2,2); pcolor(XM,YM,w_0flp); axis('off');
       ylabel(''); title('t=0'); shading interp;
       text(3.1,1.0,'+')
     subplot(2,2,4); pcolor(XM,YM,wq4flp);
       title('t=T'); shading interp
       text(3.1,1.0,'+')
     set(gcf,'DefaultTextColor','yellow')
      subplot(2,2,1); text(3.1,1.0,'+')
      subplot(2,2,2); text(3.1,1.0,'+')
      subplot(2,2,3); text(3.1,1.0,'+')
      subplot(2,2,4); text(3.1,1.0,'+')
     pause
     print -depsc trfcst.eps;
   end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Precession: examine and plot the triad precession.
%
%  Change from Triad variables to Spring variables.

Postproc=1;
if(Postproc)
  clf
           
%  Extract the modulus and phase information.
  Ampscale = gamma;
  a1_mod = abs(a1)/Ampscale;
  a2_mod = abs(a2)/Ampscale;
  a3_mod = abs(a3)/Ampscale;
  a_mod_0 = [ a1_mod(1) a2_mod(1) a3_mod(1) ];

  hold off
  plot(time,a1_mod,time,a2_mod,time,a3_mod);
  title('Component Moduli'); pause
  
  % Phases already obtained above.
  % phi_1 = phase(a1)+(sigma1)*t;
  % phi_2 = phase(a2)+(sigma2)*t;
  % phi_3 = phase(a3)+(sigma3)*t;
 
  r2d=(180/pi);
  plot(time,phi_1,time,phi_2,time,phi_3);
  title('Component Phases'); pause
  theta = (phi_1-phi_2)/2;
  plot(time,r2d*theta);
  title('Half-Phase Difference (degrees)'); pause

  muprod = mu1*mu2*mu3;
  A1 = (muprod/mu1)*a1_mod.*exp(i*phi_1);
  A2 = (muprod/mu2)*a2_mod.*exp(i*phi_2);
  A3 = (muprod/mu3)*a3_mod.*exp(i*phi_3);

  % Approximate precession angle.
  Delta_THETA = r2d*(A1(1)^2-A2(1)^2)/(A1(1)^2+A2(1)^2)

  % title('Triad Precession')
  clf; title('')
  SAmaj = abs(A1)+abs(A2);
  SAmin = abs(A1)-abs(A2);
  Xenv = SAmaj.*cos(theta);
  Yenv = SAmaj.*sin(theta);
  patch( Xenv, Yenv,'b','EdgeColor','b'); hold on;
  patch(-Xenv,-Yenv,'b','EdgeColor','b');
  axis('square'); axis('equal'); axis('off');
  pause
  print -depsc trstar.eps;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
