
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
%%%   For the Fourier Transforms, the redundant row     %%%,
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify type of Initial Conditions.
  DAYLEN=1;          %  Forecast length in days.
  NX = 61;  NY = 21;  %  Spatial resolution
  DELTA_t = 1/12;      %  Timestep in hours.

     ICtype =  0;
      
%%%%%%%%%%%%%%%%% =Pseudo-real 500 mb flow.

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
%  w is periodic in both directions. Z is not.%% 


%Pseudo-real 500mb geopotential field.

   % Define the field size.
   Z_0 = zeros(nx+1,ny+1);
   
   %    Set the seed of the random number generator
   %    to its initial value (same results each run).
            rand('state',0)
   %    Set the seed of the random number generator
   %    to yield different results each run.
        %rand('state',sum(100*clock))
   
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

for n=1:1
   
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
d2wdx2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
