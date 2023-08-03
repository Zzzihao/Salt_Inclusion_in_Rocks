% Read the Boundary of the Stress Distribution. 
% Salt Inclusion. Smear out the coefficients of Stokes euqation.

clc, clear;
tic;

% Physics
A_Salt      = 1e17;                     % Salt viscosity [Pa-s]
A_Shale     = 1e15;                     % Shale viscosity [Pa-s]
g           = 9.81;                     % Gravity [m/s2]
p_grad_h    = 0.75*22620.6;             % Minimum horizontal stress gradient [Pa/m]
p_grad_v    = 0.85*22620.6;             % Vertical stress gradient [Pa/m]
rho_Shale   = p_grad_v/g;               % Bulk Shale density [kg/m3]
rho_Salt    = 2200;                     % Salt density [kg/m3]
rho_Water   = 1020;                     % Water density [kg/m3]
D_sf        = 800;                      % Seafloor depth [m]
% Input dimensions
Lz          = 10000;                    % Domain's Z-length
aspect      = 3;                        % Aspect ratio
% Scales
rho_sc      = rho_Shale;                % Density scale
eta_sc      = A_Shale;                  % Rock viscosity scale
L_sc        = power( (eta_sc^2)/g/(rho_sc^2) , 1/3); % Length scale
PI2         = g*(L_sc^3)*(rho_sc^2)/(eta_sc^2); % Make PI2 = 1 to stabilize
Lz          = Lz/L_sc;                  % Non-dimensionalized Z-length
v_sc        = rho_sc*g*(L_sc^2)/eta_sc; % Velocity scale
p_sc        = (eta_sc^2)/(rho_sc*L_sc^2); % Pressure scale
p_water     = rho_Water*g*D_sf/p_sc;    % Seawater pressure at the top boundary
p_NMHP      = p_grad_h*Lz*L_sc/p_sc;    % Nondimensionalized minimum horizontal pressure
p_NVP       = p_grad_v*Lz*L_sc/p_sc;    % Nondimensionalized vertical pressure

% Numerics
nx          = 20;                       % Number of blocks in X-direction
nz          = 20;                       % Number of blocks in Z-direction
epsi        = 1e-8;                     % Tolerance
niter       = 5e7;                      % Maximum iterations
nout_iter   = 1000;                     % Output iterations
CFL         = 0.9/sqrt(2);              % Courant-Friedrichs-Lewy
Re          = 3/2*sqrt(10)*pi;          % Reynolds number
r           = 1.0;                      % Bulk to shear modulus ratio
err         = 2*epsi;                   % Initialize errors
err_evo1    = [];
err_evo2    = [];
iter        = 0;                        % Initialize iteration
% Pre-processing
Lx          = aspect * Lz;              % X-length
dx          = Lx/(nx-1);                % Block X-size
dz          = Lz/(nz-1);                % Block Z-size
xn          = -dx/2:dx:Lx+dx/2;         % Grid X-location
zn          = -dz/2:dz:Lz+dz/2;         % Grid Z-location
xc          = 0:dx:Lx;                  % Block center X-location
zc          = 0:dz:Lz;                  % Block center Z-location
[xc2,  zc2] = ndgrid(xc,zc);            % Block center X-Z mesh
max_lxz     = max(Lx, Lz);              % Maximum block length
Vpdt        = min(dx, dz) * CFL;        % P-wave velocity * pseudo time step
Vx          = zeros(nx+1,nz+1);         % Initialize Velocity-X
Vz          = zeros(nx+1,nz+1);         % Initialize Velocity-Z
P           = zeros(nx,nz);             % Initialize Total Pressure
txx         = ones(nx,nz);              % Initialize Deviatoric Stress XX
tzz         = ones(nx,nz);              % Initialize Deviatoric Stress ZZ
txz         = zeros(nx, nz);            % Initialize Deviatoric Stress XZ

% Salt location and density
SaltFraction= zeros(nx,nz);             % Initialize Salts Fraction
%img         = imread('SaltGeometry2D-4.JPG'); % Load the salt geometry
img         = imread('SaltGeometry2D-5.JPG'); % Load the salt geometry
bw          = imrotate(imbinarize(rgb2gray(img)),-90); % Create a binary image
bw          = cast(bw, class(SaltFraction));  % Convert binary to double
SaltFraction= imresize(1-bw, [nx, nz], 'nearest'); % Resize the image to [nx, nz]
etan        = ones(nx, nz);                   % Initialize shale viscosity
etan(SaltFraction==1.0) = A_Salt/A_Shale;     % Initialize salts viscosity
rho_Rock    = rho_Shale * ones(nx,nz);  % Initialize shale density
rho_Rock(SaltFraction==1.0) = rho_Salt; % Initialize salts density
rho_Rockz   = XZ_mean(rho_Rock);
% Smear out the coefficients of Stokes equation
for ism     = 1:10
    etan    = smooth(etan, 1.0);        % Smoothed rock viscosity
end
% Calulate the maximum rock viscosity among neighbors
etanmax         = compute_maxloc(etan);
etanmax(1, :)   = etanmax(2, :);
etanmax(end, :) = etanmax(end - 1, :);
etanmax(:, 1)   = etanmax(:, 2);
etanmax(:, end) = etanmax(:, end - 1);
% Pre-processing of numerics
dt_Rho      = zeros(nx, nz);            % Pseudo time step / rho
dt_Rho      = Vpdt*max_lxz/Re./etanmax;
Gdt         = zeros(nx, nz);            % G * pseudo time step
Gdt         = (Vpdt^2)./dt_Rho/(r+2);

% Pseudo Transient Iterations
while err > epsi && iter <= niter
    % Update Pressure
    divV    = Z_mean(diff(Vx,1,1)/dx) + X_mean(diff(Vz,1,2)/dz);
    P       = P - r*Gdt.*divV;                 
    % Strain rate and Stress
    Exx     = Z_mean(diff(Vx, 1, 1)/dx);
    Ezz     = X_mean(diff(Vz, 1, 2)/dz);
    Exz     = 0.5*( X_mean(diff(Vx,1,2)/dz) + Z_mean(diff(Vz,1,1)/dx) );
    txx     = (txx+2*Gdt.*Exx) ./ (Gdt./etan + 1);
    tzz     = (tzz+2*Gdt.*Ezz) ./ (Gdt./etan + 1);
    txz     = (txz+2*Gdt.*Exz) ./ (Gdt./etan + 1);
    % Residuals
    RX      = Z_mean(-diff(P,1,1)/dx) + X_mean(PI2*diff(txz,1,2)/dz) + Z_mean(PI2*diff(txx,1,1)/dx);
    RZ      = X_mean(-diff(P,1,2)/dz) + Z_mean(PI2*diff(txz,1,1)/dx) + X_mean(PI2*diff(tzz,1,2)/dz) - PI2*(rho_Rockz/rho_sc);
    dVx     = XZ_mean(dt_Rho) .*RX;
    dVz     = XZ_mean(dt_Rho) .*RZ;
    % Update velocity
    Vx(2:end-1, 2:end-1) = Vx(2:end-1, 2:end-1) + dVx;
    Vz(2:end-1, 2:end-1) = Vz(2:end-1, 2:end-1) + dVz;
    % BC bottom
    Vz(2:end-1,1)        = Vz(2:end-1,2) - dz*(-p_water-p_NVP+X_mean(P(:,1)))./X_mean(etan(:,1))/2;
    % BC top
    Vz(2:end-1,end)      = Vz(2:end-1,end-1) + dz*(-p_water+X_mean(P(:,end)))./X_mean(etan(:,end))/2;
    % BC left
    MM = linspace(1,0,nz+1);
    Vx(1,2:end-1)        = Vx(2,2:end-1) - dx*( -p_NMHP*MM(2:end-1)+Z_mean(P(1,:)) )/2./Z_mean(etan(1,:));
    % BC right
    Vx(end,2:end-1)      = Vx(end-1,2:end-1) + dx*( -p_NMHP*MM(2:end-1)+Z_mean(P(end,:)) )/2./Z_mean(etan(end,:));
    
    % Update iteration
    iter        = iter + 1;
    % Compare error and tolerance
    if mod(iter, nout_iter) == 0
        Vmin    = min(Vx(:));
        Vmax    = max(Vx(:));
        Pmin    = min(P(:));
        Pmax    = max(P(:));
        norm_Rx = norm(RX) / (Pmax - Pmin) * Lx / sqrt(numel(RX));
        norm_Rz = norm(RZ) / (Pmax - Pmin) * Lx / sqrt(numel(RZ));
        norm_dV = norm(divV) / (Vmax - Vmin) * Lx / sqrt(numel(divV));
        err     = max([norm_Rx, norm_Rz, norm_dV]);
        err_evo1(end+1) = err;
        err_evo2(end+1) = iter;
        fprintf('Total steps = %d, err = %.3e [norm_Rx=%.3e, norm_Ry=%.3e, norm_dV=%.3e]\n', iter, err, norm_Rx, norm_Rz, norm_dV);
    end
    % Make plots
    if mod(iter, nout_iter) == 1
        clf;
        colormap jet;
        subplot(3,3,1);
        imagesc(flipud(txx'));title('Deviatoric Stress txx')
        colorbar;
        subplot(3,3,2);
        imagesc(flipud(((tzz'))));title('Deviatoric Stress tzz')
        colorbar;
        subplot(3,3,3);
        imagesc(flipud(txz'));title('Stress xz')
        colorbar;
        subplot(3,3,4);
        imagesc(flipud(P'));title('Pressure')
        colorbar;
        subplot(3,3,5);
        imagesc(flipud(Vx'));title('Velocity x')
        colorbar;
        subplot(3,3,6);
        imagesc(flipud(Vz'));title('Velocity z')
        colorbar;
        subplot(3,3,7);
        imagesc(flipud((txx+P)'));title('stress xx')
        colorbar;
        subplot(3,3,8);
        imagesc(flipud((tzz+P)'));title('stress zz')
        colorbar;
        subplot(3,3,9);
        imagesc(flipud(etan'));
        colorbar;
        drawnow;
        set(gca,'color','white');
    end
end
