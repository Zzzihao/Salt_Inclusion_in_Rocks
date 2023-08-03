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
aspectX     = 3;                        % Aspect ratio X
aspectY     = 3;                        % Aspect ratio Y
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
nx          = 10;                       % Number of blocks in X-direction
ny          = 10;                       % Number of blocks in Y-direction
nz          = 10;                       % Number of blocks in Z-direction
epsi        = 1e-8;                     % Tolerance
niter       = 5e7;                      % Maximum iterations
nout_iter   = 1000;                     % Output iterations
nout_save   = 10000;                    % Save iterations
CFL         = 0.8/sqrt(3);              % Courant-Friedrichs-Lewy
Re          = 3/2*sqrt(10)*pi;          % Reynolds number
r           = 1.0;                      % Bulk to shear modulus ratio
err         = 2*epsi;                   % Initialize errors
err_evo1    = [];
err_evo2    = [];
iter        = 0;                        % Initialize iteration
% Pre-processing
Lx          = aspectX * Lz;             % X-length
Ly          = aspectY * Lz;             % Y-length
dx          = Lx/(nx-1);                % Block X-size
dy          = Ly/(nx-1);                % Block Y-size
dz          = Lz/(nz-1);                % Block Z-size
xn          = -dx/2:dx:Lx+dx/2;         % Grid X-location
yn          = -dy/2:dy:Ly+dy/2;         % Grid Y-location
zn          = -dz/2:dz:Lz+dz/2;         % Grid Z-location
xc          = 0:dx:Lx;                  % Block center X-location
yc          = 0:dy:Ly;                  % Block center Y-location
zc          = 0:dz:Lz;                  % Block center Z-location
[Xc,Yc,Zc]  = ndgrid(xc,yc,zc);         % Mesh location
max_lxyz    = max(max(Lx,Ly),Lz);       % Maximum block length
Vpdt        = min(min(dx,dy),dz) * CFL; % P-wave velocity * pseudo time step
Vx          = zeros(nx+1,ny+1,nz+1);    % Initialize Velocity-X
Vy          = zeros(nx+1,ny+1,nz+1);    % Initialize Velocity-Y
Vz          = zeros(nx+1,ny+1,nz+1);    % Initialize Velocity-Z
P           = zeros(nx,ny,nz);          % Initialize Total Pressure
txx         = ones(nx,ny,nz);           % Initialize Deviatoric Stress XX
tyy         = ones(nx,ny,nz);           % Initialize Deviatoric Stress YY
tzz         = ones(nx,ny,nz);           % Initialize Deviatoric Stress ZZ
txy         = zeros(nx,ny,nz);          % Initialize Deviatoric Stress XY
txz         = zeros(nx,ny,nz);          % Initialize Deviatoric Stress XZ
tyz         = zeros(nx,ny,nz);          % Initialize Deviatoric Stress YZ
P_BC        = repmat(reshape(linspace(1,0,nz+1), [1,1,nz+1]), [nx+1,ny+1,1]);

% Salt location and density
radc     = ((Xc - Lx*.5)./4).^2 +((Yc - Ly*.5)./4).^2 + ((Zc - Lz*.5)./1).^2 ;
SaltFraction= zeros(nx,ny,nz);          % Initialize Salts Fraction
SaltFraction(radc<0.2/L_sc) = 1;
etan      = ones(nx,ny,nz); 
etan(radc<0.2/L_sc) = A_Salt/A_Shale;
rho_Rock    = rho_Shale * ones(nx,ny,nz);  % Initialize shale density
rho_Rock(SaltFraction==1.0) = rho_Salt; % Initialize salts density
rho_Rockz   = XYZ_mean3D(rho_Rock);

% Smear out the coefficients of Stokes equation
for ism     = 1:10
    etan    = smooth3D(etan, 1.0);      % Smoothed rock viscosity
end
% Calulate the maximum rock viscosity among neighbors
etanmax         = compute_maxloc3D(etan);
etanmax(1, :, :)   = etanmax(2, :, :);
etanmax(end, :, :) = etanmax(end - 1, :, :);
etanmax(:, 1, :)   = etanmax(:, 2, :);
etanmax(:, end, :) = etanmax(:, end - 1, :);
etanmax(:, :, 1)   = etanmax(:, :, 2);
etanmax(:, :, end) = etanmax(:, :, end - 1);
% Pre-processing of numerics
dt_Rho      = zeros(nx,ny,nz);           % Pseudo time step / rho
dt_Rho      = Vpdt*max_lxyz/Re./etanmax;
Gdt         = zeros(nx,ny,nz);           % G * pseudo time step
Gdt         = (Vpdt^2)./dt_Rho/(r+2);

% Pseudo Transient Iterations
while err > epsi && iter <= niter
    % Update Pressure
    divV    = YZ_mean3D(diff(Vx,1,1)/dx) + XZ_mean3D(diff(Vy,1,2)/dy) + XY_mean3D(diff(Vz,1,3)/dz);
    P       = P - r*Gdt.*divV;
    % Strain rate and Stress
    Exx     = YZ_mean3D(diff(Vx, 1, 1)/dx);
    Eyy     = XZ_mean3D(diff(Vy, 1, 2)/dy);
    Ezz     = XY_mean3D(diff(Vz, 1, 3)/dz);
    Exy     = 0.5*( XZ_mean3D(diff(Vx,1,2)/dy) + YZ_mean3D(diff(Vy,1,1)/dx) );
    Exz     = 0.5*( XY_mean3D(diff(Vx,1,3)/dz) + YZ_mean3D(diff(Vz,1,1)/dx) );
    Eyz     = 0.5*( XY_mean3D(diff(Vy,1,3)/dz) + XZ_mean3D(diff(Vz,1,2)/dy) );
    txx_o   = txx;
    tyy_o   = tyy;
    tzz_o   = tzz;
    txy_o   = txy;
    txz_o   = txz;
    tyz_o   = tyz;
    txx     = (txx_o+2*Gdt.*Exx) ./ (Gdt./etan + 1);
    tyy     = (tyy_o+2*Gdt.*Eyy) ./ (Gdt./etan + 1);
    tzz     = (tzz_o+2*Gdt.*Ezz) ./ (Gdt./etan + 1);
    txy     = (txy_o+2*Gdt.*Exy) ./ (Gdt./etan + 1);
    txz     = (txz_o+2*Gdt.*Exz) ./ (Gdt./etan + 1);
    tyz     = (tyz_o+2*Gdt.*Eyz) ./ (Gdt./etan + 1);
    sigmaxx = txx + P;
    sigmayy = tyy + P;
    sigmazz = tzz + P;

    % Residuals
    RX    = YZ_mean3D(-diff(P,1,1)/dx) + XZ_mean3D(PI2*diff(txy,1,2)/dy) + XY_mean3D(PI2*diff(txz,1,3)/dz) + YZ_mean3D(PI2*diff(txx,1,1)/dx);
    RY    = XZ_mean3D(-diff(P,1,2)/dy) + YZ_mean3D(PI2*diff(txy,1,1)/dx) + XY_mean3D(PI2*diff(tyz,1,3)/dz) + XZ_mean3D(PI2*diff(tyy,1,2)/dy);
    RZ    = XY_mean3D(-diff(P,1,3)/dz) + YZ_mean3D(PI2*diff(txz,1,1)/dx) + XZ_mean3D(PI2*diff(tyz,1,2)/dy) + XY_mean3D(PI2*diff(tzz,1,3)/dz) - PI2*(rho_Rockz/rho_sc);
    dVx   = XYZ_mean3D(dt_Rho) .*RX;
    dVy   = XYZ_mean3D(dt_Rho) .*RY;
    dVz   = XYZ_mean3D(dt_Rho) .*RZ;
    
    % Update velocity
    Vx(2:end-1, 2:end-1, 2:end-1) = Vx(2:end-1, 2:end-1, 2:end-1) + dVx;
    Vy(2:end-1, 2:end-1, 2:end-1) = Vy(2:end-1, 2:end-1, 2:end-1) + dVy;
    Vz(2:end-1, 2:end-1, 2:end-1) = Vz(2:end-1, 2:end-1, 2:end-1) + dVz;
    
    % BC bottom
    Vz(2:end-1,2:end-1,1)     = Vz(2:end-1,2:end-1,2) - dz*(-p_water-p_NVP+XY_mean3D(P(:,:,1)))./XY_mean3D(etan(:,:,1))/2;
    % BC top
    Vz(2:end-1,2:end-1,end)   = Vz(2:end-1,2:end-1,end-1) + dz*(-p_water+XY_mean3D(P(:,:,end)))./XY_mean3D(etan(:,:,end))/2;
    % BC left
    Vx(1,2:end-1,2:end-1)     = Vx(2,2:end-1,2:end-1) - dx*( -p_NMHP*P_BC(1,2:end-1,2:end-1)+YZ_mean3D(P(1,:,:)) )/2./YZ_mean3D(etan(1,:,:));
    % BC right
    Vx(end,2:end-1,2:end-1)   = Vx(end-1,2:end-1,2:end-1) + dx*( -p_NMHP*P_BC(end,2:end-1,2:end-1)+YZ_mean3D(P(end,:,:)) )/2./YZ_mean3D(etan(end,:,:));
    % BC front
    Vy(2:end-1,1,2:end-1)     = Vy(2:end-1,2,2:end-1) - dy*( -p_NMHP*P_BC(2:end-1,1,2:end-1)+XZ_mean3D(P(:,1,:)) )/2./XZ_mean3D(etan(:,1,:));
    % BC back
    Vy(2:end-1,end,2:end-1)   = Vy(2:end-1,end-1,2:end-1) + dy*( -p_NMHP*P_BC(2:end-1,end,2:end-1)+XZ_mean3D(P(:,end,:)) )/2./XZ_mean3D(etan(:,end,:));
    
    % Update iteration
    iter        = iter + 1;
    
    % Compare error and tolerance
    if mod(iter, nout_iter) == 0
        Vmin    = min(Vx(:));
        Vmax    = max(Vx(:));
        Pmin    = min(P(:));
        Pmax    = max(P(:));
        norm_Rx = sqrt(sum(RX(:).^2)) / (Pmax - Pmin) * Lx / sqrt(numel(RX));
        norm_Ry = sqrt(sum(RY(:).^2)) / (Pmax - Pmin) * Ly / sqrt(numel(RY));
        norm_Rz = sqrt(sum(RZ(:).^2)) / (Pmax - Pmin) * Lx / sqrt(numel(RZ));
        norm_dV = sqrt(sum(divV(:).^2)) / (Vmax - Vmin) * Lx / sqrt(numel(divV));
        err     = max([norm_Rx, norm_Ry, norm_Rz, norm_dV]);
        err_evo1(end+1) = err;
        err_evo2(end+1) = iter;
        fprintf('Total steps = %d, err = %.3e [norm_Rx=%.3e, norm_Ry=%.3e, norm_Rz=%.3e, norm_dV=%.3e]\n', iter, err, norm_Rx, norm_Ry, norm_Rz, norm_dV);
    end
    if mod(iter, nout_save) == 1
        save('data3D.mat');
    end
end
