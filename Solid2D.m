% This code is developed based on FastICE (Thermomechanical ice
% deformation models 3D. Boundary conditions, solid rheology,
% nondimensionalized variables has been changed. New nondimensionalized
% variables and equations have been added. The purpose of the code is
% to show the stress distribution of a solid under a stress condition. 
% 
% FastICE is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% FastICE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with FastICE. If not, see <http://www.gnu.org/licenses/>.

clear
tic
% Physics
A_star  = 4e-25;                          % Coefficient1 in Chemia et al. 2009
n_pow   = 2.0;                            % Coefficient2 in Chemia et al. 2009
A       = (3^(n_pow/2+1/2)*A_star)^(-1);  % Coefficient3 in Chemia et al. 2009
rhoS    = 2200;                           % Solid density
rho_Water = 1020;                         % Seawater density
D_sf    = 600;                            % Seafloor depth
g       = 9.81;                           % Gravity
SR      = 0.5;                            % Stress ratio
% Input values
Lz      = 1000;                           % Domain's Z-length
aspect  = 5;                              % Aspect ratio
% Scales
rho_sc  = rhoS;
eta_sc  = 1e16;
L_sc    = power( (eta_sc^2)/g/(rho_sc^2) , 1/3);  % Make PI2 = 1 to stabilize
Lz      = Lz/L_sc;                        % Non-dimensionalized Z-length
PI2     = g*(L_sc^3)*(rho_sc^2)/(eta_sc^2);
v_sc    = rho_sc*g*(L_sc^2)/eta_sc;
p_sc    = (eta_sc^2)/(rho_sc*L_sc^2);
p_water = rho_Water*g*D_sf/p_sc;          % Seawater pressure at the top boundary
p_NDHSM = rhoS*g*Lz*L_sc/p_sc;            % Nondimensionalized hydrostatic maximum pressure
% Numerics
nx      = 200;
nz      = 80;
eta0    = 1e4;
epsi    = 1e-8; 
damp    = 2; 
relV    = 1/1;
relP    = 1/2;
eta_b   = 1/2;
rele    = 0.1;
niter   = 5e6;
nout_iter = 1000;
% Pre-processing
Lx          = aspect*Lz;
dx          = Lx/(nx-1); 
dz          = Lz/(nz-1);
xn          = -dx/2:dx:Lx+dx/2; 
zn          = -dz/2:dz:Lz+dz/2;
xc          = 0:dx:Lx;
zc          = 0:dz:Lz;
[xc2,  zc2] = ndgrid(xc,zc);
Vx          = zeros(nx+1,nz);
Vz          = zeros(nx,nz+1);
P           = zeros(nx,nz);
phi0        = zeros(nx,nz);
Exzn        = zeros(nx,nz);
txzn        = zeros(nx,nz);
etan        = ones(nx,nz);
txx         = ones(nx,nz);
tzz         = ones(nx,nz);
etas        = zeros(nx-1,nz-1);
dVxdt       = zeros(nx-1,nz-2);
dVzdt       = zeros(nx-2,nz-1); 
for iter = 1:niter % Pseudo-Transient cycles
    % Timesteps - GPU KERNEL 1
    etax     = 0.5*(etan(2:end  ,2:end-1) + etan(1:end-1,2:end-1));
    etaz     = 0.5*(etan(2:end-1,2:end  ) + etan(2:end-1,1:end-1));
    dtP      = relP*4.1/max(nx,nz)*etan*(1+eta_b);
    dtVx     = relV*min(dx,dz)^2./etax/(1+eta_b)/4.1;
    dtVz     = relV*min(dx,dz)^2./etaz/(1+eta_b)/4.1;
    % Pressure update and strain rates - GPU KERNEL 2
    divV     = diff(Vx,1,1)/dx + diff(Vz,1,2)/dz;
    P        = P - divV.*dtP;
    Exx      = diff(Vx,1,1)/dx - 1/2*divV;
    Ezz      = diff(Vz,1,2)/dz - 1/2*divV;
    Exz      = 0.5*(diff(Vx(2:end-1,:),1,2)/dz+diff(Vz(:,2:end-1),1,1)/dx);                  
    % Stresses - GPU KERNEL 3
    txx      = 2*etan.*( diff(Vx,1,1)./dx - 1/2*divV + eta_b*divV);
    tzz      = 2*etan.*( diff(Vz,1,2)./dz - 1/2*divV + eta_b*divV);
    txz      = 2*etas.*Exz;
    etas     = 0.25*(etan(1:end-1,1:end-1)+etan(2:end  ,1:end-1)+etan(1:end-1,2:end  ) + etan(2:end  ,2:end  ));
    txzn(2:end-1,2:end-1) = 0.25*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end  )+txz(2:end,2:end  ));
    Exzn(2:end-1,2:end-1) = 0.25*(Exz(1:end-1,1:end-1)+Exz(2:end,1:end-1)+Exz(1:end-1,2:end  )+Exz(2:end,2:end  ));
    % GPU KERNEL 4
    Exzn(:,1) = Exzn(:,2);
    % GPU KERNEL 5
    ResVx       = -diff(P(:,2:end-1),1,1)/dx + PI2*diff(txz,1,2)/dz + PI2*diff(txx(:,2:end-1),1,1)/dx;
    ResVz       = -diff(P(2:end-1,:),1,2)/dz + PI2*diff(txz,1,1)/dx + PI2*diff(tzz(2:end-1,:),1,2)/dz - PI2*(rhoS/rho_sc);
    dVxdt       = dVxdt.*(1-damp/nx) + ResVx;
    dVzdt       = dVzdt.*(1-damp/nz) + ResVz;   
    EII2        = (Exx.^2 + Ezz.^2)/2 + Exzn.^2;
    tII         = sqrt((txx.^2 + tzz.^2)/2 + txzn.^2);
    etait       = (A/eta_sc).*(rho_sc*g*L_sc.*tII).^(1-n_pow);
    % GPU KERNEL 6
    Vx(2:end-1,2:end-1) = Vx(2:end-1,2:end-1) + dVxdt.*dtVx;
    Vz(2:end-1,2:end-1) = Vz(2:end-1,2:end-1) + dVzdt.*dtVz;
    etan            = min(exp(rele*log(etait) + (1-rele)*log(etan)),eta0);
    % GPU KERNEL 7
    etan(:,[1 end]) = etan(:,[2 end-1]);
    % GPU KERNEL 8
    etan([1 end],:) = etan([2 end-1],:);
    % GPU KERNEL 9
    % Boundary conditions
    % BC left
    Vx(1,:)         = Vx(2,:) + dx*(p_water+SR*p_NDHSM*linspace(0,1,nz)-P(1,:))./etan(1,:)/2 + dx*(1/2-eta_b)*divV(1,:);
    Vz(1,2:end-1)   = Vz(2,2:end-1);
    % BC right
    Vx(end,:)       = Vx(end-1,:) - dx*(p_water+SR*p_NDHSM*linspace(0,1,nz)-P(end,:))./etan(end,:)/2 + dx*(1/2-eta_b)*divV(end,:);
    Vz(end,2:end-1) = Vz(end-1,2:end-1);
    % GPU KERNEL 10
    % BC bottom
    Vx(:,1)   = Vx(:,2);
    Vz(:,1)   = 0.0;
    % BC top             
    Vz(:,end) = Vz(:,end-1) + dz*(-p_water+P(:,end))./etan(:,end)/2 + dz*(1/2-eta_b)*divV(:,end);
    Vx(:,end) = Vx(:,end-1);
%% errors
    norm_X  = sqrt(sum(dVxdt(:).*dVxdt(:)))/(nx*nz);
    norm_Z  = sqrt(sum(dVzdt(:).*dVzdt(:)))/(nx*nz);
    norm_P  = sqrt(sum(divV(:).*divV(:)))/(nx*nz);
    % Monitoring
    if (mod(iter,100)==1) 
    iterror = max(norm_X,max(norm_Z,norm_P));
    if ((iterror < epsi) && (iter > 2)),break; end
    end
    if (mod(iter,nout_iter)==1) 
        fprintf('iter=%d, err=%1.3e \n',iter,iterror)
        figure(1)
        clf
        colormap jet
        subplot(4,2,1)
        imagesc(flipud(txx'))
        colorbar
        subplot(4,2,2)
        imagesc(flipud(((tzz'))))
        colorbar
        subplot(4,2,3)
        imagesc(flipud(txz'))
        colorbar
        subplot(4,2,4)
        imagesc(flipud(((P'))))
        colorbar
        subplot(4,2,5)
        imagesc(flipud(((etan'))))
        colorbar
        subplot(4,2,6)
        imagesc(flipud(((Vx'))))
        colorbar
        subplot(4,2,7)
        imagesc(flipud(((Vz'))))
        colorbar
        subplot(4,2,8)
        imagesc(flipud(((phi0'))))
        colorbar
        drawnow
        set(gca,'color','white')
    end
end
CPU_time=toc;
MTP_eff = (nx*nz)*iter*(2*3+2)*8/1e9/CPU_time;
fprintf('iter=%d, err=%1.3e \n',iter,iterror)
fprintf('time=%f, MTP_eff=%f \n',CPU_time,MTP_eff)

Vx_MAT   = Vx;
Vz_MAT   = Vz;
P_MAT    = P;
txx_MAT  = txx;
tzz_MAT  = tzz;
txz_MAT  = txz;
xc_MAT   = xc;
xn_MAT   = xn;
%save('data_MAT_8.mat','xc_MAT','xn_MAT','Vx_MAT','Vz_MAT','P_MAT','txx_MAT','tzz_MAT','txz_MAT')