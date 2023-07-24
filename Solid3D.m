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
A_star  = 1e-22;                          % Coefficient1 in Chemia et al. 2009
n_pow   = 2.0;                            % Coefficient2 in Chemia et al. 2009
A       = (3^(n_pow/2+1/2)*A_star)^(-1);  % Coefficient3 in Chemia et al. 2009
rhoS    = 2200;                           % Solid density
g       = 9.81;                           % Gravity
SRX     = 0.5;                            % Stress ratio X-Z
SRY     = 0.5;                            % Stress ratio Y-Z
% Input values
Lz      = 1000;                           % Z-length
aspectX = 5;
aspectY = 5;
% Scales
rho_sc  = rhoS;
eta_sc  = 1e14;
L_sc    = power( (eta_sc^2)/g/(rho_sc^2) , 1/3);  % Make PI2 = 1 to stabilize
Lz      = Lz/L_sc;                        % Non-dimensionalized Z-length
PI2     = g*(L_sc^3)*(rho_sc^2)/(eta_sc^2);
V_sc    = rho_sc*g*(L_sc^2)/eta_sc;
P_sc    = (eta_sc^2)/(rho_sc*L_sc^2);
p_NDHSM = rhoS*g*Lz*L_sc/P_sc;            % Nondimensionalized hydrostatic maximum pressure
% Numerics 
nx      = 25*2;
ny      = 25*2;
nz      = 10*2;
eta0    = 1e4;
epsi    = 1e-8; 
damp    = 2; 
relV    = 1/1;
relP    = 1/4;
eta_b   = 1/4;
rele    = 0.1;
niter   = 5e6;
nout_iter = 1000;
% Pre-processing
Lx      = aspectX*Lz;
Ly      = aspectY*Lz;
dx      = Lx/(nx-1); 
dy      = Ly/(ny-1);
dz      = Lz/(nz-1);
xc      = 0:dx:Lx; 
yc      = 0:dy:Ly;
zc      = 0:dz:Lz;
xn      = -dx/2:dx:Lx+dx/2; 
yn      = -dy/2:dy:Ly+dy/2; 
zn      = -dz/2:dz:Lz+dz/2;
[xc1,yc1,zc1] = ndgrid(xc,yc,zc);
[xn2,yc2]     = ndgrid(xn,yc);
[xc2,yn2]     = ndgrid(xc,yn);
xn      = -dx/2:dx:Lx+dx/2; 
yn      = -dy/2:dy:Ly+dy/2;
zn      = -dz/2:dz:Lz+dy/2;
Vx      = zeros(nx+1,ny,nz);
Vy      = zeros(nx,ny+1,nz);
Vz      = zeros(nx,ny,nz+1);
P       = zeros(nx,ny,nz);
phi0    = zeros(nx,ny,nz);
Exyn    = zeros(nx,ny,nz);
Exzn    = zeros(nx,ny,nz);
Eyzn    = zeros(nx,ny,nz);
txyn    = zeros(nx,ny,nz);
txzn    = zeros(nx,ny,nz);
tyzn    = zeros(nx,ny,nz);
eta_xy  = zeros(nx-1,ny-1,nz-2);
eta_xz  = zeros(nx-1,ny-2,nz-1);
eta_yz  = zeros(nx-2,ny-1,nz-1);
etan    = ones(nx,ny,nz);
dVxdt   = zeros(nx-1,ny-2,nz-2);
dVydt   = zeros(nx-2,ny-1,nz-2); 
dVzdt   = zeros(nx-2,ny-2,nz-1);
P_BCX   = repmat(reshape(linspace(SRX, 0, nz), [1, 1, nz]), [nx, ny, 1]);%*p_NDHSM;
P_BCY   = repmat(reshape(linspace(SRY, 0, nz), [1, 1, nz]), [nx, ny, 1]);%*p_NDHSM;
for iter = 1:niter % Pseudo-Transient cycles
    % Timesteps - GPU KERNEL 1
    etax     =  0.5*(etan(2:end,2:end-1,2:end-1) + etan(1:end-1,2:end-1,2:end-1));  
    etay     =  0.5*(etan(2:end-1,2:end,2:end-1) + etan(2:end-1,1:end-1,2:end-1));
    etaz     =  0.5*(etan(2:end-1,2:end-1,2:end) + etan(2:end-1,2:end-1,1:end-1));
    dtP      = relP*6.1/max(nx,max(ny,nz))*etan*(1+eta_b);
    dtVx     = relV*min(dx,min(dy,dz))^2./etax/(1+eta_b)/6.1;
    dtVy     = relV*min(dx,min(dy,dz))^2./etay/(1+eta_b)/6.1;
    dtVz     = relV*min(dx,min(dy,dz))^2./etaz/(1+eta_b)/6.1;
    % Pressure update and strain rates - GPU KERNEL 2
    divV     = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy + diff(Vz,1,3)/dz;
    P        = P - divV.*dtP;
    Exx      = diff(Vx,1,1)/dx - 1/3*divV;
    Eyy      = diff(Vy,1,2)/dy - 1/3*divV;
    Ezz      = diff(Vz,1,3)/dz - 1/3*divV;
    Exy      = 0.5*(diff(Vx(2:end-1,:,2:end-1),1,2)/dy ...
                   +diff(Vy(:,2:end-1,2:end-1),1,1)/dx); 
    Exz      = 0.5*(diff(Vx(2:end-1,2:end-1,:),1,3)/dz ...
                   +diff(Vz(:,2:end-1,2:end-1),1,1)/dx);    
    Eyz      = 0.5*(diff(Vy(2:end-1,2:end-1,:),1,3)/dz ...
                   +diff(Vz(2:end-1,:,2:end-1),1,2)/dy);   
    % Stresses - GPU KERNEL 3      
    txx      =  2*etan.*(Exx + eta_b*divV);
    tyy      =  2*etan.*(Eyy + eta_b*divV);
    tzz      =  2*etan.*(Ezz + eta_b*divV);
    txy      =  2*eta_xy.*Exy;
    txz      =  2*eta_xz.*Exz; 
    tyz      =  2*eta_yz.*Eyz; 
    eta_xy   = 0.25*(etan(1:end-1,1:end-1,2:end-1) + etan(2:end,1:end-1,2:end-1) ...
                    +etan(1:end-1,2:end,2:end-1) + etan(2:end,2:end,2:end-1));
    eta_xz   = 0.25*(etan(1:end-1,2:end-1,1:end-1) + etan(2:end,2:end-1,1:end-1) ...
                    +etan(1:end-1,2:end-1,2:end) + etan(2:end,2:end-1,2:end));                
    eta_yz   = 0.25*(etan(2:end-1,1:end-1,1:end-1) + etan(2:end-1,2:end,1:end-1) ...
                    +etan(2:end-1,1:end-1,2:end) + etan(2:end-1,2:end,2:end));         
    Exyn(2:end-1,2:end-1,2:end-1) = 0.25*(Exy(1:end-1,1:end-1,:)+Exy(2:end,1:end-1,:) ...
                                 +Exy(1:end-1,2:end,:)+Exy(2:end,2:end,:));  
    Exzn(2:end-1,2:end-1,2:end-1) = 0.25*(Exz(1:end-1,:,1:end-1)+Exz(2:end,:,1:end-1) ...
                                 +Exz(1:end-1,:,2:end)+Exz(2:end,:,2:end));
    Eyzn(2:end-1,2:end-1,2:end-1) = 0.25*(Eyz(:,1:end-1,1:end-1)+Eyz(:,2:end,1:end-1) ...
                                 +Eyz(:,1:end-1,2:end)+Eyz(:,2:end,2:end));
    txyn(2:end-1,2:end-1,2:end-1) = 0.25*(txy(1:end-1,1:end-1,:)+txy(2:end,1:end-1,:) ...
                                 +txy(1:end-1,2:end,:)+txy(2:end,2:end,:));  
    txzn(2:end-1,2:end-1,2:end-1) = 0.25*(txz(1:end-1,:,1:end-1)+txz(2:end,:,1:end-1) ...
                                 +txz(1:end-1,:,2:end)+txz(2:end,:,2:end));
    tyzn(2:end-1,2:end-1,2:end-1) = 0.25*(tyz(:,1:end-1,1:end-1)+tyz(:,2:end,1:end-1) ...
                                 +tyz(:,1:end-1,2:end)+tyz(:,2:end,2:end));                         
    % GPU KERNEL 4                  
    Exyn(:,:,1) = Exyn(:,:,2);
    Exzn(:,:,1) = Exzn(:,:,2);
    Eyzn(:,:,1) = Eyzn(:,:,2);
    % GPU KERNEL 5
    ResVx       = diff(-P(:,2:end-1,2:end-1),1,1)/dx + PI2*diff(txx(:,2:end-1,2:end-1),1,1)/dx + PI2*diff(txy,1,2)/dy + PI2*diff(txz,1,3)/dz;
    ResVy       = diff(-P(2:end-1,:,2:end-1),1,2)/dy + PI2*diff(txy,1,1)/dx + PI2*diff(tyy(2:end-1,:,2:end-1),1,2)/dy + PI2*diff(tyz,1,3)/dz;
    ResVz       = diff(-P(2:end-1,2:end-1,:),1,3)/dz + PI2*diff(txz,1,1)/dx + PI2*diff(tyz,1,2)/dy + PI2*diff(tzz(2:end-1,2:end-1,:),1,3)/dz - PI2*(rhoS/rho_sc);
    dVxdt       =  dVxdt.*(1-damp/nx) + ResVx;
    dVydt       =  dVydt.*(1-damp/ny) + ResVy;   
    dVzdt       =  dVzdt.*(1-damp/nz) + ResVz;
    EII2        = (Exx.^2 + Eyy.^2 + Ezz.^2)/2 + Exyn.^2 + Exzn.^2 + Eyzn.^2;
    tII         = sqrt((txx.^2 + tzz.^2)/2 + txzn.^2);
    etait       = (A/eta_sc).*(rho_sc*g*L_sc.*tII).^(1-n_pow);
    % GPU KERNEL 6
    Vx(2:end-1,2:end-1,2:end-1) = Vx(2:end-1,2:end-1,2:end-1) + dVxdt.*dtVx;
    Vy(2:end-1,2:end-1,2:end-1) = Vy(2:end-1,2:end-1,2:end-1) + dVydt.*dtVy;
    Vz(2:end-1,2:end-1,2:end-1) = Vz(2:end-1,2:end-1,2:end-1) + dVzdt.*dtVz;
    etan              = min(exp(rele*log(etait) + (1-rele)*log(etan)),eta0);
    % GPU KERNEL 7
    etan(:,:,[1 end]) = etan(:,:,[2 end-1]);
    % GPU KERNEL 8
    etan([1 end],:,:) = etan([2 end-1],:,:);
    % GPU KERNEL 9
    etan(:,[1 end],:) = etan(:,[2 end-1],:); 
    % GPU KERNEL 10
    % BC left
    Vx(1,2:end-1,2:end-1)   =  Vx(2,2:end-1,2:end-1) + dx*(p_NDHSM*P_BCX(1,2:end-1,2:end-1)-P(1,2:end-1,2:end-1))./etan(1,2:end-1,2:end-1)/2 + dx*(1/3-eta_b)*divV(1,2:end-1,2:end-1);
    Vy(1,2:end-1,2:end-1)   =  Vy(2,2:end-1,2:end-1);
    Vz(1,2:end-1,2:end-1)   =  Vz(2,2:end-1,2:end-1);
    % BC right
    Vx(end,2:end-1,2:end-1) =  Vx(end-1,2:end-1,2:end-1) - dx*(p_NDHSM*P_BCX(end,2:end-1,2:end-1)-P(end,2:end-1,2:end-1))./etan(end,2:end-1,2:end-1)/2 + dx*(1/3-eta_b)*divV(end,2:end-1,2:end-1);
    Vy(end,2:end-1,2:end-1) =  Vy(end-1,2:end-1,2:end-1);
    Vz(end,2:end-1,2:end-1) =  Vz(end-1,2:end-1,2:end-1);
    % GPU KERNEL 11
    % BC front
    Vx(2:end-1,1,2:end-1)   =  Vx(2:end-1,2,2:end-1);
    Vy(2:end-1,1,2:end-1)   =  Vy(2:end-1,2,2:end-1) + dy*(p_NDHSM*P_BCY(2:end-1,1,2:end-1)-P(2:end-1,1,2:end-1))./etan(2:end-1,1,2:end-1)/2 + dy*(1/3-eta_b)*divV(2:end-1,1,2:end-1);
    Vz(2:end-1,1,2:end-1)   =  Vz(2:end-1,2,2:end-1);
    % BC back
    Vx(2:end-1,end,2:end-1) =  Vx(2:end-1,end-1,2:end-1);
    Vy(2:end-1,end,2:end-1) =  Vy(2:end-1,end-1,2:end-1) - dy*(p_NDHSM*P_BCY(2:end-1,end,2:end-1)-P(2:end-1,end,2:end-1))./etan(2:end-1,end,2:end-1)/2 + dy*(1/3-eta_b)*divV(2:end-1,end,2:end-1);
    Vz(2:end-1,end,2:end-1) =  Vz(2:end-1,end-1,2:end-1);
    % GPU KERNEL 12
    % BC corner
    Vx([1 end],1,2:end-1)   =  Vx([1 end],2,2:end-1);
    Vx([1 end],end,2:end-1) =  Vx([1 end],end-1,2:end-1);
    Vy(1,[1 end],2:end-1)   =  Vy(2,[1 end],2:end-1);
    Vy(end,[1 end],2:end-1) =  Vy(end-1,[1 end],2:end-1);
    Vz(1,1,2:end-1)         =  0.5* ( Vz(1,2,2:end-1)       + Vz(2,1,2:end-1)       );
    Vz(1,end,2:end-1)       =  0.5* ( Vz(1,end-1,2:end-1)   + Vz(2,end,2:end-1)     );
    Vz(end,1,2:end-1)       =  0.5* ( Vz(end,2,2:end-1)     + Vz(end-1,1,2:end-1)   );
    Vz(end,end,2:end-1)     =  0.5* ( Vz(end,end-1,2:end-1) + Vz(end-1,end,2:end-1) );
    % BC bottom
    Vx(:,:,1)     =  Vx(:,:,2);
    Vy(:,:,1)     =  Vy(:,:,2);
    Vz(:,:,1)     =  0.0;
    % BC top             
    Vz(:,:,end)   =  Vz(:,:,end-1) + dz*P(:,:,end)./etan(:,:,end)/2 + dz*(1/3-eta_b)*divV(:,:,end); 
    Vx(:,:,end)   =  Vx(:,:,end-1);
    Vy(:,:,end)   =  Vy(:,:,end-1);
    %% errors
    norm_X  = sqrt(sum(dVxdt(:).*dVxdt(:)))/(nx*ny*nz);
    norm_Y  = sqrt(sum(dVydt(:).*dVydt(:)))/(nx*ny*nz);
    norm_Z  = sqrt(sum(dVzdt(:).*dVzdt(:)))/(nx*ny*nz);
    norm_P  = sqrt(sum(divV(:).*divV(:)))/(nx*ny*nz);
    % Monitoring
    if (mod(iter,100)==1) 
    iterror = max(max(norm_X,norm_Y),max(norm_Z,norm_P));
    if ((iterror < epsi) && (iter > 2)),break; end
    end
    if (mod(iter,nout_iter)==1) 
        fprintf('iter=%d, err=%1.3e \n',iter,iterror)
    end
end
CPU_time=toc;
MTP_eff = (nx*ny*nz)*iter*(2*4+2)*8/1e9/CPU_time;
fprintf('iter=%d, err=%1.3e \n',iter,iterror)
fprintf('time=%f, MTP_eff=%f \n',CPU_time,MTP_eff)

% Plots
imagesc(reshape(P(:,25,:),[nx,nz]))
% % % Output to Paraview
% % vtkwrite(fullfile('Output',sprintf('Pressure.vtk')),'structured_points', 'Pressure', P);
% % vtkwrite(fullfile('Output',sprintf('Viscosity.vtk')),'structured_points', 'Viscosity', etan);
% % vtkwrite(fullfile('Output',sprintf('stress_xx.vtk')),'structured_points', 'Tau_xx', txx);
% % vtkwrite(fullfile('Output',sprintf('stress_yy.vtk')),'structured_points', 'Tau_yy', tyy);
% % vtkwrite(fullfile('Output',sprintf('stress_zz.vtk')),'structured_points', 'Tau_zz', tzz);
% % vtkwrite(fullfile('Output',sprintf('stress_xy.vtk')),'structured_points', 'Tau_xy', txy);
% % vtkwrite(fullfile('Output',sprintf('stress_xz.vtk')),'structured_points', 'Tau_xz', txz);
% % vtkwrite(fullfile('Output',sprintf('stress_yz.vtk')),'structured_points', 'Tau_yz', tyz);
% % vtkwrite(fullfile('Output',sprintf('Vx.vtk')),'structured_points', 'Vx', Vx);
% % vtkwrite(fullfile('Output',sprintf('Vy.vtk')),'structured_points', 'Vy', Vy);
% % vtkwrite(fullfile('Output',sprintf('Vz.vtk')),'structured_points', 'Vz', Vz);

% xc_MAT   = xc;
% xn_MAT   = xn;
% yc_MAT   = yc;
% yn_MAT   = yn;
% zc_MAT   = zc;
% zn_MAT   = zn;
% Vx_MAT = Vx;
% Vy_MAT = Vy;
% Vz_MAT = Vz;
% P_MAT = P;
% txx_MAT  = txx;
% tyy_MAT  = tyy;
% tzz_MAT  = tzz;
% txy_MAT  = txy;
% tyz_MAT  = tyz;
% txz_MAT  = txz;
% etan_MAT = etan;
% save('data_MAT_3D_4.mat','xc_MAT','xn_MAT','yc_MAT','yn_MAT','zc_MAT','zn_MAT','Vx_MAT','Vy_MAT','Vz_MAT','P_MAT','txx_MAT','tyy_MAT','tzz_MAT','txy_MAT','tyz_MAT','txz_MAT','etan_MAT')