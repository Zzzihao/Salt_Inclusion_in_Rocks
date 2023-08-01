% This is the matlab code I wrote and translated based on the julia code by
% Ludovic Raess which is exactly the same as the original author's code

clear
tic
% Physics
lx = 10.0;          % domain extends
ly = 10.0;
mus0 = 1.0;         % matrix viscosity
musi = 1e3;        % inclusion viscosity
epsilon_bg = 1.0;   % background strain-rate
% Numerics
iterMax = 1e5;      % maximum number of pseudo-transient iterations
nout = 500;         % error checking frequency
epsilon = 1e-8;     % nonlinear absolute tolerance
CFL = 0.9/sqrt(2);
Re = 5 * pi;
r = 1.0;
nx = 99;
ny = 99;
dx = lx / nx;      % cell sizes
dy = ly / ny;
max_lxy = max(lx, ly);
Vpdt = min(dx, dy) * CFL;
xc = linspace(dx / 2, lx - dx / 2, nx); % x coordinates
yc = linspace(dy / 2, ly - dy / 2, ny); % y coordinates
yv = linspace(0, ly, ny + 1);           % y coordinates for velocity

% Array allocations
Pt = zeros(nx, ny);
divV = zeros(nx, ny);
txx = zeros(nx, ny);
tyy = zeros(nx, ny);
txy = zeros(nx - 1, ny - 1);
Rx = zeros(nx - 1, ny - 2);
Ry = zeros(nx - 2, ny - 1);
dVx = zeros(nx - 1, ny - 2);
dVy = zeros(nx - 2, ny - 1);
Mus2 = zeros(nx, ny);
Must = zeros(nx, ny);
Gdt = zeros(nx, ny);
dt_Rho = zeros(nx, ny);

% Initial conditions
Rad2 = zeros(nx, ny);
Vx = zeros(nx + 1, ny);
Vy = zeros(nx, ny + 1);
for ix = 1:nx
    for iy = 1:ny
        Rad2(ix, iy) = ((ix - 1) * dx + 0.5 * dx - 0.5 * lx)^2 + ((iy - 1) * dy + 0.5 * dy - 0.5 * ly)^2;
    end
end

for ix = 1:nx + 1
    for iy = 1:ny
        Vx(ix, iy) = -epsilon_bg * ((ix - 1) * dx - 0.5 * lx);
    end
end
for ix = 1:nx
    for iy = 1:ny + 1
        Vy(ix, iy) = epsilon_bg * ((iy - 1) * dy - 0.5 * ly);
    end
end
Mus = mus0 * ones(nx, ny);
Mus(Rad2 < 1.0) = musi;
Mus2 = Mus;

for ism = 1:10
    Mus2 = smooth(Mus, 1.0);
    Mus_temp = Mus;
    Mus = Mus2;
end

Mus_t = Mus;
Mus_t = compute_maxloc(Mus);

Mus_t(1, :) = Mus_t(2, :);
Mus_t(end, :) = Mus_t(end - 1, :);
Mus_t(:, 1) = Mus_t(:, 2);
Mus_t(:, end) = Mus_t(:, end - 1);

dt_Rho = Vpdt * max_lxy / Re./ Mus_t;
Gdt = Vpdt^2 ./ dt_Rho / (r + 2.0);

    
err = 2 * epsilon;
iter = 0;
err_evo1 = [];
err_evo2 = [];

while err > epsilon && iter <= iterMax
    if iter == 11
        tic;
    end
    divV = diff(Vx, 1, 1) / dx + diff(Vy, 1, 2) / dy;
    
    Pt = Pt - r .* Gdt .* divV;
    Exz      = 0.5*(diff(Vx(2:end-1,:),1,2)/dy+diff(Vy(:,2:end-1),1,1)/dx);  
    
    txx = (txx + 2.0 * Gdt .* diff(Vx, 1, 1) / dx) ./ (Gdt ./ Mus + 1.0);
    tyy = (tyy + 2.0 * Gdt .* diff(Vy, 1, 2) / dy) ./ (Gdt ./ Mus + 1.0);
    txy = (txy + 2.0 * XZ_mean(Gdt) .* Exz) ./ ...
             ( ( XZ_mean(Gdt) ./ XZ_mean(Mus) ) + 1.0 );
    
         
    Rx = -diff(Pt(:,2:end-1),1,1)/dx + diff(txy,1,2)/dy + diff(txx(:,2:end-1),1,1)/dx;
    Ry = -diff(Pt(2:end-1,:),1,2)/dy + diff(txy,1,1)/dx + diff(tyy(2:end-1,:),1,2)/dy;
    
    dVx = (dt_Rho(1:end-1, 2:end-1) + dt_Rho(2:end, 2:end-1)) / 2 .* Rx;
    dVy = (dt_Rho(2:end-1, 1:end-1) + dt_Rho(2:end-1, 2:end)) / 2 .* Ry;
    
    Vx(2:end-1, 2:end-1) = Vx(2:end-1, 2:end-1) + dVx;
    Vy(2:end-1, 2:end-1) = Vy(2:end-1, 2:end-1) + dVy;
    

    Vx(:, 1) = Vx(:, 2);
    Vx(:, end) = Vx(:, end - 1);

    Vy(1, :) = Vy(2, :);
    Vy(end, :) = Vy(end - 1, :);

    iter = iter + 1;
    
    if mod(iter, nout) == 0
        Vmin = min(Vx(:));
        Vmax = max(Vx(:));
        Pmin = min(Pt(:));
        Pmax = max(Pt(:));
        norm_Rx = norm(Rx) / (Pmax - Pmin) * lx / sqrt(numel(Rx));
        norm_Ry = norm(Ry) / (Pmax - Pmin) * lx / sqrt(numel(Ry));
        norm_dV = norm(divV) / (Vmax - Vmin) * lx / sqrt(numel(divV));
        err = max([norm_Rx, norm_Ry, norm_dV]);
        err_evo1(end+1) = err;
        err_evo2(end+1) = iter;
        fprintf('Total steps = %d, err = %.3e [norm_Rx=%.3e, norm_Ry=%.3e, norm_dV=%.3e]\n', iter, err, norm_Rx, norm_Ry, norm_dV);
    end

    
    
end
toc













figure(1)
        clf
        colormap jet
        subplot(3,2,1)
        imagesc(flipud(txx'))
        colorbar
        subplot(3,2,2)
        imagesc(flipud(((tyy'))))
        colorbar
        subplot(3,2,3)
        imagesc(flipud(txy'))
        colorbar
        subplot(3,2,4)
        imagesc(flipud(((Pt'))))
        colorbar
        subplot(3,2,5)
        imagesc(flipud(((Vx'))))
        colorbar
        subplot(3,2,6)
        imagesc(flipud(((Vy'))))
        colorbar
        drawnow
        set(gca,'color','white')