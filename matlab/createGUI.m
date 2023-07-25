% The purpose of this code is to create a software interface for solving 
% solid deformation problems, such as glaciers and rocks.

function createGUI()
    % Create uifigure
    fig = uifigure('Name', 'GUI', 'Position', [100, 100, 450, 450]);

    % Create an input parameter button and set a callback function
    btn = uibutton(fig, 'Text', 'Please input parameters', 'Position', ...
        [175, 20, 200, 150]);
    btn.ButtonPushedFcn = @getInput;
end


function getInput(src, event)
    % Ask for input parameters
    prompt = {'Salt viscosity [Pa-s]', 'Shale viscosity [Pa-s]', ...
        'Gravity [m/s2]', 'Horizontal stress gradient [Pa/m]', ...
        'Vertical stress gradient [Pa/m]', ...
        'Salt density [kg/m3]', 'Seawater density [kg/m3]', ...
        'Seabed depth [m]', 'Simulation domain depth [m]', ...
        'Aspect'};
    title = 'Input parameters';
    dims = [1 50]; % Size of inputs
    definput = {'1e17', '1e15', '9.81', '16965', '19228', ...
                '2200', '1020', '800', '10000', '3'}; % Default Values
    answer = inputdlg(prompt, title, dims, definput);
    
    % Inputs
    if ~isempty(answer)
        A_Salt = str2double(answer{1});
        A_Shale = str2double(answer{2});
        g = str2double(answer{3});
        p_grad_h = str2double(answer{4});
        p_grad_v = str2double(answer{5});
        rho_Salt = str2double(answer{6});
        rho_Water = str2double(answer{7});
        D_sf = str2double(answer{8});
        Lz = str2double(answer{9});
        aspect = str2double(answer{10});
        
        variables = {'A_Salt', 'A_Shale', 'g', 'p_grad_h', 'p_grad_v', ...
                     'rho_Salt', 'rho_Water', 'D_sf', 'Lz', 'aspect'};
        for i = 1:length(variables)
            assignin('base', variables{i}, str2double(answer{i}));
        end
        % Calculations
        rho_Shale = p_grad_v/g;
        rho_Shale
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
        nx          = 50;                      % Number of blocks in X-direction
        nz          = 50;                      % Number of blocks in Z-direction
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
        img         = imread('SaltGeometry2D.JPG'); % Load the salt geometry
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
            RX    = Z_mean(-diff(P,1,1)/dx) + X_mean(PI2*diff(txz,1,2)/dz) + Z_mean(PI2*diff(txx,1,1)/dx);
            RZ    = X_mean(-diff(P,1,2)/dz) + Z_mean(PI2*diff(txz,1,1)/dx) + X_mean(PI2*diff(tzz,1,2)/dz) - PI2*(rho_Rockz/rho_sc);
            dVx     = XZ_mean(dt_Rho) .*RX;
            dVz     = XZ_mean(dt_Rho) .*RZ;
            % Update velocity
            Vx(2:end-1, 2:end-1) = Vx(2:end-1, 2:end-1) + dVx;
            Vz(2:end-1, 2:end-1) = Vz(2:end-1, 2:end-1) + dVz;
            % BC bottom
            Vz(2:end-1,1)     = Vz(2:end-1,2) - dz*(-p_water-p_NVP+X_mean(P(:,1)))./X_mean(etan(:,1))/2;
            % BC top
            Vz(2:end-1,end)   = Vz(2:end-1,end-1) + dz*(-p_water+X_mean(P(:,end)))./X_mean(etan(:,end))/2;
            % BC left
            MM = linspace(1,0,nz+1);
            Vx(1,2:end-1)     = Vx(2,2:end-1) - dx*( -p_NMHP*MM(2:end-1)+Z_mean(P(1,:)) )/2./Z_mean(etan(1,:));
            % BC right
            Vx(end,2:end-1)   = Vx(end-1,2:end-1) + dx*( -p_NMHP*MM(2:end-1)+Z_mean(P(end,:)) )/2./Z_mean(etan(end,:));
    
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
                figure(1);
                clf; 
                colormap jet;
                h1 = subplot(3,3,1);
                imagesc(flipud(txx'));
                colorbar;
                subplot(3,3,2);
                imagesc(flipud(((tzz'))));
                colorbar;
                subplot(3,3,3);
                imagesc(flipud(txz'));
                colorbar;
                subplot(3,3,4);
                imagesc(flipud(P'));
                colorbar;
                subplot(3,3,5);
                imagesc(flipud(Vx'));
                colorbar;
                subplot(3,3,6);
                imagesc(flipud(Vz'));
                colorbar;
                subplot(3,3,7);
                imagesc(flipud((txx+P)'));
                colorbar;
                subplot(3,3,8);
                imagesc(flipud((tzz+P)'));
                colorbar;
                subplot(3,3,9);
                imagesc(flipud(etan'));
                colorbar;
                drawnow;
                set(gca,'color','white');
            end
        end
    else
        disp('The user did not input any parameters');
    end
end
