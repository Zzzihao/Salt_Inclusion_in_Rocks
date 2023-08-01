option = 2;
% Symbol: Tension (positive) to compression (positive)
Sxx = p_water-txx(2:end-1,2:end-1) + P(2:end-1,2:end-1);
Sxz = -0.25* (txz(1:end-1,1:end-1) + txz(2:end,1:end-1) + txz(1:end-1,2:end) + txz(2:end,2:end));
Szz = p_water-tzz(2:end-1,2:end-1) + P(2:end-1,2:end-1);
% Calculate principal stress
for i = 1:nx-2
    for j = 1:nz-2
        [V, D] = eig([Sxx(i,j), Sxz(i,j); Sxz(i,j), Szz(i,j)]);
        eigenvalues = diag(D);
        S1(i,j) = max(eigenvalues) * p_sc;
        S3(i,j) = min(eigenvalues) * p_sc;
        [~, max_index] = max(eigenvalues);
        S1_direction = V(:, max_index);
        S1_direction = S1_direction / norm(S1_direction);
        S1_directions{i,j} = S1_direction;
    end
end

U = zeros(nx-2, nz-2); 
V = zeros(nx-2, nz-2); 
for i = 1:nx-2
    for j = 1:nz-2
        U(i, j) = S1_directions{i, j}(1);
        V(i, j) = S1_directions{i, j}(2);
    end
end
figure(1);
[xc3,  zc3] = ndgrid(1:nx-2,1:nz-2);
quiver(xc3, zc3, U, V,0.5, 'color', 'k','linewidth',1); % ????
xlabel('X');
ylabel('Z');
hold on;
SaltFraction2 = SaltFraction(2:end-1,2:end-1);
saltBoundary = bwperim(SaltFraction2); 
[saltBoundaryRows, saltBoundaryCols] = find(saltBoundary);
saltBoundaryY = saltBoundaryCols;
saltBoundaryX = saltBoundaryRows;
scatter(saltBoundaryX, saltBoundaryY, 'r.');
hold off;

if option == 1
% Mohr-Column failure criterion
phi_sh    = 30;                            % Shale
C0_sh   = 1e6;
phi_sa    = 30;                            % Salt
C0_sa   = 1e7;

% Calculation
%S3(find(S3<0)) = 0;
S1_HB_sh = C0_sh + S3 * tand(phi_sh);
S1_HB_sa = C0_sa + S3 * tand(phi_sa);
SaltFraction2 = SaltFraction(2:end-1,2:end-1);
S1_HB = S1_HB_sh;
S1_HB(SaltFraction2 > 0.5) = S1_HB_sa(SaltFraction2 > 0.5);
FC = S1-S1_HB; % failure if FC > 0
figure(2);
addpath('C:\Users\Zihao\OneDrive\Documents\AAA-Postdoc\BrewerMap-master');
colors = brewermap(64, 'YlGnBu');
colormap(colors);
imagesc(flipud(FC'));colorbar;
clim = caxis;
caxis([0, clim(2)]);

hold on;
saltBoundary = bwperim(SaltFraction2);
[saltBoundaryRows, saltBoundaryCols] = find(saltBoundary);
saltBoundaryY = size(SaltFraction2, 2) - saltBoundaryCols + 1;
saltBoundaryX = saltBoundaryRows;
scatter(saltBoundaryX, saltBoundaryY, 'r.');
hold off;
set(gca,'color','white');
end

if option == 2
% Hoek-Brown failure criterion
m_sh    = 10.0;                            % Shale
s_sh    = 0.5;
C0_sh   = 6e5;
m_sa    = 5.0;                            % Salt
s_sa    = 0.3;
C0_sa   = 3e5;

% Calculation
S3(find(S3<0)) = 0;
S1_HB_sh = S3 + C0_sh*sqrt(m_sh*S3/C0_sh + s_sh);
S1_HB_sa = S3 + C0_sa*sqrt(m_sa*S3/C0_sa + s_sa);
SaltFraction2 = SaltFraction(2:end-1,2:end-1);
S1_HB = S1_HB_sh;
S1_HB(SaltFraction2 > 0.5) = S1_HB_sa(SaltFraction2 > 0.5);
FC = (S1-S1_HB)/1000000; % failure if FC > 0, [MPa]
figure(2);
addpath('C:\Users\Zihao\OneDrive\Documents\AAA-Postdoc\BrewerMap-master');
colors = brewermap(64, 'YlGnBu');
colormap(colors);
imagesc(flipud(FC'));colorbar;
clim = caxis;
%caxis([0, clim(2)]);
caxis([0, 35]);

hold on;
saltBoundary = bwperim(SaltFraction2);
[saltBoundaryRows, saltBoundaryCols] = find(saltBoundary);
saltBoundaryY = size(SaltFraction2, 2) - saltBoundaryCols + 1;
saltBoundaryX = saltBoundaryRows;
scatter(saltBoundaryX, saltBoundaryY, 'r.');
hold off;
set(gca,'color','white');

end