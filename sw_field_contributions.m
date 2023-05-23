close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
wave.f = 10 * 1e9;
dielectric.h  = 2e-3;
dielectric.er = 10;
dipole.W = 0.5 * 1e-3;
dipole.L = 15 * 1e-3;
Nrho = 200;
Nphi = 200;
Nz = 201;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% CYLINDRICAL GRID
rho = linspace(eps, 30 * 1e-3, Nrho);
phi = linspace(0, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);
z = linspace(0, dielectric.h, Nz);

%% TE1 AND TM0 PROPAGATION CONSTANTS
[krho_te, krho_tm] = find_krho(wave.k0, wave.k0 * linspace(1, sqrt(dielectric.er), 1001), ...
        'GroundSlab', dielectric.h, dielectric.er);

%% TM WAVE VECTOR COMPONENTS IN CARTESIAN
kx = krho_tm * cos(cyl_grid(:, :, 2));
ky = krho_tm * sin(cyl_grid(:, :, 2));

k_comp = NaN( [size(cyl_grid, 1, 2), 2] );
k_comp(:, :, 1) = kx;
k_comp(:, :, 2) = ky;

%% RESIDUE STRATIFIED MEDIA
[~, ~, v_tm, i_tm] = residue_stratified(wave.k0 * sqrt(dielectric.er), krho_te, krho_tm, z, ...
    'GroundSlab', dielectric.h, dielectric.er);

%% DIPOLE CURRENT FOURIER TRANSFORM
J = ft_current(wave.k0, k_comp, dipole.W, dipole.L, dielectric.er, ...
    'dipole', 'x');

%% SURFACE WAVE ELECTRIC FIELD
E = NaN( [size(cyl_grid, 1, 2), 3, length(z)] );
for z_idx = 1 : 1 : length(z)
    E(:, :, :, z_idx) = sw_fields(wave.k0 * sqrt(dielectric.er), krho_tm, v_tm(z_idx), ...
        i_tm(z_idx), J, dielectric.er, cyl_grid, 'TM');
end

%% PLOT
figure('Position', [250 250 750 400]);
plot(rho * 1e3, real(E(1, :, 1, ceil(Nz / 2))), 'LineWidth', 2.0, ...
    'DisplayName', '\Re\{E_{\rho}\}');
hold on;
plot(rho * 1e3, imag(E(1, :, 1, ceil(Nz / 2))), '--', 'LineWidth', 2.0, ...
    'DisplayName', '\Im\{E_{\rho}\}');
grid on;
xlabel('\rho / mm');
ylabel('E_{\rho}');
title(['E_{\rho} Real & Imaginary Parts @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and z = ' ...
    num2str(round(z(ceil(Nz / 2)) * 1e3, 2)) ' mm']);
saveas(gcf, 'figures\real_and_imag_variation.fig');
