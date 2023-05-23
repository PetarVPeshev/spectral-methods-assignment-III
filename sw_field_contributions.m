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
Nz = 401;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% CYLINDRICAL GRID
rho = linspace(0.2, 15, Nrho) * wave.wavelength / sqrt(dielectric.er);
phi = linspace(0, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);
z = linspace(eps, dielectric.h + 2e-3, Nz);

%% TE1 AND TM0 PROPAGATION CONSTANTS
krho_range = wave.k0 * linspace(1, sqrt(dielectric.er), 1001);
[krho_te, krho_tm] = find_krho(wave.k0, krho_range, ...
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

%% TOTAL ELECTRIC FIELD IN Z
Etotal = sqrt(abs(E(:, :, 1, :)) .^ 2 + abs(E(:, :, 2, :)) .^ 2 ...
    + abs(E(:, :, 3, :)) .^ 2);
Etotal = squeeze(Etotal);

%% PLOT
figure('Position', [250 250 750 400]);
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    real(E(1, :, 1, ceil(Nz / 4))), 'LineWidth', 2.0, ...
    'DisplayName', '\Re\{E_{\rho}\}');
hold on;
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    imag(E(1, :, 1, ceil(Nz / 4))), '--', 'LineWidth', 2.0, ...
    'DisplayName', '\Im\{E_{\rho}\}');
grid on;
xticks(0 : 3 : 15);
ylim([-7 7] * 1e5);
xlim([0.2 15]);
legend show;
legend('location', 'bestoutside');
xlabel('\rho / \lambda_{d}');
ylabel('E_{\rho}^{TM} / V/m');
title(['TM E_{\rho} Real & Imaginary Parts @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and z = ' ...
    num2str(round(z(ceil(Nz / 4)) * 1e3, 2)) ' mm']);
saveas(gcf, 'figures\real_and_imag_variation.fig');

figure('Position', [250 250 750 400]);
plot(z' * 1e3, squeeze(Etotal(1, end, :)), 'LineWidth', 2.0, ...
    'DisplayName', 'E_{total}');
hold on;
xline(dielectric.h * 1e3, '--', 'Color', [0.9290 0.6940 0.1250], ...
    'LineWidth', 2.0, 'DisplayName', 'interface');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('z / mm');
ylabel('|E_{total}^{TM}| / V/m');
title(['TM |E_{total}| Amplitude @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and \rho = 15\lambda_{d}']);
saveas(gcf, 'figures\amplitude_z.fig');

figure('Position', [250 250 750 400]);
plot(phi' * 180 / pi, squeeze(Etotal(:, end, ceil(Nz / 4))), ...
    'LineWidth', 2.0, 'DisplayName', 'E_{total}');
grid on;
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\phi / deg');
ylabel('|E_{total}^{TM}| / V/m');
title(['TM |E_{total}| Amplitude @ \rho = 15\lambda_{d}, and h = ' ...
    num2str(round(z(ceil(Nz / 4)) * 1e3, 2)) ' mm']);
saveas(gcf, 'figures\amplitude_phi.fig');
