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
Nz = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c / wave.f;
wave.k0 = 2 * pi / wave.wavelength;

%% CYLINDRICAL GRID
rho = linspace(0.2, 15, Nrho) * wave.wavelength / sqrt(dielectric.er);
phi = linspace(0, 2 * pi, Nphi);
cyl_grid = meshgrid_comb(rho, phi);
z = linspace(eps, dielectric.h + 23e-3, Nz);

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
[~, ~, v_tm, i_tm] = residue_stratified(wave.k0, krho_te, krho_tm, z, ...
    'GroundSlab', dielectric.h, dielectric.er);

%% DIPOLE CURRENT FOURIER TRANSFORM
J = ft_current(wave.k0, k_comp, dipole.W, dipole.L, dielectric.er, ...
    'dipole', 'x');

E = NaN( [size(cyl_grid, 1, 2), 3, length(z)] );
Etotal = NaN( [size(cyl_grid, 1, 2), length(z)] );
for z_idx = 1 : 1 : length(z)
    %% SURFACE WAVE ELECTRIC FIELD
    % Continuity across the interface issue
%     if z(z_idx) <= dielectric.h
        E(:, :, :, z_idx) = sw_fields(wave.k0, krho_tm, v_tm(z_idx), ...
            i_tm(z_idx), J, dielectric.er, cyl_grid, 'TM');
%     else
%         E(:, :, :, z_idx) = sw_fields(wave.k0, krho_tm, v_tm(z_idx), ...
%             i_tm(z_idx), J, 1, cyl_grid, 'TM');
%     end

    %% TOTAL ELECTRIC FIELD IN Z
    Etotal(:, :, z_idx) = sqrt(abs(E(:, :, 1, z_idx)) .^ 2 ...
        + abs(E(:, :, 2, z_idx)) .^ 2 + abs(E(:, :, 3, z_idx)) .^ 2);
end

% Etotal = sqrt(abs(E(:, :, 1, :)) .^ 2 + abs(E(:, :, 2, :)) .^ 2 ...
%     + abs(E(:, :, 3, :)) .^ 2);
% Etotal = squeeze(Etotal);

%% PLOT
cyl_envelope = sqrt(krho_tm / (2 * pi)) ./ sqrt(rho);
% Real and imaginary in radial distance
figure('Position', [250 250 750 500]);
subplot(2, 1, 1);
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    real(E(1, :, 1, 41)), 'LineWidth', 2.0, ...
    'DisplayName', '\Re\{E_{\rho}\}');
hold on;
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    imag(E(1, :, 1, 41)), '--', 'LineWidth', 2.0, ...
    'DisplayName', '\Im\{E_{\rho}\}');
hold on;
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    420.9176 * cyl_envelope, '-.', 'Color', [0.4660 0.6740 0.1880], ...
    'LineWidth', 2.0, 'DisplayName', 'cyl envelope')
grid on;
xticks(0 : 3 : 15);
ylim([-5 5] * 1e4);
xlim([0.2 15]);
legend show;
legend('location', 'bestoutside');
ylabel('E_{\rho}^{TM} / V/m');
subplot(2, 1, 2);
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    real(E(1, :, 3, 41)), 'LineWidth', 2.0, ...
    'DisplayName', '\Re\{E_{z}\}');
hold on;
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    imag(E(1, :, 3, 41)), '--', 'LineWidth', 2.0, ...
    'DisplayName', '\Im\{E_{z}\}');
hold on;
plot(rho * sqrt(dielectric.er) / wave.wavelength, ...
    242.7672 * cyl_envelope, '-.', 'Color', [0.4660 0.6740 0.1880], ...
    'LineWidth', 2.0, 'DisplayName', 'cyl envelope')
grid on;
xticks(0 : 3 : 15);
ylim([-5 5] * 1e4);
xlim([0.2 15]);
legend show;
legend('location', 'bestoutside');
xlabel('\rho / \lambda_{d}');
ylabel('E_{z}^{TM} / V/m');
sgtitle(['TM Real & Imaginary Parts @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and z = ' ...
    num2str(round(z(41) * 1e3, 2)) ' mm'], ...
    'FontWeight', 'bold', 'FontSize', 11);
saveas(gcf, 'figures\E_rho_variation.fig');

% Total E in z
figure('Position', [250 250 750 400]);
plot(z' * 1e3, squeeze(Etotal(1, end, :)), 'LineWidth', 2.0, ...
    'DisplayName', 'E_{total}');
hold on;
xline(dielectric.h * 1e3, '--', 'Color', [0.4940 0.1840 0.5560], ...
    'LineWidth', 2.0, 'DisplayName', 'interface');
grid on;
xlim([min(z) max(z)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('z / mm');
ylabel('|E^{TM}| / V/m');
title(['TM |E| Amplitude @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and \rho = 15\lambda_{d}']);
saveas(gcf, 'figures\Et_z_variation.fig');

% Total in z
figure('Position', [250 250 750 500]);
subplot(2, 1, 1);
plot(z' * 1e3, squeeze(abs(E(1, end, 1, :))), 'LineWidth', 2.0, ...
    'DisplayName', '|E_{\rho}|');
hold on;
xline(dielectric.h * 1e3, '--', 'Color', [0.4940 0.1840 0.5560], ...
    'LineWidth', 2.0, 'DisplayName', 'interface');
grid on;
xlim([min(z) max(z)] * 1e3);
legend show;
legend('location', 'bestoutside');
ylabel('|E_{\rho}^{TM}| / V/m');
subplot(2, 1, 2);
plot(z' * 1e3, squeeze(abs(E(1, end, 2, :))), 'LineWidth', 2.0, ...
    'DisplayName', '|E_{z}|');
hold on;
xline(dielectric.h * 1e3, '--', 'Color', [0.4940 0.1840 0.5560], ...
    'LineWidth', 2.0, 'DisplayName', 'interface');
grid on;
xlim([min(z) max(z)] * 1e3);
legend show;
legend('location', 'bestoutside');
xlabel('z / mm');
ylabel('|E_{z}^{TM}| / V/m');
sgtitle(['TM Real & Imaginary Parts @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and \rho = ' ...
    num2str(rho(end) * sqrt(dielectric.er) / wave.wavelength) ...
    ' \lambda_{d}'], 'FontWeight', 'bold', 'FontSize', 11);
saveas(gcf, 'figures\Eabs_z_variation.fig');

% Real and imaginary in z
figure('Position', [250 250 750 500]);
subplot(2, 1, 1);
plot(z' * 1e3, squeeze(real(E(1, end, 1, :))), 'LineWidth', 2.0, ...
    'DisplayName', '\Re\{E_{z}\}');
hold on;
plot(z' * 1e3, squeeze(imag(E(1, end, 1, :))), '--', 'LineWidth', 2.0, ...
    'DisplayName', '\Im\{E_{z}\}');
hold on;
xline(dielectric.h * 1e3, '--', 'Color', [0.4940 0.1840 0.5560], ...
    'LineWidth', 2.0, 'DisplayName', 'interface');
grid on;
legend show;
legend('location', 'bestoutside');
ylabel('E_{\rho}^{TM} / V/m');
subplot(2, 1, 2);
plot(z' * 1e3, squeeze(real(E(1, end, 2, :))), 'LineWidth', 2.0, ...
    'DisplayName', '\Re\{E_{z}\}');
hold on;
plot(z' * 1e3, squeeze(imag(E(1, end, 2, :))), '--', 'LineWidth', 2.0, ...
    'DisplayName', '\Im\{E_{z}\}');
hold on;
xline(dielectric.h * 1e3, '--', 'Color', [0.4940 0.1840 0.5560], ...
    'LineWidth', 2.0, 'DisplayName', 'interface');
grid on;
legend show;
legend('location', 'bestoutside');
xlabel('z / mm');
ylabel('E_{z}^{TM} / V/m');
sgtitle(['TM Real & Imaginary Parts @ \phi = ' ...
    num2str(phi(1) * 180 / pi) ' deg, and \rho = ' ...
    num2str(rho(end) * sqrt(dielectric.er) / wave.wavelength) ...
    ' \lambda_{d}'], 'FontWeight', 'bold', 'FontSize', 11);
saveas(gcf, 'figures\E_z_variation.fig');

% Total E in phi
figure('Position', [250 250 750 400]);
plot(phi' * 180 / pi, squeeze(Etotal(:, end, 41)), ...
    'LineWidth', 2.0, 'DisplayName', 'E_{total}');
grid on;
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\phi / deg');
ylabel('|E_{total}^{TM}| / V/m');
title(['TM |E_{total}| Amplitude @ \rho = 15\lambda_{d}, and h = ' ...
    num2str(round(z(41) * 1e3, 2)) ' mm']);
saveas(gcf, 'figures\Et_phi_variation.fig');

% Total in phi
figure('Position', [250 250 750 500]);
subplot(2, 1, 1);
plot(phi' * 180 / pi, squeeze(abs(E(:, end, 1, 41))), ...
    'LineWidth', 2.0, 'DisplayName', '|E_{\rho}|');
grid on;
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
ylabel('|E_{\rho}^{TM}| / V/m');
subplot(2, 1, 2);
plot(phi' * 180 / pi, squeeze(abs(E(:, end, 2, 41))), ...
    'LineWidth', 2.0, 'DisplayName', '|E_{z}|');
grid on;
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('\phi / deg');
ylabel('|E_{z}^{TM}| / V/m');
sgtitle(['TM Real & Imaginary Parts @ \rho = ' num2str(rho(end) * ...
    sqrt(dielectric.er) / wave.wavelength) ' \lambda_{d}, and h = ' ...
    num2str(round(z(41) * 1e3, 2)) ' mm'], ...
    'FontWeight', 'bold', 'FontSize', 11);
saveas(gcf, 'figures\Eabs_phi_variation.fig');

% Real and imaginary in phi
figure('Position', [250 250 750 500]);
subplot(2, 1, 1);
plot(phi' * 180 / pi, squeeze(real(E(:, end, 1, 41))), ...
    'LineWidth', 2.0, 'DisplayName', '\Re\{E_{z}\}');
hold on;
plot(phi' * 180 / pi, squeeze(imag(E(:, end, 1, 41))), '--', ...
    'LineWidth', 2.0, 'DisplayName', '\Im\{E_{z}\}');
grid on;
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
ylabel('E_{\rho}^{TM} / V/m');
subplot(2, 1, 2);
plot(phi' * 180 / pi, squeeze(real(E(:, end, 2, 41))), ...
    'LineWidth', 2.0, 'DisplayName', '\Re\{E_{z}\}');
hold on;
plot(phi' * 180 / pi, squeeze(imag(E(:, end, 2, 41))), '--', ...
    'LineWidth', 2.0, 'DisplayName', '\Im\{E_{z}\}');
grid on;
xlim([0 360]);
legend show;
legend('location', 'bestoutside');
xlabel('z / mm');
ylabel('E_{z}^{TM} / V/m');
sgtitle(['TM Real & Imaginary Parts @ \rho = ' num2str(rho(end) * ...
    sqrt(dielectric.er) / wave.wavelength) ' \lambda_{d}, and h = ' ...
    num2str(round(z(41) * 1e3, 2)) ' mm'], ...
    'FontWeight', 'bold', 'FontSize', 11);
saveas(gcf, 'figures\E_phi_variation.fig');
