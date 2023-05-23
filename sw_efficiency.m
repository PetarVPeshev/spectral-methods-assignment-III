close all;
clear;
clc;

if ~exist([pwd() '\figures'], 'dir')
    mkdir('figures');
end

addpath('../spectral-methods-library');
c = physconst('LightSpeed');

%% PARAMETERS
dielectric.h  = 2e-3;
dielectric.er = 10;
dipole.W = 0.5e-3;
dipole.L = 5.3e-3;
N = 101;
Ntheta = 200;
Nphi = 800;
R = 1;

%% DEPENDENT PARAMETERS
wave.f  = linspace(1, 15, N) * 1e9;
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;

%% SPHERICAL COORDINATES
theta = linspace(eps, pi / 2 - 0.1 * pi / 180, Ntheta);
phi = linspace(0, 2 * pi, Nphi);
sph_grid = meshgrid_comb(theta, phi);

dph = phi(2) - phi(1);

%% FAR-FIELD Z
z = R * cos(sph_grid(:, :, 1));

%% ELEMENTARY CURRENT
J_elem = zeros( [size(sph_grid, 1, 2), 2] );
J_elem(:, :, 1) = 1;

fs_power_elem = NaN(1, length(wave.f));
rad_power_elem = NaN(1, length(wave.f));
Psw_elem = NaN(1, length(wave.f));
for f_idx = 1 : 1 : length(wave.f)
    [k_comp, k] = wave_vector(1, wave.k0(f_idx), sph_grid);
    krho = sqrt(k_comp(:, :, 1) .^ 2 + k_comp(:, :, 2) .^ 2);
    
    %% RADIATED POWER BY ELEMENTARY DIPOLE IN FREE SPACE
    [fs_vte, fs_ite, fs_vtm, fs_itm] = tx_fs(wave.k0(f_idx), krho, z, dielectric.h);
    FS_SGF = spectral_gf(1, wave.k0(f_idx), k_comp(:, :, 1), k_comp(:, :, 2), fs_vtm, fs_vte, fs_itm, fs_ite, 'E', 'J');
    FS_E = farfield(wave.k0(f_idx), R, sph_grid, k_comp(:, :, 3), z, FS_SGF, J_elem, dielectric.h);
    [~, ~, fs_power_elem(f_idx)] = directivity(1, FS_E, sph_grid, R);
    fs_power_elem(f_idx) = 2 * fs_power_elem(f_idx);
    
    %% RADIATED POWER BY ELEMENTARY DIPOLE
    [vte, ite, vtm, itm] = stratified_media(wave.k0(f_idx), krho, z, 'GroundSlab', dielectric.h, dielectric.er);
    SGF = spectral_gf(1, wave.k0(f_idx), k_comp(:, :, 1), k_comp(:, :, 2), vtm, vte, itm, ite, 'E', 'J');
    E = farfield(wave.k0(f_idx), R, sph_grid, k_comp(:, :, 3), z, SGF, J_elem, dielectric.h);
    [~, ~, rad_power_elem(f_idx)] = directivity(1, E, sph_grid, R);

    %% SURFACE WAVE POWER FOR ELEMENTARY DIPOLE
    krho_range = wave.k0(f_idx) * linspace(1, sqrt(dielectric.er), 1001);
    [~, krho_tm] = find_krho(wave.k0(f_idx), krho_range, ...
        'GroundSlab', dielectric.h, dielectric.er);
    Psw_elem(f_idx) = sw_power_elem(wave.k0(f_idx), dielectric.er, dielectric.h, krho_tm, 'TM');
end

%% DIPOLE CURRENT
rad_power_dipole = NaN(1, length(wave.f));
Psw_dipole = NaN(1, length(wave.f));
for f_idx = 1 : 1 : length(wave.f)
    [k_comp, k] = wave_vector(1, wave.k0(f_idx), sph_grid);
    krho = sqrt(k_comp(:, :, 1) .^ 2 + k_comp(:, :, 2) .^ 2);
    
    %% RADIATED POWER BY ELEMENTARY DIPOLE
    [vte, ite, vtm, itm] = stratified_media(wave.k0(f_idx), krho, z, 'GroundSlab', dielectric.h, dielectric.er);
    SGF = spectral_gf(1, wave.k0(f_idx), k_comp(:, :, 1), k_comp(:, :, 2), vtm, vte, itm, ite, 'E', 'J');
    J = ft_current(wave.k0(f_idx), k_comp, dipole.W, dipole.L, dielectric.er, 'dipole', 'x');
    E = farfield(wave.k0(f_idx), R, sph_grid, k_comp(:, :, 3), z, SGF, J, dielectric.h);
    [~, ~, rad_power_dipole(f_idx)] = directivity(1, E, sph_grid, R);

    %% SURFACE WAVE POWER FOR ELEMENTARY DIPOLE
    krho_range = wave.k0(f_idx) * linspace(1, sqrt(dielectric.er), 1001);
    [~, krho_tm] = find_krho(wave.k0(f_idx), krho_range, ...
        'GroundSlab', dielectric.h, dielectric.er);
    Psw_dipole(f_idx) = sw_power_elem(wave.k0(f_idx), dielectric.er, dielectric.h, krho_tm, 'TM');
    
    kx_sw = krho_tm * cos(phi);
    ky_sw = krho_tm * sin(phi);
    k_comp_sw = NaN( [1, length(phi)] );
    k_comp_sw(:, :, 1) = kx_sw;
    k_comp_sw(:, :, 2) = ky_sw;
    Jsw = ft_current(wave.k0(f_idx), k_comp_sw, dipole.W, dipole.L, dielectric.er, 'dipole', 'x');
    Iphi = sum( (abs(Jsw(:, :, 1)) .^ 2) .* (cos(phi) .^ 2) ) * dph;

    Psw_dipole(f_idx) = Psw_dipole(f_idx) * Iphi / pi;
end

%% UNIFORM CURRENT
rad_power_uniform = NaN(1, length(wave.f));
Psw_uniform = NaN(1, length(wave.f));
for f_idx = 1 : 1 : length(wave.f)
    [k_comp, k] = wave_vector(1, wave.k0(f_idx), sph_grid);
    krho = sqrt(k_comp(:, :, 1) .^ 2 + k_comp(:, :, 2) .^ 2);
    
    %% RADIATED POWER BY UNIFORM CURRENT
    [vte, ite, vtm, itm] = stratified_media(wave.k0(f_idx), krho, z, 'GroundSlab', dielectric.h, dielectric.er);
    SGF = spectral_gf(1, wave.k0(f_idx), k_comp(:, :, 1), k_comp(:, :, 2), vtm, vte, itm, ite, 'E', 'J');
    J = ft_current(wave.k0(f_idx), k_comp, dipole.W, dipole.L, dielectric.er, 'uniform', 'x');
    E = farfield(wave.k0(f_idx), R, sph_grid, k_comp(:, :, 3), z, SGF, J, dielectric.h);
    [~, ~, rad_power_uniform(f_idx)] = directivity(1, E, sph_grid, R);

    %% SURFACE WAVE POWER FOR UNIFORM CURRENT
    krho_range = wave.k0(f_idx) * linspace(1, sqrt(dielectric.er), 1001);
    [~, krho_tm] = find_krho(wave.k0(f_idx), krho_range, ...
        'GroundSlab', dielectric.h, dielectric.er);
    Psw_uniform(f_idx) = sw_power_elem(wave.k0(f_idx), dielectric.er, dielectric.h, krho_tm, 'TM');
    
    kx_sw = krho_tm * cos(phi);
    ky_sw = krho_tm * sin(phi);
    k_comp_sw = NaN( [1, length(phi)] );
    k_comp_sw(:, :, 1) = kx_sw;
    k_comp_sw(:, :, 2) = ky_sw;
    Jsw = ft_current(wave.k0(f_idx), k_comp_sw, dipole.W, dipole.L, 'uniform', 'x');
    Iphi = sum( (abs(Jsw(:, :, 1)) .^ 2) .* (cos(phi) .^ 2) ) * dph;

    Psw_uniform(f_idx) = Psw_uniform(f_idx) * Iphi / pi;
end

%% EFFICIENCY
eta_elem = rad_power_elem ./ (rad_power_elem + Psw_elem);
eta_dipole = rad_power_dipole ./ (Psw_dipole + rad_power_dipole);
eta_uniform = rad_power_uniform ./ (Psw_uniform + rad_power_uniform);

%% PLOT
figure('Position', [250 250 750 400]);
plot(wave.f * 1e-9, rad_power_elem ./ fs_power_elem, 'LineWidth', 2.0, ...
    'DisplayName', 'radiated');
hold on;
plot(wave.f * 1e-9, Psw_elem ./ fs_power_elem, 'LineWidth', 2.0, ...
    'DisplayName', 'surface wave')
grid on;
xticks(min(wave.f * 1e-9) : 2 : max(wave.f * 1e-9));
xlim([min(wave.f * 1e-9) max(wave.f * 1e-9)]);
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('|P| / |P_{rad,FS}|');
title(['Normalized P_{rad} & P_{SW} @ elementary current, h = ' ...
    num2str(dielectric.h * 1e3) ' mm, and \epsilon_{r} = ' ...
    num2str(dielectric.er)]);
saveas(gcf, 'figures\rad_and_sw_elementary.fig');

figure('Position', [250 250 750 400]);
plot(wave.f * 1e-9, eta_elem * 100, 'LineWidth', 2.0, ...
    'DisplayName', '\eta');
hold on;
plot(wave.f * 1e-9, eta_dipole * 100, 'LineWidth', 2.0, ...
    'DisplayName', '\eta, dipole');
hold on;
plot(wave.f * 1e-9, eta_uniform * 100, 'LineWidth', 2.0, ...
    'DisplayName', '\eta, uniform');
grid on;
xticks(min(wave.f * 1e-9) : 2 : max(wave.f * 1e-9));
xlim([min(wave.f * 1e-9) max(wave.f * 1e-9)]);
legend show;
legend('location', 'bestoutside');
xlabel('f / GHz');
ylabel('\eta / %');
title(['Efficiency @ h = ' num2str(dielectric.h * 1e3) ' mm, ' ...
    'and \epsilon_{r} = ' num2str(dielectric.er)]);
saveas(gcf, 'figures\efficiency.fig');
