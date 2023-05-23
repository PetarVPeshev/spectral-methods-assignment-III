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
top_medium.er = 1;
wave.f  = linspace(1, 20, 1001) * 1e9;
N = 1001;

%% DEPENDENT PARAMETERS
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;
krho_norm = linspace(1, sqrt(dielectric.er), N);

krho_te = NaN(1, length(wave.f));
krho_tm = NaN(1, length(wave.f));
for idx = length(wave.f) : -1 : 1
    krho = krho_norm * wave.k0(idx);

    [krho_te(idx), krho_tm(idx)] = find_krho(wave.k0(idx), krho, ...
        'GroundSlab', dielectric.h, dielectric.er);
end

%% PLOT PROPAGATION CONSTANT
wavelength_d = wave.wavelength / sqrt(dielectric.er);
figure('Position', [250 250 750 400]);
plot(dielectric.h ./ wavelength_d, krho_te ./ wave.k0, ...
    'LineWidth', 2.0, 'DisplayName', 'TE1');
hold on;
plot(dielectric.h ./ wavelength_d, krho_tm ./ wave.k0, ...
    'LineWidth', 2.0, 'DisplayName', 'TM0');
grid on;
xlim([min(dielectric.h ./ wavelength_d) max(dielectric.h ./ wavelength_d)]);
legend show;
legend('location', 'bestoutside');
xlabel('h / \lambda_{d}');
ylabel('k_{\rho}/k_{0}');
title(['Normalized k_{\rho} @ h = ' num2str(dielectric.h * 1e3) ...
    ' mm, and \epsilon_{r} = ' num2str(dielectric.er)]);

%% PLOT DENOMINATOR ZEROS
[D_te, D_tm] = dispersion_eqn(wave.k0(end), krho_norm * wave.k0(end), ...
    'GroundSlab', dielectric.h, dielectric.er);

figure('Position', [250 250 750 400]);
yyaxis left;
plot(krho / wave.k0(end), abs(1 ./ D_te), 'LineWidth', 2.0, ...
    'DisplayName', 'TE');
hold on;
ylim([0 3]);
ylabel('1 / D^{TE}');
yyaxis right;
plot(krho / wave.k0(end), abs(1 ./ D_tm), 'LineWidth', 2.0, ...
    'DisplayName', 'TM');
grid on;
ylabel('1 / D^{TM}');
ylim([0 0.1]);
xlim([min(krho / wave.k0(end)) max(krho / wave.k0(end))]);
legend show;
legend('location', 'bestoutside');
xlabel('k_{\rho} / k_{0}');
title(['1 / D = 1 / (Z_{0}+Z_{sd}) @ f = ' num2str(wave.f(end) * 1e-9) ...
    ' GHz, h = ' num2str(dielectric.h * 1e3) ' mm, and ' ...
    '\epsilon_{r} = ' num2str(dielectric.er)]);
