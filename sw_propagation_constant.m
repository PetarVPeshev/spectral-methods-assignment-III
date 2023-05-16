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

for idx = length(wave.f) : -1 : length(wave.f)
    krho = krho_norm * wave.k0(idx);

    [D_te, D_tm] = dispersion_eqn(wave.k0(idx), krho, ...
        'GroundSlab', dielectric.h, dielectric.er);

    % Peaks
    [~, peak_te] = findpeaks( abs(1 ./ D_te) );
    [~, peak_tm] = findpeaks( abs(1 ./ D_tm) );

    % Guess krho
    krho_p_te = krho(peak_te);
    krho_p_tm = krho(peak_tm);

    
end

%% PLOT DENOMINATOR ZEROS
figure('Position', [250 250 750 400]);
yyaxis left;
plot(krho / wave.k0(idx), abs(1 ./ D_te), 'LineWidth', 2.0, ...
    'DisplayName', 'TE');
hold on;
ylim([0 3]);
ylabel('1 / D^{TE}');
yyaxis right;
plot(krho / wave.k0(idx), abs(1 ./ D_tm), 'LineWidth', 2.0, ...
    'DisplayName', 'TM');
grid on;
ylabel('1 / D^{TM}');
ylim([0 0.1]);
xlim([min(krho / wave.k0(idx)) max(krho / wave.k0(idx))]);
legend show;
legend('location', 'bestoutside');
xlabel('k_{\rho} / k_{0}');
title(['1 / D = 1 / (Z_{0}+Z_{sd}) @ f = ' num2str(wave.f(end) * 1e-9) ...
    ' GHz, h = ' num2str(dielectric.h * 1e3) ' mm, and ' ...
    '\epsilon_{r} = ' num2str(dielectric.er)]);
