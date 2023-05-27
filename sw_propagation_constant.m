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
dielectric.er = [10 5];
N = 1001;

%% DEPENDENT PARAMETERS
wave.f  = linspace(1, 20, N) * 1e9;
wave.wavelength = c ./ wave.f;
wave.k0 = 2 * pi ./ wave.wavelength;

for er_idx = 1 : 1 : length(dielectric.er)
    krho_norm = linspace(1, sqrt(dielectric.er(er_idx)), N);
    
    krho_te = NaN(1, length(wave.f));
    krho_tm = NaN(1, length(wave.f));
    for idx = length(wave.f) : -1 : 1
        krho = krho_norm * wave.k0(idx);
    
        [krho_te(idx), krho_tm(idx)] = find_krho(wave.k0(idx), krho, ...
            'GroundSlab', dielectric.h, dielectric.er(er_idx));
    end
    
    %% PLOT PROPAGATION CONSTANT
    wavelength_d = wave.wavelength / sqrt(dielectric.er(er_idx));
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
        ' mm, and \epsilon_{r} = ' num2str(dielectric.er(er_idx))]);
    saveas(gcf, ['figures\propagation_constant_er_' ...
        num2str(dielectric.er(er_idx)) '.fig']);
end
