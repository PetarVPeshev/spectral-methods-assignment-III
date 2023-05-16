function [D_te, D_tm] = dispersion_eqn(k0, krho, varargin)
%DISPERSION_EQN Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;
    if strcmp(varargin{1}, 'GroundSlab')
        slab_length = varargin{2};
        dielectric_er = varargin{3};

        dielectric_k = k0 * sqrt(dielectric_er);
        dielectric_impedance = wave_impedance / sqrt(dielectric_er);
    
        air_kz = -1j * sqrt( - k0 ^ 2 + krho .^ 2 );
        dielectric_kz = - 1j * sqrt( - dielectric_k ^ 2 + krho .^ 2 );
        
        % Air impedance
        Zair_te = wave_impedance * k0 ./ air_kz;
        Zair_tm = wave_impedance * air_kz / k0;
    
        % Substrate impedance
        Zs_te = dielectric_impedance * dielectric_k ./ dielectric_kz;
        Zs_tm = dielectric_impedance * dielectric_kz / dielectric_k;
    
        % Substrate input impedance
        Zin_te = 1j * Zs_te .* tan(dielectric_kz * slab_length);
        Zin_tm = 1j * Zs_tm .* tan(dielectric_kz * slab_length);
    
        % Denominator
        D_te = Zair_te + Zin_te;
        D_tm = Zair_tm + Zin_tm;
    end
end

