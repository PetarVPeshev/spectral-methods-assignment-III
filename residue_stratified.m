function varargout = residue_stratified(k0, krho_te, krho_tm, z, varargin)
%RESIDUES_STRATIFICATION Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;
    if strcmp(varargin{1}, 'GroundSlab')
        slab_length = varargin{2};
        dielectric_er = varargin{3};

        % Air propagation vector
        air_kz_te = - 1j * sqrt( - k0 ^ 2 + krho_te .^ 2 );
        air_kz_tm = - 1j * sqrt( - k0 ^ 2 + krho_tm .^ 2 );
        
        % Dielectric propagation vector
        dielectric_k = k0 * sqrt(dielectric_er);
        dielectric_kz_te = - 1j * sqrt( - dielectric_k ^ 2 + krho_te .^ 2 );
        dielectric_kz_tm = - 1j * sqrt( - dielectric_k ^ 2 + krho_tm .^ 2 );
                
        % Air impedance
        Zair_te = wave_impedance * k0 ./ air_kz_te;
        Zair_tm = wave_impedance * air_kz_tm / k0;
        
        % Substrate impedance
        Zs = wave_impedance / sqrt(dielectric_er);
        Zs_te = Zs * dielectric_k ./ dielectric_kz_te;
        Zs_tm = Zs * dielectric_kz_tm / dielectric_k;
        
        % Substrate input impedance
        Zin_te = 1j * Zs_te .* tan(dielectric_kz_te * slab_length);
        Zin_tm = 1j * Zs_tm .* tan(dielectric_kz_tm * slab_length);
        
        % Derivative of dispersion function
        delta_k = k0 / 500;
        
        [D_te_1, ~] = dispersion_eqn(k0, krho_te + delta_k / 2, ...
            'GroundSlab', slab_length, dielectric_er);
        [D_te_2, ~] = dispersion_eqn(k0, krho_te - delta_k / 2, ...
            'GroundSlab', slab_length, dielectric_er);
        D_te_prime = (D_te_1 - D_te_2) / delta_k;
                
        [~, D_tm_1] = dispersion_eqn(k0, krho_tm + delta_k / 2, ...
            'GroundSlab', slab_length, dielectric_er);
        [~, D_tm_2] = dispersion_eqn(k0, krho_tm - delta_k / 2, ...
            'GroundSlab', slab_length, dielectric_er);
        D_tm_prime = (D_tm_1 - D_tm_2) / delta_k;
        
        % Parallel impedance
        Z_te = Zair_te .* Zin_te ./ D_te_prime;
        Z_tm = Zair_tm .* Zin_tm ./ D_tm_prime;
        
        % Dielectric medium voltage and current
        v_te = Z_te .* sin(dielectric_kz_te .* z) / sin(dielectric_kz_te * slab_length);
        i_te = 1j * (Z_te / Zs) * cos(dielectric_kz_te .* z) / sin(dielectric_kz_te * slab_length);
        v_tm = Z_tm .* sin(dielectric_kz_tm .* z) / sin(dielectric_kz_tm * slab_length);
        i_tm = 1j * (Z_tm / Zs) * cos(dielectric_kz_tm .* z) / sin(dielectric_kz_tm * slab_length);

        varargout{1} = v_te;
        varargout{2} = i_te;
        varargout{3} = v_tm;
        varargout{4} = i_tm;
    else
        error('Error. Invalid argument.');
    end
end

