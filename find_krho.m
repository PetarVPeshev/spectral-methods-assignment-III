function [krho_te, krho_tm] = find_krho(k, krho, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(varargin{1}, 'GroundSlab')
        slab_length = varargin{2};
        dielectric_er = varargin{3};

        [D_te, D_tm] = dispersion_eqn(k, krho, ...
            'GroundSlab', slab_length, dielectric_er);
    
        % Peaks
        [~, peak_te] = findpeaks( abs(1 ./ D_te) );
        [~, peak_tm] = findpeaks( abs(1 ./ D_tm) );
        
        % Guess krho
        krho_te = krho(peak_te);
        krho_tm = krho(peak_tm);
        
        % Step k
        delta_k = k / 500;
        
        % Newton's method for TE krho
        krho_te_prev = 0;
        if isempty(krho_te)
            krho_te = k;
            krho_te_prev = k;
        end

        while abs(krho_te - krho_te_prev) > 0.00001
            [D_te, ~] = dispersion_eqn(k, krho_te, ...
                'GroundSlab', slab_length, dielectric_er);
        
            [D_te_1, ~] = dispersion_eqn(k, krho_te + delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
            [D_te_2, ~] = dispersion_eqn(k, krho_te - delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
        
            D_te_prime = (D_te_1 - D_te_2) / delta_k;
            krho_te_prev = krho_te;
            krho_te = krho_te - D_te / D_te_prime;
            
            if imag(krho_te) ~= 0
                krho_te = k;
                break;
            end
        end
        
        % Newton's method for TM krho
        krho_tm_prev = 0;
        if isempty(krho_te)
            krho_tm = k;
            krho_tm_prev = k;
        end
        
        while abs(krho_tm - krho_tm_prev) > 0.00001 && ~isempty(krho_tm)
            [~, D_tm] = dispersion_eqn(k, krho_tm, ...
                'GroundSlab', slab_length, dielectric_er);
        
            [~, D_tm_1] = dispersion_eqn(k, krho_tm + delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
            [~, D_tm_2] = dispersion_eqn(k, krho_tm - delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
        
            D_tm_prime = (D_tm_1 - D_tm_2) / delta_k;
            krho_tm_prev = krho_tm;
            krho_tm = krho_tm - D_tm / D_tm_prime;
            
            if imag(krho_tm) ~= 0
                krho_tm = k;
                break;
            end
        end
    end
end