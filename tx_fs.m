function varargout = tx_fs(k0, krho, z, slab_length)
%STRATIFIED_MEDIA Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;

    % Air propagation vector
    air_kz = - 1j * sqrt( - k0 ^ 2 + krho .^ 2 );
    
    % Air impedance
    Zair_te = wave_impedance * k0 ./ air_kz;
    Zair_tm = wave_impedance * air_kz / k0;

    % Substrate input impedance
    Zin_te = Zair_te;
    Zin_tm = Zair_tm;

    % Parallel impedance
    Z_te = Zair_te .* Zin_te ./ (Zair_te + Zin_te);
    Z_tm = Zair_tm .* Zin_tm ./ (Zair_tm + Zin_tm);

    % Air medium voltage and current
    v_te = Z_te .* exp(1j * air_kz .* (slab_length - z));
    i_te = v_te ./ Zair_te;
    v_tm = Z_tm .* exp(1j * air_kz .* (slab_length - z));
    i_tm = v_tm ./ Zair_tm;

    varargout{1} = v_te;
    varargout{2} = i_te;
    varargout{3} = v_tm;
    varargout{4} = i_tm;
end

