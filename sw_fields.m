function E = sw_fields(k0, krho, v, i, J, er, cyl_grid, component)
%SW_FIELDS Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;
    rho = cyl_grid(:, :, 1);
    phi = cyl_grid(:, :, 2);

    if strcmp(component, 'TE')
        error('Error. Not implemented.');
    elseif strcmp(component, 'TM')
        dielectric_impedance = wave_impedance ./ sqrt(er);
        k = k0 * sqrt(er);

        const = 1j * sqrt(krho / (2 * pi)) * exp(1j * pi / 4);
        const = const .* cos(phi) .* exp(- 1j * krho * rho) ./ sqrt(rho);

        Erho = v .* J(:, :, 1) .* const;
        Ez = - dielectric_impedance * (krho / k) * i .* J(:, :, 1) .* const;

        E = zeros( [size(J, 1, 2), 3] );
        E(:, :, 1) = Erho;
        E(:, :, 3) = Ez;
    else
        error('Error. Invalid argument.');
    end
end

