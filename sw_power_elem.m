function Psw = sw_power_elem(k0, er, h, ksw, component)
%SW_POWER Summary of this function goes here
%   Detailed explanation goes here
    wave_impedance = 376.730313668;
    if strcmp(component, 'TM')
        Iphi = pi;
    
        % Numerical solution for surface wave power @ dielectric
        zd = linspace(0, h, 2001);
        dzd = zd(2) - zd(1);
        [~, ~, ~, id_tm] = residue_stratified(k0, 0, ksw, zd, ...
            'GroundSlab', h, er);
        Iz_d = sum(abs(id_tm) .^ 2) * dzd / er;
        % Numerical solution for surface wave power @ dielectric surface
        z = linspace(h + h / 1000, 3 * h, 2001);
        dz = z(2) - z(1);
        [~, ~, ~, i_tm] = residue_stratified(k0, 0, ksw, z, ...
            'GroundSlab', h, er);
        Iz_s = sum(abs(i_tm) .^ 2) * dz;
        Iz = wave_impedance * (Iz_d + Iz_s) / k0;
    
        Psw = 0.5 * (ksw ^ 2) * Iz * Iphi / (2 * pi);
    else
        erro('Error. Invalid argument.');
    end
end

