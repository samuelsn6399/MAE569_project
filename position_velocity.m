function [r,v] = position_velocity(a,e,nu,mu,Omega,omega,i,isdeg)
%position_velocity outputs position and velocity vectors for given orbital
%elements
    if isdeg; nu = deg2rad(nu); Omega = deg2rad(Omega); 
        omega = deg2rad(omega); i = deg2rad(i); end
    p = a*(1-e^2);
    r_mag = p/(1+e*cos(nu));
    r_peri = [r_mag*cos(nu), r_mag*sin(nu), 0]';
    v_peri = sqrt(mu/p)*[-sin(nu), e+cos(nu), 0]';
    R = peri2geo(Omega,omega,i);
    r = R*r_peri;
    v = R*v_peri;
end

