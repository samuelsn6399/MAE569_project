function [a_mag,e_mag,i,Omega,omega,nu] = orbitalElements(r,v,mu)
%orbitalElements obtains the orbital elements of a satelite given its
%   velocity and distance relative to the center of the earth.
%
%   Input the distance and velocity vectors to obtain the 6-orbital elemets
%   to fully define a satelite's orbital motion

% if gravitational parameter mu is not specified default to Earth's
if isempty(mu)
    mu = 3.986e5; % [km^3/s^2]
end
r_mag = sqrt(sum(r.^2)); % magnitude of r vector
v_mag = sqrt(sum(v.^2)); % magnitude of v vector

% 3-fundamental vectors
h = cross(r,v); % orbital momentum vector
h_mag = sqrt(sum(h.^2)); % magnitude of h vector
n = cross([0 0 1], h); %
n_mag = sqrt(sum(n.^2)); % magnitude of n vector
e = (1/mu)*((v_mag^2-mu/r_mag)*r - (dot(r,v)*v)); % essentricity vector
e_mag = sqrt(sum(e.^2)); % magnitude of e vector
e_unit = e./e_mag; % directional unit vector of e

% orbital elements
p_mag = h_mag^2/mu; % semi-minor axis
p = p_mag.*e_unit; % semi-minor axis vector
a_mag = p_mag/(1-e_mag^2); % semi-major axis
a = a_mag.*(-e_unit); % semi-major axis vector
i = rad2deg(acos(h(3)/h_mag)); % angle on inclination
Omega = rad2deg(acos(n(1)/n_mag)); % angle of
if n(2)<0
    Omega = 360-Omega;
end
omega = rad2deg(acos(dot(n,e)/(n_mag*e_mag))); % angle of
if e(3)<0
    omega = 360-Omega;
end
nu = rad2deg(acos(dot(e,r)/(e_mag*r_mag))); % true anomaly
if dot(r,v)<0
    nu = 360-nu;
end
end