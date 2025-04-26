function [r, v] = rvFromOrbitalEle(a, e, i, Omega, omega, theta)

mu = 1; % helio

% calc perifocal r vector
p = a*(1-e^2);
r_mag = p/(1+e*cosd(theta));
r_peri = [r_mag*cosd(theta);r_mag*sind(theta);0];

% calc perifocal v vector
A = sqrt(mu/p);
v_peri = A.*[-sind(theta);(e+cosd(theta));0];

% convert perifocal to heliocentric
% R_peri2helio = peri2geo(Omega, omega, i);
% r_helio = R_peri2helio*r_peri;
% v_helio = R_peri2helio*v_peri;
[r_helio,v_helio] = perifocal2Equatorial(r_peri,v_peri,i,Omega,omega);
% output
r = r_helio;
v = v_helio;