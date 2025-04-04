function [r_ijk, v_ijk] = perifocal2Equatorial(r_pkw,v_pkw,i,Omega,omega)
%perifocal2Equatorial obtains r and v in the earth centric equatorial plane
%   from r and v in the satelite or perifocal plane (2d)
%
%   Detailed explanation goes here

R = [cosd(Omega)*cosd(omega)-sind(Omega)*sind(omega)*cosd(i) -cosd(Omega)*sind(omega)-sind(Omega)*cosd(omega)*cosd(i) sind(Omega)*sind(i);
     sind(Omega)*cosd(omega)+cosd(Omega)*sind(omega)*cosd(i) -sind(Omega)*sind(omega)+cosd(Omega)*cosd(omega)*cosd(i) -cosd(Omega)*sind(i);
     sind(omega)*sind(i) cosd(omega)*sind(i) cosd(i)];

r_ijk = R*r_pkw;
v_ijk = R*v_pkw;
end