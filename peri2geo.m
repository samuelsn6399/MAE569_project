function R = peri2geo(Omega,omega,inc)
%Peri2geo gives the R matrix needed to convert a vector from perifocal
%coordinates to geocentric-equatorial coordinates given the three angles.
R = [cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(inc),-cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(inc),sin(Omega)*sin(inc);
    sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(inc),-sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(inc),-cos(Omega)*sin(inc);
    sin(omega)*sin(inc),cos(omega)*sin(inc),cos(inc)];
end

