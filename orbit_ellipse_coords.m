function transform_ellipse = orbit_ellipse_coords(a,e,Omega,omega,inc,isdeg)
%orbit_ellipse_coords plots an ellipse in perifocal coordinates, then
%transforms it to GeoCentric Equatorial coordinates. 
if isdeg; Omega = deg2rad(Omega); omega = deg2rad(omega); 
    inc = deg2rad(inc); end

R = peri2geo(Omega,omega,inc);
c = a*e;
b = sqrt(a^2-c^2);

ellipse_x = a*cos(linspace(0,2*pi)) - c;
ellipse_y = b*sin(linspace(0,2*pi));
ellipse_z = zeros(size(ellipse_x));

ellipse_coords = [ellipse_x; ellipse_y; ellipse_z];
transform_ellipse = R * ellipse_coords;
end

