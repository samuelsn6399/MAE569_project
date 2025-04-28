function E = eccentric_anomoly(e,nu,isdeg)
%eccentric_anomoly Calclulates the eccentric anomoly from e and nu
% If nu is in degrees, set isdeg to true. Returns E in radians.
if nargin < 3; isdeg = false; end
if isdeg; nu = deg2rad(nu); end
E = acos((e + cos(nu))/(1 + e*cos(nu)));
if nu > pi()
    E = 2*pi()-E;
end
end

