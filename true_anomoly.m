function nu = true_anomoly(e,E,isdeg)
%true_anomoly Calclulates the true anomoly from E and e
% If E is in degrees, set isdeg to true. Returns nu in same unit as E.
if nargin < 3; isdeg = false; end
if isdeg; E = deg2rad(E); end

nu = acos((e - cos(E))/(e*cos(E) - 1));
if isdeg; nu = rad2deg(nu); end
end

