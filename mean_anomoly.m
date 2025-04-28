function M = mean_anomoly(e,E,isdeg)
%mean_anomoly Calclulates the mean anomoly from E and e
% If E is in degrees, set isdeg to true. Returns M in same unit as E.
if nargin < 3; isdeg = false; end
if isdeg; E = deg2rad(E); end
M = E - e*sin(E);
if isdeg; M = rad2deg(M); end
end

