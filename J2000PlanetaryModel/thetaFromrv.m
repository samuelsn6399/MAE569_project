function theta = thetaFromrv(r,v)

% heliocentric
mu = 1;

% magnitudes
r_mag = sqrt(dot(r,r));
v_mag = sqrt(dot(v,v));

% calc eccentricity vector
e = 1/mu*((v_mag^2-mu/r_mag)*r - dot(r,v)*v);
e_mag = sqrt(dot(e,e));
if dot(r,v)>0
    theta = rad2deg(acos((dot(e,r))/(e_mag*r_mag)));
else
    theta = 360-rad2deg(acos((dot(e,r))/(e_mag*r_mag)));
end

