data = readtable('planets.csv');
data.i = data.i*2*pi/180;
data.Omega = data.Omega*2*pi/180;
data.omega = data.omega*2*pi/180;
data.theta = data.theta*2*pi/180;

mu = 1; % i'm guessing - Noa

for i = 1:height(data)
    p = data.a(i)*(1-data.eccentricity(i)^2);
    rmag = p/(1+data.eccentricity(i)^2*cos(data.theta(i)));
    rpqw = [rmag*cos(data.theta(i)), rmag*sin(data.theta(i)), 0];
    vpqw = [sqrt(mu/p)*sin(-data.theta(i)), sqrt(mu/p)*(data.eccentricity(i)+cos(data.theta(i))), 0];
    Omega = data.Omega(i);
    omega = data.omega(i);
    inc = data.i(i);
    Rtilde = peri2geo(Omega, omega, inc);
    rijk = (Rtilde*rpqw')';
    vijk = (Rtilde*vpqw')';
    data.ri0(i) = rijk(1);
    data.rj0(i) = rijk(2);
    data.rk0(i) = rijk(3);
    data.vi0(i) = vijk(1);
    data.vj0(i) = vijk(2);
    data.vk0(i) = vijk(3);
end

% random testing
r0 = [0 1 0];
v0 = [0 0 1];
h0 = cross(r0, v0);
ecc = 1./mu .* ((sqrt(sum(v0.^2)).^2-mu./sqrt(sum(r0.^2))).*r0-sum(r0.*v0).*v0);
semimajor = sqrt(sum(h0.^2)).^2./mu./(1-sqrt(sum(ecc.^2)).^2);

universal_ToF_x(r0, v0, semimajor, 0.5*3.14159265, 1)