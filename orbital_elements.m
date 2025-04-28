function [e,i,Omega,omega,nu,a,p,h_mag] = orbital_elements(r,v,mu,isdeg)
    %Calculates Orbital Elements from position and velocity vectors.
    %Returns [e,i,Omega,omega,nu,a,p]
    if nargin < 4; isdeg = false; end
    if size(r) == [1,3]; r = r'; end
    if size(v) == [1,3]; v = v'; end

    h = cross(r,v);
    h_mag = norm(h);
    n = cross([0,0,1]',h);
    n_mag = norm(n);
    r_mag = norm(r);
    v_mag = norm(v);
    e_vec = 1/mu*((v_mag^2-(mu/r_mag))*r-(dot(r,v)*v));
    e = norm(e_vec);
    i = acos(h(3)/h_mag);
    if n(2) > 0
        Omega = acos(n(1)/n_mag);
    else
        Omega = 2*pi - acos(n(1)/n_mag);
    end
    if e_vec(3) > 0
        omega = acos(dot(n,e_vec)/(n_mag*e));
    else
        omega = 2*pi - acos(dot(n,e_vec)/(n_mag*e));
    end
    if dot(r,v) > 0
        nu = acos(dot(e_vec,r)/(e*r_mag));
    else
        nu = 2*pi - acos(dot(e_vec,r)/(e*r_mag));
    end
    p = h_mag^2/mu;
    a = p/(1-e^2);

    if isdeg; i = rad2deg(i); nu = rad2deg(nu);
        Omega = rad2deg(Omega); omega = rad2deg(omega); end

end

