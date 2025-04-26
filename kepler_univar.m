function [r1, v1] = kepler_univar(r0, v0, deltat, mu, verbose, soft)
    if nargin < 5
        verbose = false;
    end
    if nargin < 6
        soft = 1;
    end
    h0 = cross(r0, v0);
    ecc = 1./mu .* ((mag(v0).^2-mu./mag(r0)).*r0-sum(r0.*v0).*v0);
    semimajor = mag(h0).^2./mu./(1-mag(ecc).^2);
    r0dotv0 = sum(r0.*v0);
    magr0 = mag(r0);
    tol = 1e-8;
    ncnt = 0;
    maxn = 20;
    % x_n = 0.*sqrt(mu).*deltat./semimajor.*ones(1, len); % guess idk
    x_n = sqrt(mu).*deltat./semimajor;
    t_n = 0;
    
    % define S and C!!
    funcz = @(x, a) x.^2./a;
    z0 = funcz(x_n, semimajor);
    if verbose
        fprintf("\nKepler's universal variable method\nOrbit is: ")
    end
    outstr = "";
    if z0 > tol
        S = @(z) (sqrt(z) - sin(sqrt(z)))./sqrt(z.^3);
        C = @(z) (1 - cos(sqrt(z)))/z;
        outstr = sprintf("elliptic\n");
    elseif z0 < -tol
        S = @(z) (sinh(sqrt(-z)) - sqrt(-z))./sqrt(-z.^3);
        C = @(z) (1 - cosh(sqrt(-z)))/z;
        outstr = sprintf("hyperbolic\n");
        x_n = 0;
    else
        C = @(z) (1/2 - z/factorial(4) + z.^2/factorial(6));
        S = @(z) (1/6 - z/factorial(5) + z.^2/factorial(7));
        outstr = sprintf("parabolic\n");
        x_n = sqrt(deltat);
    end
    if verbose
        fprintf(outstr)
        fprintf("n\tx_n\t  delta t  \t(dt/dx)_n\n")
    end
    z_n = z0;
    r_n = 0;
    while ncnt < maxn && abs(deltat-t_n) > tol
        ncnt = ncnt + 1;
        z_n = funcz(x_n, semimajor);
        t_n = (x_n.^3.*S(z_n) + r0dotv0./sqrt(mu).*x_n.^2.*C(z_n) + magr0.*x_n.*(1-z_n.*S(z_n)))./sqrt(mu);
        r_n = x_n.^2.*C(z_n) + r0dotv0./sqrt(mu).*x_n.*(1-z_n.*S(z_n)) + magr0.*(1-z_n.*C(z_n));
        dtdx_n = r_n./sqrt(mu);
        if verbose
            fprintf("%d\t%4.4f\t  %4.4f  \t%4.4f\n", ncnt, x_n, deltat-t_n, dtdx_n);
        end
        x_n = x_n + soft*(deltat - t_n)./dtdx_n;
    end
    f = 1 - x_n.^2./magr0.*C(z_n);
    g = t_n - x_n.^3./sqrt(mu).*S(z_n);
    fdot = x_n.*sqrt(mu)./magr0./mag(r_n).*(z_n.*S(z_n)-1);
    gdot = 1 - x_n.^2./r_n.*C(z_n);
    r1 = f.*r0 + g.*v0;
    v1 = fdot.*r0 + gdot.*v0;
    if abs(f*gdot - g*fdot - 1) > 0.001 && verbose
        disp("UniVar: f-g-cross failed")
    end
end