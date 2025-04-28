function [E, M, dMdE, MMn] = Newton_ToF(E_init,M_init,delta_t,n,M_0,e,errtol)
%NEWTON_TOF calculates E and M using Newton's method for an ellipse.
%   delta_t and n should have units such that n*delta_t is unitless.
if nargin < 7; errtol = 1e-4; end

iter = 1;
max_iter = 100;

E = E_init;
M = M_init;

M_1 = n*delta_t + M_0;
k = floor(M_1/(2*pi));
M_1 = M_1 - k*(2*pi);
while iter <= max_iter
    dMdE(iter) = 1 - e*cos(E(iter));
    MMn(iter) = (M_1-M(iter));
    E(iter + 1) = E(iter) + MMn(iter)/dMdE(iter);
    M(iter + 1) = (E(iter+1) - e*sin(E(iter+1)));
    if abs(MMn(iter)) < errtol
        break
    end
    iter = iter + 1;
end
end

