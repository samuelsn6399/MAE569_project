function [r,v,x,t,z,S,C,dtdx,t_t] = universal_ToF_z(r0,v0,delta_t,mu,x_guess,errtol,maxiter,new_plot,show_iter)
%universal_ToF given a radius, velocity, and time, find the value of x and
%r
if nargin < 6; errtol = 0.0001; end
if nargin < 7; maxiter = 100; end
e = norm(1/mu*((norm(v0)^2-mu/norm(r0))*r0 - dot(r0,v0)*v0));
h = norm(cross(r0,v0));
a = h^2/(mu*(1-e^2));
if nargin < 5 || ~x_guess; x_guess = sqrt(mu)*delta_t/a; end
x = x_guess;
z = x^2/a;
S = S_ToF(z);
C = C_ToF(z);
t = (x^3*S + dot(r0,v0)/sqrt(mu)*x^2*C + norm(r0)*x*(1 - z*S))/sqrt(mu);
iter = 1;
while iter <= maxiter
    dtdx(iter) = 1/sqrt(mu)*(x(iter)^2*C(iter) + dot(r0,v0)/sqrt(mu)*x(iter)*(1 - z(iter)*S(iter)) + norm(r0)*(1-z(iter)*C(iter)));
    t_t(iter) = delta_t - t(iter);
    x(iter + 1) = x(iter) + t_t(iter)/dtdx(iter);
    z(iter + 1) = x(iter + 1)^2/a;
    S(iter + 1) = S_ToF(z(iter + 1));
    C(iter + 1) = C_ToF(z(iter + 1));
    t(iter + 1) = (x(iter + 1)^3*S(iter + 1) + dot(r0,v0)/sqrt(mu)*x(iter + 1)^2*C(iter + 1) + norm(r0)*x(iter + 1)*(1 - z(iter + 1)*S(iter + 1)))/sqrt(mu);
    if abs(t_t(iter)) < errtol
        break
    end
    iter = iter + 1;
end
mag_r = x(end)^2*C(end) + dot(r0,v0)/sqrt(mu)*x(end)*(1 - z(end)*S(end)) + norm(r0)*(1 - z(end)*S(end));
f = 1-x(end)^2/norm(r0)*C(end);
g = t(end) - x(end)^3/sqrt(mu)*S(end);
f_dot = sqrt(mu)/(norm(r0)*mag_r)*x(end)*(z(end)*S(end) - 1);
g_dot = 1 - x(end)^2/mag_r*C(end);
r = f*r0 + g*v0;
v = f_dot*r0 + g_dot*v0;

if new_plot
    figure
    hold on
    minx = floor(min(x))-2;
    maxx = ceil(max(x))+2;
    x_plot = linspace(minx,maxx);
    z_plot = x_plot.^2./a;
    t_plot = (x_plot.^3.*S_ToF(z_plot) + dot(r0,v0)/sqrt(mu)*x_plot.^2.*C_ToF(z_plot) + norm(r0)*x_plot.*(1 - z_plot.*S_ToF(z_plot)))/sqrt(mu);

    plot(x_plot,t_plot,[minx,maxx],[delta_t,delta_t],'k--',LineWidth=2)
    xlabel x
    ylabel t(x)


    plot(x,t,'r--*',LineWidth=2,MarkerSize=10);
    title("Newton's Method Graph")
    legend('t(x)','x','Newton Iterations',Location='best')
    xlim 'tight'
    hold off
end
if show_iter
    for iter = 1:size(x,2)-1
        if iter == 1
            fprintf(' i          x_n        Δt-tn        dt/dx        x_n+1\n')
            fprintf('__ ____________ ____________ ____________ ____________\n')
        end
        fprintf('%2d %12.5g %12.5g %12.5g %12.5g\n',iter,x(iter),t_t(iter),dtdx(iter),x(iter+1))
    end
end
end


