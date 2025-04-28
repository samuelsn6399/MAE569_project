function [r,v,x,t,dtdx,t_t] = universal_ToF_x(r0,v0,a,delta_t,mu,x_guess,errtol,maxiter,new_plot,show_iter)
%universal_ToF given a radius, velocity, and time, find the value of x and
%r
if nargin <= 5 
    x_guess = sqrt(mu)*delta_t/a; 
end

if nargin < 7 
    errtol = 0.0001; 
end
if nargin < 8 
    maxiter = 100;
    new_plot = false;
    show_iter = false;
end
x = x_guess;
t = a/sqrt(mu)*(x-sqrt(a)*sin(x/sqrt(a))) + dot(r0,v0)/mu*a*(1-cos(x/sqrt(a))) + norm(r0)*sqrt(a/mu)*sin(x/sqrt(a));
iter = 1;
while iter <= maxiter
    dtdx(iter) = 1/sqrt(mu)*(a + a*(dot(r0,v0)/sqrt(mu*a)*sin(x(iter)/sqrt(a)) + (norm(r0)/a - 1)*cos(x(iter)/sqrt(a))));
    t_t(iter) = delta_t - t(iter);
    x(iter + 1) = x(iter) + t_t(iter)/dtdx(iter);
    t(iter + 1) = a/sqrt(mu)*(x(iter + 1)-sqrt(a)*sin(x(iter + 1)/sqrt(a))) + dot(r0,v0)/mu*a*(1-cos(x(iter + 1)/sqrt(a))) + norm(r0)*sqrt(a/mu)*sin(x(iter + 1)/sqrt(a));
    if abs(t_t(iter)) < errtol
        break
    end
    iter = iter + 1;
end
mag_r = a + a*(dot(r0,v0)/sqrt(mu*a)*sin(x(end)/sqrt(a)) + (norm(r0)/a - 1)*cos(x(end)/sqrt(a)));
f = 1 - a/norm(r0)*(1-cos(x(end)/sqrt(a)));
g = a^2/sqrt(mu*a)*(dot(r0,v0)/sqrt(mu*a)*(1-cos(x(end)/sqrt(a))) + norm(r0)/a*sin(x(end)/sqrt(a)));
f_dot = -sqrt(mu*a)/(norm(r0)*mag_r)*sin(x(end)/sqrt(a));
g_dot = 1 - a/mag_r + a/mag_r*cos(x(end)/sqrt(a));
r = f*r0 + g*v0;
v = f_dot*r0 + g_dot*v0;

if new_plot
    figure
    hold on
    minx = floor(min(x))-2;
    maxx = ceil(max(x))+2;
    x_plot = linspace(minx,maxx,(maxx-minx)*15);
    mu = 1;
    t_plot = a/sqrt(mu)*(x_plot-sqrt(a)*sin(x_plot/sqrt(a))) + dot(r0,v0)/mu*a*(1-cos(x_plot/sqrt(a))) + norm(r0)*sqrt(a/mu)*sin(x_plot/sqrt(a));

    plot(x_plot,t_plot,[minx,maxx],[delta_t,delta_t],'k--',LineWidth=2)
    xlabel x
    ylabel t(x)


    plot(x,t,'r--*',LineWidth=2,MarkerSize=10);
    title("Newton's Method Graph")
    legend('t(x)','Δt','Newton Iterations',Location='best')
    xlim 'tight'
    hold off
end
if show_iter
    fprintf(' i          x_n        Δt-tn        dt/dx        x_n+1\n')
    fprintf('__ ____________ ____________ ____________ ____________\n')
    for iter = 1:size(x,2)-1
        fprintf('%2d %12.5g %12.5g %12.5g %12.5g\n',iter,x(iter),t_t(iter),dtdx(iter),x(iter+1))
    end
end
end

