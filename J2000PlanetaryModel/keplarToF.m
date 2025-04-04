function [r, v, theta] = keplarToF(a,e,i,Omega,omega,r_0,v_0,theta_0,ToF)

r_0mag = sqrt(dot(r_0,r_0));
v_0mag = sqrt(dot(v_0,v_0));

mu = 1; % helio
n = sqrt(mu./(a.^3)); % mean motion

% define universal keplar functions
z_func = @(x,a) x.^2./a;
S1_func = @(z) (sqrt(z)-sin(sqrt(z)))./sqrt(z.^3);
S2_func = @(z) (sinh(sqrt(-z))-sqrt(-z))./sqrt((-z).^3);
S3_func = @(z) 1./factorial(3)-z./factorial(5)+z.^2./factorial(7);
C1_func = @(z) (1-cos(sqrt(z)))./z;
C2_func = @(z) (1-cosh(sqrt(-z)))./z;
C3_func = @(z) 1./factorial(2)-z./factorial(4)+z.^2./factorial(6);
r_func = @(x,z,S,C,r,r_mag,v) x.^2.*C ...
                            + dot(r,v)./sqrt(mu).*x.*(1-z.*S) ...
                            + r_mag.*(1-z.*C); % position calculation
t_func = @(x,z,S,C,r,r_mag,v) x.^3.*S ...
                            + dot(r,v)./sqrt(mu).*x.^2.*C ...
                            + r_mag.*x.*(1-z.*S); % time calculation
dtdx_func = @(r) r./sqrt(mu); % derivative of time w/ respect to x variable
MPRAE = @(curr, prev) 100.*abs(curr-prev)./curr; % magnitude percent relative absolute error
f = @(x,C,r_0mag) 1-x.^2./r_0mag.*C;
f_dot = @(x,z,S,r_0mag,r_mag) sqrt(mu).*x./(r_0mag.*r_mag).*(z.*S-1);
g = @(x,t,S) t-x.^3./sqrt(mu).*S;
g_dot = @(x,C,r_mag) 1-x.^2./r_mag.*C;

% newton method for iterative solve setup
iter_max = 100;
tol = 0.001; % MPRAE tolerance

% preallocate variables
r_n = zeros(iter_max,1);
x_n = zeros(iter_max,1);
z_n = zeros(iter_max,1);
S_n = zeros(iter_max,1);
C_n = zeros(iter_max,1);
dtdx_n = zeros(iter_max,1);
t_n = zeros(iter_max,1);
ea = zeros(iter_max,1);
% initialize with initial guess
x_n(1) = sqrt(1)*ToF/a; % INITIAL GUESS: ELLIPTICAL
% x_n(1) = 1; % INITIAL GUESS: GENERAL
z_n(1) = z_func(x_n(1),a);
if z_n(1)>0.001 % elliptical case
    S_n(1) = S1_func(z_n(1));
    C_n(1) = C1_func(z_n(1));
elseif z_n(1)<-0.001 % hyperbolic case
    S_n(1) = S2_func(z_n(1));
    C_n(1) = C2_func(z_n(1));
else % parabolic case
    S_n(1) = S3_func(z_n(1));
    C_n(1) = C3_func(z_n(1));
end

r_n(1) = r_func(x_n(1),z_n(1),S_n(1),C_n(1),r_0,r_0mag,v_0);
dtdx_n(1) = dtdx_func(r_n(1));
t_n(1) = t_func(x_n(1),z_n(1),S_n(1),C_n(1),r_0,r_0mag,v_0);
ea(1) = MPRAE(ToF,t_n(1));

% iterative approach: solve for x
for i = 1:iter_max
    x_n(i+1) = x_n(i) + (ToF-t_n(i))./dtdx_n(i);
    z_n(i+1) = z_func(x_n(i+1),a);
    if z_n(i+1)>0.001 % elliptical case
        S_n(i+1) = S1_func(z_n(i+1));
        C_n(i+1) = C1_func(z_n(i+1));
    elseif z_n(i+1)<-0.001 % hyperbolic case
        S_n(i+1) = S2_func(z_n(i+1));
        C_n(i+1) = C2_func(z_n(i+1));
    else % parabolic case
        S_n(i+1) = S3_func(z_n(i+1));
        C_n(i+1) = C3_func(z_n(i+1));
    end
    r_n(i+1) = r_func(x_n(i+1),z_n(i+1),S_n(i+1),C_n(i+1),r_0,r_0mag,v_0);
    dtdx_n(i+1) = dtdx_func(r_n(i+1));
    t_n(i+1) = t_func(x_n(i+1),z_n(i+1),S_n(i+1),C_n(i+1),r_0,r_0mag,v_0);
    ea(i+1) = MPRAE(ToF,t_n(i+1));

    if ea(i+1)<tol
        break
    end
end

% index through each r and v vector entry (number of orbits given)

x_1 = x_n(1:find(x_n(:,1),1,'last'));
r_1 = r_n(1:length(x_1));
z_1 = z_n(1:length(x_1));
S_1 = S_n(1:length(x_1));
C_1 = C_n(1:length(x_1));
dtdx_1 = dtdx_n(1:length(x_1));
t_1 = t_n(1:length(x_1));
ea_1 = ea(1:length(x_1));
% solve f,g,and r
f_1 = f(x_1(end),C_1(end),r_0mag);
g_1 = g(x_1(end),t_1(end),S_1(end));
% solve f_dot,g_dot, and v
f_dot1 = f_dot(x_1(end),z_1(end),S_1(end),r_0mag,r_1(end));
g_dot1 = g_dot(x_1(end),C_1(end),r_1(end));
% check 1=(f)(g_dot)-(f_dot)(g)
check = f_1.*g_dot1-f_dot1.*g_1; % expect ones
r_1vec = f_1.*r_0(:)+g_1.*v_0(:);
v_1vec = f_dot1.*r_0(:)+g_dot1.*v_0(:);

% calc theta
theta_1 = thetaFromrv(r_1vec,v_1vec);

% output
r = r_1vec;
v = v_1vec;
theta = theta_1;