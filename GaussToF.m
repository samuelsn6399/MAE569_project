function [v1,v2] = GaussToF(r1,r2,ToF, short)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    short = true;
end
mu = 1;
iter_max = 100;
tol = 1e-6; % MPRAE tolerance
r1_mag = sqrt(dot(r1,r1));
r2_mag = sqrt(dot(r2,r2));

%% define gauss universal variable functions
S1_func = @(z) (sqrt(z)-sin(sqrt(z)))./sqrt(z.^3);
S2_func = @(z) (sinh(sqrt(-z))-sqrt(-z))./sqrt((-z).^3);
S3_func = @(z) 1./factorial(3)-z./factorial(5)+z.^2./factorial(7);
C1_func = @(z) (1-cos(sqrt(z)))./z;
C2_func = @(z) (1-cosh(sqrt(-z)))./z;
C3_func = @(z) 1./factorial(2)-z./factorial(4)+z.^2./factorial(6);
y_func = @(z,A,r1,r2,S,C) r1+r2-A*(1-z*S)/sqrt(C);
x_func = @(y,C) sqrt(y/C);
t_func = @(x,y,A,S) (x.^3.*S ...
                            + A*sqrt(y))/sqrt(mu); 
dSdz_func = @(z,S,C) (C-3*S)/(2*z);
dCdz_func = @(z,S,C) (1-z*S-2*C)/(2*z);
dtdz_func = @(x,y,A,S,C,dS,dC) (x^3*(dS-((3*S*dC)/(2*C))) ...
                            + A/8*(3*S*sqrt(y)/C+A/x))/sqrt(mu); 
MPRAE = @(curr, prev) 100.*abs(curr-prev)./curr;
f = @(y,r1) 1-y/r1;
f_dot = @(x,z,S,r1,r2) -sqrt(mu)*x/(r1*r2)*(1-z*S);
g = @(y,A) A*sqrt(y/mu);
g_dot = @(y,r2) 1-y/r2;

%% Shortway
dtheta_s = acos(dot(r1,r2)/(r1_mag*r2_mag));
if ~short
    dtheta_s = 2*pi - dtheta_s;
end
A_s = (sqrt(r1_mag*r2_mag)*sin(dtheta_s))/(sqrt(1-cos(dtheta_s)));

%% preallocate variables
z_n = zeros(iter_max,1);
S_n = zeros(iter_max,1);
C_n = zeros(iter_max,1);
dSdz_n = zeros(iter_max,1);
dCdz_n = zeros(iter_max,1);
dtdz_n = zeros(iter_max,1);
x_n = zeros(iter_max,1);
y_n = zeros(iter_max,1);
t_n = zeros(iter_max,1);
ea = zeros(iter_max,1);

% initialize with initial guess
z_n(1) = 1; % INITIAL GUESS
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
y_n(1) = y_func(z_n(1),A_s,r1_mag,r2_mag,S_n(1),C_n(1));
x_n(1) = x_func(y_n(1),C_n(1));
t_n(1) = t_func(x_n(1),y_n(1),A_s,S_n(1));
dSdz_n(1) = dSdz_func(z_n(1),S_n(1),C_n(1));
dCdz_n(1) = dCdz_func(z_n(1),S_n(1),C_n(1));
dtdz_n(1) = dtdz_func(x_n(1),y_n(1),A_s,S_n(1),C_n(1),dSdz_n(1),dCdz_n(1));
ea(1) = MPRAE(ToF,t_n(1));

% iterative approach: solve for x
for i = 1:iter_max
    z_n(i+1) = z_n(i) + (ToF-t_n(i))/dtdz_n(i);
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
    y_n(i+1) = y_func(z_n(i+1),A_s,r1_mag,r2_mag,S_n(i+1),C_n(i+1));
    x_n(i+1) = x_func(y_n(i+1),C_n(i+1));
    t_n(i+1) = t_func(x_n(i+1),y_n(i+1),A_s,S_n(i+1));
    dSdz_n(i+1) = dSdz_func(z_n(i+1),S_n(i+1),C_n(i+1));
    dCdz_n(i+1) = dCdz_func(z_n(i+1),S_n(i+1),C_n(i+1));
    dtdz_n(i+1) = dtdz_func(x_n(i+1),y_n(i+1),A_s,S_n(i+1),C_n(i+1),dSdz_n(i+1),dCdz_n(i+1));
    ea(i+1) = MPRAE(ToF,t_n(i+1));
    if ea(i+1)<tol
        break
    end
end

% index through each r and v vector entry (number of orbits given)

z_1 = z_n(1:find(x_n(:,1),1,'last'));
x_1 = x_n(1:length(z_1));
y_1 = y_n(1:length(z_1));
S_1 = S_n(1:length(z_1));

% solve f,g
f_1 = f(y_1(end),r1_mag);
g_1 = g(y_1(end),A_s);
% solve f_dot,g_dot
f_dot1 = f_dot(x_1(end),z_1(end),S_1(end),r1_mag,r2_mag);
g_dot1 = g_dot(y_1(end),r2_mag);


% check 1=(f)(g_dot)-(f_dot)(g)
check = f_1*g_dot1-f_dot1*g_1; % expect one
% solve v1,v2
v_1vec = (r2-f_1*r1)/g_1;
v_2vec = (g_dot1*r2-r1)/g_1;

% output
v1 = v_1vec;
v2 = v_2vec;
end