%% Index J2000.mat
% load the variable into workspace first!

% j2000.mat is a structured array; meaning it uses a combination of dot
% indexing and numeric indexing. [ex. J2000(1).planet] returns the first
% planet's name, which is 'mercury'


% index to a specific planet using J2000(i); where i is:
% 1 - Mercury
% 2 - Venus
% 3 - Earth
% 4 - Mars
% 5 - Jupiter
% 6 - Saturn
% 7 - Uranus
% 8 - Neptune
% 9 - Pluto

% the dot index names are as follows
% J2000.planet
% J2000.a
% J2000.e
% J2000.i
% J2000.Omega
% J2000.omega
% J2000.date
% J2000.T
% J2000.rvFromOrbitalElev
% J2000.theta

%% examples:
% index to earths r,v,and theta values on December 25 2025
% Note: The time 11:58 must be specified
clear;clc;
%                                        (   Y, M, D, H,MI, S)
date_index = find(J2000(3).date==datetime(2025,12,25,11,58,00))
r = J2000(3).r(date_index,:) % return the entire r vector
v = J2000(3).r(date_index,:) % return the entire v vector
theta = J2000(3).theta(date_index) % return theta
T = J2000(3).T(date_index) % return the time in days since J2000