%% Load J2000 Planetary Data, Initialize Vectors, Setup Unit Conversions
clc;clear;
% set time of flight (time span of planetary model starting from year 2000)
ToF = 365*80+1; % [days, including the day initial day (start at day 1)]
dT = 1; % [days]
T_day = 1:dT:ToF; % Day 1: January 1, 2000 at 11:58 am UTC
date = datetime('2000-01-01 11:58:00.00')+T_day'-1;

% define unit conversions
TU2s = 5.023e6; % [s/TU]
s2TU = 1/TU2s; % [TU/s]
min2TU = s2TU*60; % [TU/min]
hr2TU = min2TU*60; % [TU/hr]
day2TU = hr2TU*24; % [TU/day]

% convert to Astronomical Units for calculations
T_TU = T_day*day2TU;

% initialize J2000 elements as a struct array
J2000 = struct('planet',{'Mercury','Venus  ','Earth  ','Mars   ', ...
                         'Jupiter','Saturn ','Uranus ','Neptune', ...
                         'Pluto  '}, ...
               'a', {0.387099,0.723332,1.000000,1.523662,5.203363, ...
                     9.537070,19.19126,30.06896,39.48169}, ...
               'e', {0.205631,0.006773,0.016710,0.093412,0.048393, ...
                     0.054151,0.047168,0.008586,0.248808}, ...
               'i', {7.004870,3.394710,0.000050,1.850610,1.305300, ...
                     2.484460,0.769860,1.769170,17.14175}, ...
               'Omega', {48.331670,76.680690,-11.26064,49.578520,100.55615, ...
                         113.71504,74.229880,131.72169,110.30347}, ...
               'omega', {29.124780,54.852290,114.20783,286.46230,-85.80230, ...
                         -21.28310,96.734360,-86.75034,110.30347}, ...
          ...  % initialize time dependent elements  
               'date',date, ...
               'T',zeros(length(T_day),1), ...
               'r',zeros(length(T_day),3), ...
               'v',zeros(length(T_day),3), ...
               'theta', {[174.7944;zeros(length(T_day)-1,1)], ...
                         [50.44675;zeros(length(T_day)-1,1)], ...
                         [-2.48284;zeros(length(T_day)-1,1)], ...
                         [19.41248;zeros(length(T_day)-1,1)], ...
                         [19.55053;zeros(length(T_day)-1,1)], ...
                         [-42.4876;zeros(length(T_day)-1,1)], ...
                         [142.2679;zeros(length(T_day)-1,1)], ...
                         [259.9087;zeros(length(T_day)-1,1)], ...
                         [14.86205;zeros(length(T_day)-1,1)]});

%% Run Time of Flight Calculation For the Specified Time Interval
% iterate for each of the 9 planets 
% iterate for each day
for p = 1:9
    for i = 1:length(T_day)
        % initial position
        if i == 1
            J2000(p).T(1) = T_day(1);
            [r_0,v_0] = rvFromOrbitalEle(J2000(p).a,J2000(p).e,J2000(p).i, ...
                                     J2000(p).Omega,J2000(p).omega, ...
                                    [J2000(p).theta(1)]);
            J2000(p).r(1,:) = r_0;
            J2000(p).v(1,:) = v_0;
        else
        J2000(p).T(i) = T_day(i);
        [r,v,theta] = keplarToF(J2000(p).a,J2000(p).e,J2000(p).i, ...
                                J2000(p).Omega,J2000(p).omega, ...
                                J2000(p).r(1,:),J2000(p).v(1,:), ...
                                J2000(p).theta(1),T_TU(i));
        J2000(p).r(i,:) = r;
        J2000(p).v(i,:) = v;
        J2000(p).theta(i) = theta;
        % update progress
        fprintf('planet: %5.1f/9 Day: %5.1f/%5.1f\n', p,i*dT+1,ToF);
        end
    end
end

%% Save To .mat File
save('J200.mat','J2000')

%% Display Solar System
clf;
figure(1)
hold on
for plot=1:9
    plot3(J2000(plot).r(:,1),J2000(plot).r(:,2),J2000(plot).r(:,3))
end
plot3(0,0,0,'.','Color','r','MarkerSize',20)
hold off
legend('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus', ...
        'Neptune','Pluto','Sun','Location','bestoutside')
title('Solar System Model')
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])
grid on
view(0,0)

%% Find Planetary Positions at December 25, 2025 at 8:37 am UTC.
DoI0 = datetime('2000-01-01 11:58:00.00');
DoI1 = datetime('2025-12-25 08:37:00.00');
DoI_diff = seconds(DoI1-DoI0);
DoI_TU = DoI_diff*s2TU;

% Initialize table data storage
roI_all = zeros(9,3);
voI_all = zeros(9,3);
thetaoI_all = zeros(9,1);

% calc for each planet
for p1 = 1:9
    [roI,voI,thetaoI] = keplarToF(J2000(p1).a,J2000(p1).e,J2000(p1).i, ...
        J2000(p1).Omega,J2000(p1).omega, ...
        J2000(p1).r(1,:),J2000(p1).v(1,:), ...
        J2000(p1).theta(1),DoI_TU);
    roI_all(p1, :) = roI;
    voI_all(p1, :) = voI;
    thetaoI_all(p1) = thetaoI;
end
planet = {J2000.planet};
% Create table
T = table(planet', roI_all(:,1), roI_all(:,2), roI_all(:,3), ...
          voI_all(:,1), voI_all(:,2), voI_all(:,3), thetaoI_all, ...
          'VariableNames', {'Planet', 'r_x', 'r_y', 'r_z', ...
                            'v_x', 'v_y', 'v_z', 'theta'});

% Display table
disp(T);
% save table
save('2025_12_25_08_37_00.mat','T')