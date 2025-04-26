%% Optimize Jupiter 2 Earth Arrival + Depart Dates
clear;clc;
%% Constants and Conversions
AUTU2kms = 29.784852; % [(km/s)/(AU/TU)]
TU2day = 58.132821; % [day/TU]
TU2s = 5.023e6; % [s/TU]
s2TU = 1/TU2s; % [TU/s]
min2TU = s2TU*60; % [TU/min]
hr2TU = min2TU*60; % [TU/hr]
day2TU = hr2TU*24; % [TU/day]


%% Initialize Iterative Variables
% Iterative variables must be same size and whole days
jupiter_dDate1 = datetime(2050,01,01,11,58,00);
dDate1 = 0:10:1000;
jupiter_dDate1_idx = jupiter_dDate1+dDate1;
ToF12_dayidx = 500:10:1500; % [day]
deltaV = zeros(length(jupiter_dDate1_idx),length(ToF12_dayidx));
deltaV1 = zeros(length(jupiter_dDate1_idx),length(ToF12_dayidx));
deltaV2 = zeros(length(jupiter_dDate1_idx),length(ToF12_dayidx));
deltaV3 = zeros(length(jupiter_dDate1_idx),length(ToF12_dayidx));

%  Iterate on ToF To Flyby W/ Hohmman Transfer post fly-by
for j=1:length(jupiter_dDate1_idx)
    for i=1:length(ToF12_dayidx)
        [deltaV(j,i),deltaV1(j,i),deltaV2(j,i),deltaV3(j,i)] = Jupiter2EarthV3_func(jupiter_dDate1_idx(j), ToF12_dayidx(i));
        fprintf('Depart: %3d/%3d | ToF: %3d/%3d |\n',j,length(jupiter_dDate1_idx),i,length(ToF12_dayidx))
    end
end

figure(4);clf(4);
X = meshgrid(jupiter_dDate1_idx)';
Y = meshgrid(ToF12_dayidx);
Z = deltaV;
surf(X,Y,Z)
xlabel('Date of Departure [month-year]','Interpreter','latex')
ylabel('ToF [days]','Interpreter','latex')
zlabel('$\delta$ V [km/s]','Interpreter','latex')
title('Optimize: Time of Flight and Departture Date','Interpreter','latex')
grid on