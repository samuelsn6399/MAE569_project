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

%% Iterate On Depart Date and Flyby Date
jupiter_dDate1 = datetime(2050,04,12,11,58,00);
% earth_aDate2 = datetime(2052,06,25,11,58,00);

%% Initialize Iterative Variables
% Vary ToF; fix depart date
ToF12_dayidx = 500:1800; % [day]
earth_aDate2idx = jupiter_dDate1+zeros(1,length(ToF12_dayidx)); % [date time] create a vector of datetime
earth_aDate2idx = ToF12_dayidx+earth_aDate2idx;


deltaV = zeros(1,length(ToF12_dayidx));

%% load J2000 Data
load("J2000.mat")
jdate_index = find(J2000(5).date==jupiter_dDate1);
r_jd1 = J2000(5).r(jdate_index,:)'; % return the entire r vector
v_jd1 = J2000(5).v(jdate_index,:)'; % return the entire v vector
T_jd1 = J2000(5).T(jdate_index); % return the time in days since J2000

% Build earth arrival vectors
r_ea2 = zeros(3,length(ToF12_dayidx));
v_ea2 = zeros(3,length(ToF12_dayidx));
T_ea2 = zeros(3,length(ToF12_dayidx));
edate_index = zeros(1,length(ToF12_dayidx));
for i = 1:length(ToF12_dayidx)
    edate_index(1,i) = find(J2000(3).date==earth_aDate2idx(1,i));
    r_ea2(:,i) = J2000(3).r(edate_index(1,i),:); % return the entire r vector
    v_ea2(:,i) = J2000(3).v(edate_index(1,i),:); % return the entire v vector
    T_ea2(:,i) = J2000(3).T(edate_index(1,i)); % return the time in days since J2000
end


%% Part 1: Jupiter 2 Earth
r1 = r_jd1;
r2 = r_ea2;
ToF12_day = T_ea2-T_jd1;
ToF12 = ToF12_day.*day2TU;
v1 = zeros(3,length(ToF12_dayidx));
v2 = zeros(3,length(ToF12_dayidx));
v2_unit = zeros(3,length(ToF12_dayidx));
v2_mag = zeros(1,length(ToF12_dayidx));
theta_s2 = zeros(1,length(ToF12_dayidx));
a_s1 = zeros(1,length(ToF12_dayidx));
e_s1 = zeros(1,length(ToF12_dayidx));
i_s1 = zeros(1,length(ToF12_dayidx));
O_s1 = zeros(1,length(ToF12_dayidx));
o_s1 = zeros(1,length(ToF12_dayidx));
theta_s1 = zeros(1,length(ToF12_dayidx));
for i = 1:length(ToF12_dayidx)
    [v1(:,i), v2(:,i)] = GaussToF(r1',r2(:,i)',ToF12_dayidx(1,i)*day2TU);
    v2_mag(1,i) = sqrt(dot(v2(:,i),v2(:,i))); % [AU/TU]
    v2_unit(:,i) = v2(:,i)./v2_mag(1,i);
    [a_s1(1,i), e_s1(1,i), i_s1(1,i), O_s1(1,i), o_s1(1,i), theta_s1(1,i)] = orbitalElements(r1',v1(:,i)',1);
    theta_s2(1,i) = thetaFromrv(r2(:,i),v2(:,i));
end

%% Part 2: Earth Flyby
v_infea1AU = v2-v_ea2; % [AU/TU]
v_infea1_ijk = v_infea1AU.*AUTU2kms; % [km/s] convert to km/s
v_infea1_ijk_mag = sqrt(dot(v_infea1_ijk,v_infea1_ijk)); % [km/s]
mu_e = 3.986e5; % [km^3/s^2]
E_hyper1 = v_infea1_ijk_mag.^2./2;
a_hyper1 = -mu_e./(2.*E_hyper1);
% recomend between 1000km and GEO altitude (No satelites to run into)
rp_hyper1 = 1000+6378; % [km]
e_hyper1 = 1-rp_hyper1./a_hyper1;
delta = -2.*asin(1./e_hyper1); % rad


infa2infd = zeros(3,3,length(ToF12_dayidx));
r_soie = zeros(3,length(ToF12_dayidx));
v_infed1_ijk = zeros(3,length(ToF12_dayidx));
v_infed1_mag = zeros(1,length(ToF12_dayidx));
v3 = zeros(3,length(ToF12_dayidx));
v3_mag = zeros(1,length(ToF12_dayidx));
for i = 1:length(ToF12_dayidx)
    infa2infd(:,:,i) = [cos(delta(1,i)) -sin(delta(1,i)) 0;
                        sin(delta(1,i))  cos(delta(1,i)) 0;
                        0                0               1];

    % convert to perifocal to calc flyby: already have e and a above, so find h
    % to get i, O, and o
    R_soie = J2000(3).a*(150.47e6)*(5.972e24/1.989e30)^(2/5); % [km]
    r_soie = R_soie.*(-v2_unit(:,i)); % [km]

    [a_flyby,e_flyby,i_flyby,O_flyby,o_flyby,nu_flyby] = orbitalElements(r_soie(:,i)',v_infea1_ijk(:,i)',mu_e);
    [~,~,R_per2equ] = perifocal2Equatorial([0;0;0],[0;0;0],i_flyby,O_flyby,O_flyby);

    v_infea1_pkw = R_per2equ\v_infea1_ijk(:,i);
    v_infea1_pkw_mag = sqrt(dot(v_infea1_pkw,v_infea1_pkw));
    v_infed1_pkw = infa2infd(:,:,i)*v_infea1_pkw; % [km/s]
    v_infed1_pkw_mag = sqrt(dot(v_infed1_pkw,v_infed1_pkw));
    v_infed1_ijk(:,i) = R_per2equ*v_infed1_pkw; % [km/s]


    v_infed1_mag(:,i) = sqrt(dot(v_infed1_ijk(:,i),v_infed1_ijk(:,i))); % [km/s]
    v_infed1AU = v_infed1_ijk(:,i).*(1/AUTU2kms); % [AU/TU]
    v3(:,i) = v_infed1AU+v_ea2; % [AU/TU]
    v3_mag(:,i) = sqrt(dot(v3(:,i),v3(:,i))); % [AU/TU]
end

r3 = r2;

%% Part 3: Deep Space Manuevers

theta_s3 = zeros(1,length(ToF12_dayidx));
a_s3 = zeros(1,length(ToF12_dayidx));
e_s3 = zeros(1,length(ToF12_dayidx));
i_s3 = zeros(1,length(ToF12_dayidx));
O_s3 = zeros(1,length(ToF12_dayidx));
o_s3 = zeros(1,length(ToF12_dayidx));

for i = 1:length(ToF12_dayidx)
    [a_s3(1,i), e_s3(1,i), i_s3(1,i), O_s3(1,i), o_s3(1,i), theta_s3(1,i)] = orbitalElements(r3(:,i),v3(:,i),1);

    E_s3 = acos((e_s3+cosd(theta_s3))/(1+e_s3*cosd(theta_s3))); % eccentric anomaly post fly-by
    E_s4 = pi; % [rad] eccentric anomaly at apogee
    theta_s4_desired = 180;
    if theta_s3(1,i)>180
        k=1;
        E_s3 = 2*pi-E_s3;
    elseif theta_s3<=180
        k=0;
    else
        disp("your code is broken")
    end

    ToF34 = sqrt(a_s3^3)*(2*pi*k+(E_s4-e_s3*sin(E_s4))-(E_s3-e_s3*sin(E_s3))); % TU from fly-by to apogee of post-flyby
    ToF34_day = ToF34*TU2day;

    [r4,v4_deepSpaceTraj,theta_s4] = keplarToF(a_s3,e_s3,i_s3,O_s3,o_s3,r3,v3,theta_s3,getTU(floor(ToF34_day),'d'));
    r4=r4';
end


%% Calc Delta V1
r_cj = 10000+69911; % [km]
mu_j = 126.687e6; % [km^3/s^2]
v_cj = sqrt(mu_j/r_cj); % [km/s]
v_infj = (v_jd1-v1)*AUTU2kms; % [km/s]
v_infj_mag = sqrt(dot(v_infj,v_infj)); % [km/s]
v_hyperj = sqrt((v_infj_mag^2+2*mu_j/r_cj)); % [km/s]
deltaV1 = v_hyperj-v_cj; % [km/s]

%% Calc Delta V2 For Hohman
ToF45_day_min = 300; % [days] avoid imaginary solutions
ToF45_day_max = 500; % [days] avoid multiple laps
ToF45_day_idx = ToF45_day_min:ToF45_day_max;
% initialize variables for storage
deltaV2_iter = zeros(length(ToF45_day_idx),1);
deltaV3_iter = zeros(length(ToF45_day_idx),1);
for idx = 1:length(ToF45_day_idx)

    ToF45_day = ToF45_day_idx(idx); % days from apogee of post-flyby orbit to earth arrival
    ToF45 = getTU(ToF45_day,'d'); % [TU]

    ToF35 = ToF34+ToF45;
    ToF35_day = floor(ToF34_day+ToF45_day); % nearest whole day

    earth_aDate5 = daysadd(earth_aDate2,ToF35_day); % date of second earth arrival post-flyby
    edate_index = find(J2000(3).date==earth_aDate5);
    r_ea5 = J2000(3).r(edate_index,:); % return the entire r vector
    v_ea5 = J2000(3).v(edate_index,:); % return the entire v vector
    T_ea5 = J2000(3).T(edate_index); % return the time in days since J2000

    r5 = r_ea5;

    [v4,v5] = GaussToF(r4,r5,ToF45); % this v4 is the required v4 to acheive desired orbit
    % need v4_deepspaceTraj to find deltaV
    [a_s4, e_s4, i_s4, O_s4, o_s4, theta_s4prime] = orbitalElements(r4,v4,1);

    deltaV2_AUTUvector = v4_deepSpaceTraj'-v4; % [AU/TU]
    deltaV2_AUTUmag = sqrt(dot(deltaV2_AUTUvector,deltaV2_AUTUvector)); % [AU/TU]
    deltaV2_iter(idx) = deltaV2_AUTUmag*AUTU2kms; % [km/s]

    r_ce = 2000+6378; % [km]
    v_ce = sqrt(mu_e/r_ce); % [km/s]
    v_infea2 = (v_ea5-v5)*AUTU2kms; % [km/s]
    v_infea2_mag = sqrt(dot(v_infea2,v_infea2)); % [km/s]
    v_hyperea2 = sqrt((v_infea2_mag^2+2*mu_e/r_ce)); % [km/s]
    deltaV3_iter(idx) = v_hyperea2-v_ce; % [km/s]

end

%% Calc Delta V3
[deltaV3, idx_min] = min(deltaV3_iter);
deltaV2 = deltaV2_iter(idx_min);
ToF45_day = ToF45_day_idx(idx_min);

ToF45 = getTU(ToF45_day,'d'); % [TU]

ToF35 = ToF34+ToF45;
ToF35_day = floor(ToF34_day+ToF45_day); % nearest whole day

earth_aDate5 = daysadd(earth_aDate2,ToF35_day); % date of second earth arrival post-flyby
edate_index = find(J2000(3).date==earth_aDate5);
r_ea5 = J2000(3).r(edate_index,:); % return the entire r vector
v_ea5 = J2000(3).v(edate_index,:); % return the entire v vector
T_ea5 = J2000(3).T(edate_index); % return the time in days since J2000
T_ed3 = T_ea2;

r5 = r_ea5;

[v4,v5] = GaussToF(r4,r5,ToF45); % this v4 is the required v4 to acheive desired orbit
% need v4_deepspaceTraj to find deltaV
[a_s4, e_s4, i_s4, O_s4, o_s4, theta_s4prime] = orbitalElements(r4,v4,1);

T_deepspaceTraj = sqrt(4*pi^2*a_s3^3); % [TU] period of post-flyby helio orbit
T_deepspaceTraj_day = T_deepspaceTraj*TU2day;;
T_deepspaceTraj_dayidx = 1:1:T_deepspaceTraj*TU2day;

T_deepspaceMan = sqrt(4*pi^2*a_s4^3); % [TU] period of post-flyby helio orbit
T_deepspaceMan_day = T_deepspaceMan*TU2day;
T_deepspaceMan_dayidx = 1:1:T_deepspaceMan*TU2day;

T_e2e = T_ed3:1:T_ea5;
T_e24 = T_ed3:1:T_ed3+floor(ToF34_day);

deltaV1
deltaV2
deltaV3