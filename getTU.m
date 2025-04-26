function TU = getTU(t, from)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% define unit conversions
TU2s = 5.023e6; % [s/TU]
s2TU = 1/TU2s; % [TU/s]
min2TU = s2TU*60; % [TU/min]
hr2TU = min2TU*60; % [TU/hr]
day2TU = hr2TU*24; % [TU/day]

if isempty(from)
    TU = t*day2TU;
elseif from == 'hr'
    TU = t*hr2TU;
elseif from == 'min'
    TU = t*min2TU;
elseif from == 's'
    TU = t*s2TU;
elseif from == 'd'
    TU = t*day2TU;
else
    TU = t*day2TU;
end