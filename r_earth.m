function r = r_earth(units)
%Returns the radius of Earth
    if nargin == 0
        r = 6378.136; %km
    elseif units == "english"
        r = 3963.190; %mi
    else
        r = 6378.136; %km
    end

end

