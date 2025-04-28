function r = r_GEO(units)
%Returns the radius of Earth
    if nargin == 0
        r = 42164; %km
    elseif units == "english"
        r = 26199; %mi
    else
        r = 42164; %km
    end

end

