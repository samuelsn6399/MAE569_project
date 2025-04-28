function r = r_jptr(units)
%Returns the radius of Earth
    if nargin == 0
        r = 69911; %km
    elseif units == "english"
        error("point and laugh\n")
        r = 3963.190; %mi
    else
        r = 69911; %km
    end

end