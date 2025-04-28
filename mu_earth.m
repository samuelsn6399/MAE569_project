function r = mu_earth(units)
%Returns the gravitational parameter of Earth
    if nargin == 0
        r = 3.986004418e5; %km^3/s^2
    elseif units == "english"
        r = 1.407646882e16; %ft^3/s^2
    else
        r = 3.986004418e5; %km^3/s^2
    end

end

