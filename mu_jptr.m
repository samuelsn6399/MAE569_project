function r = mu_jptr(units)
%Returns the gravitational parameter of Earth
    if nargin == 0
        r = 1.26686534e8; %km^3/s^2
    elseif units == "english"
        error("somebody wants english units. lol.\n")
        r = 1.407646882e16; %ft^3/s^2
    else
        r = 1.26686534e8; %km^3/s^2
    end

end