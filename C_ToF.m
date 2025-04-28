function C = C_ToF(z)
%C_ToF calucluates C given Z and a for the universal time of flight
%formulation
C = zeros(size(z));
for i = 1:length(z)
    if z(i) > 0.01
        C(i) = (1-cos(sqrt(z(i))))/z(i);
    elseif z(i) < -0.01
        C(i) = (1-cosh(sqrt(-z(i))))/z(i);
    else
        C(i) = 1/factorial(2);
        for k = 1:50
            C(i) = C(i) + ((-z(i))^k)/factorial(2*k + 2);
        end
    end
end
end

