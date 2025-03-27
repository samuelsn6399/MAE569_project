function S = S_ToF(z)
%S_TOF calucluates S given Z and a for the universal time of flight
%formulation
S = zeros(size(z));
for i = 1:length(z)
    if z(i) > 0.01
        S(i) = (sqrt(z(i))-sin(sqrt(z(i))))/sqrt(z(i)^3);
    elseif z(i) < -0.01
        S(i) = (sinh(sqrt(-z(i))) - sqrt(-z(i)))/sqrt((-z(i))^3);
    else
        S(i) = 1/factorial(3);
        for k = 1:50
            S(i) = S(i) + ((-z(i))^k)/factorial(2*k + 3);
        end
    end
end
end



