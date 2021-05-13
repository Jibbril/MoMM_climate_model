clc, clf, clear

global CO2Emissions CO2ConcRCP45

i = 1:5;
A = [0.113 0.213 0.258 0.273 0.143];
tau = [2 12.2 50.4 243.3 Inf];
M_0 = 600;
M_0s = [0 140 560 1680];
len = length(CO2Emissions);
conversion_factor = 0.469;


Ms = zeros([len 1]);
for i = 0:len-1
    Ms(i+1) = M(i, M_0, A, tau, CO2Emissions);
end

hold on
plot(0:len-1, Ms*conversion_factor);
plot(0:len-1, CO2ConcRCP45)
legend(["Calculated"  "From data"])


function res = tau_i(tau_0,t,CO2Emissions)
    k = 3.06*10^(-3);
    U = 0;
    for i=0:t-1
        U = U + CO2Emissions(i+1);
    end
res = tau_0 * (1 + k * U);
end

function res = I(t,t_tilde, A, tau, CO2Emissions)
    sum = 0;
    for i=1:5
        sum = sum +  A(i)*exp(-t/tau_i(tau(i), t_tilde, CO2Emissions));
    end
    
    res =  sum;
end

function res = M(t, M_0, A, tau, CO2Emissions)
    sum = M_0;
    for t_tilde=0:t
        sum = sum + I(t-t_tilde,t, A, tau, CO2Emissions) * CO2Emissions(t_tilde + 1);
    end
    res = sum;
end
