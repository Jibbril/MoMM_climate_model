clc, clf, clear

i = 1:5;
A = [0.113 0.213 0.258 0.273 0.143];
tau = [2 12.2 50.4 243.3 Inf];
M_0 = 600;
M_0s = [0 140 560 1680];


Is = zeros([1 500]);
hold on
for i=1:length(M_0s)
    for j = 1:length(Is)
        Is(j) = I(j, A, tau, M_0s(i));
    end
    plot(0:length(Is)-1, Is);
end
legend(["0 GtC", "140 GtC", "560 GtC", "1680 GtC"])



function res = tau_i(tau_0, M_0)
    k = 3.06*10^(-3);
    res = tau_0 * (1 + k * M_0);
end

function res = I(t, A, tau, M_0)
    sum = 0;
    for i=1:5
        sum = sum +  A(i)*exp(-t/tau_i(tau(i), M_0));
    end
    
    res =  sum;
end