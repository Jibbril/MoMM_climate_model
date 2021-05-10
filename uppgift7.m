clc,clf, clear

global CO2Emissions CO2ConcRCP45

NPP_0 = 60;
conversion_factor = 0.469;
A = [0.113 0.213 0.258 0.273 0.143];
tau = [2 12.2 50.4 243.3 Inf];
M_0 = 600;
M_0s = [0 140 560 1680];
len = length(CO2Emissions);
k = 3.06*10^(-3);
beta = 0.5;

BF = [600 600 1500];
F = [
    0 60 0;
    15 0 45;
    45 0 0
];

alpha = [
    F(1,1)/BF(1) F(1,2)/BF(1) F(1,3)/BF(1);
    F(2,1)/BF(2) F(2,2)/BF(2) F(2,3)/BF(2);
    F(3,1)/BF(3) F(3,2)/BF(3) F(3,3)/BF(3);
];



U = @(t) CO2Emissions(t);

f = @(t,B) [
    alpha(3,1)*B(3,t) + alpha(2,1)*B(2,t) - NPP(B(1,t), beta) + U(t); 
    NPP(B(1,t), beta) - alpha(2,3)*B(2,t) - alpha(2,1)*B(2,t); 
    alpha(2,3)*B(2,t) - alpha(3,1)*B(3,t)
];
 


B = zeros([3 len]);
B(:,1) = BF;
dB = zeros([3 len]);
span = 1:len;
Tau = zeros([1 5]);
I_ = zeros([1 5]);
B_water = zeros([1 len]);

for t = 1:len-1
    dB(:,t) = f(t,B);
    IU = 0;
    for i=1:t
        for j=1:5
            Tau(j) = tau(j) * (1 + k * sum(dB(1,1:t-1)));
            I_(j) = A(j) * exp(-(t-i)/Tau(j));
        end
        IU = IU + sum(I_) * dB(1,i);
    end
    
    B(:, t+1) = dB(:, t) + [BF(1) B(2,t) B(3,t)]' + [IU 0 0]';
    
    B_water(t+1) = sum(CO2Emissions(1:t)) - (sum(B(:,t+1)) - sum(BF));
end

B(1,:) = B(1,:) - BF(1);
B(2,:) = B(2,:) - BF(2);
B(3,:) = B(3,:) - BF(3);
deltaCO2ConcRCP45 = CO2ConcRCP45 - CO2ConcRCP45(1);


subplot(1,2,1)
axis([0 800 0 400])
hold on
plot(span,B(1,:) * conversion_factor);
plot(span,B(2,:) * conversion_factor);
plot(span,B(3,:) * conversion_factor);
plot(span,B_water * conversion_factor);
legend(["B_{air}" "B_{plants}" "B_{ground}" "B_{water}"], "Location", "northwest")
hold off


subplot(1,2,2)
axis([0 800 0 400])
hold on
plot(span, B(1,:) * conversion_factor);
plot(span, deltaCO2ConcRCP45)
legend(["Calculated" "From data"])
hold off


function res = NPP(B1, beta)
    NPP_0 = 60;
    BF = 600;
    res = NPP_0 * (1 + beta * log(B1/BF));
end