clc,clf, clear

global CO2Emissions CO2ConcRCP45

beta = 0.35;
NPP_0 = 60; % Fråga om denna!
len = 736;
conversion_factor = 0.469;


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

NPP = @(B1) NPP_0 * (1 + beta * log(B1/BF(1)));

U = @(t) CO2Emissions(t);

f = @(t,B) [
    alpha(3,1)*B(3) + alpha(2,1)*B(2) - NPP(B(1)) + U(t); 
    NPP(B(1)) - alpha(2,3)*B(2) - alpha(2,1)*B(2); 
    alpha(2,3)*B(2) - alpha(3,1)*B(3);
];


B = zeros([3 len]);
B(:,1) = BF;
T = len;
span = 1:len;

for i = 1:T-1
    B(:,i+1) = B(:,i) + f(i,B(:,i));
end

figure(1)
hold on
plot(span,B(1,:));
plot(span,B(2,:));
plot(span,B(3,:));
hold off

legend(["B1" "B2" "B3"])

figure(2)
hold on
plot(span, B(1,:) * conversion_factor);
plot(span, CO2ConcRCP45)
legend(["Calculated" "From data"])
hold off



