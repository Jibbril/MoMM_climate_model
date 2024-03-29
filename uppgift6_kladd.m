clc,clf, clear

global CO2Emissions CO2ConcRCP45

beta = 0.35;
NPP_0 = 60;
conversion_factor = 0.469;
A = [0.113 0.213 0.258 0.273 0.143];
tau = [2 12.2 50.4 243.3 Inf];
M_0 = 600;
M_0s = [0 140 560 1680];
len = length(CO2Emissions);


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
    I(1,t, A, tau, CO2Emissions, alpha, B, beta)*(alpha(3,1)*B(3,t) + alpha(2,1)*B(2,t) - NPP(B(1,t), beta) + U(t)); 
    NPP(B(1,t), beta) - alpha(2,3)*B(2,t) - alpha(2,1)*B(2,t); 
    alpha(2,3)*B(2,t) - alpha(3,1)*B(3,t)
];
% Han tyckte att man skulle l�sa det med tre ekvationer, inte skapa en 
% fj�rde box. F� till koppling mellan M och B. T�nk inte vattnet f�rst
% och sedan land, utan t�nk att l�sa systemet p� land givet utsl�pp, sen
% tar havet upp lite. Upptag mellan natur och atmosf�r g�r snabbare i 
% verkligeheten �n upptag mellan atmosf�r och vatten. Man kanske kan g�ra
% det i loopen? Typ l�sa systemet f�r en viss iteration, sen dra av det som
% ska till vattnet p� r�tt st�lle? 

% F�rmodligen r�tt nu, men gl�mde �ndra utsl�ppen i Tau funktionen ocks�!
% Har �tg�rdat en del men blir fortfarande inte bra. Debugga processen, kan
% vara n�got indexfel i tau funktionen. Annars finns den fungerande koden i
% f�reg�ende commit. 


B = zeros([3 len]);
B(:,1) = BF;
span = 1:len;

for i = 1:len-1
    B(:,i+1) = B(:,i) + f(i,B);
end


subplot(1,2,1)
axis([0 800 600 2200])
hold on
plot(span,B(1,:));
plot(span,B(2,:));
plot(span,B(3,:));
legend(["B1" "B2" "B3"])
hold off

subplot(1,2,2)
hold on
plot(span, B(1,:) * conversion_factor);
plot(span, CO2ConcRCP45)
legend(["Calculated" "From data"])
hold off


function res = NPP(B1, beta)
    NPP_0 = 60;
    BF = 600;
    res = NPP_0 * (1 + beta * log(B1/BF));
end


% Funktioner fr�n uppgift 4
function res = tau_i(tau_0,t,CO2Emissions, alpha, B, beta)
    k = 3.06*10^(-3);
    U = 0;
    for i=0:t-1
        U = U + (alpha(3,1)*B(3,i+1) + alpha(2,1)*B(2,i+1) - NPP(B(1,i+1), beta) + CO2Emissions(i+1));
        %U = U + CO2Emissions(i+1);
    end
res = tau_0 * (1 + k * U);
end

function res = I(t,t_tilde, A, tau, CO2Emissions, alpha, B, beta)
    sum = 0;
    for i=1:5
        sum = sum +  A(i)*exp(-t/tau_i(tau(i), t_tilde, CO2Emissions, alpha, B, beta));
    end
    
    res =  sum;
end

% Funktionen nedan �r inte l�ngre relevant utifr�n �ndringar som gjorts
% ovan. Ska den anv�ndas m�ste den formateras om och ges fler argument.
function res = M(t, M_0, A, tau, CO2Emissions)
    sum = M_0;
    for t_tilde=0:t
        sum = sum + I(t-t_tilde,t, A, tau, CO2Emissions) * CO2Emissions(t_tilde + 1);
    end
    res = sum;
end