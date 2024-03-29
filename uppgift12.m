%% a)
clc, clf, clear
plot_scenarios();

%% b)
clc, clf, clear
global TAnomali

% Ska utg� ifr�n f�rindustriell tid, resonera kring l�mligt val av

TMax = 436;
T = min(TMax,TMax);
nasa_data_start = 115;
conversion_factor = 0.469;
lambda = 0.8; % 0.5-1.3 K/Wm^2
k = 0.5; % 0.2-1 W/Km^2
s = 1;
reference_period = [85 135]; %1850-1900

CO2_scenario_i = create_co2_scenario_i();
CO2_scenario_ii = create_co2_scenario_ii();
CO2_scenario_iii = create_co2_scenario_iii();

delta_T1 = run_scenario(TMax, CO2_scenario_i,lambda,k,s,reference_period);
delta_T2 = run_scenario(TMax, CO2_scenario_ii,lambda,k,s,reference_period);
delta_T3 = run_scenario(TMax, CO2_scenario_iii,lambda,k,s,reference_period);

hold on
plot(1:T,delta_T1(1:T))
plot(1:T,delta_T2(1:T))
plot(1:T,delta_T3(1:T))
plot(nasa_data_start:nasa_data_start+length(TAnomali)-1,TAnomali)

legend(["Scenario i" "Scenario ii" "Scenario iii" "Nasa data"], "Location", "NorthWest")


%% c)
clc, clf, clear
global TAnomali

TMax = 336;
T = min(TMax,TMax);
nasa_data_start = 115;
conversion_factor = 0.469;
lambda = 0.8; % 0.5-1.3 K/Wm^2
k = 0.5; % 0.2-1 W/Km^2
s = 1;
reference_period = [85 135]; %1850-1900

CO2_scenario_i = create_co2_scenario_i();
CO2_scenario_ii = create_co2_scenario_ii();
CO2_scenario_iii = create_co2_scenario_iii();

delta_T1 = run_scenario(TMax, CO2_scenario_i,lambda,k,s,reference_period);
delta_T2 = run_scenario(TMax, CO2_scenario_ii,lambda,k,s,reference_period);
delta_T3 = run_scenario(TMax, CO2_scenario_iii,lambda,k,s,reference_period);

hold on
axis([0 350 -1.5 4])
plot(1:T,delta_T1(1:T))
plot(1:T,delta_T2(1:T))
plot(1:T,delta_T3(1:T))
plot(nasa_data_start:nasa_data_start+length(TAnomali)-1,TAnomali)

legend(["Scenario i" "Scenario ii" "Scenario iii" "Nasa data"], "Location", "NorthWest")

end_val_i = delta_T1(end);
end_val_ii = delta_T2(end);
end_val_iii = delta_T3(end);

disp("End value i: " + num2str(end_val_i))
disp("End value ii: " + num2str(end_val_ii))
disp("End value iii: " + num2str(end_val_iii))
disp("Temperature span: " + num2str(end_val_iii-end_val_i))




function plot_scenarios()
    CO2_scenario_i = create_co2_scenario_i();
    CO2_scenario_ii = create_co2_scenario_ii();
    CO2_scenario_iii = create_co2_scenario_iii();
    len = length(CO2_scenario_i);
    
    hold on 
    plot(1:len,CO2_scenario_i);
    plot(1:len,CO2_scenario_ii);
    plot(1:len,CO2_scenario_iii);
    
    legend(["Scenario i" "Scenario ii" "Scenario iii"], "Location", "NorthWest")
end

function res = run_scenario(TMax, CO2_scenario,lambda,k,s, reference_period)
    BF = [600 600 1500];
    %TStart = min(length(TAnomali), max(nasa_data_start, 500));

    co2 = co2_model(BF, CO2_scenario);
    RF = calculate_RF(co2(1,:),BF(1));
    RF = add_aerosol_and_other_RF(RF,s);
    delta_T = heat_model(RF,TMax,lambda,k);

    avg = mean(delta_T(reference_period(1):reference_period(2)));
    delta_T = delta_T - avg;
    res = delta_T;
end

function res = heat_model(RF,T, lambda,k)
c = 4186; %J/kgK
rho = 1020; %kg/m^3
h = 50; %m
d = 2000; %m
C_1 = c*h*rho / (60*60*24*365); % W*yr/Km^2
C_2 = c*d*rho / (60*60*24*365); % W*yr/


f = @(t,delta_T, RF) [
    1/C_1 * (RF - delta_T(1,t)/lambda - k * (delta_T(1,t) - delta_T(2,t)));
    1/C_2 * k * (delta_T(1,t) - delta_T(2,t)) 
];


delta_T = zeros([2 T]);

for i = 1:T-1
    delta_T(:,i+1) = delta_T(:,i) + f(i,delta_T, RF(i));
end

res = delta_T(1,:);
end

function res = co2_model(BF,CO2_scenario)
A = [0.113 0.213 0.258 0.273 0.143];
tau = [2 12.2 50.4 243.3 Inf];
len = length(CO2_scenario);
k = 3.06*10^(-3);
beta = 0.5;

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



U = @(t) CO2_scenario(t);

f = @(t,B) [
    alpha(3,1)*B(3,t) + alpha(2,1)*B(2,t) - NPP(B(1,t), beta) + U(t); 
    NPP(B(1,t), beta) - alpha(2,3)*B(2,t) - alpha(2,1)*B(2,t); 
    alpha(2,3)*B(2,t) - alpha(3,1)*B(3,t)
];
 


B = zeros([3 len]);
B(:,1) = BF;
dB = zeros([3 len]);
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
    
    B_water(t+1) = sum(CO2_scenario(1:t)) - (sum(B(:,t+1)) - sum(BF));
end

res = B;

end

function res = add_aerosol_and_other_RF(RF,s)
global totRadForcAerosols totRadForcExclCO2AndAerosols

RFAerosols = s * totRadForcAerosols;
res = RF + RFAerosols + totRadForcExclCO2AndAerosols;
end

function res = create_co2_scenario_i()
    global CO2Emissions

    CO2 = CO2Emissions;

    start_index = 256;
    end_index = start_index + 49;
    start_val = CO2Emissions(start_index);

    linear_decrease = @(t) start_val*(1 - 1/49*(t-start_index));

    for i = start_index:end_index
        CO2(i) = linear_decrease(i);
    end

    CO2(end_index:end) = 0;
    res = CO2;
end

function res = create_co2_scenario_ii()
    global CO2Emissions

    CO2 = CO2Emissions;

    start_index = 256;
    start_val = CO2Emissions(start_index);
    
    CO2(start_index:end) = start_val;
    res = CO2;
end

function res = create_co2_scenario_iii()
    global CO2Emissions
    CO2 = CO2Emissions;

    start_index = 256;
    end_index = start_index + 79;
    start_val = CO2Emissions(start_index);
    
    % Ska vara g�nger 2.5 inte 1.5
    
    linear_increase = @(t) start_val *(1  + 1.5*(t-start_index)/79);

    for i = start_index:end_index
        CO2(i) = linear_increase(i);
    end

    CO2(end_index:end) = CO2(end_index);
    
    res = CO2;
    
end

function res = calculate_RF(concentrations,baseline)
    res = 5.35 * log(concentrations/baseline);
end

function res = NPP(B1, beta)
    NPP_0 = 60;
    BF = 600;
    res = NPP_0 * (1 + beta * log(B1/BF));
end