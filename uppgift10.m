%% a)
clc, clf, clear
lambda = 0.8; % 0.5-1.3 K/Wm^2
k = 0.5; % 0.2-1 W/Km^2
RF = 1;
c = 4186; %J/kgK
rho = 1020; %kg/m^3
h = 50; %m
d = 2000; %m
C_1 = c*h*rho / (60*60*24*365); % W*yr/Km^2
C_2 = c*d*rho / (60*60*24*365); % W*yr/
T = 5000;


f = @(t,delta_T) [
    1/C_1 * (RF - delta_T(1,t)/lambda - k * (delta_T(1,t) - delta_T(2,t)));
    1/C_2 * k * (delta_T(1,t) - delta_T(2,t)) 
];


delta_T = zeros([2 T]);
span = 1:T;

for i = 1:T-1
    delta_T(:,i+1) = delta_T(:,i) + f(i,delta_T);
end

subplot(1,2,1)
hold on
axis([0 T 0 1.5])
plot(span, delta_T(1,:))
plot(span, delta_T(2,:))
legend(["$\Delta T_1$" "$\Delta T_2$"], "interpreter", "latex")



e_folding = 1-exp(-1);
delta_T0 = delta_T - e_folding * delta_T(1, T);

subplot(1,2,2)
hold on 
axis([0 T -0.5 1.5])
plot(span, delta_T0(1,:))
plot(span, delta_T0(2,:))
legend(["$\Delta T_1$" "$\Delta T_2$"], "interpreter", "latex")

t1 = 0;
t2 = 0;
for i=1:T
    if delta_T0(1,i) > 0 && t1 == 0
        t1 = i;
    end
    if delta_T0(2,i) > 0 && t2 == 0
        t2 = i;
        break
    end
end

disp("E-folding for top box: " + num2str(t1) + " years")
disp("E-folding for bottom box: " + num2str(t2) + " years")


final_radiation = delta_T(1,T)/lambda;
disp("Final radiation: " + num2str(final_radiation))
disp("Radiative forcing: " + num2str(RF))


%% b)
clc, clf, clear
lambda = 0.8; % 0.5-1.3 K/Wm^2
k = 0.5; % 0.2-1 W/Km^2
RF = 1;
c = 4186; %J/kgK
rho = 1020; %kg/m^3
h = 50; %m
d = 2000; %m
C_1 = c*h*rho / (60*60*24*365); % W*yr/Km^2
C_2 = c*d*rho / (60*60*24*365); % W*yr/
T = 6000;


f = @(t,delta_T) [
    1/C_1 * (RF - delta_T(1,t)/lambda - k * (delta_T(1,t) - delta_T(2,t)));
    1/C_2 * k * (delta_T(1,t) - delta_T(2,t)) 
];


delta_T = zeros([2 T]);
span = 1:T;

for i = 1:T-1
    delta_T(:,i+1) = delta_T(:,i) + f(i,delta_T);
end

subplot(1,2,1)
hold on
axis([0 T 0 1.5])
plot(span, delta_T(1,:))
plot(span, delta_T(2,:))
legend(["$\Delta T_1$" "$\Delta T_2$"], "interpreter", "latex")



e_folding = 1-exp(-1);
delta_T0 = delta_T - e_folding * delta_T(1, T);

subplot(1,2,2)
hold on 
axis([0 T -0.5 1.5])
plot(span, delta_T0(1,:))
plot(span, delta_T0(2,:))
legend(["$\Delta T_1$" "$\Delta T_2$"], "interpreter", "latex")

t1 = 0;
t2 = 0;
for i=1:T
    if delta_T0(1,i) > 0 && t1 == 0
        t1 = i;
    end
    if delta_T0(2,i) > 0 && t2 == 0
        t2 = i;
        break
    end
end

disp("E-folding for top box: " + num2str(t1) + " years")
disp("E-folding for bottom box: " + num2str(t2) + " years")


final_radiation = delta_T(1,T)/lambda;
disp("Final radiation: " + num2str(final_radiation))
disp("Radiative forcing: " + num2str(RF))

%% c)
clc, clf, clear
lambda = 0.8; % 0.5-1.3 K/Wm^2
k = 0.5; % 0.2-1 W/Km^2
RF = 1;
c = 4186; %J/kgK
rho = 1020; %kg/m^3
h = 50; %m
d = 2000; %m
C_1 = c*h*rho / (60*60*24*365); % W*yr/Km^2
C_2 = c*d*rho / (60*60*24*365); % W*yr/
T = 200;


f = @(t,delta_T) [
    1/C_1 * (RF - delta_T(1,t)/lambda - k * (delta_T(1,t) - delta_T(2,t)));
    1/C_2 * k * (delta_T(1,t) - delta_T(2,t)) 
];


delta_T = zeros([2 T]);
span = 1:T;

for i = 1:T-1
    delta_T(:,i+1) = delta_T(:,i) + f(i,delta_T);
end

subplot(1,2,1)
hold on
axis([0 T 0 1.5])
plot(span, delta_T(1,:))
plot(span, delta_T(2,:))
legend(["$\Delta T_1$" "$\Delta T_2$"], "interpreter", "latex")



e_folding = 1-exp(-1);
delta_T0 = delta_T - e_folding * delta_T(1, T);

subplot(1,2,2)
hold on 
axis([0 T -0.5 1.5])
plot(span, delta_T0(1,:))
plot(span, delta_T0(2,:))
legend(["$\Delta T_1$" "$\Delta T_2$"], "interpreter", "latex")

t1 = 0;
t2 = 0;
for i=1:T
    if delta_T0(1,i) > 0 && t1 == 0
        t1 = i;
    end
    if delta_T0(2,i) > 0 && t2 == 0
        t2 = i;
        break
    end
end

disp("E-folding for top box: " + num2str(t1) + " years")
disp("E-folding for bottom box: " + num2str(t2) + " years")


final_radiation = delta_T(1,T)/lambda;
disp("Final radiation: " + num2str(final_radiation))
disp("Radiative forcing: " + num2str(RF))





