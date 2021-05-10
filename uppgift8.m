clc, clf, clear

global CO2ConcRCP45 CO2RadForc

len = length(CO2ConcRCP45);
span = linspace(1,len,len);
pCO2_0 = CO2ConcRCP45(1);

RFs = zeros([1 len]);

for i=1:len 
    RFs(i) = 5.35 * log(CO2ConcRCP45(i)/pCO2_0);
end

hold on
plot(span, RFs)
plot(span,CO2RadForc)
legend(["Calculated" "CO2RadForc"])
