clc, clf, clear

global totRadForcAerosols totRadForcExclCO2AndAerosols

s = 1;
RFAerosols = s * totRadForcAerosols;

len = length(totRadForcAerosols);
span = linspace(1,len,len);



hold on
plot(span, RFAerosols)
plot(span,totRadForcExclCO2AndAerosols)
plot(span,RFAerosols + totRadForcExclCO2AndAerosols)
legend(["RF Aerosols" "RF Other" "RF Aerosols and Other"])
