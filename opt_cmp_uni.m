load('CanadianLynxData');
y = log(lynx);
y = y-mean(y);
fs = 1;
MAX_OSC = 6;
MAX_AR = 20;
tic
[osc_param1,osc_AIC1] = osc_decomp(y,fs,MAX_OSC,MAX_AR,'quasi-newton',false);
time1=toc
tic
[osc_param2,osc_AIC2] = osc_decomp(y,fs,MAX_OSC,MAX_AR,'quasi-newton',true);
time2=toc
tic
[osc_param3,osc_AIC3] = osc_decomp(y,fs,MAX_OSC,MAX_AR,'trust-region',true);
time3=toc
