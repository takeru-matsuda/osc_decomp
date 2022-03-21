load('NorthSouthSunspotData');
Y = log(dat+1);
Y = Y-mean(Y,2)*ones(1,size(Y,2));
fs = 12;
MAX_OSC = 6;
MAX_VAR = 20;
tic
[osc_param1,osc_AIC1] = osc_decomp(Y,fs,MAX_OSC,MAX_VAR,'quasi-newton',false);
time1=toc
tic
[osc_param2,osc_AIC2] = osc_decomp(Y,fs,MAX_OSC,MAX_VAR,'quasi-newton',true);
time2=toc
tic
[osc_param3,osc_AIC3] = osc_decomp(Y,fs,MAX_OSC,MAX_VAR,'trust-region',true);
time3=toc
