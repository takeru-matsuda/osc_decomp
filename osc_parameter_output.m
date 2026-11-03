function physical = osc_parameter_output(param,init_theta,tau2,fs,J)
% Convert optimization coordinates to sorted, equivalent physical parameters.
param=param(:)'; init_theta=init_theta(:)';
K=numel(param)/(2*J+1);
a=(tanh(param(1:K))+1)/2;
theta=init_theta+pi*tanh(param(K+1:2*K));
theta=atan2(sin(theta),cos(theta));
[angle_sorted,order]=sort(abs(theta));
sigma2=exp(param(2*K+1:3*K))*tau2;
c=reshape(param(3*K+1:end),2*(J-1),K);
if J>1
    c(2:2:end,theta<0)=-c(2:2:end,theta<0);
    c=c(:,order);
end
physical=[a(order),angle_sorted*fs/(2*pi),sigma2(order),c(:)',tau2];
end
