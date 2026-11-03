function S = osc_component_spectrum(a,theta,H,freq)
% Spectral density per radian/sample, with unit innovation variance.
freq=freq(:);
A=abs(1-a*exp(1i*(theta-freq))).^(-2);
B=abs(1-a*exp(1i*(theta+freq))).^(-2);
J=size(H,1); S=zeros(J,J,numel(freq));
for i=1:numel(freq)
    S(:,:,i)=H*[A(i)+B(i),1i*(A(i)-B(i)); ...
        -1i*(A(i)-B(i)),A(i)+B(i)]*H'/(4*pi);
end
end
