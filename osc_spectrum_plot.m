function []=osc_spectrum_plot(Y,fs,osc_a,osc_f,osc_sigma2,osc_tau2,osc_c)
    K = length(osc_a);
	J = size(Y,1);
    if J == 1
        osc_c = zeros(1,2*K);
        osc_c(1:2:end) = 1;
    end
	T = size(Y,2);
    H = zeros(J,2*K);
    H(1,1:2:2*K) = 1;
    kk = 1;
    for k=1:K
        for j=2:J
            H(j,2*k-1:2*k) = osc_c(kk:kk+1);
            kk = kk+2;
        end
    end
    for j=1:J
        peri = zeros(1,T/2);
        spect = zeros(K,T/2);
        for i=1:T/2
            peri(i) = abs(Y(j,:)*exp(-1i*2*pi*i/T*[1:T])')^2/2/pi/T;
            for k=1:K
                a = osc_a(k);
                theta = 2*pi*osc_f(k)/fs;
                A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
                b = (A-2*a*cos(theta)+sqrt((A-2*a*cos(theta))^2-4))/2;
                spect(k,i) = -norm(H(j,2*k-1:2*k))^2*osc_sigma2(k)*a*cos(theta)/b*abs(1+b*exp(-1i*2*pi*i/T))^2/abs(1-2*a*cos(theta)*exp(-1i*2*pi*i/T)+a^2*exp(-1i*4*pi*i/T))^2/2/pi;
            end
        end
        noise = osc_tau2/2/pi*ones(1,T/2);
        figure,hold on
        p1 = plot([1:T/2]/T*fs,log10(peri),'k--');
        p2 = plot([1:T/2]/T*fs,log10(spect(1,:)),'r-');
        for k=1:K
            plot([1:T/2]/T*fs,log10(spect(k,:)),'r-');
        end
        p3 = plot([1:T/2]/T*fs,log10(noise),'g+-');
        p4 = plot([1:T/2]/T*fs,log10(sum(spect)+noise),'b*-');
        if J == 1
            legend([p1 p2 p3 p4],{'periodogram','oscillator','noise','sum of oscillator & noise'});
        else
            legend([p1 p2 p3 p4],{sprintf('periodogram of y%d',j),'oscillator','noise','sum of oscillator & noise'});
        end
        xlabel('frequency','FontSize',20);
        ylabel('log10 power','FontSize',20);
        set(gca,'FontSize',16);
    end
end
