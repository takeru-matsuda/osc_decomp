function osc_spectrum_plot(Y,fs,osc_a,osc_f,osc_sigma2,osc_tau2,osc_c)
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
    
    m = ceil(T/2);
    freq_idx = (1:m)';
    freq_rad = 2*pi*freq_idx/T;
    freq_hz = freq_idx/T*fs;
    exp_idx = exp(-1i * freq_rad * (1:T));
    
    for j=1:J
        peri = abs(Y(j,:) * exp_idx').^2 / (2*pi*T);
        spect = zeros(K, m);
        for k=1:K
            a = osc_a(k);
            theta = 2*pi*osc_f(k)/fs;
            A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
            b = (A-2*a*cos(theta)+sqrt((A-2*a*cos(theta))^2-4))/2;
            
            H_norm2 = norm(H(j,2*k-1:2*k))^2;
            numer = abs(1 + b * exp(-1i * freq_rad)).^2;
            denom = abs(1 - 2*a*cos(theta)*exp(-1i * freq_rad) + a^2 * exp(-1i * 2 * freq_rad)).^2;
            spect(k,:) = -H_norm2 * osc_sigma2(k) * a * cos(theta) / b .* numer' ./ denom' / (2*pi);
        end
        noise = (osc_tau2 / (2*pi)) * ones(1, m);
        
        figure, hold on
        p1 = plot(freq_hz, log10(peri), 'k--');
        p2 = plot(freq_hz, log10(spect(1,:)), 'r-');
        for k=2:K
            plot(freq_hz, log10(spect(k,:)), 'r-');
        end
        p3 = plot(freq_hz, log10(noise), 'g+-');
        p4 = plot(freq_hz, log10(sum(spect, 1) + noise), 'b*-');
        
        if J == 1
            legend([p1 p2 p3 p4], {'periodogram','oscillator','noise','sum of oscillator & noise'});
        else
            legend([p1 p2 p3 p4], {sprintf('periodogram of y%d',j),'oscillator','noise','sum of oscillator & noise'});
        end
        xlabel('frequency', 'FontSize', 20);
        xlim([0 fs/2]);
        ylabel('log10 power', 'FontSize', 20);
        set(gca, 'FontSize', 16);
    end
end
