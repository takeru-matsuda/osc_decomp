function osc_spectrum_plot(Y,fs,osc_a,osc_f,osc_sigma2,osc_tau2,osc_c)
    K = length(osc_a);
    J = size(Y,1);
    if J == 1
        osc_c = zeros(1,2*K);
        osc_c(1:2:end) = 1;
    end
    T = size(Y,2);
    H = osc_observation_matrix(osc_c,J,K);

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
            S = osc_component_spectrum(a,theta,H(j,2*k-1:2*k),freq_rad);
            spect(k,:) = osc_sigma2(k)*real(reshape(S,1,[]));
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

