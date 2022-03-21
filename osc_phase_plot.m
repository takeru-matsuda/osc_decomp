function []=osc_phase_plot(osc_phase,osc_mean,osc_cov,fs,K,conf,nmc)
    arguments
        osc_phase
        osc_mean
        osc_cov
        fs (1,1)
        K (1,1) {mustBeInteger(K),mustBePositive(K)}
        conf (1,1) = 2*normcdf(1)-1 % one sigma
        nmc (1,1) = 10^3;
    end
    T = size(osc_mean,2);
    phase1 = zeros(1,T);
    phase2 = zeros(1,T);
    seeds = randn(2,nmc);
    figure
    for k=1:K
        for t=1:T
            tmp = angle_conf_MC(osc_mean(2*k-1:2*k,t,K),osc_cov(2*k-1:2*k,2*k-1:2*k,t,K),conf,seeds);
            phase1(t) = tmp(1);
            phase2(t) = tmp(2);
        end
        subplot(K,1,k),hold on
        plot_phase_area((1:T)/fs,osc_phase(k,1:T,K),phase1,phase2,[.8,.8,.8]);
        plot_phase_nocross((1:T)/fs,osc_phase(k,1:T,K),'k-',2);
        xlim([1 T]/fs);
        set(gca,'FontSize',12);
        ax = gca;
        set(ax,'YTick',[-pi 0 pi]);
        set(ax,'YTickLabel',{'-3.14', '0', '3.14'});
    end
end
