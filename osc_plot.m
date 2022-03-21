function []=osc_plot(osc_mean,osc_cov,fs,K,conf)
    arguments
        osc_mean
        osc_cov
        fs (1,1)
        K (1,1) {mustBeInteger(K),mustBePositive(K)}
        conf (1,1) = 2*normcdf(1)-1 % one sigma
    end
    figure
	T = size(osc_mean,2);
	for k=1:K
    	subplot(K,1,k),hold on
    	std_bar = reshape(sqrt(osc_cov(2*k-1,2*k-1,1:T,K)),1,T);
	    xx = [(1:T)/fs (T:-1:1)/fs]';
    	yy = [osc_mean(2*k-1,1:T,K)-norminv(conf+(1-conf)/2)*std_bar osc_mean(2*k-1,T:-1:1,K)+norminv(conf+(1-conf)/2)*std_bar(T:-1:1)]';
    	fill(xx,yy,[.8,.8,.8],'EdgeColor','none');
    	plot((1:T)/fs,osc_mean(2*k-1,1:T,K),'k-');
	    xlim([1 T]/fs);
    	set(gca,'FontSize',12);
	end
end
