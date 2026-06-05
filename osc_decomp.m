function [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp(Y,fs,MAX_OSC,algorithm,obj_func,ss_interval)
    arguments
        Y
        fs (1,1)
        MAX_OSC (1,1) {mustBeInteger(MAX_OSC),mustBePositive(MAX_OSC)} = 5
        algorithm = 'quasi-newton_grad1'
        obj_func = 'whittle'
        ss_interval = 1
    end
%    warning('off','MATLAB:nearlySingularMatrix')
%    warning('off','MATLAB:illConditionedMatrix')
    if size(Y,1) == 1
        [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_uni(Y,fs,MAX_OSC,algorithm,obj_func,ss_interval);
    else
        [osc_param,osc_AIC,osc_mean,osc_cov,osc_phase] = osc_decomp_multi(Y,fs,MAX_OSC,algorithm,obj_func,ss_interval);    
    end
    warning('on','MATLAB:nearlySingularMatrix')
    warning('on','MATLAB:illConditionedMatrix')
end


