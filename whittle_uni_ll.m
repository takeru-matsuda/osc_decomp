function [mll,tau2hat,grad] = whittle_uni_ll(param,init_theta,freq,p)    
    param = param(:);
    K = length(param) / 3;
    
    u = param(1:K)';
    v = param(K+1:2*K)';
    w = param(2*K+1:3*K)';
    
    a = (tanh(u) + 1) / 2;
    theta = init_theta(:)' + tanh(v) * pi;
    sigma2 = exp(w);

    freq_col = freq(:);
    p_col = p(:);
    m = length(freq_col);

    ThetaMinus = theta - freq_col;
    ThetaPlus  = theta + freq_col;

    den1 = 1 - 2 * a .* cos(ThetaMinus) + a.^2;
    den2 = 1 - 2 * a .* cos(ThetaPlus)  + a.^2;

    sp_vec = sum((sigma2 / 2) .* (1 ./ den1 + 1 ./ den2), 2) + 1;

    tau2hat = mean(p_col ./ sp_vec);
    logspsum = sum(log(sp_vec));

    ll = -m * (log(tau2hat) + 1) - logspsum;
    mll = -ll;

    if nargout > 2
        weight = 1 ./ sp_vec - p_col ./ (tau2hat * sp_vec.^2);
        
        dsp_da = - sigma2 .* ( (a - cos(ThetaMinus)) ./ (den1.^2) + (a - cos(ThetaPlus)) ./ (den2.^2) );
        da_du = 2 * a .* (1 - a);
        dsp_du = dsp_da .* da_du;
        grad_u = sum(weight .* dsp_du, 1)';
        
        dsp_dtheta = - (a .* sigma2) .* ( sin(ThetaMinus) ./ (den1.^2) + sin(ThetaPlus) ./ (den2.^2) );
        dtheta_dv = pi * (1 - tanh(v).^2);
        dsp_dv = dsp_dtheta .* dtheta_dv;
        grad_v = sum(weight .* dsp_dv, 1)';
        
        dsp_dsigma2 = 0.5 * (1 ./ den1 + 1 ./ den2);
        dsigma2_dw = sigma2;
        dsp_dw = dsp_dsigma2 .* dsigma2_dw;
        grad_w = sum(weight .* dsp_dw, 1)';
        
        grad = [grad_u; grad_v; grad_w];
        grad = -grad;
    end
end