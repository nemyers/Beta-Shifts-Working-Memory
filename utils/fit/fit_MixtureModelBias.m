function [fit,cfg] = fit_MixtureModelBias(x,data,bias,cfg)
%        [fit,cfg] = fit_MixtureModelBias(x,data,bias,cfg)
%
% fit mixture model by minimising MSE

    if nargin < 4
        cfg = struct;
    end
    if ~isfield(cfg,'fitmean')
        cfg.fitmean = false;
    end
    if ~isfield(cfg,'offset')
        if isfield(cfg,'meancenter')
            if cfg.meancenter
                cfg.offset = circ_mean(data);
            else
                cfg.offset = 0;
            end
        else
            cfg.meancenter = 'false';
            cfg.offset = 0;
        end
    elseif ~isfield(cfg,'meancenter')
        if cfg.offset == circ_mean(data)
            cfg.meancenter = true;
        else
            cfg.meancenter = false;
        end
    end
    if ~isfield(cfg,'fitfunction')
        cfg.fitfunction = 'vonmises';
    end
    if ~isfield(cfg,'minparam')
        cfg.minparam   = [0 0 0 0];
    end
    switch cfg.fitfunction
        case 'vonmises'
            if ~isfield(cfg,'maxparam')
                cfg.maxparam   = [5 1 1 1];
            end
            if ~isfield(cfg,'startparam')
                cfg.startparam = [1 0.3 0.6 0.15];
            end
        case 'gaussian'
            if ~isfield(cfg,'maxparam')
                cfg.maxparam   = [3 1 1 1];
            end
            if ~isfield(cfg,'startparam')
                cfg.startparam = [0.5 0.3 0.5 0.5];
            end
    end
    if cfg.fitmean
        offs = [];
    else
        offs = cfg.offset;
    end
    startparam = cfg.startparam;
    minparam   = cfg.minparam;
    maxparam   = cfg.maxparam;
    
    algopt     = optimset('fmincon');
    %alternatives: 'active-set' or 'sqp' or 'interior-point' for fmincon
    algopt     = optimset(algopt,'Algorithm','interior-point','Display','off');
    switch cfg.fitfunction
        case 'vonmises'
            [Pvec fval] = fmincon(@(P)mle_mixpdfVM(x,data,bias,[P offs]),startparam,[],[],[],[],minparam,maxparam,[],algopt);
        case 'gaussian'
            [Pvec fval] = fmincon(@(P)mle_mixpdfGauss(x,data,bias,[P offs]),startparam,[],[],[],[],minparam,maxparam,[],algopt);
    end
    fit = [Pvec fval];
    %von Mises mixture model
    function ll = mle_mixpdfVM(x, y, bias, param)
        K  = param(1);
        K2 = param(2)*K;
        Pm = param(3);
        Pg = param(4);
        if length(param)<5, mu = 0; else, mu = param(5); end
        p  = circ_vmpdf(x,mu,K)*Pm + size(x,2)/(2*pi)*Pg - circ_vmpdf(x,bias,K2)*(1-Pm-Pg);
        %p  = reshape(p,size(x));
        p  = p./abs(sum(p,2));
        %ll = -sum(log(p));
        ll = sum(sum((p-y).^2,1),2);
    end
    %Gaussian mixture model (low sigma only!)
    function ll = mle_mixpdfGauss(x, y, bias, param)
        sigma  = param(1);
        sigma  = param(2)*sigma;
        Pm     = param(3);
        Pg     = param(4);
        if length(param)<5, mu =0; else, mu = param(5); end
        p      = normpdf(x,mu,sigma)*Pm + size(x,2)/(2*pi)*Pg - normpdf(x,bias,sigma2)*(1-Pm-Pg);
        %p  = reshape(p,size(x));
        p = p./abs(sum(p,2));
        %ll     = -sum(log(p));
        ll = sum(sum((p-y).^2,1),2);
    end
end