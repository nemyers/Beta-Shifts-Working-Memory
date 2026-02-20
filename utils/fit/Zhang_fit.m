function [fit,cfg] = Zhang_fit(data,cfg)
%        [fit,cfg] = Zhang_fit(data,cfg)
%
% fit maximum-likelihood mixture model as described 
% by Zhang & Luck, Nature 2008

    if nargin < 2
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
        cfg.minparam   = [0   0];
    end
    switch cfg.fitfunction
        case 'vonmises'
            if ~isfield(cfg,'maxparam')
                cfg.maxparam   = [100 1];
            end
            if ~isfield(cfg,'startparam')
                cfg.startparam = [5 0.5];
            end
        case 'gaussian'
            if ~isfield(cfg,'maxparam')
                cfg.maxparam   = [3 1];
            end
            if ~isfield(cfg,'startparam')
                cfg.startparam = [0.5 0.5];
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
            [Pvec fval] = fmincon(@(P)mle_mixpdfVM(data,[P offs]),startparam,[],[],[],[],minparam,maxparam,[],algopt);
        case 'gaussian'
            [Pvec fval] = fmincon(@(P)mle_mixpdfGauss(data,[P offs]),startparam,[],[],[],[],minparam,maxparam,[],algopt);
    end
    fit = [Pvec fval];
    %von Mises mixture model
    function ll = mle_mixpdfVM(x, param)
        K  = param(1);
        Pm = param(2);
        if length(param)<3, mu =0; else, mu = param(3); end
        p  = circ_vmpdf(x,mu,K)*Pm + 1/(2*pi)*(1-Pm);
        ll = -sum(log(p));
    end
    %Gaussian mixture model (low sigma only!)
    function ll = mle_mixpdfGauss(x, param)
        sigma  = param(1);
        Pm     = param(2);
        if length(param)<3, mu =0; else, mu = param(3); end
        p      = normpdf(x,mu,sigma)*Pm + 1/(2*pi)*(1-Pm);
        ll     = -sum(log(p));
    end
end