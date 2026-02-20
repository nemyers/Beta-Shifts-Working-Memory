function herr = plotBarsWithDots(plotdat,colormat,plot_sig,plot_sig_diff,plot_exact,maxdat,maxdat2,showBaseline)
    if nargin<8
        showBaseline = true;
    end
    if nargin<7
        estim_maxdat2 = true;
    else
        estim_maxdat2 = false;
    end
    if nargin<6
        estim_maxdat = true;
    else
        estim_maxdat = false;
    end
    if nargin<5
        plot_exact=true;
    end
    if nargin<4
        plot_sig_diff=true;
    end
    if nargin<3
        plot_sig=true;
    end
    hold on
    linestyle = '-';
    
    nt = size(plotdat,2);
    mubar = mean(plotdat,1);
    sdbar = std(bsxfun(@minus,plotdat,mean(plotdat,2)),[],1)./sqrt(size(plotdat,1));
    for ibar = 1:nt
        hbar(ibar) = bar(ibar,mubar(1,ibar),'facecolor',colormat(ibar,:),'EdgeColor', 'none', 'BarWidth', 0.7);
        if ~showBaseline, hbar(ibar).ShowBaseLine='off'; end
    end
    herr = ploterr(1:nt, mubar, [], sdbar, 'k.', 'abshhxy', 0);
    set(herr(1), 'marker', 'none'); % remove marker
    
    if size(plotdat,1)>1
        for isub = 1:size(plotdat,1)
            plot([1.2 1.8],plotdat(isub,:),['k' linestyle],'markersize',15,'color',ones(1,3)*0.7);
        end
    end
    if mean(mean(plotdat))>0
        if estim_maxdat,  maxdat = max(plotdat(:))+range(plotdat(:))*0.05;  end
        if estim_maxdat2, maxdat2 = max(plotdat(:))+range(plotdat(:))*0.15; end
    else
        if estim_maxdat,  maxdat = min(plotdat(:))-range(plotdat(:))*0.05;  end
        if estim_maxdat2, maxdat2 = min(plotdat(:))-range(plotdat(:))*0.15; end
    end
    if plot_sig
        pval = ttestfast(plotdat);
        for ibar = 1:nt    
            mysigstar(gca, [-0.25 0.25]+ibar, maxdat, pval(ibar),'k',plot_exact); 
        end;
    end
    if plot_sig_diff
        pval = ttestfast(diff(plotdat,[],2));
        mysigstar(gca, [1 2], maxdat2, pval,'k',plot_exact); 
    end
end