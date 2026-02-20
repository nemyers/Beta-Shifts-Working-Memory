function herr = plot4BarsWithDots(plotdat,colormat,plot_sig,plot_sig_diff,plot_exact)
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
    
    nt = size(plotdat,2);
    nx = size(plotdat,1);
    mubar = mean(plotdat,1);
    sdbar = std(bsxfun(@minus,plotdat,mean(plotdat,2)),[],1)./sqrt(size(plotdat,1));
    %sdbar = std(plotdat,[],1)./sqrt(size(plotdat,1));
    for ibar = 1:nt
        bar(ibar+0.15,mubar(1,ibar),'facecolor',colormat(ibar,:),'EdgeColor', 'none', 'BarWidth', 0.3);
    end
    herr = ploterr([1:nt]+0.15, mubar, [], sdbar, 'k.', 'abshhxy', 0);
    set(herr(1), 'marker', 'none'); % remove marker
    for ibar = 1:nt
        plot(ibar+(rand(nx,1)-0.5)*0.1-0.15,plotdat(:,ibar),'k.','markersize',15,'color',colormat(ibar,:));
    end
    if mean(mean(plotdat))>0
        maxdat = max(plotdat(:))+range(plotdat(:))*0.05;
        maxdat2 = max(plotdat(:))+range(plotdat(:))*0.15;
        maxdat3 = max(plotdat(:))+range(plotdat(:))*0.25;
        maxdat4 = max(plotdat(:))+range(plotdat(:))*0.35;
    else
        maxdat = min(plotdat(:))-range(plotdat(:))*0.05;
        maxdat2 = min(plotdat(:))-range(plotdat(:))*0.15;
        maxdat3 = min(plotdat(:))-range(plotdat(:))*0.25;
        maxdat4 = min(plotdat(:))-range(plotdat(:))*0.35;
    end
    if plot_sig
        pval_vs0 = ttestfast(plotdat);
        for ibar = 1:nt    
            mysigstar(gca, [-0.25 0.25]+ibar, maxdat, pval_vs0(ibar),'k',plot_exact); 
        end
    end
    if plot_sig_diff
        pval_diff = ttestfast(diff(plotdat(:,[1 2]),[],2));
       % pval_diff = [pval_diff ttestfast(diff(plotdat(:,[3 4]),[],2))];
        for ipval = 1:length(pval_diff)
            mysigstar(gca, [1 2]+2*(ipval-1), maxdat2, pval_diff(ipval),'k',plot_exact); end

        pval_diff = ttestfast(diff(plotdat(:,[1 3]),[],2));
%        pval_diff = [pval_diff ttestfast(diff(plotdat(:,[2 4]),[],2))];
        mysigstar(gca, [1 3], maxdat3, pval_diff(1),'k',plot_exact);
 %       mysigstar(gca, [2 4], maxdat4, pval_diff(2),'k',plot_exact); 
    end
end