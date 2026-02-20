function herr = plot2BarsWithDots(plotdat,colormat,plot_sig,plot_sig_diff,plot_exact,plot_dots_on_top,xshift)
    if nargin<7
        xshift = 0;
    end
    if nargin<6
        plot_dots_on_top = true;
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
    
    nt = size(plotdat,2);
    nx = size(plotdat,1);
    mubar = mean(plotdat,1);
    sdbar = std(bsxfun(@minus,plotdat,mean(plotdat,2)),[],1)./sqrt(size(plotdat,1));

    xlocations = [1:nt];
    xlocations(1:2:end)  = xlocations(1:2:end)+xshift;
    xlocations(2:2:end)  = xlocations(1:2:end)-xshift;

    for ibar = 1:nt
        bar(xlocations(ibar),mubar(1,ibar),'facecolor',colormat(ibar,:),'EdgeColor', 'none', 'BarWidth', 0.7);
    end

    if plot_dots_on_top
        for ibar = 1:nt
            plot(xlocations(ibar)+(rand(nx,1)-0.5)*0.1,plotdat(:,ibar),'k.','markersize',15,'color',ones(1,3)*0.7);
        end
    end
    herr = ploterr(xlocations(1:nt), mubar, [], sdbar, 'k.', 'abshhxy', 0);
    set(herr(1), 'marker', 'none'); % remove marker
    
    if mean(mean(plotdat))>0
        maxdat = max(plotdat(:))+range(plotdat(:))*0.05;
        maxdat2 = max(plotdat(:))+range(plotdat(:))*0.15;
    else
        maxdat = min(plotdat(:))-range(plotdat(:))*0.05;
        maxdat2 = min(plotdat(:))-range(plotdat(:))*0.15;
    end
    if plot_sig
        pval_vs0 = ttestfast(plotdat);
        for ibar = 1:nt    
            mysigstar(gca, [-0.25 0.25]+xlocations(ibar), maxdat, pval_vs0(ibar),'k',plot_exact); 
        end
    end
    if plot_sig_diff
        for itest = 1:floor(nt/2)
            pval_diff = ttestfast(diff(plotdat(:,[1 2]+2*(itest-1)),[],2));
            mysigstar(gca, xlocations([1 2]+2*(itest-1)), maxdat2, pval_diff,'k',plot_exact); 
        end
    end
end