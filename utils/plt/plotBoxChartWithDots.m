function herr = plotBoxChartWithDots(plotdat,xdat,colormat,plot_lines,plot_sig,plot_sig_diff,plot_exact,maxdat,maxdat2)
    if nargin<9
        estim_maxdat2 = true;
    else
        estim_maxdat2 = false;
    end
    if nargin<8
        estim_maxdat = true;
    else
        estim_maxdat = false;
    end
    if nargin<7
        plot_exact=true;
    end
    if nargin<6
        plot_sig_diff=true;
    end
    if nargin<5
        plot_sig=true;
    end
    if nargin<4
        plot_lines = false;
    end
    hold on
    linestyle = '-';
    
    nt = size(plotdat,2);
    np = size(plotdat,1);
    mubar = mean(plotdat,1);
    offsets = (rand(np,nt)-0.5)*0.2;

    if size(colormat,1)==1
        colormat = repmat(colormat,[nt 1]); end
    
    if plot_lines
    if size(plotdat,1)>1
        for isub = 1:size(plotdat,1)
            plot(xdat+offsets(isub,:),plotdat(isub,:),['k' linestyle],'markersize',15,'color',ones(1,3)*0.7);
        end
    end
    end
    
    for ibar = 1:nt
        boxchart(xdat(ibar)*ones(np,1),plotdat(:,ibar),'BoxFaceColor',colormat(ibar,:),'Notch','off','MarkerColor',[0 0 0]);
        plot(xdat(ibar)+[-0.25 +0.25],ones(1,2)*mubar(1,ibar),'k-','linewidth',0.5);
        for isub = 1:np
            hp = scatter(xdat(ibar)+offsets(isub,ibar),plotdat(isub,ibar),25,'markeredgecolor',[0 0 0],'markerfacecolor',colormat(ibar,:),'MarkerFaceAlpha',0.5);
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