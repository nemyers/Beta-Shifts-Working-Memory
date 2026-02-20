function [stats, cfg] = plot1Dcluster(data,cfg)
% cluster correct and plot significance
cfgfields = fieldnames(cfg);
if ~ismember(cfgfields,'nperm')
    cfg.nperm = 1000;
end
if ~ismember(cfgfields,'minp')
    cfg.minp = 1/cfg.nperm;
end
if ~ismember(cfgfields,'pthresh')
    cfg.pthresh = 0.05;
end
if ~ismember(cfgfields,'pclusterformingthresh')
    cfg.pclusterformingthresh = 0.05;
end
if ~ismember(cfgfields,'ylims')
    cfg.ylims = get(gca,'ylim');
end
if ~ismember(cfgfields,'linespacing')
    cfg.linespacing = 0.05;
end

stats   = struct;
stats.pthresh = cfg.pthresh;
stats.nperm   = cfg.nperm;
stats.sigtime = cfg.xvals;
stats.pvals   = [];
stats.barhandles = [];

nsig = size(data,3);
plotcolors = cfg.plotcolor;
for isig = 1:nsig
    [praw] = ttestfast(data(:,:,isig));
    if cfg.nperm
        pcorr   = ClusterCorrection2(data(:,:,isig),cfg.nperm,cfg.pclusterformingthresh);
    else
        pcorr = praw;
    end
    stats.pvals(isig).praw  = praw;
    stats.pvals(isig).pcorr = pcorr;
    
    cfg.lineheight = min(cfg.ylims) + range(cfg.ylims)*cfg.linespacing*(isig-0.50);
    cfg.pheight    = min(cfg.ylims) + range(cfg.ylims)*cfg.linespacing*(isig+0.25);
    cfg.plotcolor  = plotcolors(isig,:);
    [barhandle, cfg] = plotsigbar(pcorr,cfg);
    stats.barhandles(isig).handle = barhandle;
end