function [barhandle, cfg] = plotsigbar(pcorr,cfg)
% plot sig bars %

nt = length(pcorr);

cfgfields = fieldnames(cfg);
if ~ismember(cfgfields,'plotcolor')
    cfg.plotcolor = [0.0 0.0 0.0];
end
if ~ismember(cfgfields,'linewidth')
    cfg.linewidth = 4;
end
if ~ismember(cfgfields,'pthresh')
    cfg.pthresh = 0.05;
end
if ~ismember(cfgfields,'lineheight')
    ylims = get(gca,'ylim');
    cfg.lineheight = min(ylims) + range(ylims)*0.05;
end
if ~ismember(cfgfields,'xvals')
    cfg.xvals = 1:nt;
end
if ~ismember(cfgfields,'stepsize')
    cfg.stepsize = 1/4;
end
if ~ismember(cfgfields,'verbose')
    cfg.verbose = 'true';
end

if ~ismember(cfgfields,'plotpval')
    cfg.plotpval = true;
end
if ~ismember(cfgfields,'pheight')
    ylims = get(gca,'ylim');
    cfg.pheight = min(ylims) + range(ylims)*0.10;
end
if ~ismember(cfgfields,'plotpexact')
    cfg.plotpexact = true;
end
if ~ismember(cfgfields,'pname')
    cfg.pname = '';
end
if ~ismember(cfgfields,'minp')
    cfg.minp = 1/1000;
end
if ~ismember(cfgfields,'plotline')
    cfg.plotline = false;
end
if ~ismember(cfgfields,'numonly')
    cfg.numonly = false;
end
if ~ismember(cfgfields,'minifz')
    cfg.minifz = 11;
end
if ~ismember(cfgfields,'fz')
    cfg.fz = 16;
end
if ~ismember(cfgfields,'numsigdig')
    cfg.numsigdig = 3;
end

stepsize = median(diff(cfg.xvals))*cfg.stepsize;

b        = bwconncomp(pcorr<cfg.pthresh);

barhandle = [];
if b.NumObjects>0
    if cfg.verbose
        if b.NumObjects==1
            fprintf('\n 1 cluster found.\n');
        else
            fprintf('\n %d clusters found.\n',b.NumObjects);
        end
    end
    for ib = 1:b.NumObjects
        currid = b.PixelIdxList{ib};
        if nnz(currid)>0
            xvec = [cfg.xvals(currid(1))-stepsize cfg.xvals(currid(end))+stepsize];
            barhandle(ib) = plot(xvec,ones(length(xvec),1)*cfg.lineheight,'k-','linewidth',cfg.linewidth,'color',cfg.plotcolor);
            if cfg.plotpval
                plotpval(gca,xvec,[1 1]*cfg.pheight,nanmean(pcorr(1,currid)),cfg.plotcolor,cfg.plotpexact,cfg.pname,cfg.minp,cfg.plotline,cfg.numonly,cfg.minifz,cfg.fz,cfg.numsigdig);
            end
        end
        if cfg.verbose
            fprintf('Cluster %d: p=%g, %1.3f-%1.3f\n',ib,pcorr(currid(1)),cfg.xvals(currid([1 end])));
        end
    end
else
    if cfg.verbose
        fprintf('\nNo clusters found.\nMinimum corrected p-value: p=%g\n',min(pcorr));
    end
end
end