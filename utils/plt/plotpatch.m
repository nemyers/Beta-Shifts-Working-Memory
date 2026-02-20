function [linehandle patchhandle dothandle cfg] = plotpatch(y,timevec,patchcolor,cfg)
%[linehandle patchhandle dothandle cfg] = plotpatch(y,timevec,patchcolor,cfg)

nsubs = size(y,1);
ntime = size(y,2);
ncond = size(y,3);

if ndims(y) > 3
    error('plotpatch: Input Data Matrix has too many dimensions!');
end

if nargin < 4
   cfg           = [];
end

if ~isfield(cfg,'shading')
    cfg.shading = true;
end
if ~isfield(cfg,'errorline')
    cfg.errorline = false;
end
if ~isfield(cfg,'linewidth')
    cfg.linewidth = 4;
end
if ~isfield(cfg,'errorlinewidth')
    cfg.errorlinewidth = 1;
end
if ~isfield(cfg,'alpha')
    cfg.alpha = 0.5;
end
if ~isfield(cfg,'plotdots')
    cfg.plotdots = false;
end
if ~isfield(cfg,'dotsize')
    cfg.dotsize = 12;
end
if ~isfield(cfg,'plotline')
    cfg.plotline = true;
end
if ~isfield(cfg,'linestyle')
    cfg.linestyle = '-';
end
if ~isfield(cfg,'circavg')
    cfg.circavg = false;
end

if nsubs == 1,
    cfg.shading = false;
end

if nargin < 3 || isempty(patchcolor)
    if ncond <= 12
        patchcolor = (linspecer(ncond,'qualitative'));
    else
        patchcolor = (linspecer(ncond,'sequential'));
    end
end
cfg.patchcolor = patchcolor;

if nargin < 2 | isempty(timevec)
   timevec = 1:ntime; 
end

xp  = [timevec fliplr(timevec)];

patchhandle = [];
linehandle  = [];
dothandle   = [];
%hold on
for icondition = 1:ncond
    if ~cfg.circavg
        ymu = squeeze(nanmean(y(:,:,icondition),1));
        ysd = squeeze(nanstd(y(:,:,icondition),[],1))./sqrt(size(y,1));
        %ysd = squeeze(nanstd(bsxfun(@minus,y(:,:,icondition),nanmean(y,3)),[],1))./sqrt(size(y,1));
        alpha = 0.05;
        t_crit = tinv(1-alpha/2,nsubs-1);
        errorbound = t_crit * (ysd);
        %ysd = errorbound;
    else
        ymu = squeeze(circ_mean(y(:,:,icondition),[],1));
        ysd = squeeze(circ_std(y(:,:,icondition),[],1))./sqrt(size(y,1));    
    end
    yp  = [ymu+ysd fliplr(ymu-ysd)];
    cp  = patchcolor(icondition,:);
    hold on
    if cfg.shading
        patchhandle(icondition) = patch(xp,yp,cp,'EdgeColor','none','FaceAlpha',cfg.alpha);
    end
    if cfg.errorline
        patchhandle(icondition,1) = plot(timevec,ymu+ysd,['k' cfg.linestyle],'linewidth',cfg.errorlinewidth,'color',cp); 
        patchhandle(icondition,2) = plot(timevec,ymu-ysd,['k' cfg.linestyle],'linewidth',cfg.errorlinewidth,'color',cp); 
    end
    if cfg.plotline
        linehandle(icondition)  = plot(timevec,ymu,['k' cfg.linestyle],'linewidth',cfg.linewidth,'color',cp);
    end
    if cfg.plotdots
       dothandle(icondition) = plot(timevec,ymu,'ko','Markerfacecolor',cp,'markeredgecolor',[0 0 0],'markersize',cfg.dotsize,'linewidth',1);
    end
end

end