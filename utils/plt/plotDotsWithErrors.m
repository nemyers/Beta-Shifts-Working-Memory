function [linehandle dothandle cfg] = plotDotsWithErrors(y,x,dotcolor,cfg)
%[linehandle dothandle cfg] = plotDotsWithErrors(y,dotcolor,cfg)

nsubs = size(y,1);
ncond = size(y,2);

if ndims(y) > 2
    error('plotDotsWithErrors: Input Data Matrix has too many dimensions!');
end

if nargin < 4
   cfg           = [];
end

if ~isfield(cfg,'linewidth')
    cfg.linewidth = 0.5;
end
if ~isfield(cfg,'dotsize')
    cfg.dotsize = 20;
end
if ~isfield(cfg,'circavg')
    cfg.circavg = false;
end

if nargin < 3 || isempty(dotcolor)
    if ncond <= 12
        dotcolor = (linspecer(ncond,'qualitative'));
    else
        dotcolor = (linspecer(ncond,'sequential'));
    end
end
if size(dotcolor,1)==1
    dotcolor = repmat(dotcolor,[ncond,1]); end
cfg.dotcolor = dotcolor;
if nargin < 2
    x = 1:ncond;
end

% get average and s.e.m.
if ~cfg.circavg
    ymu  = squeeze(nanmean(y,1));
    ysd = squeeze(nanstd(y,[],1))./sqrt(size(y,1));
else
    ymu  = squeeze(circ_mean(y,[],1));
    ysd = squeeze(circ_std(y,[],1))./sqrt(size(y,1));
end

% loop over conditions
h = [];
linehandle  = [];
dothandle   = [];
hold on
for ih = 1:size(ymu,2)
    h = ploterr(x(ih),ymu(1,ih), ...
        {[ymu(1,ih) - ysd(1,ih)].*0 [ymu(1,ih) + ysd(1,ih)].*0}, ...
        {ymu(1,ih) - ysd(1,ih) ymu(1,ih) + ysd(1,ih)}, ...
        '.', 'abshhxy', 0);
    set(h(1),'color', cfg.dotcolor(ih, :), 'markersize', cfg.dotsize);
    set(h(2),'color', cfg.dotcolor(ih, :), 'linewidth', cfg.linewidth);
    set(h(3),'color', cfg.dotcolor(ih, :), 'linewidth', cfg.linewidth);
    dothandle(ih,1) = h(1);
    linehandle(ih,1) = h(2);
    linehandle(ih,2) = h(3);
end
end