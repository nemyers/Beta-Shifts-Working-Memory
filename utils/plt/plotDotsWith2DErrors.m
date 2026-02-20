function [dothandle linehandle cfg] = plotDotsWith2DErrors(y,dotcolor,cfg)
%[dothandle linehandle cfg] = plotDotsWithErrors(y,dotcolor,cfg)

nsubs = size(y,1);
ncond = size(y,3);

if ndims(y) > 3
    error('plotDotsWithErrors: Input Data Matrix has too many dimensions!');
end

if nargin < 3
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

if nargin < 2 || isempty(dotcolor)
    if ncond <= 12
        dotcolor = (linspecer(ncond,'qualitative'));
    else
        dotcolor = (linspecer(ncond,'sequential'));
    end
elseif size(dotcolor,1) == 1
    dotcolor = repmat(dotcolor,[ncond,1]);
end
cfg.dotcolor = dotcolor;

% get average and s.e.m.
if ~cfg.circavg
    ymu = permute(nanmean(y,1),[2 3 1]);
    ysd = permute(nanstd(y,[],1),[2 3 1])./sqrt(size(y,1));
    %ysd = permute(nanstd(y,[],1),[2 3 1]);
else
    ymu = permute(circ_mean(y,[],1),[2 3 1]);
    ysd = permute(circ_std(y,[],1),[2 3 1])./sqrt(size(y,1));
end

% loop over conditions
h = [];
linehandle  = [];
dothandle   = [];
hold on
for ih = 1:size(ymu,2)
    h = ploterr(ymu(1,ih),ymu(2,ih), ...
        {ymu(1,ih) - ysd(1,ih) ymu(1,ih) + ysd(1,ih)}, ...
        {ymu(2,ih) - ysd(2,ih) ymu(2,ih) + ysd(2,ih)}, ...
        '.', 'abshhxy', 0);
    set(h(1),'color', cfg.dotcolor(ih, :), 'markersize', cfg.dotsize);
    set(h(2),'color', cfg.dotcolor(ih, :), 'linewidth', cfg.linewidth);
    set(h(3),'color', cfg.dotcolor(ih, :), 'linewidth', cfg.linewidth);
    dothandle(ih,1) = h(1);
    linehandle(ih,1) = h(2);
    linehandle(ih,2) = h(3);
end
end