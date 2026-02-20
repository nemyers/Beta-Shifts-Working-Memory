function h = plotpval(ax, xpos, ypos, pval, color, pexact,pname,minp,plotline,numonly,minifz,fz,numsigdig)
% h = plotpval(ax, xpos, ypos, pval, color, pexact,pname,minp,plotline,numonly,minifz,fz,numsigdig)
% add pval next to significance bar in a figure
% adapted from mysigstar.m, written by Anne Urai 
% (in turn adapted from sigstar)

if ~exist('minp','var'); minp = 1/1000; end %assume 1/1000 as the lowest possible non-zero p-value
if ~exist('pname','var'); pname = ''; end %add name of test
if ~exist('pexact','var'); pexact = false; end %print exact pvalue instead of stars?
if ~exist('ax', 'var'); ax = gca; end
if ~exist('color', 'var'); color = 'k'; end
if ~exist('plotline', 'var'); plotline = false; end
if ~exist('numonly', 'var'); numonly = false; end
if ~exist('minifz', 'var'); minifz = 11; end
if ~exist('fz', 'var'); fz = 16; end
if ~exist('numsigdig','var'); numsigdig = 3; end;

% draw line
if numel(xpos) > 1 & plotline,
    hold on
    % plot the horizontal line
    p = plot(ax, [xpos(1), xpos(2)], ...
        [ypos ypos], '-', 'LineWidth', 0.5, 'color', color);
    
    % use white background
    txtBg = 'w';
    %txtBg = 'none';
else
    txtBg = 'none';
end

fontweight = 'bold';
if ~pexact
    % get sig stars
if pval < 1e-3
    txt = '***';
elseif pval < 1e-2
    txt = '**';
elseif pval < 0.05
    txt = '*';
elseif ~isnan(pval),
    % this should be smaller
    txt = 'n.s.';
    %txt = '';
    fz = minifz; fontweight = 'normal';
else
    return
end

else
    %get exact p value
    if pval == 0
        nzeros = ceil(log10(1/minp));%estimate leading zeros
        if numonly
            txt = sprintf(['%0.' num2str(nzeros) 'f'],minp);
        else
            txt = sprintf(['p<%0.' num2str(nzeros) 'f'],minp);
        end
        fz  = minifz; fontweight = 'normal';
    elseif ~isnan(pval)
        if numonly
            txt = sprintf(['%1.' sprintf('%d',numsigdig) 'f'],pval);
            %txt = sprintf(['%g'],pval);
        else
            %txt = sprintf('p=%1.3f',pval);
            txt = sprintf(['p=%1.' sprintf('%d',numsigdig) 'f'],pval);
            %txt = sprintf('p=%g',pval);
        end
        fz  = minifz; fontweight = 'normal';
    end
end

if ~isempty(pname), txt = [txt ' (' pname ')']; end
% draw the stars in the bar
h = text(double(mean(xpos)), double(mean(ypos)), txt, ...
    'horizontalalignment', 'center', 'backgroundcolor', ...
    txtBg, 'margin', 1, 'fontsize', fz, 'fontweight', fontweight, 'color', color, 'Parent', ax);
end
