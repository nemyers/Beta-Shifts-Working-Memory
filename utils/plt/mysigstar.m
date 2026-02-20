function h = mysigstar(ax, xpos, ypos, pval, color, pexact,pname,fz)
% replaces sigstar, which doesnt work anymore in matlab 2014b

if ~exist('fz','var'); fz = 16; end %add name of test
if ~exist('pname','var'); pname = ''; end %add name of test
if ~exist('pexact','var'); pexact = false; end %print exact pvalue instead of stars?
if ~exist('ax', 'var'); ax = gca; end
if ~exist('color', 'var'); color = 'k'; end

if numel(ypos) > 1,
    assert(ypos(1) == ypos(2), 'line wont be straight!');
    ypos = ypos(1);
end

% draw line
hold on;
if numel(xpos) > 1,
    % plot the horizontal line
    p = plot(ax, [xpos(1), xpos(2)], ...
        [ypos ypos], '-', 'LineWidth', 0.5, 'color', color);
    
    % use white background
    txtBg = 'w';
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
    elseif pval < 0.10
        txt = '+';
        fontweight = 'normal';
        fz = fz*14/16;
    elseif ~isnan(pval),
        % this should be smaller
        txt = 'n.s.';
        %txt = '';
        fz = fz*11/16; fontweight = 'normal';
    else
        return
    end
else
    %get exact p value
    if pval == 0
        txt = 'p<.001';
        fz  = fz*11/16; fontweight = 'normal';
    elseif ~isnan(pval)
        txt = sprintf('p=%g',pval);
        fz  = fz*11/16; fontweight = 'normal';
    end
end

if ~isempty(pname), txt = [txt ' (' pname ')']; end
% draw the stars in the bar
h = text(double(mean(xpos)), double(mean(ypos)), txt, ...
    'horizontalalignment', 'center', 'backgroundcolor', ...
    txtBg, 'margin', 1, 'fontsize', fz, 'fontweight', fontweight, 'color', color, 'Parent', ax);
end
