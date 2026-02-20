function [y] = randdv_sigmoid(dim,dvmu,fixdvmu)
%  RANDDV_SIGMOID  Random decision value generator using sigmoid PDF
%  [y] = randdv_sigmoid(dim,dvmu,[fixdvmu])

if nargin < 3
    fixdvmu = false;
end
if nargin < 2
    error('Missing input argument(s).');
end
if size(dim) ~= 2
    error('dim should be [n(trials) n(samples)].');
end

dvmu = min(max(dvmu,-0.5),+0.5);

if abs(dvmu) == 0.5
    beta = sign(dvmu)*inf;%set slope to Inf - Heavyside step function
else
    beta = fzero(@(b)getdvmu(b)-dvmu,1);%get slope parameter for sig. pdf
end

x = linspace(-1,+1,1000);
p = 1./(1+exp(-beta*x));

if fixdvmu
    y = [];
    while size(y,1) < dim(1)
        ycur = rngits(dim,x,p);
        icur = abs(mean(ycur,2)-dvmu) < 0.001;
        y = cat(1,y,ycur(icur,:));
    end
    y = y(1:dim(1),:);
else
    y = rngits(dim,x,p);
end

end

function [dvmu] = getdvmu(beta)
x = linspace(-1,+1,1000);
dvmu = trapz(x,x./(1+exp(-beta*x)));
end