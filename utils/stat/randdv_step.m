function [y] = randdv_step(dim,dvmu,fixdvmu)
%  RANDDV_STEP  Random decision value generator using step PDF
%  [y] = randdv_step(dim,dvmu,[fixdvmu])

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

x = -1:0.001:+1;
p = 0.5*ones(size(x))+dvmu*sign(x);

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