function [x,iter] = randexp(min,avg,max,siz)

if nargin < 4
    siz = [1,1];
end
if nargin < 3
    error('Wrong input argument(s).');
end


x = min-log(rand(siz))*(avg-min);
iter = 0;

while true
    iter = iter+1;
    isup = find(x(:) > max);
    nsup = numel(isup);
    if nsup == 0
        break
    end
    x(isup) = min-log(rand(nsup,1))*(avg-min);
end

end