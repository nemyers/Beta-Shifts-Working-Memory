function y = randwithcorr(x,rho,muy,sigmay)
% randwithcorr - random vector y with a given correlation coefficient against vector x
% usage: y = randwithcorr(x,rho,muy,sigmay)
% 
% arguments: (input)
%  x - vector (or array) of numbers, with any general distribution
%      (must have numel(x) >= 3, and x may not be constant)
%      There is no requirement that x be normally distributed.
%
%  rho - scalar, numeric - defines the correlation coefficient of
%      the new vector y that will be correlated to x. Required:
%        -1 <= rho <= 1
%
%  muy - (optional) scalar, numeric - defines the mean of y
%      Default value = 0
%
%  sigmay - (optional) scalar, numeric - defines the standard deviation of y
%      Default value = 1
%
% arguments: (output)
%  y - vector of array of numbers of the same size and shape as x.
%
%      corr(x(:),y(:)) == rho
%
%      Note that the correlation coefficient here will be exactly as given.

% check that x has more than two elements, AND that it is not constant.
if (numel(x) <= 2) || all(diff(x(:)) == 0)
  error('x must have at least three elements, not all of which can be the same.')
end

% default for sigmay = 1
if (nargin < 4) || isempty(sigmay)
  sigmay = 1;
elseif ~isscalar(sigmay) || (sigmay <= 0)
  error('if provided, sigmay must be scalar, positive numeric')
end

% default for muy = 0
if (nargin < 3) || isempty(muy)
  muy = 0;
elseif ~isscalar(muy)
  error('if provided, muy must be scalar')
end

% rho has no default
if ~isscalar(rho) || (rho <-1) || (rho > 1)
  error('rho must be scalar, -1 <= rho, rho <= 1')
end

xsize = size(x);
x = x(:);

% special cases for -1 and 1
if rho == 1
  y = x;
  return
elseif rho == -1
  y = -x;
  return
elseif rho < 0
  % being mentally lazy here. just worry about a positive correlation.
  rho = -rho;
  x = -x;
end

% to make things simple, shift and scale x to have mean == 0, std == 1.
% this does not require anything about the distribution of x, merely that
% x is not a constant vector.
x = (x - mean(x))/std(x);

% generate a random noise vector that is orthogonal to x,
% but also has a zero sum. The simple solution is to start with
% a general random vecor, then subtract off the unwanted components,
% thus projecting into the null space. Whenever remains will be orthogonal
% to the indicated vectors.
y0 = randn(xsize);
M = orth([ones(xsize),x]); % orth is efficient for 2 columns, even if they are long
dp = y0'*M;
y0 = y0 - M*dp';
y0 = y0/std(y0);
% y0 is a random vector, with...
%   mean zero,
%   unit standard deviation,
%   orthogonal to x.

% now we can do the special case for 0 correlation too
if rho == 0
  y = y0;
  % shift and scale y now to have the desired mean and std
    y = muy + y*sigmay;

    % finally, reshape the variable to be the same shape as x
    y = reshape(y,xsize);
  return
end

% we now have two vectors, x and y, that are orthogonal to each other,
% so they have a correlation of zero. Some linear combination of them
% will have the desired correlation though. And since x and y0 both have
% a unit standard deviation, we can not think of the correlation like a
% slope. So choose beta such that
%
%   y = (1-beta)*y0 + beta*x = y0 + beta*(x-y0)
fun = @(beta) rho - corr(x,y0 + beta*(x-y0));
beta = fzero(fun,[0,1]);

y = y0 + beta*(x-y0);

% shift and scale y now to have the desired mean and std
y = muy + y*sigmay;

% finally, reshape the variable to be the same shape as x
y = reshape(y,xsize);

