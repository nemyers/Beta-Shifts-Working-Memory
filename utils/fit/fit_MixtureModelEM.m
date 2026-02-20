function [B, LL, RW, iter, outvals] = fit_MixtureModelEM (X, T, NT, B_start)
% --> www.paulbays.com, modified N. Myers 2023

if (nargin<2 || size(X,2)>1 || size(T,2)>1 || size(X,1)~=size(T,1) || nargin>2 && ~isempty(NT) && (size(NT,1)~=size(X,1) || size(NT,1)~=size(T,1)))     
    error('Input is not correctly dimensioned');
    return; 
end

if (nargin>3 && (B_start(1)<0 || any(B_start(2:4)<0) || any(B_start(2:4)>1) || abs(sum(B_start(2:4))-1) > 10^-6))
    error('Invalid model parameters');
    return;
end

MaxIter = 10^4; MaxdLL = 10^-4;

n = size(X,1); 

if (nargin<3) 
    NT = zeros(n,0); nn = 0;
else
    nn = size(NT,2);
end

% Default starting parameters
if (nargin<4)    
    K = 5; Pt = 0.5; K2= 5;
    if (nn>0) Pn = 0.3; else Pn = 0; end
    Pu = 1-Pt-Pn;
else
    K = B_start(1); 
    Pt = B_start(2); Pn = B_start(3); Pu = B_start(4);
    K2 = B_start(5);
end

E  = X-T; E = mod(E + pi, 2*pi) - pi;
NE = repmat(X,1,nn)-NT; NE = mod(NE + pi, 2*pi) - pi;

LL = nan; dLL = nan; iter = 0;
outvals = nan(MaxIter,4);
while (1)
    iter = iter + 1;
    
    Wt = Pt * vonmisespdf(E,0,K);
    Wg = Pu * ones(n,1)/(2*pi);

    if nn==0
        Wn = zeros(size(NE));
    else
        Wn = Pn/nn * vonmisespdf(NE,0,K2);
    end
    
    W = sum([Wt Wn Wg],2);%sum of pdfs for each trial
    
    dLL = LL-sum(log(W));
    LL = sum(log(W));
    if (abs(dLL) < MaxdLL | iter > MaxIter), break; end
    
    Pt = sum(Wt./W)/n;%rel weght of target across experiment
    Pn = sum(sum(Wn,2)./W)/n; 
    Pu = sum(Wg./W)/n;
            
    outvals(iter,:) = [Pt Pn Pu LL];
    
    %estimate kappa for targets
    rw = [(Wt./W) ]; %rel weight of target per trial
    S = [sin(E) sin(NE)]; C = [cos(E) cos(NE)];
    r = [sum(sum(S.*rw)) sum(sum(C.*rw))]; %weighted estimation of kappa
    
    if sum(sum(rw))==0%Uniform distr
        K = 0;
    else
        R = sqrt(sum(r.^2))/sum(sum(rw));        
        K = A1inv(R);
    end
    
    if nn==0
        K2 = 0;
    else
        %estimate kappa2 for non-targets
        rw = [(Wn./repmat(W,1,nn))]; %rel weight of non-target per trial
        S = [sin(E) sin(NE)]; C = [cos(E) cos(NE)];
        r = [sum(sum(S.*rw)) sum(sum(C.*rw))]; %weighted estimation of kappa

        if sum(sum(rw))==0%Uniform distr
            K2 = 0;
        else
            R = sqrt(sum(r.^2))/sum(sum(rw));        
            K2 = A1inv(R);
        end
    end
    
    if n<=15
        if K2<2
            K2 = max(K2-2/(n*K2), 0);
        else
            K2 = K2 * (n-1)^3/(n^3+n);
        end
        
        if nn>0
            if K2<2
                K2 = max(K2-2/(n*K2), 0);
            else
                K2 = K2 * (n-1)^3/(n^3+n);
            end
        end
    end
end

RW = [(Wt./W) (Wn./repmat(W,1,nn))];

if iter>MaxIter
    warning('JV10_function:MaxIter','Maximum iteration limit exceeded.');
    B = [NaN NaN NaN NaN NaN]; LL = NaN; RW = NaN;
else  
    B = [K Pt Pn Pu K2];
end

%%

function K = A1inv(R)

if (0 <= R & R < 0.53)
    K = 2 * R + R^3 + (5 * R^5)/6;
elseif (R < 0.85)
    K = -0.4 + 1.39 * R + 0.43/(1 - R);
else
    K = 1/(R^3 - 4 * R^2 + 3 * R);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2010 Paul Bays. This program is free software: you can     %
%   redistribute it and/or modify it under the terms of the GNU General  %
%   Public License as published by the Free Software Foundation.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%