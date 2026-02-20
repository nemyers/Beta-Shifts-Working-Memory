  function [efs,F,cdfs,p]=repanova2D(d,D,fn,gg,alpha,print_output)	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [efs,F,cdfs,p,eps,dfs,b,y2,sig]=repanova(d,D,fn,gg,alpha)
%
% General N-way (OxPxQxR...) repeated measures ANOVAs (no nonrepeated factors)
%
% Input:
%
% d = data	A matrix with rows = replications (eg subjects) and
%		              columns = conditions 
%
% D = factors	A vector with as many entries as factors, each entry being
%		the number of levels for that factor
%
%		Data matrix d must have as many columns (conditions) as
%		the product of the elements of the factor matrix D
%
%		First factor rotates slowest; last factor fastest
% 
% 	Eg, in a D=[2 3] design: factor A with 2 levels; factor B with 3:
%	    data matrix d must be organised:
%
%		A1B1	A1B2	A1B3	A2B1	A2B2	A2B3
% 	rep1
%	rep2
%	...
%	
% Output:
%
% efs 	= effect, eg [1 2] = interaction between factor 1 and factor 2
% F   	= F value
% cdfs 	= corrected df's (using Greenhouse-Geisser)
% p     = p-value
% eps   = epsilon
% dfs   = original dfs
% b     = betas
% y2    = cell array of means for each level in a specific ANOVA effect
% sig   = cell array of significant effects (uncorrected and Bonferroni corrected)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
   print_output = true; 
end

if nargin<5
	alpha=0.05;
end

if nargin<4
	gg=1;		% G-G correction
end

if nargin<3		% No naming of factors provided
   for f=1:length(D)
	fn{f}=sprintf('%d',f);
   end
end

Nf = length(D);		% Number of factors
Nd = prod(D);		% Number of conditions
Ne = 2^Nf - 1;		% Number of effects
Nr = size(d,1);		% Number of replications (eg subjects)
nt = size(d,3);     % Number of timepoints
sig=cell(2,1);

if size(d,2) ~= Nd
	error(sprintf('data has %d conditions; design only %d',size(d,2),Nd))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc = cell(Nf,2);	% create main effect/interaction component contrasts
for f = 1 : Nf
	sc{f,1} = ones(D(f),1);
	sc{f,2} = detrend(eye(D(f)),0);
end 

sy = cell(Nf,2);	% create main effect/interaction components for means
for f = 1 : Nf
	sy{f,1} = ones(D(f),1)/D(f);
	sy{f,2} = eye(D(f));
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1 : Ne		% Go through each effect
	cw = num2binvec(e,Nf)+1;
	c  = sc{1,cw(Nf)};	% create full contrasts
	for f = 2 : Nf
		c = kron(c,sc{f,cw(Nf-f+1)});
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	df1 = rank(c);
	df2 = df1*(Nr-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for it = 1:nt
        y         = d(:,:,it) * c;% project data to contrast sub-space        
        b{e,it}   = mean(y);			
        ss(1,it)  =  sum(y*b{e,it}');
        mse(1,it) = (sum(diag(y'*y)) - ss(1,it))/df2;
    end
	mss       =  ss/df1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	efs{e}    = Nf+1-find(cw==2);			% codes which effect 
	F(e,:)    = mss./mse;
	dfs(e,:)  = [df1 df2];
	cdfs(e,:) = dfs(e,:);
	p(e,:)    = 1-spm_Fcdf(F(e,:),cdfs(e,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-function to code all main effects/interactions

function b = num2binvec(d,p)

if nargin<2
	p = 0;		% p is left-padding with zeros option
end

d=abs(round(d));

if(d==0)
	b = 0;
else
	b=[];
 	while d>0
		b=[rem(d,2) b];
		d=floor(d/2);
 	end
end

b=[zeros(1,p-length(b)) b];

