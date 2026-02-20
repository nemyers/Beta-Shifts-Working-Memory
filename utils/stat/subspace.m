function theta = subspace(A,B)
%SUBSPACE Angle between subspaces.
%   SUBSPACE(A,B) finds the angle between two
%   subspaces specified by the columns of A and B.
%
%   If the angle is small, the two spaces are nearly linearly
%   dependent.  In a physical experiment described by some
%   observations A, and a second realization of the experiment
%   described by B, SUBSPACE(A,B) gives a measure of the amount
%   of new information afforded by the second experiment not
%   associated with statistical errors of fluctuations.
%
%   Class support for inputs A, B:
%      float: double, single

%   The algorithm used here ensures that both small and 
%   large (i.e. close to pi/2) angles are
%   computed accurately, and it allows subspaces of different
%   dimensions following the definition in [2]. The first issue
%   is crucial.  The second issue is not so important; but
%   since the definition from [2] coinsides with the standard
%   definition when the dimensions are equal, there should be
%   no confusion - and subspaces with different dimensions may
%   arise in problems where the dimension is computed as the
%   numerical rank of some inaccurate matrix.
%
%   MATLAB's subspace Rev. 5.5-5.8 fails to provide correct answers
%   for angles smaller than e-8. 
%   MATLAB's subspace Rev. 5.9-5.10.4.3 fails to provide correct answers
%   for angles close to pi/2, makes unnecessary data copying 
%   and has a for loop that should have been vectorized
%
%   For a more general code, which computes all angles
%   and corresponding canonical vectors in a general scalar product, see
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=55&objectType=file

%   Copyright (c) 2000 Andrew Knyazev knyazev@na-net.ornl.gov
%   License: BSD (free software)  
%   This revision is a fix for MATLAB's 
%   $Revision: 5.10.4.3  $Date: 2007/09/18 02:15:38 $
%
%   References:
%   [1] A. Bjorck & G. Golub, Numerical methods for computing
%       angles between linear subspaces, Math. Comp. 27 (1973),
%       pp. 579-594.
%   [2] P.-A. Wedin, On angles between subspaces of a finite
%       dimensional inner product space, in B. Kagstrom &
%       A. Ruhe (Eds.), Matrix Pencils, Lecture Notes in
%       Mathematics 973, Springer, 1983, pp. 263-285.
%   [3] A. V. Knyazev and M. E. Argentati, Principal Angles between Subspaces
%       in an A-Based Scalar Product: Algorithms and Perturbation Estimates. 
%       SIAM Journal on Scientific Computing, 23 (2002), no. 6, 2009-2041.
%       http://epubs.siam.org:80/sam-bin/dbq/article/37733

if issparse(A) | issparse(B)
   error('The code does not work and is not intended for sparse data.')
end

if size(A,1) ~= size(B,1)
   error('Row dimensions of A and B must be the same.')
end

% Compute orthonormal bases, using SVD in "orth" to avoid problems
% when A and/or B is nearly rank deficient.
% Expensive, but very stable.
A = orth(A);
B = orth(B);

threshold=sqrt(2)/2; % Define the threshold for determining when an angle is small

% Compute the most accurate way, according to [1,3].
s = svd(A'*B);
costheta = min(s); %This is the cosine of the angle, accurate for large angles
if costheta < threshold % Check for small angles 
   theta = acos(min(1,costheta));
else
   if size(A,2)<size(B,2) %ignore angles due to different sizes as in [2,3]
   sintheta = norm(A' - (A'*B)*B');
   else
   sintheta = norm(B' - (B'*A)*A');
   end
   theta = asin(min(1,sintheta)); % recompute using sine
end
