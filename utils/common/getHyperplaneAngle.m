function theta = getHyperplaneAngle(A,B)
%theta = getHyperplaneAngle(A,B)
%each column vector is one vector defining the plane (i.e., 1 PC)

% Compute the projection
B = B - A*(A'*B);
% Make sure its magnitude is less than 1.
theta = asin(norm(B));
end

