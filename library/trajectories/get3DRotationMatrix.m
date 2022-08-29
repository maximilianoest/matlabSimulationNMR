function rotationMatrix = get3DRotationMatrix(angleInRadians,rotationAxis)
%
% creates a rotation matrix that rotates a vector around a line connecting
% the origin an the given point.
%
% example1:
% rotate around a random direction a random amount and then back
% the result should be an Identity matrix
% r = rand(4,1);
% rotationmat3D(r(1),[r(2),r(3),r(4)]) * get3DRotationmatrix(-r(1) ...
% ,[r(2),r(3),r(4)])
%
% Bileschi 2009

lengthOfAxis = norm(rotationAxis);
if (lengthOfAxis < eps)
   error('axis direction must be non-zero vector');
end

rotationMatrix = nan(3);
rotationAxis = rotationAxis/lengthOfAxis;

axisX = rotationAxis(1);
axisY = rotationAxis(2);
axisZ = rotationAxis(3);

axisXSquared = axisX^2;
axisYSquared = axisY^2;
axisZSquared = axisZ^2;

cosine = cos(angleInRadians);
sinus = sin(angleInRadians);

rotationMatrix(1,1) = axisXSquared+(axisYSquared+axisZSquared)*cosine;
rotationMatrix(1,2) = axisX*axisY*(1-cosine)-axisZ*sinus;
rotationMatrix(1,3) = axisX*axisZ*(1-cosine)+axisY*sinus;
rotationMatrix(2,1) = axisX*axisY*(1-cosine)+axisZ*sinus;
rotationMatrix(2,2) = axisYSquared+(axisXSquared+axisZSquared)*cosine;
rotationMatrix(2,3) = axisY*axisZ*(1-cosine)-axisX*sinus;
rotationMatrix(3,1) = axisX*axisZ*(1-cosine)-axisY*sinus;
rotationMatrix(3,2) = axisY*axisZ*(1-cosine)+axisX*sinus;
rotationMatrix(3,3) = axisZSquared+(axisXSquared+axisYSquared)*cosine;