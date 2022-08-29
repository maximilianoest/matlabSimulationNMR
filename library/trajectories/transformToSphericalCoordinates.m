function [theta,phi] = transformToSphericalCoordinates( ...
    nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ)
% This function calculates the transformation from cartesian coordiates
% to the angles theta and phi.

hypotuseXY = hypot(nearestNeighboursX,nearestNeighboursY);
theta = pi/2-atan2(nearestNeighboursZ,hypotuseXY);
phi = atan2(nearestNeighboursY,nearestNeighboursX);

end
