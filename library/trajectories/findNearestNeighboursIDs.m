function [nearestNeighbourIndex,nearestNeighbourDistancesPow3] ...
    = findNearestNeighboursIDs(numberOfNearestNeighbours,relativeXPositions...
    ,relativeYPositions,relativeZPositions)
% The indices which are return by this function are ordered with increasing
% distance to the considered atom.

numberOfNearestNeighbours = numberOfNearestNeighbours + 1;

distances = sqrt(relativeXPositions.^2+relativeYPositions.^2 ...
    +relativeZPositions.^2);
inverseMeanDistances = 1./mean(distances,2);
inverseDistances = [];
sumIndex = 1;

for i=1:numberOfNearestNeighbours
  smallestDistance = max(inverseMeanDistances);
  closestNeighbourId = inverseMeanDistances == smallestDistance;
  numberOfClosestElements = sum(closestNeighbourId);
  inverseDistances(sumIndex:sumIndex+numberOfClosestElements-1) = ...
      smallestDistance;
  nearestNeighbourIndex(sumIndex:sumIndex+numberOfClosestElements-1) = ...
      find(inverseMeanDistances == smallestDistance); 
  sumIndex = sumIndex + numberOfClosestElements; 
  inverseMeanDistances(closestNeighbourId) = min(inverseMeanDistances)-1;  
end
nearestNeighbourIndex = nearestNeighbourIndex(1:numberOfNearestNeighbours);
nearestNeighbourIndex(1) = [];

nearestNeighbourDistancesPow3 = distances(nearestNeighbourIndex,:).^3;  

end
