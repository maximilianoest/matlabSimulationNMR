function [trajX,trajY,trajZ] = ...
    calculateRelativePositions(trajX,trajY,trajZ,atomNumber)

trajX = bsxfun(@minus,trajX,trajX(atomNumber,:));
trajY = bsxfun(@minus,trajY,trajY(atomNumber,:));
trajZ = bsxfun(@minus,trajZ,trajZ(atomNumber,:));

end