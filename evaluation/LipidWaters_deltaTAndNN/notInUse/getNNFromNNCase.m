function nearestNeighbours = getNNFromNNCase(whichNNCase)
tmp = strsplit(whichNNCase{1},'nearestNeighbours');
nearestNeighbours = str2double(tmp{2});
end
