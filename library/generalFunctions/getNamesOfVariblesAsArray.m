function variableNamesArray = getNamesOfVariblesAsArray(varargin)
variableNamesArray = string.empty();

for index = 1:length(varargin)
    variableNamesArray(end+1) = inputname(index); %#ok<AGROW>
end

end