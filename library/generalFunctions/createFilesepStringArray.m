function fileSeparatorStringArray = createFilesepStringArray( ...
    numberOfFileSeparatorsInARow)
fileSeparatorStringArray = string([]);
for fileSepNr = 1:numberOfFileSeparatorsInARow
    fileSeparatorStringArray(end+1) = string(filesep); %#ok<AGROW>
end
end