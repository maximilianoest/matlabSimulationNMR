function matFileIndices = getMatFileIndicesFromIndexFile(indexFilePath ...
    ,groupsToSearchIn,unassignedGroup)

if length(string(unassignedGroup)) > 1
   error("Only implementation for single unassigned group implemented."); 
end

groupNamesOfInterest = [groupsToSearchIn, unassignedGroup];

indexFileID = fopen(indexFilePath);
line = fgetl(indexFileID);
writeLines = false;
indexMatrix = cell(length(groupNamesOfInterest),1);
while line ~= -1
    indicesInLine = str2num(line); %#ok<ST2NM>
    if isempty(indicesInLine)
        groupName = erase(line,{'[ ',' ]'});
        if any(groupName == groupNamesOfInterest)
            writeLines = true;
            positionInIndexMatrix = find(groupName == groupNamesOfInterest);
        else
            writeLines = false;
        end
        indicesInLine = str2num(fgetl(indexFileID)); %#ok<ST2NM>
    end
    
    if writeLines
        indexMatrix{positionInIndexMatrix} ...
            = [indexMatrix{positionInIndexMatrix} indicesInLine];
    end
    
    line = fgetl(indexFileID);
end

matFileIndices = cell(length(groupsToSearchIn),1);
unassignedIndices = indexMatrix{groupNamesOfInterest == unassignedGroup};
for indexNr = 1 : length(unassignedIndices)
    index = unassignedIndices(indexNr);
    foundAGroup = false;
    for assignedGroupNr = 1 : length(indexMatrix) - 1
        isIndexInGroup = any(indexMatrix{assignedGroupNr} == index);
        if isIndexInGroup
            if foundAGroup
                error("Index assigned to more than one group.");
            end
            matFileIndices{assignedGroupNr}(end+1) = indexNr;
            foundAGroup = true;
        end
    end
end
end