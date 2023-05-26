function [wmFieldStrengthArray,r1WMCollectionWithCounterLimit ...
    ,gmFieldStrengthArray,r1GMCollectionWithCounterLimit] ...
    = getCollectionsFromTable(relaxationTimesTable ...
    ,dataPointsCountLowerLimit,fieldStrengthsToIgnore)

fieldStrengthArray = [];
lowerLimitFieldStrengthsToIgnore = fieldStrengthsToIgnore - 0.0001;
upperLimitFieldStrengthsToIgnore = fieldStrengthsToIgnore + 0.0001;

for lineNr = 1:size(relaxationTimesTable,1)
    fieldStrength = relaxationTimesTable.fieldStrength(lineNr);
    if ~any((lowerLimitFieldStrengthsToIgnore < fieldStrength) ...
            .*(upperLimitFieldStrengthsToIgnore > fieldStrength))
         if ~any(fieldStrengthArray == fieldStrength)
             fieldStrengthArray(end+1) ...
                 = relaxationTimesTable.fieldStrength(lineNr); %#ok<AGROW>
            
         end
    end
end

fieldStrengthArray = sort(fieldStrengthArray);
lowerLimitFieldStrengthsToConsider = fieldStrengthArray - 0.0001;
upperLimitFieldStrengthsToConsider = fieldStrengthArray + 0.0001;

r1WMCollection = cell(1,length(fieldStrengthArray));
r1GMCollection = cell(1,length(fieldStrengthArray));

for lineNr = 1:size(relaxationTimesTable,1)
    fieldStrength = relaxationTimesTable.fieldStrength(lineNr);
    if ~any((lowerLimitFieldStrengthsToConsider < fieldStrength) ...
            .*(upperLimitFieldStrengthsToConsider > fieldStrength))
        continue;
    end
    t1WM = relaxationTimesTable.T1_WM(lineNr);
    if ~isnan(t1WM)
        r1WMCollection{fieldStrengthArray == fieldStrength}(end+1) ...
            = 1/t1WM;
    end   
    
    
    t1GM = relaxationTimesTable.T1_GM(lineNr);
    if ~isnan(t1GM)
        r1GMCollection{fieldStrengthArray == fieldStrength}(end+1) ...
            = 1/t1GM;
    end
    
end

r1WMCollectionWithCounterLimit = {};
wmFieldStrengthArray = [];

r1GMCollectionWithCounterLimit = {};
gmFieldStrengthArray = [];

for collectionElementNr = 1:size(r1WMCollection,2)
    if length(r1WMCollection{collectionElementNr}) ...
            > dataPointsCountLowerLimit
        r1WMCollectionWithCounterLimit{end+1} ...
            = r1WMCollection{collectionElementNr}; %#ok<AGROW>
        wmFieldStrengthArray(end+1) = fieldStrengthArray( ...
            collectionElementNr); %#ok<AGROW>
    elseif dataPointsCountLowerLimit == 1
        fieldStrength = fieldStrengthArray(collectionElementNr);
        fprintf("Only 1 data point was found for a field " ...
            + "strength of %.3f T in WM --> searching for STD in table. \n" ...
            ,fieldStrength);
        indexInTable = find( ...
            relaxationTimesTable.fieldStrength > (fieldStrength - 0.001) ...
            & relaxationTimesTable.fieldStrength < (fieldStrength + 0.001));
        standardDev = relaxationTimesTable.T1_WMSTD(indexInTable);
        if ~isnan(standardDev)
            fprintf("For %s at WM a STD of %.4f was found. \n" ...
                ,relaxationTimesTable.Reference{indexInTable},standardDev);
            r1WMCollectionWithCounterLimit{end+1} ...
                = r1WMCollection{collectionElementNr}; %#ok<AGROW>
            wmFieldStrengthArray(end+1) = fieldStrengthArray( ...
                collectionElementNr); %#ok<AGROW>
        end
    end
end

for collectionElementNr = 1:size(r1GMCollection,2)
    if length(r1GMCollection{collectionElementNr}) ...
            > dataPointsCountLowerLimit
        r1GMCollectionWithCounterLimit{end+1} ...
            = r1GMCollection{collectionElementNr}; %#ok<AGROW>
        gmFieldStrengthArray(end+1) = fieldStrengthArray( ...
            collectionElementNr); %#ok<AGROW>
    elseif dataPointsCountLowerLimit == 1
        fieldStrength = fieldStrengthArray(collectionElementNr);
        fprintf("Only 1 data point was found for field " ...
            + "strength of %.3f T in GM --> searching for STD in table. \n" ...
            ,fieldStrength);
        indexInTable = find( ...
            relaxationTimesTable.fieldStrength > (fieldStrength - 0.001) ...
            & relaxationTimesTable.fieldStrength < (fieldStrength + 0.001));
        standardDev = relaxationTimesTable.T1_GMSTD(indexInTable);
        
        if ~isnan(standardDev)
            fprintf("For %s at GM a STD of %.4f was found. \n" ...
                ,relaxationTimesTable.Reference{indexInTable},standardDev);
            r1GMCollectionWithCounterLimit{end+1} ...
                = r1WMCollection{collectionElementNr}; %#ok<AGROW>
            gmFieldStrengthArray(end+1) = fieldStrengthArray( ...
                collectionElementNr); %#ok<AGROW>
        end
        
    end
end

end