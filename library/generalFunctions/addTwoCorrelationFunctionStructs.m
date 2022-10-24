function corrFuncSaverStruct = addTwoCorrelationFunctionStructs( ...
    corrFuncSaverStruct,corrFuncCalculatedStruct,thetaNumber,phiNumber)

saverStructFieldNames = fieldnames(corrFuncSaverStruct);
calculationStructFieldNames = fieldnames(corrFuncCalculatedStruct);

differentFieldNames = setdiff(saverStructFieldNames ...
    ,calculationStructFieldNames);
if ~isempty(differentFieldNames)
    error('addTwoStructsWithSameFieldNames:notTheSameFieldNames' ...
        ,'structs do not contain the same field names');
end

fieldNamesCount = length(saverStructFieldNames);

for fieldNr = 1:fieldNamesCount
    fieldName = saverStructFieldNames{fieldNr};
    corrFuncSaverThetaPhi = squeeze(corrFuncSaverStruct ...
        .(fieldName)(thetaNumber,phiNumber,:))';
    corrFuncCalculated = corrFuncCalculatedStruct.(fieldName);
    if length(corrFuncSaverThetaPhi) ...
            ~= length(corrFuncCalculatedStruct.(fieldName))
        error('addTwoStructsWithSameFieldNames:arraysDoNotHaveTheSameLength' ...
            ,'The arrays within the structs do not have the same lengths');
    end
    corrFuncSaverStruct.(fieldName)(thetaNumber,phiNumber,:) ...
        = corrFuncSaverThetaPhi + corrFuncCalculated;
end

end