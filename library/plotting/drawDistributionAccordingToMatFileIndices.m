function drawDistributionAccordingToMatFileIndices(meanPositions ...
    ,groupsToSearchIn,matFileIndices)

dimensionNames = ["X","Y","Z"];
   
fig = initializeFigure();
legendEntries = {};
for groupNr = 1 : length(groupsToSearchIn)
    groupName = groupsToSearchIn(groupNr);
    for dimension = 1 : 3
        initializeSubplot(fig,2,2,dimension);
        histogram(meanPositions(matFileIndices{groupNr},dimension) ...
            ,'BinWidth',0.1);
        title(sprintf("Distribution along %s",dimensionNames(dimension)));
        ylabel("$\#$ H atoms");
        if dimension > 1
            lgd = legend();
            lgd.Visible = 'off';
        end
    end
    legendEntries{end+1} = groupName; %#ok<AGROW>
end
initializeSubplot(fig,2,2,1);
legend(legendEntries);

end