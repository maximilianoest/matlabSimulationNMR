function writeIntoPlotWindow(figureHandle,subplotIndices,textToWrite)
if length(subplotIndices) ~= 3
    error("Subplot indices are not given correctly");
end

ax = initializeSubplot(figureHandle,subplotIndices(1),subplotIndices(2) ...
    ,subplotIndices(3));
text(0.05,0.5,textToWrite,'interpreter','latex');
set(ax,'visible','off');
lgd = legend();
lgd.Visible = 'Off';

end