function [colors] = getColorsFromColorConfigurationFile()

colorConfiguration = readConfigurationFile('colors.txt');
configColors = colorConfiguration.colors;
configColors = split(configColors,";");
colors = [];

explicitFTColor = '#6699ff';
explicitFTColor = sscanf(explicitFTColor(2:end),'%2x%2x%2x',[1 3])/255;

for colorNr = 1:length(configColors)
    color = configColors{colorNr};
    colors(end+1,1:3) = sscanf(color(2:end),'%2x%2x%2x',[1 3])/255;
end

end
