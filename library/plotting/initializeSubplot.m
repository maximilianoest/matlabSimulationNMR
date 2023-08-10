function currentSubFigure = initializeSubplot( ...
    parentFigure,lines,columns,position)
%  INITIALIZESUBPLOT(parentFigure,lines,columns,position) takes the parent
%  figure and plots a subplot at the position index with a structure given
%  by lines and columns. Return value is the subfigure (= axes objec). The
%  properties given by the parent figure are used to initialize the subplot


% if ~isempty(varargin)
%     if mod(length(varargin),2) ~= 0
%         error("Wrong key value pair");
%     end
%     for elementNr = 1:length(varargin)
%         
%     end
% end

xMinorGrid = parentFigure.CurrentAxes.XMinorGrid;
yMinorGrid = parentFigure.CurrentAxes.YMinorGrid;
zMinorGrid = parentFigure.CurrentAxes.ZMinorGrid;


axisFontSize = get(parentFigure.CurrentAxes,'FontSize');
titleFontSize =  get(parentFigure.CurrentAxes.Title,'FontSize');
if ~isempty(parentFigure.CurrentAxes.Legend)
    legendFontSize = parentFigure.CurrentAxes.Legend.FontSize;
else
    legendFontSize = 12;
end
lineWidth = get(parentFigure,'DefaultLineLineWidth');

subplot(lines,columns,position)
set(gca,'DefaultLineLineWidth',lineWidth);
ax = gca;
set(gca,'FontSize',axisFontSize);
set(gca,'TickLabelInterpreter','latex')
ax.FontSize = axisFontSize;
ax.Title.FontSize = titleFontSize;
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';

set(ax,'XMinorGrid',xMinorGrid,'YMinorGrid',yMinorGrid...
    ,'ZMinorGrid',zMinorGrid);

if isprop(parentFigure.CurrentAxes,'Legend')
    lgd = legend();
    set(lgd,'Interpreter','latex')
    lgd.FontSize = legendFontSize;
end
currentSubFigure = gca;
hold on
end