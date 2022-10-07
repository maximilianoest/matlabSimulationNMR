function currentSubFigure = initializeSubplot( ...
    parentFigure,lines,columns,position)
%  INITIALIZESUBPLOT(parentFigure,lines,columns,position) takes the parent
%  figure and plots a subplot at the position index with a structure given
%  by lines and columns. Return value is the subfigure (= axes objec). The
%  properties given by the parent figure are used to initialize the subplot



xMinorGrid = parentFigure.CurrentAxes.XMinorGrid;
yMinorGrid = parentFigure.CurrentAxes.YMinorGrid;
zMinorGrid = parentFigure.CurrentAxes.ZMinorGrid;

subplot(lines,columns,position)
set(gcf,'DefaultLineLineWidth',get(parentFigure,'DefaultLineLineWidth'));
ax = gca;
set(gca,'FontSize',get(parentFigure.CurrentAxes,'FontSize'))
set(gca,'TickLabelInterpreter','latex')
ax.FontSize = get(parentFigure.CurrentAxes,'FontSize');
ax.Title.FontSize = get(parentFigure.CurrentAxes.Title,'FontSize');
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';

set(ax,'XMinorGrid',xMinorGrid,'YMinorGrid',yMinorGrid...
    ,'ZMinorGrid',zMinorGrid);

if isprop(parentFigure.CurrentAxes,'Legend')
    lgd = legend();
    set(lgd,'Interpreter','latex')
    lgd.FontSize = parentFigure.CurrentAxes.Legend.FontSize;
end
currentSubFigure = gca;
hold on
end