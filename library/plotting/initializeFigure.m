function currentFigure = initializeFigure(varargin)
% INITIALIZEFIGURE:
% - KEY                 DEFAULT VALUE
% - posAndSize:         [50 50 900 600]
% - overallFontSize     16
% - titleFontSize       16 (NOT overwritten by "overallFontSize" if given)
% - axisFontSize        16 (NOT overwritten by "overallFontSize" if given)
% - lineWidth           1.5
% - legend              true
% - legendFontSize      16
%
% Return value is the figure handler of the created figure which can be
% used with: set(0,'CurrentFigure',currentFigure);
%
% Syntax examples:
% initializeFigure('posAndSize', [70 70 600 400] ...
%   ,'titleFontSize', 20,'axisFontSize',15,'lineWidth',5)
% plot(sin(0:0.1:2*pi))
% ylabel('testY $\varphi$')
% xlabel('testX $\theta$')
% title('testTitle $\pi$')
%
% initializeFigure('posAndSize', [70 70 600 400] ...
%   ,'titleFontSize', 20,'axisFontSize',15,'lineWidth',5)
% plot(sin(0:0.1:2*pi))
% ylabel(sprintf('testY $\\varphi$: %d',12))
% xlabel('testX $\theta$')
% title('testTitle $\pi$')

% 
values.posAndSize = [50 50 900 600];
values.overallFontSize = 16;
values.titleFontSize = values.overallFontSize;
values.axisFontSize = values.overallFontSize;
values.lineWidth = 1.5;
values.grid = 'minor';
values.legend = true;
values.legendFontSize = values.overallFontSize;
titleFontSizeAlreadyWritten = false;
axisFontSizeAlreadyWritten = false;
legendFontSizeAlreadyWritten = false;

for inputElementNr = 1:2:size(varargin,2)
    key = varargin{inputElementNr};
    value = varargin{inputElementNr+1};
    switch key
        case 'posAndSize'
            values.posAndSize = value;
        case 'overallFontSize'
            values.overallFontSize = value;
            if ~titleFontSizeAlreadyWritten
                values.titleFontSize = value;
            end
            if ~axisFontSizeAlreadyWritten
                values.axisFontSize = value;
            end
            if ~legendFontSizeAlreadyWritten
                values.legendFontSize = value;
            end
        case 'titleFontSize'
            values.titleFontSize = value;
            titleFontSizeAlreadyWritten = true;
        case 'axisFontSize'
            values.axisFontSize = value;
            axisFontSizeAlreadyWritten = true;
        case 'lineWidth'
            values.lineWidth = value;
        case 'grid'
            if sum(strcmp(value,{'on','off','minor'}))
                values.grid = value;
            else
                wrongValueErrorMessage(key,value);
            end
        case 'legend'
            if value == true || value == false
                values.legend = value;
            else
                wrongValueErrorMessage(key,value);
            end
        case 'legendFontSize'
            values.legendFontSize = value;
        otherwise
            wrongKeyErrorMessage(key);
    end 
end


figure('Position',values.posAndSize);
set(gcf,'DefaultLineLineWidth',values.lineWidth)
ax = gca;
set(gca,'FontSize',values.overallFontSize)
set(gca,'TickLabelInterpreter','latex')
ax.FontSize = values.axisFontSize;
ax.Title.FontSize = values.titleFontSize;
ax.Title.Interpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.Interpreter = 'latex';
grid(ax,values.grid);
if values.legend
    lgd = legend();
    set(lgd,'Interpreter','latex')
    lgd.FontSize = values.legendFontSize;
end
currentFigure = gcf;
hold on
end

function wrongKeyErrorMessage(key)
error('initializeFigure:unknownKey',['The key "%s"' ...
    ' is not implemented but could be by you.'] ...
    ,key);
end


function wrongValueErrorMessage(key,value)
error('initializeFigure:wrongValue',['The value: "%s"' ...
    ' for the key: "%s" is wrong.'],value,key);
end

