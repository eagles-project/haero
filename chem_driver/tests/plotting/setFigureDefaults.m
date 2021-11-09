%===============================================================================
%===============================================================================
%
% This function sets the default figure settings for the current session of
% Matlab to some reasonable settings.
%
%===============================================================================
%
% Input:
%   fig (optional): if provided, it will set the default settings for a
%       particular figure. Otherwise, it will set for all new figures in the
%       current session
% 
% Outputs:
%   figSettings: a struct containing the default figure settings that are
%                set by this function
%         color: 7 x 3 array containing the color order for the default
%                colormap
%        rgbcmy: struct containing some good-looking basic colors 
%                (red, green, blue, cyan, magenta, yellow)
%       xxxList: cell arrays containing a list of markers/lines/dashed lines
%                to cycle through
% 
%===============================================================================
%===============================================================================

function [figSettings, color, rgbcmy, mList, lineList, dashList] =...
                                                          setFigureDefaults(fig)

% must close all open figures, since settings will only apply to new
% figures
close all

if nargin == 0
%     set defaults to graphics root (i.e., all figures for current session)
    fig = groot;
elseif nargin > 1
    error('SetFigureDefaults may only be called with 1 or 0 arguments')
end

% get line color order
color = lines(7);

% default line width and marker size
lWidth = 3;
mSize = 10;

% vector of where the figure will be positioned by default--would recommend
% changing this based on your screen setup
% [left, bottom, width, height]
loc = [2000 50 1500 900]; % this is for Mike's 2 monitor setup
% loc = [200, 200, 1500, 900]; % this is for a single monitor
% loc = [50, 50, 1000, 600]; % this is for laptop screen

% font size for legend and axes labels
fSize = 24;

% text interpreter to enable latex syntax on labels and axes
textInterp = 'latex';

set(fig, 'defaultLineLineWidth', lWidth);
set(fig, 'defaultLineMarkerSize', mSize);

set(fig, 'defaultScatterLineWidth', lWidth);

set(fig, 'defaultFigurePosition', loc);
set(fig, 'defaultAxesFontSize', fSize);
set(fig, 'defaultAxesTickLabelInterpreter', textInterp);
set(fig, 'defaultTextInterpreter', textInterp);

set(fig, 'defaultLegendInterpreter', textInterp);
set(fig, 'defaultLegendFontSize', fSize);

set(fig, 'defaultTextFontSize', fSize);
set(fig, 'defaultTextFontSize', fSize);

% scatter(nan, nan, '+', 'markerEdgeColor', 'k');
% [~, icons] = legend('a');
% 
% % marker size in legend
% set(icons(2), 'defaultGroupMarkerSize', mSize);
% set(icons(2), 'defaultGroupLineWidth', lWidth);

figSettings = get(groot, 'default');

rgbcmy         = struct;
rgbcmy.red     = [254 0 0] / 255;
rgbcmy.green   = [0 179 0] / 255;
rgbcmy.blue    = [0 76 232] / 255;
rgbcmy.cyan    = [0 215 255] / 255;
rgbcmy.magenta = [212 38 255] / 255;
rgbcmy.yellow  = [255 167 0] / 255;

mList = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '.'};
lineList = {'-+', '-o', '-*', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p',...
            '-h', '-.'};
dashList = {'--+', '--o', '--*', '--x', '--s', '--d', '--^', '--v', '-->',...
            '--<', '--p', '--h', '--.'};

end

