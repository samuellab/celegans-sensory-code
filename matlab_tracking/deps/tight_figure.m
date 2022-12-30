function h = tight_figure(fignum)
% y = TIGHT_FIGURE(fignum)
%
%   Create a figure with tight axes.
%   https://www.mathworks.com/matlabcentral/answers/352024-programmatically-performing-expand-axes-to-fill-figure
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

if nargin == 0
    h = figure;
else
    h = figure(fignum);
end

AxesH = axes;
drawnow;
InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
