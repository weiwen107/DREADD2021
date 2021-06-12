function ax_h = supersubplot(fig,m,n,p)
%supersubplot- Organize axes across mulitple figures
%
%  Aax_h = supersubplot(fig,m,n,p)
%
%  Returns a set of axes that can be arranged across multiple figures.
%
%  Inputs:
%       fig- figure number of the first figure in series
%       m- the number of rows of plots in each figure
%       n- the number of columns of plots in each figure
%       p- the plot number to use, where the numbers go from left to right,
%       top to bottom, and then continue across figures
%
%  Outputs:
%       ax_h- the axes handle
%  Example:
%       fig = figure;
%       ax = supersubplot(fig,3,3,12);
%       will return a subplot in position 3 on a second figure

if isempty(fig), fig = figure;
end
if ~ishandle(fig), fig = figure;
end
ud = get(fig,'userdata');
if isempty(ud), ud = fig;
end

figure_number = ceil(p/(m*n));
if length(ud)<figure_number
    ud(figure_number)=figure;
else
    figure(ud(figure_number));
end
plotnumber = p-m*n*(figure_number-1);
set(fig,'userdata',ud);
ax_h = subplot(m,n,plotnumber);
