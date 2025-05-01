

function CP = xyPlot(obj, timeAxisName, xName, yName, CA1, CA2, plotArg_)

% CA1 is color coded by time; 
% CA2 is color coded by velocity;
if nargin < 7
    plotArg_ = {};
end


i = 0;
if ~isempty(CA1)
    i = i+1;
    CP(i) = xyPlot1(obj, timeAxisName, xName, yName, CA1, plotArg_);
end
if ~isempty(CA2)
    i = i+1;
    CP(i) = xyPlot2(obj, timeAxisName, xName, yName, CA2, plotArg_);
end


end



function CP = xyPlot1(obj, timeAxisName, xName, yName, CA1, plotArg_)

% CA1 is color coded by time; 
% CA2 is color coded by velocity;

time = obj.read(timeAxisName);
x = obj.read(xName);
y = obj.read(yName);

CM = colormap('jet');
if isempty(plotArg_)   % not handling plotArg_ correctly.
    CP = patch(CA1, [x(:); nan], [y(:); nan], [time(:); nan], 'EdgeColor','interp');
else
    CP = scatter(CA1, [x(:); nan], [y(:); nan], 0.5, [time(:); nan]);
end
colormap(CA1, CM);
CB = colorbar('peer', CA1);
CB.Location = 'northoutside';
CB.Position(2) = sum(CA1.Position([2,4])); drawnow;
CB.Position(2) = sum(CA1.Position([2,4]))+ 0.02*CA1.Position(4);
CB.Position(4) = 0.02*CA1.Position(4);
CB.Label.String = 'Color coded by time';

plot(CA1, [x(1)], [y(1)], 'o', 'color', CM(1,:),'linewidth', 6);
plot(CA1, [x(end)], [y(end)], '^', 'color', CM(end,:),'linewidth', 2);

CA1.XAxis.Label.String = xName;
CA1.YAxis.Label.String = yName;
% CA1.Title.String = 'Trace';


end


function CP = xyPlot2(obj, timeAxisName, xName, yName, CA2, plotArg_)

% CA1 is color coded by time; 
% CA2 is color coded by velocity;

time = obj.read(timeAxisName);
x = obj.read(xName);
y = obj.read(yName);
v = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2)./(time(2:end)-time(1:end-1));
v(v > 3*std(v)) = 3*std(v);    %%%%%%%%%%%%%%%%%%%%%%%%%%
x = (x(2:end)+x(1:end-1))/2;
y = (y(2:end)+y(1:end-1))/2;

CM = colormap('jet');
CB = colorbar('peer', CA2);
CB.Location = 'northoutside';
CB.Position(2) = sum(CA2.Position([2,4])); drawnow;
CB.Position(2) = sum(CA2.Position([2,4]))+ 0.02*CA2.Position(4);
CB.Position(4) = 0.02*CA2.Position(4);
CB.Label.String = 'Color coded by speed';

if ~isempty(plotArg_)
    CP = scatter(CA2, [x(:); nan], [y(:); nan], 0.5, [v(:); nan]);
else
    CP = patch(CA2, [x(:); nan], [y(:); nan], [v(:); nan], 'EdgeColor','interp');
end
colormap(CA2, CM);
plot(CA2, [x(1)], [y(1)], 'o', 'color', CM(1,:),'linewidth', 6);
plot(CA2, [x(end)], [y(end)], '^', 'color', CM(end,:),'linewidth', 2);

CA2.XAxis.Label.String = xName;
CA2.YAxis.Label.String = yName;
% CA2.Title.String = 'Trace';


end






