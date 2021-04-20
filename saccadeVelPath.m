function [v, s, t] = saccadeVelPath(amp, dt, toPlot)

if ~exist('toPlot', 'var')
    toPlot = false;
end

[dur,pkv] = saccadeProfile(amp);
[xData, yData] = deal([0 dur/2 dur]',[0 pkv 0]');
ft = fittype( 'poly2' );
[fitresult, ~] = fit( xData, yData, ft , 'Normalize', 'on');

N = floor(dur/(dt*1000));
t = (0:N)*dt;
v = feval(fitresult, t*1000);
if v(end) ~= 0
    v(end+1) = 0;
    t(end+1) = t(end) + dt;
end
%v(v < 1) = 0;
s = cumtrapz(t,v);

if toPlot
    % Plot fit with data.
    figure( 'Name', 'Saccade Profile' );
    h = plot( fitresult, xData, yData );
    xlabel('Time (ms)');
    ylabel('Velocity (deg/s)');
    grid on
    hold on
    plot(t*1000,v,'k*');
end
end