function [isi, ma] = interSaccInterval(sacc, win, type, Fs, toPlot)

if ~exist('toPlot','var')
    toPlot = true;
end

lim = 10000;    % limit equal 10 seconds

if ~exist('win', 'var')
    win = 10;
end

if ~exist('type', 'var')
    type = 'pkn';
end

if ~exist('Fs', 'var')
    Fs = 250;
end

isi = [diff([sacc.(type)])*(1/Fs)*1000]';
isi(isi<0) = [];    % remove negative interval
isi(isi>lim) = [];
if win ~= 0
    ma = tsmovavg(isi,'e',win,1);
else
    ma = [];
end

if toPlot
    figure
    subplot(1,2,1)
    plot(isi,'-*k');
    if win ~= 0
        hold on
        plot(ma,'-r');
    end
    h = subplot(1,2,2);
    histnc(h,isi,20,'k','median');
    
end

end