%% The Parameters of Eye movements
% http://www.liv.ac.uk/~pcknox/teaching/Eymovs/params.htm

%% Main Sequence
% There is a fairly fixed relationship between peak velocity and amplitude
%  D=2.2A+21  (see Fig. 4.3, pg72,   in Carpenter, 1988).

A = 2:.1:20;
D = 2.2*A+21;
figure('units','pixels','Position',[0 0 1024 768]);
plot(A,D,'k');
xlabel('Saccade Amplitude (deg)','FontName','Arial','FontWeight','Bold','FontSize',12);
ylabel('Saccade Duration (ms)','FontName','Arial','FontWeight','Bold','FontSize',12);
axis([0 21 0 75]);
set(gca,'TickDir','Out')
a = findobj(gcf); % get the handles associated with the current figure
allaxes = findall(a,'Type','axes');
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
set(allaxes,'FontName','Arial','FontWeight','Bold','FontSize',12);