function [sacc, states] = saccadeDetectVT(w, th, ax, Fs, typ, time)
% SACCADEDETECTVT perform iVT detection on velocity data
% Usage:
%   sacc = saccadeDetectVT(w, th, ax, Fs, typ);
%
% Inputs:
%       w - angular velocity
%      th - velocity threshold
%      ax - axis to be used
%      Fs - Sample Frequency
%     typ - latency type
%

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Jun 2016; Version: v2; Last revision: 2016-06-09
%
% Changelog:
%   v2 - bugfix locsUP size

%------------------------------- BEGIN CODE -------------------------------
if ~exist('ax', 'var')
    ax = 't';
end
if ~exist('typ', 'var')
    typ = 'pkn';
end
if ~exist('time', 'var')
    doPkt = false;
else
    doPkt = true;
    [~,n] = size(time);
    if n > 1
        time = time(:,1);
    end
end

%   prepare sctruct
sacc = struct('pkv',nan,'pkn',nan,'pkt',nan,'ini',nan,'fin',nan,...
    'dur',nan,'lat',nan,'amp',nan,'ampT',nan,'ampA',nan,'ampE',nan,...
    'ampR',nan,'mainc','n','saxes','n','over',false);

T = 1/Fs;

velA = rad2deg(w(:,4));
velE = rad2deg(w(:,3));
velR = rad2deg(w(:,2));
velT = rad2deg(w(:,5));

switch ax
    case 'a'
        vel = abs(velA);
    case 'e'
        vel = abs(velE);
    case 'r'
        vel = abs(velR);
    case 't'
        vel = velT;   % always positive!
end


L = length(vel);
states = ones( L,1);

[~, locsUP] = crossdet(vel, th);
locsUP = locsUP(locsUP < L);
for ii = 1:length(locsUP)
    sacc(ii).ini = locsUP(ii);
    id = find(vel(locsUP(ii):end) < th, 1, 'first');
    if isempty(id)
        if all(vel(locsUP(ii):end) > th)
            [~, id] = min(vel(locsUP(ii):end));
            id = id + 1;
        end
    end
    if ~isempty(id)
        sacc(ii).fin = locsUP(ii) + id - 2;
        if length(sacc(ii).ini:sacc(ii).fin) >= 1
            states(sacc(ii).ini:sacc(ii).fin) = 2;
        end
        sacc(ii).dur = (sacc(ii).fin - sacc(ii).ini) * T * 1000;
        [Y,Im] = max(vel(sacc(ii).ini:sacc(ii).fin));
        sacc(ii).pkv = Y;
        sacc(ii).pkn = sacc(ii).ini + Im - 1;
        ids = sacc(ii).ini:sacc(ii).fin;
        sacc(ii).ampT = sum(velT(ids))*T;
        sacc(ii).ampA = sum(velA(ids))*T;
        sacc(ii).ampE = sum(velE(ids))*T;
        sacc(ii).ampR = sum(velR(ids))*T;
        switch ax
            case 't'
                sacc(ii).amp = sacc(ii).ampT;
            case 'a'
                sacc(ii).amp = sacc(ii).ampA;
            case 'e'
                sacc(ii).amp = sacc(ii).ampE;
            case 'r'
                sacc(ii).amp = sacc(ii).ampR;
        end
        maxA = max(abs(velA(ids)));
        maxE = max(abs(velE(ids)));
        maxR = max(abs(velR(ids)));
        [~,I] = max([maxA maxE maxR]);
        if isempty(I), I = 4; end
        switch I
            case 1
                sacc(ii).mainc = 'a';
            case 2
                sacc(ii).mainc = 'e';
            case 3
                sacc(ii).mainc = 'r';
            otherwise
                sacc(ii).mainc = 'n';
        end
        
        if doPkt
            sacc(ii).pkt = time(sacc(ii).pkn);
        end
        if ii > 1
            sacc(ii).lat = (sacc(ii).(typ) - sacc(ii-1).(typ))*T*1000;
        else
            sacc(ii).lat = [];
        end
        
        sacc(ii).saxes = ax;
        sacc(ii).over = false;
    end
end

% figure, ah(1) = subplot(2,1,1); plot(velT, 'm'); hold on; line([1 length(velT)],[1 1]*th); ylim([0 1500]); ah(2) = subplot(2,1,2); stem(z, 'k'); linkaxes(ah, 'x');
% figure, th=30; ah(1) = subplot(3,1,1); plot(velT, 'm'); hold on; line([1 length(velT)],[1 1]*th); ylim([0 1500]); ah(2) = subplot(3,1,2); stem(z, 'k'); ah(3) = subplot(3,1,3); plot(rad2deg(w(:,5)), 'm'); ylim([0 1500]); linkaxes(ah, 'x');

end
%-------------------------------- END CODE --------------------------------