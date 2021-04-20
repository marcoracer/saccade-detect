function sacc = saccadeDetectStates(states, w, ax, Fs, typ, time)
% SACCADEDETECTSTATES perform convertion from states to saccade
% Usage:
%   sacc = saccadeDetectStates(states, w, ax, Fs, typ, time)
%

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% May 2016; Version: v1; Last revision: 2016-05-30
%
% Changelog:

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

velT = rad2deg( w(:, 5) );

dz = diff(states);
up = find(dz > 0) + 1;
dn = find(dz < 0) + 1;

upN = length(up);
dnN = length(dn);
if upN ~= dnN
    warning('saccadeDetectHMM: up and dn should be same lenght!');
    [~,I] = min([upN dnN]);
    if I == 1
        if dn(1) < up(1)
            dn = dn(2:end);
        else
            dn = dn(1:end-1);
        end
    else
        if up(1) < dn(1)
            up = up(1:end-1);
        else
            up = up(2:end);
        end
    end
    
end

for ii = 1:length(up)
    sacc(ii).ini = up(ii);
    sacc(ii).fin = dn(ii) - 1;
    sacc(ii).dur = (sacc(ii).fin - sacc(ii).ini) * T * 1000;
    [Y,Im] = max(velT(sacc(ii).ini:sacc(ii).fin));
    sacc(ii).pkv = Y;
    sacc(ii).pkn = sacc(ii).ini + Im - 1;
    ids = sacc(ii).ini:sacc(ii).fin;
    sacc(ii).amp = sum(velT(ids)) * T;
    sacc(ii).mainc = 'n';
    % fixing cut shifting
    sacc(ii).ini = up(ii);
    sacc(ii).fin = dn(ii) - 1;
    sacc(ii).pkn = sacc(ii).ini + Im - 1;
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
%-------------------------------- END CODE --------------------------------