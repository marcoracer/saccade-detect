function [sacc, states] = saccadeDetectHMM(model, w, ax, Fs, typ, time)
% SACCADEDETECTHMM perform iHMM detection on velocity data
% Usage:
%   sacc = saccadeDetectHMM(model, w, ax, Fs, typ)
%
% Inputs:
%   model - dml.hmm model
%       w - angular velocity
%      ax - axis to be used
%      Fs - Sample Frequency
%     typ - latency type
%

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Abr 2016; Version: v4; Last revision: 2016-06-7
%
% Changelog:
%   v2 - add pkt & doc
%   v3 - bugfix sacc(ii).amp & abs(vel) & output states
%   v4 - add fixStates

%------------------------------- BEGIN CODE -------------------------------
cut = 1;       % do not use first points!
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

if ~exist('model', 'var')
    doModel = true;
else
    doModel = false;
end

if ischar(model)
    doModel = true;
end

if doModel
    model = dml.hmm('verbose', true);
    model.nhidden = 2;
    model.transmat = [0.9987 0.0013; 0.0083 0.9917]; %[0.9991 0.0009; 0.0085 0.9915];
    model.mu = [11.3827 610.3481]; %[5.1997 232.8332];
    model.Sigma(1,1,1:2) = [41.2451 3.6538e+06]; %[9.3084 3.1618e+04];
    model.prior = [1 0]';
    model.mixmat = [1 1]';
end

velA = rad2deg(w(cut:end,4));
velE = rad2deg(w(cut:end,3));
velR = rad2deg(w(cut:end,2));
velT = rad2deg(w(cut:end,5));

switch ax
    case 'a'
        x = abs(velA);
    case 'e'
        x = abs(velE);
    case 'r'
        x = abs(velR);
    case 't'
        x = velT;   % always positive!
end

X = reshape(x',[1 size(x,2) size(x,1)]);

%tic
Z = model.test(X);
%toc

states = squeeze(Z);
states = fixStates(states);

dz = diff(states);
up = find(dz > 0) + 1;
dn = find(dz < 0) + 1;

if isempty(up) || isempty(dn)
    warning('saccadeDetectHMM: no saccade detected!');
    return;
end

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
else
    if dn(1) < up(1)
        dn = dn(2:end);
        up = up(1:end-1);
    end
end

for ii = 1:length(up)
    sacc(ii).ini = up(ii);
    sacc(ii).fin = dn(ii) - 1;
    sacc(ii).dur = (sacc(ii).fin - sacc(ii).ini) * T * 1000;
    [Y,Im] = max(x(sacc(ii).ini:sacc(ii).fin));
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
    % fixing cut shifting
    sacc(ii).ini = up(ii) + cut - 1;
    sacc(ii).fin = dn(ii) + cut - 2;
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

% figure, ah(1) = subplot(2,1,1); plot(velT, 'm'); hold on; line([1 length(velT)],[1 1]*th); ylim([0 1500]); ah(2) = subplot(2,1,2); stem(z, 'k'); linkaxes(ah, 'x');
% figure, th=30; ah(1) = subplot(3,1,1); plot(velT, 'm'); hold on; line([1 length(velT)],[1 1]*th); ylim([0 1500]); ah(2) = subplot(3,1,2); stem(z, 'k'); ah(3) = subplot(3,1,3); plot(rad2deg(w(:,5)), 'm'); ylim([0 1500]); linkaxes(ah, 'x');

end
%-------------------------------- END CODE --------------------------------
function states = fixStates(states)
%bkp = states;
if numel( strfind( reshape(states,1,[]), ...
        [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2] ) ) > 1
    dz = diff(states);
    if numel(find(dz == 0)) < numel(find(dz ~= 0))
        states(dz == 0) = 2;
        states(dz ~= 0) = 1;
    else
        states(dz ~= 0) = 2;
        states(dz == 0) = 1;
    end
    states(end) = states(end-1);
elseif numel(find(states==1)) < numel(find(states==2))
    states(states == 1) = 3;
    states(states == 2) = 1;
    states(states == 3) = 2;
end
end