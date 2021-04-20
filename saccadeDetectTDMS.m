function [sacc, dett] = saccadeDetectTDMS(data, w, ax, th, lim, mTail, ...
    Fs, typ, toSave, name, mdist, toPrint)
% SACCADEDETECTTDMS perform iVT detection on TDMS data
% Usage:
%   [sacc, dett] = saccadeDetectTDMS(data,w,ax,th,lim,mTail,Fs,typ,toSave,name,mdist)
%
% Inputs:
%    data - TDMS raw data
%       w - angular velocity
%      ax - axis to be used
%      th - threshold
%     lim - maximum distance of consecutive saccades
%   mTail - maximum window
%      Fs - Sample Frequency
%     typ - latency type
%  toSave - boolean to save detections
%    name - filename
%   mdist - minimum distance of consecutive saccades
% toPrint - to print

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Aug 2014; Version: v14; Last revision: 2016-07-26
%
% Changelog:
%   v? - integrate estimate velocity instead of diffe
%   v? - select used axis
%   v4 - add peak minimum distance
%   v5 - fix pkDel absolute value & prepare struct
%   v6 - add lat and dett
%   v7 - add sig before each procc and find -vel error in last
%   v8 - change use of mTail
%   v9 - use
%  v10 - add mdist
%  v11 - bugfix save and name
%  v12 - add typ & add pkt
%  v13 - change searchIDS & toPrint
%  v14 - bugfix peak and dep at first and last point

%   TODO
%   what is sacc(scc).over?

%------------------------------- BEGIN CODE -------------------------------
T = 1/Fs;
[N, L] = size(data);
if L > 1
    data = data(:,1);
end
toPlot = false;

if ~exist('mdist','var')
    mdist = 120;
end

if ~exist('toSave', 'var')
    toSave = false;
end

if ~exist('toPrint', 'var')
    toPrint = false;
end

if ~exist('name', 'var')
    name = sprintf('%sDetect',datestr(now,'yyyymmdd'));
else
    if strcmp(name,'')
        name = sprintf('%sDetect',datestr(now,'yyyymmdd'));
    end
end

if ~exist('typ', 'var')
    typ = 'pkn';
end

if toSave
    if ischar(toSave)
        [saveDir,name,~] = fileparts(toSave);
        toSave = true;
    elseif islogical(toSave) && exist('name', 'var')
        [saveDir,name,~] = fileparts(name);
    else
        if ispc
            userDir = winqueryreg('HKEY_CURRENT_USER',...
                ['Software\Microsoft\Windows\CurrentVersion\' ...
                'Explorer\Shell Folders'],'Personal');
        else
            userDir = char(java.lang.System.getProperty('user.home'));
        end
        saveDir = uigetdir(userDir);
    end
    oldDir = pwd;
    cd(saveDir);
end

stra = {'Azimuth','Elevation','Roll','Axis'};

pkDel = 0;
velA = rad2deg(w(:,4));
velE = rad2deg(w(:,3));
velR = rad2deg(w(:,2));
velT = rad2deg(w(:,5));
[pksA,depA,pidxA,didxA] = peakdet(velA, th, 'zero', mdist);
[pksE,depE,pidxE,didxE] = peakdet(velE, th, 'zero', mdist);
[pksR,depR,pidxR,didxR] = peakdet(velR, th, 'zero', mdist);
[pksT,~,pidxT,didxT] = peakdet(velT, th, 'threshold', mdist);
if ~isempty(didxT), fprintf('There is depressions in velT!!!'); end

idr = [];

ids = find(pksA < lim);
pkDel = pkDel + abs(length(pksA)-length(ids));
pksA = pksA(ids); pidxA = pidxA(ids);
ids = find(depA > -lim);
pkDel = pkDel + abs(length(depA)-length(ids));
depA = depA(ids); didxA = didxA(ids);
ids = find(pksE < lim);
pkDel = pkDel + abs(length(pksE)-length(ids));
pksE = pksE(ids); pidxE = pidxE(ids);
ids = find(depE > -lim);
pkDel = pkDel + abs(length(depE)-length(ids));
depE = depE(ids); didxE = didxE(ids);
ids = find(pksR < lim);
pkDel = pkDel + abs(length(pksR)-length(ids));
pksR = pksR(ids); pidxR = pidxR(ids);
ids = find(depR > -lim);
pkDel = pkDel + abs(length(depR)-length(ids));
depR = depR(ids); didxR = didxR(ids);
ids = find(pksT < lim);
pkDel = pkDel + abs(length(pksT)-length(ids));
pksT = pksT(ids); pidxT = pidxT(ids);

if toPrint
    fprintf('There were %d peaks/dep excluded\n',pkDel);
end

switch ax
    case 'a'
        IDS = 1;
    case 'e'
        IDS = 2;
    case 'r'
        IDS = 3;
    case 't'
        IDS = 4;
    otherwise % all
        IDS = 1:4;
        
end

scc = 1;    % saccade count

for jj = IDS
    if jj == 1
        vel = velA;
        pidx = pidxA;
        pks = pksA;
        dep = depA;
        didx = didxA;
        ax = 'a';
        ID = 5;
    elseif jj == 2
        vel = velE;
        pidx = pidxE;
        pks = pksE;
        dep = depE;
        didx = didxE;
        ax = 'e';
        ID = 6;
    elseif jj == 3
        vel = velR;
        pidx = pidxR;
        pks = pksR;
        dep = depR;
        didx = didxR;
        ax = 'r';
        ID = 7;
    else
        vel = velT;
        pidx = pidxT;
        pks = pksT;
        dep = [];
        didx = [];
        ax = 't';
        ID = 8;
    end
    
    pks = pks(pidx > 1);    % remove peak at the first point
    pidx = pidx(pidx > 1);
    pks = pks(pidx < N);    % remove peak at the last point
    pidx = pidx(pidx < N);
    dep = dep(didx > 1);    % remove dep at the first point
    didx = didx(didx > 1);
    dep = dep(didx < N);    % remove dep at the last point
    didx = didx(didx < N);
    
    if toPrint
        fprintf('...Processing %s with %d sacc + and %d sacc -\n', ...
            stra{jj},numel(pks),numel(dep));
    end
    
    %   prepare sctruct
    sacc = struct('pkv',nan,'pkn',nan,'pkt',nan,'ini',nan,'fin',nan,...
        'dur',nan,'lat',nan,'amp',nan,'ampT',nan,'ampA',nan,'ampE',nan,...
        'ampR',nan,'mainc','n','saxes','n','over',false);
    dett = struct('sigp',cell(1),'meanp',cell(1),'varp',cell(1),...
        'xp',cell(1),'nxp',[],'idsp',[],'sigv',cell(1),'meanv',cell(1),...
        'varv',cell(1),'xv',cell(1),'nxv',[],'idsv',[]);
    
    for ii = 1:numel(pks)
        sacc(scc).pkv = pks(ii);
        sacc(scc).pkn = pidx(ii);
        sacc(scc).pkt = data(sacc(scc).pkn);
        
        % find ini
        ids = searchIDS(pidx(ii), 'UP', vel, mTail, th, N);
        sig = vel(ids);
        [nUP, ~, Tx, gam, s2, A0, x] = knowDCUnknowJumpTimeDetector(sig, th);
        dett(scc).sigp = sig;
        dett(scc).meanp = A0;
        dett(scc).varp = s2;
        dett(scc).xp = x;
        dett(scc).nxp = nUP;
        dett(scc).idsp = ids;
        if ~isempty(nUP)
            if toPlot
                plotTxVsGamma(sig,th,Tx,gam,nUP);
                if toSave
                    saveas(gcf, fullfile(saveDir, sprintf('%s%cPkUp%04d',...
                        name, upper(ax), scc)),'fig');
                else
                    pause();
                end
            elseif toSave
                plotTxVsGamma(sig,th,Tx,gam,nUP);
                %saveas(gcf,fullfile(saveDir,sprintf('%sPkUp%04d',name,scc)),'fig');
                print(gcf,'-depsc', fullfile(saveDir, sprintf('%s%cPkUp%04d',...
                    name, upper(ax), scc)));
                close;
            end
            sacc(scc).ini = ids(nUP);
        else
            sacc(scc).ini = [];
        end
        % find fin
        ids = searchIDS(pidx(ii), 'DN', vel, mTail, th, N);
        sig = flip(vel(ids));
        [nDN, ~, Tx, gam, s2, A0, x] = knowDCUnknowJumpTimeDetector(sig,th);
        dett(scc).sigv = sig;
        dett(scc).meanv = A0;
        dett(scc).varv = s2;
        dett(scc).xv = x;
        dett(scc).nxv = nDN;
        dett(scc).idsv = ids;
        if ~isempty(nDN)
            if toPlot
                plotTxVsGamma(sig,th,Tx,gam,nDN);
                if toSave
                    saveas(gcf, fullfile(saveDir, sprintf('%s%cPkDn%04d',...
                        name, upper(ax), scc)),'fig');
                else
                    pause();
                end
            elseif toSave
                plotTxVsGamma(sig,th,Tx,gam,nDN);
                %saveas(gcf,fullfile(saveDir,sprintf('%sPkDn%04d',name,scc)),'fig');
                print(gcf,'-depsc', fullfile(saveDir, sprintf('%s%cPkDn%04d',...
                    name, upper(ax), scc) ) );
                close;
            end
            sacc(scc).fin = ids(end - nDN+1);
            sacc(scc).dur = (sacc(scc).fin - sacc(scc).ini)*T*1000;
        else
            sacc(scc).fin = [];
            sacc(scc).dur = [];
        end
        
        if ~isempty(sacc(scc).ini) && ~isempty(sacc(scc).fin)
            if ~isnan(sacc(scc).ini) && ~isnan(sacc(scc).fin)
                %sacc(scc).amp = sum(diffe([data(sacc(scc).ini:sacc(scc).fin,ID)]));
                sacc(scc).amp = sum(vel(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampT = sum(velT(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampA = sum(velA(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampE = sum(velE(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampR = sum(velR(sacc(scc).ini:sacc(scc).fin))*T;
                maxA = max(abs(velA(sacc(scc).ini:sacc(scc).fin)));
                maxE = max(abs(velE(sacc(scc).ini:sacc(scc).fin)));
                maxR = max(abs(velR(sacc(scc).ini:sacc(scc).fin)));
                [~,I] = max([maxA maxE maxR]);
                if isempty(I), I = 4; end
                switch I
                    case 1
                        sacc(scc).mainc = 'a';
                    case 2
                        sacc(scc).mainc = 'e';
                    case 3
                        sacc(scc).mainc = 'r';
                    otherwise
                        sacc(scc).mainc = 'n';
                end
                if scc > 1
                    scl = 1;
                    while isempty(sacc(scc - scl).fin) && scl < scc -1
                        scl = scl + 1;
                    end
                    if scl < scc
                        sacc(scc).lat = (sacc(scc).(typ) - ...
                            sacc(scc - scl).(typ)) * T * 1000;
                    else
                        sacc(scc).lat = [];
                    end
                end
            end
        end
        %max(data(sacc(scc).ini:sacc(scc).fin,ID)) - ...
        %    min(data(sacc(scc).ini:sacc(scc).fin,ID));
        sacc(scc).saxes = ax;
        sacc(scc).over = false;
        scc = scc +1;
    end
    
    for ii = 1:numel(dep)
        sacc(scc).pkv = dep(ii);
        sacc(scc).pkn = didx(ii);
        sacc(scc).pkt = data(sacc(scc).pkn);
        
        % find ini
        ids = searchIDS(didx(ii), 'UP', -vel, mTail, th, N);
        sig = -vel(ids);
        [nUP, ~, Tx, gam, s2, A0, x] = knowDCUnknowJumpTimeDetector(sig, th);
        dett(scc).sigp = sig;
        dett(scc).meanp = A0;
        dett(scc).varp = s2;
        dett(scc).xp = x;
        dett(scc).nxp = ids(nUP);
        dett(scc).idsp = ids;
        if ~isempty(nUP)
            if toPlot
                plotTxVsGamma(sig,th,Tx,gam,nUP);
                if toSave
                    saveas(gcf,fullfile(saveDir, sprintf('%s%cDpUp%04d',...
                        name, upper(ax), scc)),'fig');
                else
                    pause();
                end
            elseif toSave
                plotTxVsGamma(sig,th,Tx,gam,nUP);
                %saveas(gcf,fullfile(saveDir,sprintf('%sDepUp%04d',name,scc)),'fig');
                print(gcf,'-depsc', sprintf('%s%cDpUp%04d',...
                    name, upper(ax), scc));
                close;
            end
            sacc(scc).ini = ids(nUP);
        end
        % find fin
        ids = searchIDS(didx(ii), 'DN', -vel, mTail, th, N);
        sig = flip(-vel(ids));
        [nDN, ~, Tx, gam, s2, A0, x] = knowDCUnknowJumpTimeDetector(sig,th);
        dett(scc).sigv = sig;
        dett(scc).meanv = A0;
        dett(scc).varv = s2;
        dett(scc).xv = x;
        dett(scc).nxv = ids(length(ids) - nDN+1);
        dett(scc).idsv = ids;
        if ~isempty(nDN)
            if toPlot
                plotTxVsGamma(sig,th,Tx,gam,nDN);
                if toSave
                    saveas(gcf,fullfile(saveDir, sprintf('%s%cDpDn%04d',...
                        name, upper(ax), scc) ),'fig');
                else
                    pause();
                end
            elseif toSave
                plotTxVsGamma(sig,th,Tx,gam,nDN);
                %saveas(gcf,fullfile(saveDir,sprintf('%sDpDn%04d',name,scc)),'fig');
                print(gcf,'-depsc', sprintf('%s%cDpDn%04d',...
                    name, upper(ax), scc) );
                close;
            end
            sacc(scc).fin = ids(end-nDN+1);
            sacc(scc).dur = (sacc(scc).fin - sacc(scc).ini)*T*1000;
        end
        
        if ~isempty(sacc(scc).ini) && ~isempty(sacc(scc).fin)
            if ~isnan(sacc(scc).ini) && ~isnan(sacc(scc).fin)
                %sacc(scc).amp = sum(diffe([data(sacc(scc).ini:sacc(scc).fin,ID)]));
                sacc(scc).amp = sum(vel(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampT = sum(velT(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampA = sum(velA(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampE = sum(velE(sacc(scc).ini:sacc(scc).fin))*T;
                sacc(scc).ampR = sum(velR(sacc(scc).ini:sacc(scc).fin))*T;
                %max(data(sacc(scc).ini:sacc(scc).fin,ID)) - ...
                %    min(data(sacc(scc).ini:sacc(scc).fin,ID));
                maxA = max(abs(velA(sacc(scc).ini:sacc(scc).fin)));
                maxE = max(abs(velE(sacc(scc).ini:sacc(scc).fin)));
                maxR = max(abs(velR(sacc(scc).ini:sacc(scc).fin)));
                [~,I] = max([maxA maxE maxR]);
                if isempty(I), I = 4; end
                switch I
                    case 1
                        sacc(scc).mainc = 'a';
                    case 2
                        sacc(scc).mainc = 'e';
                    case 3
                        sacc(scc).mainc = 'r';
                    otherwise
                        sacc(scc).mainc = 'n';
                end
                if scc > 1
                    scl = 1;
                    while isempty(sacc(scc - scl).fin) && scl < scc -1
                        scl = scl + 1;
                    end
                    if scl < scc
                        sacc(scc).lat = (sacc(scc).(typ) - ...
                            sacc(scc - scl).(typ)) * T * 1000;
                    else
                        sacc(scc).lat = [];
                    end
                end
            end
        end
        sacc(scc).saxes = ax;
        sacc(scc).over = false;
        scc = scc +1;
    end
end

if toSave
    cd(oldDir);
end

end
%-------------------------------- END CODE --------------------------------
function ids = searchIDS(idi, typ, vel, mTail, th, N)
%ids = searchIDS(pidx(ii), 'UP', vel, mTail, th, N);
if strcmp(typ,'UP')
    lc = find(vel(1:idi) < th ,1, 'last');
    lc = lc + 2;
    if isempty(lc)
        lc = idi;
    end
    ids = lc-mTail:lc;
    ids = ids(ids > 0 & ids < N);
    lc = find(vel(ids) < th, 1, 'first');
    if ~isempty(lc)
        ids = ids(lc:end);
    end
elseif strcmp(typ,'DN')
    lc = find(vel(idi:end) < th ,1, 'first');
    lc = lc + idi -3;
    if isempty(lc)
        lc = idi;
    end
    ids = lc:lc+mTail;
    ids = ids(ids > 0 & ids < N);
    lc = find(vel(ids) < th, 1, 'last');
    if ~isempty(lc)
        ids = ids(1:lc);
    end
end
%lc = ids(lc);

idx = find(ids < 2,1,'last');
if ~isempty(idx)
    ids = ids(idx:end);
end
idx = find(ids > N-2,1,'first');
if ~isempty(idx)
    ids = ids(1:idx);
end
end
