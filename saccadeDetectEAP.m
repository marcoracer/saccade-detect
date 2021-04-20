function [sacc, dett] = saccadeDetectEAP(vel, Fs, th, mdist, mTail, typ)
%saccadeDetectEAP perform iVTD detection on horizontal and vertical
% projection data from TDMS
%
% Inputs:
%    out - cell with EMA saccade/fixation data
%     Fs - Sampling frequency
%     th - threshold
%  mdist - minimal distance between saccades (default 100 ms)
%  mTail - minimum window for search of the onset/offset
%    typ - latency type
%
% Outputs:
% . sacc - 
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: saccadeDetectHMM, saccadeDetectTDMS;

% Author: Marco Borges, Data Scientist/Biomedical Engineer
% Independent Researcher
% email address: marcoafborges@gmail.com
% Website: https://github.com/marcoracer
% Aug 2018; Version: v1; Last revision: 2018-08-31
% Changelog:

%------------------------------- BEGIN CODE -------------------------------

if ~exist('mdist','var')
    mdist = 10;
end

if ~exist('mTail','var')
    mTail = 100;
end

if ~exist('Fs','var')
    Fs = 1000;
end

if ~exist('typ', 'var')
    typ = 'pkn';
end

dt = 1/Fs;

%   prepare sctruct
sacc = struct('pkv',nan,'pkn',nan,'ini',nan,'fin',nan,'dur',nan,...
                'lat',nan,'amp',nan);
dett = struct('sigp',cell(1),'meanp',cell(1),'varp',cell(1),...
        'xp',cell(1),'nxp',[],'idsp',[],'sigv',cell(1),'meanv',cell(1),...
        'varv',cell(1),'xv',cell(1),'nxv',[],'idsv',[]);

scc = 1;    % saccade count

[pks, dep, pidx, didx] = peakdet(vel, th, 'th', mdist);

N = length(vel);

pks = pks(pidx > 1);    % remove peak at the first point
pidx = pidx(pidx > 1);
pks = pks(pidx < N);    % remove peak at the last point
pidx = pidx(pidx < N);
dep = dep(didx > 1);    % remove dep at the first point
didx = didx(didx > 1);
dep = dep(didx < N);    % remove dep at the last point
didx = didx(didx < N);

for ii = 1:numel(pks)
    sacc(scc).pkv = pks(ii);
    sacc(scc).pkn = pidx(ii);

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
        sacc(scc).fin = ids(end - nDN+1);
        sacc(scc).dur = (sacc(scc).fin - sacc(scc).ini) * dt * 1000;
    else
        sacc(scc).fin = [];
        sacc(scc).dur = [];
    end

    % gathering saccade related parameters
    if ~isempty(sacc(scc).ini) && ~isempty(sacc(scc).fin)
        if ~isnan(sacc(scc).ini) && ~isnan(sacc(scc).fin)
            sacc(scc).amp = sum(vel(sacc(scc).ini:sacc(scc).fin)) * dt;
            if scc > 1
                scl = 1;
                while isempty(sacc(scc - scl).fin) && scl < scc -1
                    scl = scl + 1;
                end
                if scl < scc
                    sacc(scc).lat = (sacc(scc).(typ) - ...
                        sacc(scc - scl).(typ)) * dt * 1000;
                else
                    sacc(scc).lat = [];
                end
            end
        end
    end 
    scc = scc +1;
end
end
%-------------------------------- END CODE --------------------------------
function ids = searchIDS(idi, typ, vel, mTail, th, N)
    %ids = searchIDS(pidx(ii), 'UP', vel, mTail, th, N);
    if strcmp(typ, 'UP')
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
    elseif strcmp(typ, 'DN')
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

    idx = find(ids < 2,1, 'last');
    if ~isempty(idx)
        ids = ids(idx:end);
    end
    idx = find(ids > N-2,1, 'first');
    if ~isempty(idx)
        ids = ids(1:idx);
    end
end