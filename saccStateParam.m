function sacc = saccStateParam(sccStates,w,Fs)
% SACCSTATEPARAM stract saccade parameters from sequence of states
%                usually originated from HMM detection
% sacc = saccStateParam(sccStates,w,Fs)

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Sep 2015; Version: v1; Last revision: 2015-09-30

T = 1/Fs;
velA = rad2deg(w(:,4));
velE = rad2deg(w(:,3));
velR = rad2deg(w(:,2));
velT = rad2deg(w(:,5));

[~,ids] = findpeaks(sccStates);

%   prepare sctruct
sacc = struct('pkv',nan,'pkn',nan,'ini',nan,'fin',nan,'dur',nan,...
    'lat',nan,'amp',nan,'ampT',nan,'ampA',nan,'ampE',nan,'ampR',nan,...
    'mainc','n','saxes',[],'over',[]);

for ii = 1:length(ids)
    idd = find(sccStates(ids(ii):end) == 1 ,1);
    idd = ids(ii)+idd-1;
    [pkv, idm] = max(velT(ids(ii):idd));
    sacc(ii).pkv = pkv;
    sacc(ii).pkn = ids(ii)+idm-1;
    sacc(ii).ini = ids(ii);
    sacc(ii).fin = idd;
    sacc(ii).dur = (sacc(ii).fin - sacc(ii).ini)*T*1000;
    if ii > 1
        sacc(ii).lat = (sacc(ii).pkn - sacc(ii-1).pkn)*T*1000;
    else
        sacc(ii).lat = [];
    end
    sacc(ii).amp = sum(velT(sacc(ii).ini:sacc(ii).fin))*T;
    sacc(ii).ampT = sum(velT(sacc(ii).ini:sacc(ii).fin))*T;
    sacc(ii).ampA = sum(velA(sacc(ii).ini:sacc(ii).fin))*T;
    sacc(ii).ampE = sum(velE(sacc(ii).ini:sacc(ii).fin))*T;
    sacc(ii).ampR = sum(velR(sacc(ii).ini:sacc(ii).fin))*T;
    maxA = max(abs(velA(sacc(ii).ini:sacc(ii).fin)));
    maxE = max(abs(velE(sacc(ii).ini:sacc(ii).fin)));
    maxR = max(abs(velR(sacc(ii).ini:sacc(ii).fin)));
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
end

end