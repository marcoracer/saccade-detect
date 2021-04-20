function [comp hd] = saccadeDiff(sacc1,sacc2, method, absolute, frac, out)
% SACCADEDIFF compares two saccade-parameters vectors
%     comp = saccadeDiff(sacc1,sacc2, method, absolute)

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Sep 2015; Version: v2; Last revision: 2015-10-05
% v2 : add fraction

toDist = false;

if nargin >= 2
    switch nargin
        case 2
            method = 'one-by-one';
            absolute = true;
            frac = false;
        case 3
            absolute = true;
            frac = false;
        case 4
            frac = false;
        case 6
            toDist = true;
    end
else
    error('saccadeDiff: must have two saccade vectors as input!');
end

N = min([length(sacc1) length(sacc2)]);

pkv1 = [sacc1.pkv];
pkv2 = [sacc2.pkv];
pkn1 = [sacc1.pkn];
pkn2 = [sacc2.pkn];
ini1 = [sacc1.ini];
ini2 = [sacc2.ini];
fin1 = [sacc1.fin];
fin2 = [sacc2.fin];
dur1 = [sacc1.dur];
dur2 = [sacc2.dur];
amp1 = [sacc1.amp];
amp2 = [sacc2.amp];

switch method
    case 'one-by-one'
        for ii = 1:N
            if frac
                comp(ii).ini = ini2(ii) / ini1(ii);
                comp(ii).fin = fin2(ii) / fin1(ii);
                comp(ii).dur = dur2(ii) / dur1(ii);
                comp(ii).amp = amp2(ii) / amp1(ii);
                comp(ii).pkn = pkn2(ii) / pkn1(ii);
                comp(ii).pkv = pkv2(ii) / pkv1(ii);
            else
                comp(ii).ini = ini2(ii) - ini1(ii);
                comp(ii).fin = fin2(ii) - fin1(ii);
                comp(ii).dur = dur2(ii) - dur1(ii);
                comp(ii).amp = amp2(ii) - amp1(ii);
                comp(ii).pkn = pkn2(ii) - pkn1(ii);
                comp(ii).pkv = pkv2(ii) - pkv1(ii);
            end
            if absolute
                comp(ii).ini = abs(comp(ii).ini);
                comp(ii).fin = abs(comp(ii).fin);
                comp(ii).dur = abs(comp(ii).dur);
                comp(ii).amp = abs(comp(ii).amp);
                comp(ii).pkn = abs(comp(ii).pkn);
                comp(ii).pkv = abs(comp(ii).pkv);
            end
            if toDist
                h = TDMS2DQ([out(ini1(ii),:); out(ini2(ii),:);....
                    out(fin1(ii),:); out(fin2(ii),:);...
                    out(pkn1(ii),:); out(pkn2(ii),:);], EulerSequence.ZYX, ...
                    DataFormatType.DOUBLE_POSITION_QUATERNION, false);
                dini = h(1)'*h(2);
                dfin = h(3)'*h(4);
                dpkn = h(5)'*h(6);
                comp(ii).dini = norm(rad2deg(dini.rotation_angle));
                comp(ii).dfin = norm(rad2deg(dfin.rotation_angle));
                comp(ii).dpkn = norm(rad2deg(dpkn.rotation_angle));
                hd(ii).dini = dini;
                hd(ii).dfin = dfin;
                hd(ii).dpkn = dpkn;
            end
        end
    case 'nearest-peak'
        for ii = 1:N
            tmp = pkv2-pkv1(ii);
            [~,idx] = min(abs(tmp));
            if frac
                comp(ii).ini = ini2(idx) / ini1(ii);
                comp(ii).fin = fin2(idx) / fin1(ii);
                comp(ii).dur = dur2(idx) / dur1(ii);
                comp(ii).amp = amp2(idx) / amp1(ii);
                comp(ii).pkn = pkn2(idx) / pkn1(ii);
                comp(ii).pkv = pkv2(idx) / pkv1(ii);
            else
                comp(ii).ini = ini2(idx) - ini1(ii);
                comp(ii).fin = fin2(idx) - fin1(ii);
                comp(ii).dur = dur2(idx) - dur1(ii);
                comp(ii).amp = amp2(idx) - amp1(ii);
                comp(ii).pkn = pkn2(idx) - pkn1(ii);
                comp(ii).pkv = pkv2(idx) - pkv1(ii);
            end
            if absolute
                comp(ii).ini = abs(comp(ii).ini);
                comp(ii).fin = abs(comp(ii).fin);
                comp(ii).dur = abs(comp(ii).dur);
                comp(ii).amp = abs(comp(ii).amp);
                comp(ii).pkn = abs(comp(ii).pkn);
                comp(ii).pkv = abs(comp(ii).pkv);
            end
            if toDist
                h = TDMS2DQ([out(ini1(ii),:); out(ini2(idx),:);....
                    out(fin1(ii),:); out(fin2(idx),:);...
                    out(pkn1(ii),:); out(pkn2(idx),:);], EulerSequence.ZYX, ...
                    DataFormatType.DOUBLE_POSITION_QUATERNION, false);
                dini = h(1)'*h(2);
                dfin = h(3)'*h(4);
                dpkn = h(5)'*h(6);
                comp(ii).dini = norm(rad2deg(dini.rotation_angle));
                comp(ii).dfin = norm(rad2deg(dfin.rotation_angle));
                comp(ii).dpkn = norm(rad2deg(dpkn.rotation_angle));
                hd(ii).dini = dini;
                hd(ii).dfin = dfin;
                hd(ii).dpkn = dpkn;
            end
        end
end

end