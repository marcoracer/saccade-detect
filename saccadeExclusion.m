function [sacc, n] = saccadeExclusion(sacc)
% SACCADEEXCLUSION remove all saccade itens with empty ini, fin or dur
% Usage:
%   [sacc, n] = saccadeExclusion(sacc)

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% Sep 2015; Version: v1; Last revision: 2015-09-30

N = length(sacc);
n = 0; aux = 1;
while aux <= N
    if isempty(sacc(aux).ini) || isempty(sacc(aux).fin) || isempty(sacc(aux).dur)
        sacc(aux) = [];
        n = n +1;
        N = N -1;
    else
        aux = aux +1;
    end
end