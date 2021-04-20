function [saccAll, toRem] = selectSaccades(saccAll, type, varargin)
% saccAll = selectSaccades(saccAll, type, varargin)

% Author: Marco Borges, Ph.D. Student, Computer/Biomedical Engineer
% UFMG, PPGEE, Neurodinamica Lab, Brazil
% email address: marcoafborges@gmail.com
% Website: http://www.cpdee.ufmg.br/
% ??? 2016; Version: v4; Last revision: 2016-06-11
%
% Changelog:
%   v2 - bugfix when is empty
%   v3 - add type pkv
%   v4 - add type dur

%------------------------------- BEGIN CODE -------------------------------
%sacc = struct('pkv',nan,'pkn',nan,'ini',nan,'fin',nan,'dur',nan,...
%    'lat',nan,'amp',nan,'ampT',nan,'ampA',nan,'ampE',nan,'ampR',nan,...
%    'mainc','n','saxes','n','over',false);

N = length(saccAll);

switch type
    case 'amp'
        if isnumeric(varargin{1})
            toRem = [];
            num = varargin{1};
            if length(num) ~= 2
                error('selectSaccades: amp parameter must be 2 values!');
            end
            for ii = 1:N
                if saccAll(ii).amp < num(1) || saccAll(ii).amp > num(2)
                    toRem = [toRem ii];
                end
            end
            saccAll(toRem) = [];
        end
    case 'pkv'
        if isnumeric(varargin{1})
            toRem = [];
            num = varargin{1};
            if length(num) ~= 2
                error('selectSaccades: pkv parameter must be 2 values!');
            end
            for ii = 1:N
                if saccAll(ii).pkv < num(1) || saccAll(ii).pkv > num(2)
                    toRem = [toRem ii];
                end
            end
            saccAll(toRem) = [];
        end
    case 'dur'
        if isnumeric(varargin{1})
            toRem = [];
            num = varargin{1};
            if length(num) ~= 2
                error('selectSaccades: dur parameter must be 2 values!');
            end
            for ii = 1:N
                if saccAll(ii).dur < num(1) || saccAll(ii).dur > num(2)
                    toRem = [toRem ii];
                end
            end
            saccAll(toRem) = [];
        end
    case 'blk'
        ids = [];
        for ii = 1:N
            if isempty(saccAll(ii).amp) || isempty(saccAll(ii).ini) || ...
                    isempty(saccAll(ii).fin) || isempty(saccAll(ii).dur)
                ids = [ids, ii];
            elseif isnan(saccAll(ii).amp) || isnan(saccAll(ii).ini) || ...
                    isnan(saccAll(ii).fin) || isnan(saccAll(ii).dur)
                ids = [ids, ii];
            end
        end
        saccAll(unique(ids)) = [];
end

end
%-------------------------------- END CODE --------------------------------