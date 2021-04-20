function [time,angPos,angVel,nScc,states] = simSaccades(samples, dt, MSPS, ...
    MA, minT, distT, doLimit, doRange)
%SIMSACCADES generates saccadic angular displacement and angular velocity
% [time,angPos,angVel] = simSaccades(samples, dt, MSPM, MA, minT, distT, doLimit, doRange)

% v2 - add limit angle range
% v3 - add range
% v4 - add minT : minimum Delay Time & MA : max amplitude
% v5 - add distT

if ~exist('dt', 'var')
    dt = 1/500;
end

if ~exist('MSPS', 'var') % mean saccade per second
    MSPS = 6;
end

if ~exist('MA', 'var') % mean saccade amplitude
    MA = 50;
end

if ~exist('minT', 'var')
    minT = 50;  % minimum time in ms
end

if ~exist('distT', 'var')
    distT = 'normal';  % default : normal distribution
end

if ~exist('doLimit', 'var')
    doLimit = true;
end

if ~exist('doRange', 'var')
    doRange = true;
end

time = (0:samples-1)*dt;
angVel = zeros(samples,1);
states = ones(samples,1);
nScc = round(abs(MSPS + 2*randn) * time(end));

tScc = randi(samples,nScc,1);
if doLimit
    sig = rand(nScc,1);
end

for ii = 1:nScc
    switch distT
        case 'normal'
            amp = abs(MA + 30*randn);
        case 'gamma'
            amp = gamrnd(20,1);
    end
    
    [v,~,~] = saccadeVelPath(amp, dt);
    n = length(v);
    npo = round((minT/1000)/dt);
    if npo < 1
        error('simSaccades : npn less than 1');
    end
    if tScc(ii)+n-1 < samples - npo
        ids = tScc(ii) - npo:tScc(ii) + n-1 + npo;
        ids = ids(ids > 0);
        if all(angVel(ids) == 0) % OLD NOTE+1 to add one space
            if doLimit
                if sig(ii) < .5
                    v = -v;
                end
            end
            angVel(tScc(ii):tScc(ii)+n-1) = v;
            states(tScc(ii):tScc(ii)+n-1) = 2;
        else
            nScc = nScc -1;
        end
    else
        nScc = nScc -1;
    end
end

angPos = cumtrapz(time, angVel);
if doRange
    while(min(angPos) < -180 || max(angPos) > 180)
        for ii = 1:length(angPos)
            if angPos(ii) < -180
                angPos(ii) = 360 + angPos(ii);
            elseif angPos(ii) > 180
                angPos(ii) = -360 + angPos(ii);
            end
        end
    end
end
%angPos2 = cumsum(angVel)*dt;

%if any(angPos ~= angPos2)
%    warning('take a look at integration methods!');
%end

end