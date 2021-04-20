function [dur,pkv] = saccadeProfile(amp)
% SACCADEPROFILE returns the duration for a saccade amplitude
% [dur,pkv] = saccadeProfile(amp)

dur = 2.2*amp + 21;
pkv = amp*1.72/(dur/1000);

end

%vertex = -b/(2*a);