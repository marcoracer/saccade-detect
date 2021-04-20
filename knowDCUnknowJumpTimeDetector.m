function [n0, xn, Tx, gam, s2, A0, x] = knowDCUnknowJumpTimeDetector(sig,th,a0,pFA,doPlot)
%[n0, xn, Tx, gam] = knowDCUnknowJumpTimeDetector(sig,th,a0,pFA)

%   v1 - add var(1:ii)
%   v2 - return s2, x and doPlot
%   v3 - bugfix ids
%   v4 - change doAverage 20160711

[m,n] = size(sig);

if ~exist('pFA', 'var')
    pFA = 1E-3;
end

if ~exist('a0', 'var')
    a0 = mean(sig(sig < th));
    s2_ = var(sig(sig < th));
    doAverage = false;
elseif isnan(a0)
    doAverage = true;
end

if ~exist('doPlot', 'var')
    doPlot = false;
end

Tx = zeros(m-1,1);
gam = zeros(m-1,1);
PFA = Qinv(pFA);
s2 = cell(1,m-1);
x = cell(1,m-1);
A0 = cell(1,m-1);

for ii = 2:m
    if doAverage
        a0 = mean(sig(1:ii));
        s2_ = var(sig(1:ii));
    end
    A0{ii-1} = a0;
    sig_ = sig(ii:m);
    x{ii-1} = sig_;
    s2{ii-1} = s2_;
    Tx(ii-1) = sum(sig_ -a0 -th/2);
    gam(ii-1) = sqrt( s2_ / (m-ii) ) * PFA;
end
ids = find(Tx > gam);
[~,n0] = max(Tx(ids));
n0 = ids(n0);
xn = sig(n0);

if doPlot
    plotTxVsGamma(sig,th,Tx,gam,n0);
end

end

function y = Qinv(x)
%  This program computes the inverse Q function or the value
%  which is exceeded by a N(0,1) random variable with a
%  probability of x.
y = sqrt(2) * erfinv(1 - 2*x);
end