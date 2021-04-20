function sccstates = saccadeStates(data, sacc)

fixationState = 1;
saccadeState = 2;
L = length(data);
sccstates = ones(L,1) * fixationState;
N = length(sacc);
for ii = 1:N
    sccstates(sacc(ii).ini:sacc(ii).fin) = saccadeState;
end

end