function wSim = simTauchen(ResTauchen, nSim)

wIdx    = nan(nSim, 1);               % Simulations: in state units
wIdx(1) = ceil((ResTauchen.in.N)/2);  % Start at middle of grid

Pi      = ResTauchen.Pi;
PiCumul = cumsum(Pi, 2);

for j=2:nSim
    wIdx(j) =  find(PiCumul(wIdx(j-1), :) >= rand(), 1, 'first');    
end

wSim = ResTauchen.wGrid(wIdx);

end