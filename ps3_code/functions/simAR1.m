function wSim = simAR1(rho, wBar, vareps, nSim)

muw  = wBar/(1-rho);  % Unconditional mean


wSim = nan(nSim,1);
wSim(1) = muw;



for j=2:nSim
    wSim(j) = wBar + rho * wSim(j-1) + normrnd(0, sqrt(vareps));
end

end