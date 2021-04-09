function out = aInterp(aVec, aPol, aStar)


[~, Ny] = size(aPol);

out = nan(Ny,1);

% for j = 1:Ny
%     
%     aPol_j = aPol(:, j);
%     
%     jHi = find(aVec >= aStar, 1, 'first');
%     jLo = find(aVec <= aStar, 1, 'last');
%     
%     if jHi == jLo
%         out(j) = aPol_j(jHi);
%     else  % Taylor approximation
%         out(j) = aPol_j(jLo) + (aStar - aVec(jLo)) * (aPol_j(jHi) - aPol_j(jLo)) / (aVec(jHi) - aVec(jLo));
%     end
% end

for j = 1:Ny
   out(j) = interp1(aVec, aPol(:, j), aStar, 'linear', 'extrap');
    
end



end