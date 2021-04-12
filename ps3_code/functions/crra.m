function out = crra(c, gamma)

if gamma == 1
    out       = log(c);
    out(c<=0) = -inf;
    
else
    out       = c .^ (1-gamma) / (1-gamma);  
    out(c<=0) = -inf;
end

end