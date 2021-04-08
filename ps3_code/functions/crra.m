function out = crra(c, gamma)

if gamma == 1
    out = log(c);
else
    out = c .^ (1-gamma) / (1-gamma);
end

end