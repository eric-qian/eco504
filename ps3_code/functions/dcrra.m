function out = dcrra(c, gamma)


out       = c .^ (-gamma);
out(c<=0) = -inf;


end