function c = invdcrra(u, gamma)


c = u.^(-1/gamma);


c(c<0) = -inf;

end