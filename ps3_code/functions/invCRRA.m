function c = invCRRA(u, gamma)

if gamma == 1
    c = exp(u);
else
    c = ((1-gamma)*u).^(1/(1-gamma));
end
end