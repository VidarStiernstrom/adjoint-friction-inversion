function  swF_s = swF_s(s, taur, taup, D_c)
    swF_s = (s<D_c).*(-(taup -taur) ./ D_c);

