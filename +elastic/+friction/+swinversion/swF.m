function  swF = swF(s, taur, taup, D_c)
    if s < D_c
        swF = taup - (taup -taur) * s / D_c;
    else
        swF = taur;
    end
