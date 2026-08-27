function swF_p = getPartialDerivativeFun(param)
    switch param
        case 'taup'
            swF_p = @elastic.friction.swinversion.swF_taup;
        case 'taur'
            swF_p = @elastic.friction.swinversion.swF_taur;
        case 'D_c'
            swF_p = @elastic.friction.swinversion.swF_D_c;
        case 'tau0'
            swF_p = @elastic.friction.swinversion.swF_tau0;
        otherwise   
            error('Parameter %s not supported.', param);
    end