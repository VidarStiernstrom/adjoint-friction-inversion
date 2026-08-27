function  swF_tau0 = swF_tau0(s, taur, taup, D_c, tau_str_L)
    % Nov 10 2025: add tau_str_L
    default_arg('tau_str_L', NaN);

    swF_tau0 = - ones(size(s));
    
    if ~isnan(tau_str_L)
        taustr = (((s<D_c)&(s>0)).*(taup-(taup-taur).*s./D_c) + (s>=D_c).*taur + (s<=0).*taup);
        swF_tau0 = (taustr <= tau_str_L) .* swF_tau0;
    end
