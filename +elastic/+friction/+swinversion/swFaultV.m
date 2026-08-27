function Vstar = swFaultV(s, eta, tau_l, taup, taur, D_c, tau_str_L)
    % Oct 16 2025: add tau_str_L
    default_arg('tau_str_L', NaN);

    taustr = (((s<D_c)&(s>0)).*(taup-(taup-taur).*s./D_c) + (s>=D_c).*taur + (s<=0).*taup);
    if ~isnan(tau_str_L)
        taustr = min(taustr, tau_str_L);
    end

    % Note that taul is negative!
    Vstar = (-tau_l>taustr).*(-tau_l- taustr)./eta;
    taustar = (-tau_l-eta.*Vstar);