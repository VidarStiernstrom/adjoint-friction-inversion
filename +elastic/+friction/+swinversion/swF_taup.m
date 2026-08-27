function  swF_taup = swF_taup(s, taur, taup, D_c, tau_str_L)
    % Nov 10 2025: add tau_str_L
    default_arg('tau_str_L', NaN);

    swF_taup = (s < D_c) .* (1 - s ./ D_c);

    if ~isnan(tau_str_L)
        taustr = (((s<D_c)&(s>0)).*(taup-(taup-taur).*s./D_c) + (s>=D_c).*taur + (s<=0).*taup);
        swF_taup = (taustr <= tau_str_L) .* swF_taup;
    end

    % if s < D_c
    %     swF_taup = 1 - s / D_c;
    % else
    %     swF_taup = 0;
    % end
