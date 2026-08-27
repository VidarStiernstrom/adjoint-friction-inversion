function  swF_D_c = swF_D_c(s, taur, taup, D_c, tau_str_L)
    % Nov 10 2025: add tau_str_L
    default_arg('tau_str_L', NaN);

    swF_D_c = (s < D_c) .* (taup - taur) .* s ./ (D_c.^2) ;

    if ~isnan(tau_str_L)
        taustr = (((s<D_c)&(s>0)).*(taup-(taup-taur).*s./D_c) + (s>=D_c).*taur + (s<=0).*taup);
        swF_D_c = (taustr <= tau_str_L) .* swF_D_c;
    end

    % if s < D_c
    %     swF_D_c = (taup -taur) * s / (D_c^2);
    % else
    %     swF_D_c = 0;
    % end
