function  swF_taur = swF_taur(s, taur, taup, D_c, tau_str_L)
    % Nov 10 2025: add tau_str_L
    default_arg('tau_str_L', NaN);

    swF_taur = (s < D_c) .* s ./ D_c + (s >= D_c) .* 1;

    if ~isnan(tau_str_L)
        taustr = (((s<D_c)&(s>0)).*(taup-(taup-taur).*s./D_c) + (s>=D_c).*taur + (s<=0).*taup);
        swF_taur = (taustr <= tau_str_L) .* swF_taur;
    end

    % if s < D_c
    %     swF_taur = s / D_c;
    % else
    %     swF_taur = 1;
    % end
