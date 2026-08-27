function taustar = swFaulttauadj(defrupt, taustr_adj, eta, tau_l)
    % Note that taul and taustr is negative!
    Vstar = defrupt.*(-tau_l- taustr_adj)./eta;
    taustar = (-tau_l-eta.*Vstar);

    % Vstar = - tau_l ./ eta;
    % taustar = (-tau_l-eta.*Vstar);

    % taustar = zeros(size(tau_l));
    % taustar =  -tau_l ./ eta;

    % taustar = defrupt.*(- taustr_adj) + (1 - defrupt) .* (- tau_l);
