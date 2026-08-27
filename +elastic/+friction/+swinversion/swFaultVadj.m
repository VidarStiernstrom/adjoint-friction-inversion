function Vstar = swFaultVadj(defrupt, taustr_adj, eta, tau_l)
    % Note that taul and taustr is negative! <- is taustr really negative?
    Vstar = defrupt.*(-tau_l- taustr_adj)./eta;
    taustar = (-tau_l-eta.*Vstar);

    % Vstar = defrupt.* (- tau_l) ./ eta;
    % Vstar = - tau_l ./ eta;
    % taustar = (-tau_l-eta.*Vstar);

    % Vstar = -tau_l ./ eta;
    % Vstar = zeros(size(tau_l));

    % Vstar = defrupt.*(- tau_l - taustr_adj)/eta;

