% V threshold definition: might not work for traction free fault
function rupt = defrupt(V)
    rupt = (V>0) .* 1;

% function rupt = defrupt(tau, s, taup, taur, D_c)
%     taustr = (((s<D_c)&(s>0)).*(taup-(taup-taur).*s./D_c) + (s>=D_c).*taur + (s<=0).*taup);
%     rupt = (tau>=taustr) .* 1;

  
% old definition

    % rupt = (abs(V) > 1e-10) .* 1;

    % rupt = ones(size(V));
    % rupt = zeros(size(V));

    % rupt = not(abs(V) > 1e-10) .* 1;

    %rupt = (abs(V) < 1e-10) .* 1; %  the opposite of the original definition

    % rupt = (V == 0).*1; % This is the opposite of the original definition
    
% function rupt = defrupt(s)
%     rupt = (s>0);
