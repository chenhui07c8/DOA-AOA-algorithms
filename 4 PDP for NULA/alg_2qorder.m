function [r1,r2] = alg_2qorder(psi, psi_ind, d, lambda)
    % d = d_all;
    % psi = p_test;
    ds = abs(mean(diff(d(psi_ind))));  % difference set
    ps = mean(wrapping(diff(psi(psi_ind))));    % small phase difference
    r1 = real(asin(ps*lambda/2/pi/ds)*180/pi);  
    [a, b] = sort(d, 'descend');
    a_max = a(1);
    b_max = b(1);
    a_min = a(end-1);
    b_min = b(end-1);
    % small distance
    a = a_min;
    b = b_min;
    phi_large = (a)/ds*ps;          % large phase difference
    psi_large = wrapping(psi(b));
    phi_cycle = [-1 0 1] + round((phi_large-psi_large)/2/pi);
    phi_candidate = psi_large + phi_cycle*2*pi;
    [~,b] = min(abs(phi_candidate - phi_large));
    phi_large = phi_candidate(b);
%     r2 = real(asin(phi_large*lambda/2/pi/a)*180/pi);
    phi_large_ref = phi_large*a_max/a_min;
    % large distance
    a = a_max;
    b = b_max;
%     phi_large = (a)/ds*ps;          % large phase difference
    phi_large = phi_large_ref;
    psi_large = wrapping(psi(b));
    phi_cycle = [-1 0 1] + round((phi_large-psi_large)/2/pi);
    phi_candidate = psi_large + phi_cycle*2*pi;
    [~,b] = min(abs(phi_candidate - phi_large_ref));
    phi_large = phi_candidate(b);
    r2 = real(asin(phi_large*lambda/2/pi/a)*180/pi);
end