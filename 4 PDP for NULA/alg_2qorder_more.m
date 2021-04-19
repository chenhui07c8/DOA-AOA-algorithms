function [r1,r2] = alg_2qorder_more(psi, psi_ind, d, lambda)
    ds = abs(mean(diff(d(psi_ind))));  % difference set
    ps = mean(wrapping(diff(psi(psi_ind))));    % small phase difference
    r1 = real(asin(ps*lambda/2/pi/ds)*180/pi);  
    [a, b]= max(d);
    [a, b] = sort(d, 'descend');
    a = a(2);
    b = b(2);
    phi_large = (a)/ds*ps;          % large phase difference
    psi_large = wrapping(psi(b));
    phi_cycle = [-1 0 1] + round((phi_large-psi_large)/2/pi);
    phi_candidate = psi_large + phi_cycle*2*pi;
    [~,b] = min(abs(phi_candidate - phi_large));
    phi_large = phi_candidate(b);
    r2 = real(asin(phi_large*lambda/2/pi/a)*180/pi);
end