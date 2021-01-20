function [r1,r2] = alg_2qorder_more(psi, psi_ind, d, lambda)

    % small distance
    ds = abs(mean(diff(d(psi_ind))));
    ps = mean(wrapping(diff(psi(psi_ind))));
    r1 = real(asin(ps*lambda/2/pi/ds)*180/pi);
%     error1(i) = r1-angle0;
    % large distance
    [a, b]= max(d);
    phi_large = (a)/ds*ps;
    psi_large = wrapping(psi(b));
    phi_cycle = [-1 0 1] + round((phi_large-psi_large)/2/pi);
    phi_candidate = psi_large + phi_cycle*2*pi;
    [~,b] = min(abs(phi_candidate - phi_large));
    phi_large = phi_candidate(b);
    r2 = real(asin(phi_large*lambda/2/pi/a)*180/pi);
%     error2(i) = r2-angle0;


end