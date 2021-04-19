function [r1,r2] = alg_2qorder_triplet(psi, d, lambda)
    
    % small distance
    ds = d(1)-d(3);
    ps = wrapping(psi(1)-psi(3));
    r1 = real(asin(ps*lambda/2/pi/ds)*180/pi);
    % large distance
    phi_large = (d(1)+d(3))/ds*ps;
    psi_large = wrapping(psi(2));
    phi_cycle = [-1 0 1] + round((phi_large-psi_large)/2/pi);
    phi_candidate = psi_large + phi_cycle*2*pi;
    [~,b] = min(abs(phi_candidate - phi_large));
    phi_large = phi_candidate(b);
    r2 = real(asin(phi_large*lambda/2/pi/d(2))*180/pi);



end