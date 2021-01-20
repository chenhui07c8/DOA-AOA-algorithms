function [Lmin, Lv] = alg_indicator_estimation(doa, d, P, U, lambda, thetam)

angle0 = doa;
phi_max = 2*pi*sin(thetam/180*pi)*d/lambda;
phi = 2*pi*sin(angle0/180*pi)*d/lambda;
psi = wrapping(phi);
rx_p = project_nd(d, psi);
[m, ~] = size(P);
if(m == 1)
    dis_p = norm(P - rx_p);
    Lmin_doa = sqrt(length(d))*pi;
else
    dis_p = vecnorm(P - rx_p, 2, 2);
    [B,I] = sort(dis_p);
    for ind = 2:length(I)
        hat_target = psi + U(I(ind),:);
        if(abs(hat_target(end)) < phi_max(end))
            Lmin_doa = B(ind);
            break;
        end
    end
end

Lmin_mirror = ones(size(psi))*Lmin_doa;
for i = 1:length(psi)
    psi_temp = psi;
    if(psi(i)+pi < Lmin_doa/2)
        psi_temp(i) = psi_temp(i) + 2*pi;
    elseif(pi-psi(i)< Lmin_doa/2)
        psi_temp(i) = psi_temp(i) - 2*pi;
    end
    rx_p = project_nd(d, psi_temp);

    dis_p = vecnorm(P - rx_p, 2, 2);
    [B,I] = sort(dis_p);
    if(B(1) < 10e-9)
        Lmin_mirror(i) = B(2);
    else
        Lmin_mirror(i) = B(1);
    end
%     Lmin_doa = B(ind);
end
Lmin_mirror;
Lmin_doa;
Lmin = min([Lmin_mirror Lmin_doa]);
Lv = norm(d/lambda);


end