% this function provides some disccussions on the PDP algorithm ot improve
% the performance.
function [output] = alg_projection_constrain_mirror(target, P, U, freq, d, thetam)
    global v;
    if(length(d)==1)
        ratio = norm(d.*2*pi*freq/v); 
    %******************* normal projection algorithm **********************
        rx_p = project_nd(freq/1000, target);
        [m, n] = size(P);
        dis_p = zeros(1,m);
        if(m == 1)
            dis_p = norm(P - rx_p);
        else
            for i = 1:length(P)
                dis_p(i) = norm(P(i,:)-rx_p);
            end    
        end  
        [~, a2] = min(dis_p);
        real_target = target + U(a2,:);
        real_proj = project_nd(freq,real_target);
        dt = sign(mean(real_target))*norm(real_proj - real_target)/ratio;
        out_est = real(asin(dt)*180/pi);
        output = out_est;
    
    else
        % original PDP
        phi_max = 2*pi*freq*sin(thetam/180*pi)*d/v;
        rx_p = project_nd(d, target);
        [m, n] = size(P);
        dis_p = zeros(1,m);
        if(m == 1)
            dis_p = norm(P - rx_p);
        else
            for i = 1:length(P)
                dis_p(i) = norm(P(i,:)-rx_p);
            end    
        end  
        [B, I] = sort(dis_p);
%         real_target = d/norm(d)*dis_to_pro + U(I(1),:) + P(I(1),:); 
%         output = real(asin(real_target(1)*v/2/pi/freq/d(1))*180/pi);

        % constrained PDP within the candidate area
        for ind = 1:length(I)
            dis_to_pro = d*target'/norm(d);
            real_target = d/norm(d)*dis_to_pro + U(I(ind),:) + P(I(ind),:);
            if(abs(real_target(1)) < phi_max(1))
                Lc = B(ind);
                break;
            end
        end
        output = real(asin(real_target(1)*v/2/pi/freq/d(1))*180/pi);
        
        % mirrored PDP to remove outliers
        Lmin_doa = 0.5;
        for mirror_ind = 1:length(target)
            mirror = 0;
            psi_temp = target;
            if(psi_temp(mirror_ind)+pi < Lmin_doa/2)
                psi_temp(mirror_ind) = psi_temp(mirror_ind) + 2*pi;
                mirror = 1;
            elseif(pi-psi_temp(mirror_ind)< Lmin_doa/2)
                psi_temp(mirror_ind) = psi_temp(mirror_ind) - 2*pi;
                mirror = 1;
            end
            if(mirror == 1)
                rx_p = project_nd(d, psi_temp);
                dis_p = vecnorm(P - rx_p, 2, 2);
                [B,I] = sort(dis_p);
                if(B(1) < Lc)
                    Lc = B(1);
                    for ind = 1:length(I)
                        dis_to_pro = d*psi_temp'/norm(d);
                        real_target = d/norm(d)*dis_to_pro + U(I(ind),:) + P(I(ind),:);
                        if(abs(real_target(1)) < phi_max(1))
                            Lc = B(ind);
                            break;
                        end
                    end
                    output = real(asin(real_target(1)*v/2/pi/freq/d(1))*180/pi);
                end
            end
        end
        
        
    end
    
    if(output>70) output = 70; end
    if(output<-70) output = -70; end
end