function [output] = alg_projection_constrain(target, P, U, freq, d, thetam)
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
        
%         [~, a2] = min(dis_p);
        [B, I] = sort(dis_p);
        for ind = 1:length(I)
            dis_to_pro = d*target'/norm(d);
            real_target = d/norm(d)*dis_to_pro + U(I(ind),:) + P(I(ind),:);
            if(abs(real_target(1)) < phi_max(1))
                break;
            end
        end
        output = real(asin(real_target(1)*v/2/pi/freq/d(1))*180/pi);
    end
    if(output>70) output = 70; end
    if(output<-70) output = -70; end
end