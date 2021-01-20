function [output] = alg_projection_constrained(target, search_area)
global d freq v projections changes;
ratio = norm(d.*2*pi*freq/v); 
%************************ projection algorithm **********************
    [m,n] = size(projections);
    rx_p = project_nd(freq/1000, target);
    dis_p = zeros(1,m);
    for i = 1:m
        dis_p(i) = norm(projections(i,:)-rx_p);
    end    
    [dis_est, a2] = min(dis_p);
    real_target = target + changes(a2,:);
    real_proj = project_nd(freq,real_target);
    dt = sign(mean(real_target))*norm(real_proj - real_target)/ratio;
    out_est = real(asin(dt)*180/pi);

%******************* constrained projection algorithm **********************
    if(length(search_area) == 2)
    % stable tracking
        [~,dis_sorted] = sort(dis_p);
        for cand_i = 1:min(3,(length(dis_sorted)-1))
            if((out_est>=search_area(1) && out_est<=search_area(2)))
                break;
            else
                real_target = target + changes(dis_sorted(cand_i+1),:);
                real_proj = project_nd(freq,real_target);
                dt = sign(mean(real_target))*norm(real_proj - real_target)/ratio;
                out_est = real(asin(dt)*180/pi);  
            end
        end
    end
    
    if(out_est<= search_area(2) && out_est>=search_area(1))
        output = out_est;
    elseif(out_est<search_area(1))
        output = search_area(1);
    else
        output = search_area(2);
    end
end