% Random ferns aoa algorithm
% Author: hui.chen@kaust.edu.sa
% Revised on 28-Mar-2017
% Cleaned redundant part on 2018-09-19

function output = alg_aoa_fern(phase_test,r1,r2,fern_mat,phase_all)
    area = 90;
    fern_ind = unique([r1 r2]);
    error = zeros(1,180) + 500;

    for row = 1:length(fern_ind)
        error_temp = wrapping(phase_all(fern_ind(row),:)-phase_test);
        error(fern_ind(row)) = sum(abs(error_temp));
    end
    value = sign(error(r1)-error(r2)+0.001);
    v2 = [zeros(1,90-area) value*fern_mat(90-area+1:90+area,:)'];
    candidate = find(v2>=length(fern_mat(1,:)));
    if(isempty(candidate))
        candidate = find(v2>=length(fern_mat(1,:))-2);
        if(isempty(candidate))
            candidate = 90;
        end
    end

    for row = 1:length(candidate)
        current = candidate(row);
        error_temp = wrapping(phase_all(current,:) - phase_test); 
        error(candidate(row)) = sum(abs(error_temp'));
    end

    [~, a2] = min(error);
    if(a2~=1 && a2 ~= length(phase_all)) 
        error_temp = wrapping(phase_all(a2-1:a2+1,:) - repmat(phase_test,3,1));
        error(a2-1:a2+1) = sum(abs(error_temp'));
        xc = [a2-1 error(a2-1);a2 error(a2);a2+1 error(a2+1)];
        xc2 = 0.5*(xc(3,2)-xc(1,2))/(2*xc(2,2)-xc(3,2)-xc(1,2));
        output = (a2)-90 + xc2;
    end
end
