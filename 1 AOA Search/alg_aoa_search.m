% Search based aoa algorithm
% Author: hui.chen@kaust.edu.sa
% Revised on 28-Mar-2017
% Cleaned redundant part on 2018-09-19

function output = alg_aoa_search(phi,phase_all)
% raw search
    ang_range = -89:90;       % search area from -89 to 90.
    error = zeros(1,length(ang_range));
    for i = 1:length(ang_range)
        error_temp = wrapping(phi-phase_all(i,:));
        error(i) = sum(abs(error_temp));
    end
    % figure;plot(error,'-^')
    [~, a2] = min(error);
    output = ang_range(a2);

% fine search
    if(a2~=1 && a2 ~= length(ang_range)) 
        xc = [a2-1 error(a2-1);a2 error(a2);a2+1 error(a2+1)];
        xc2 = 0.5*(xc(3,2)-xc(1,2))/(2*xc(2,2)-xc(3,2)-xc(1,2));
        output = output + xc2;
    end

end
