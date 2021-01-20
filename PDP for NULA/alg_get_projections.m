% Get all the projection points.
% hui.chen@kaust.edu.sa
% I will try to make this code clean later.....
function [P, U] = alg_get_projections(theta, d, v, freq)

% single frequency case / >=3 receivers
if(length(freq)==1)
    max_pro = d*sin(theta/180*pi)/v*freq*2*pi;   % max phase difference, project on the figure;
    min_pro = d*sin(-theta/180*pi)/v*freq*2*pi; 
    cur = wrapping(min_pro);    % search projection points
    cur_real = min_pro;   % real 
    U = [];
    change = cur_real - cur;
    stop = [];
    start = [];
    max_pi = pi*ones(size(d));
    P = [];
    U = [];
    j = 1;
    while(sum(cur_real < max_pro)==length(d))
        start = [start; cur];
        U = [U; change];
        j = j+1;
        [~,m2] = min((max_pi-cur)./d);
        cur = (pi - cur(m2)).*d/d(m2) + cur;
        cur(m2) = cur(m2) - 2*pi;
        change(m2) = change(m2) + 2*pi;
        cur_real = cur + change;
    end
    
    P = zeros(size(start));
    for i = 1:length(start(:,1))
        P(i,:) = project_nd(d, start(i,:));
    end
% >=2 frequency case / 2 receivers
elseif(length(freq)>=2)        
    max_pro = d*sin(theta/180*pi)/v*freq*2*pi;   % max phase difference, project on the figure;
    cur = zeros(size(freq));    % search projection points
    cur_real = zeros(size(freq));   % real 
    change = zeros(size(freq));
    stop = [];
    start = [];
    max_pi = pi*ones(size(freq));
    P = [];
    U = [];
    j = 1;
    while(sum(cur_real < max_pro)==length(freq))
        start = [start; cur];
        U = [U; change];
        j = j+1;
        [~,m2] = min(1000*(max_pi-cur)./freq);
        cur = (pi - cur(m2)).*freq/freq(m2) + cur;
        cur(m2) = cur(m2) - 2*pi;
        change(m2) = change(m2) + 2*pi;
        cur_real = cur + change;
    end
    
    p_boarder = [start(1:j-1,:); -start(2:j-1,:)];
    U = [U(1:end,:); -U(2:end,:)];
    P = zeros(size(p_boarder));
    for i = 1:length(p_boarder(:,1))
        P(i,:) = project_nd(freq, p_boarder(i,:));
    end
% hold on;plot(projections(:,1),projections(:,2),'ro')
end




