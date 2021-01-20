% Get all the projection points.
% hui.chen@kaust.edu.sa
% I will try to make this code clean later.....
function [projections,changes,pro_start, pro_stop] = get_projections(theta, d, v, freq)

% 1 frequency case
if(length(d)==2)
    max_pro = d*sin(theta/180*pi)/(v/freq)*2*pi;
    cur = [0 0];
    cur_real = [0 0];
    change = [0 0];
    stop = [];
    start = [0 0];
    while((cur_real(1)<max_pro(1))&&(cur_real(2)<max_pro(2)))
        if(diff(([pi pi]-cur)./d)>0) % reach right
            temp = (pi-cur(1))*d(2)/d(1)+cur(2);
            stop_temp = [pi temp];
            stop = [stop; stop_temp];
            start_temp = stop_temp-[2*pi 0];
            start = [start; start_temp];
            change_temp = change(end,:)+[2*pi 0];
            change = [change; change_temp];
            cur_real = start(end,:)+change(end,:);
            cur = start(end,:);
        elseif(diff(([pi pi]-cur)./d)<0) % reach top
            temp = (pi-cur(2))*d(1)/d(2)+cur(1);
            stop_temp = [temp pi];
            stop = [stop; stop_temp];
            start_temp = stop_temp-[0 2*pi];
            start = [start; start_temp];
            change_temp = change(end,:)+[0 2*pi];
            change = [change; change_temp];
            cur_real = start(end,:)+change(end,:);
            cur = start(end,:);
        else                            % reach the corner
            stop_temp = [pi pi];
            stop = [stop; stop_temp];
            start_temp = stop_temp-[2*pi 2*pi];
            start = [start; start_temp];
            change_temp = change(end,:)+[2*pi 2*pi];
            change = [change; change_temp];
            cur_real = start+change(end,:);
            cur = start(end,:);
        end
    end
    % process start and stop points
    pro1 = start(1:end-1,:);
    pro2 = stop;
    pro2 = [pro2; -pro1(2:end,:)];
    pro1 = [pro1; -stop(2:end,:)];
    pro1(1,:) = -pro2(1,:);
    changes = [change(1:end-1,:); -change(2:end-1,:)];
    pro_start = pro1;
    pro_stop = pro2;
    projections = zeros(size(pro1));
    for i = 1:length(pro1(:,1))
        projections(i,:) = project_nd(d, pro1(i,:));
    end
% 2 frequency case    
elseif(length(freq)>=2)        
    max_pro = d*sin(theta/180*pi)/v*freq*2*pi;   % max phase difference, project on the figure;
    cur = zeros(size(freq));    % search projection points
    cur_real = zeros(size(freq));   % real 
    change = zeros(size(freq));
    stop = [];
    start = zeros(size(freq));
    max_pi = pi*ones(size(freq));
    while(sum(cur_real < max_pro)==length(freq))
        [~,m2] = min(1000*(max_pi-cur)./freq);
        temp = (pi - cur(m2)).*freq/freq(m2) + cur;
        stop_temp = temp; 
        stop = [stop; stop_temp];
        change_bit = zeros(1,length(freq));
        change_bit(m2) = 2*pi;
        temp = temp-change_bit;
        start_temp = temp;    
        start = [start; start_temp];
        change_temp = change(end,:)+change_bit;
        change = [change; change_temp];
        cur_real = start(end,:)+change(end,:);
        cur = start(end,:);
    end
    
    % process start and stop points
    pro1 = start(1:end-1,:);
    pro2 = stop;
    pro2 = [pro2; -pro1(2:end,:)];
    pro1 = [pro1; -stop(2:end,:)];
    pro1(1,:) = -pro2(1,:);
    changes = [change(1:end-1,:); -change(2:end-1,:)];
    pro_start = pro1;
    pro_stop = pro2;
    projections = zeros(size(pro1));
    for i = 1:length(pro1(:,1))
        projections(i,:) = project_nd(freq, pro1(i,:));
    end
% hold on;plot(projections(:,1),projections(:,2),'ro')
end




