% Sub-Program: Locating gesture area
% find the gesture area for a sequence of gestures 
% cc: cross correlation value between the received signal and the reference
% start point of each gesture and their length in blocks
% Created on 2017-02-17
% Created by Chen Hui (hui.chen@kaust.edu.sa)
% Revised on 2017-03-24

function [cc,gesture_area,peak_pos] = Sub_Locating_Active_Area(y1,tx_ref,thres)

b_len = length(tx_ref);
cc = conv(y1,fliplr(tx_ref));
cc = cc(length(tx_ref):end);
% figure;plot(abs(cc));
peak = zeros(1,floor(length(y1)/b_len));
peak_pos = peak;

for i = 1:length(peak)
    [a,b] = max(abs(cc(1+(i-1)*b_len:i*b_len)));
    peak(i) = a;
    peak_pos(i) = b + (i-1)*b_len;
end
% figure;plot(peak);
% figure;plot(peak_pos)
% figure;plot(y1(peak_pos(50):peak_pos(50)+1000))
peak = [peak 0];
cc_flag = zeros(size(peak));
cc_flag((abs(peak)>thres))=1;
% cc_flag = [cc_flag 0];

% figure;plot(cc_flag,'*');axis([0 1400 -0.2 1.2])
counter_active = 1;
active = 0; % gesture active
start = [];
stop = [];
for i = 1:length(peak)
    if(cc_flag(i) && ~active)
        start(counter_active) = i;
        active = 1;
    end
    if(~cc_flag(i) && active)
        stop(counter_active) = i;
        counter_active = counter_active + 1;
        active = 0;
    end
end
% locating gesture area
% gesture_area = zeros(1,2);

% if(~isempty((start(2:end)-stop(1:end-1))<5))
if(length(start)>length(stop))
    start = start(1:end-1);
end

    temp2 = find((start(2:end)-stop(1:end-1))<5);
% end
if(~isempty(temp2)) 
    for i = 1:length(temp2)  
        stop(temp2(end-i+1))=stop(temp2(end-i+1)+1);
        stop(temp2(end-i+1)+1) = 0;
    end
end

temp1 = find((stop-start)>10);

gesture_area = [start(temp1) ; stop(temp1)]';
% gesture_area(1,:)
% active_area = gesture_area(:,2) - gesture_area(:,1);
%% find the bases of each gesture
% base = zeros(1,length(gesture_area(:,1)));
% for i = 1:length(active_area)
%     start_temp(i) = (gesture_area(i,1))*b_len;
%     [~,base(i)] = max(abs(cc(start_temp:start_temp+b_len)));
% end
% base = base + start_temp - b_len;
% active_area = active_area-1;
end