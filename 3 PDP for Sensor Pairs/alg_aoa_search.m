% add candidates for unwrapping_search.
% Revised on 28-Mar-2017
%
%
function angle_out = alg_aoa_search(g1,f1,search_range,d)

% step1: load data
ang = search_range(1):1:search_range(2);
% d = 2.5/sqrt(2);
diff_range = -sin(ang/180*pi)*d;
% [mm, nn] = size(f1);
v = 343;
phi = g1;

% step2: generate all the possible phase difference (no need to do this every time, just load data)
diff_fai = zeros(length(diff_range),length(f1));% angle and frequency
for i = 1:length(diff_range)
    for j = 1:length(f1)
        diff_fai(i,j) = -diff_range(i)/(v/f1(j))*2*pi;
        while(diff_fai(i,j)>=pi || diff_fai(i,j)<-pi)
            if(diff_fai(i,j) > pi)
                diff_fai(i,j) = diff_fai(i,j)-2*pi;
            end
            if(diff_fai(i,j) <= -pi)
                diff_fai(i,j) = diff_fai(i,j)+2*pi;
            end
        end
    end
end

% step3: calculate error and find the min
error = zeros(1,length(diff_range));
    for i = 1:length(diff_range)
        error_temp = (phi-diff_fai(i,:));
        % noise may cause error. (eg: 0 and 1.99*pi are the same)
        extra_l = find(error_temp>1.5*pi);
        error_temp(extra_l) = error_temp(extra_l)-2*pi;
        extra_s = find(error_temp<-1.5*pi);
        error_temp(extra_s) = error_temp(extra_s)+2*pi;
        error(i) = sum(abs(error_temp));
    end
% figure;plot(error)
[a1, a2] = min(error);
% step4: do 2rd order curve fitting
if(a2==1) 
    a2 = a2+1;
end
if(a2 == length(ang)) 
    a2 = a2-1;
end
angle_out = ang(a2);

end
