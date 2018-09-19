% Alg-2: Accelerated AOA algorithm using random ferns
% Author: hui.chen@kaust.edu.sa
% Last revised on 2018-09-19
close all;
clear all;
clc;
%% 1d phase difference calculation.
s1 = [0.8839 0];    % sensor 1 location in a 2-d plane
s2 = [-0.8839 0];
d = norm(s1-s2);    % distance between the two sensors
distance = 1000;    % assuming a far field model
v_speed = 343;      % speed of the sound
fs = 192000;        % sampling frequency
nfft = 1024;
F = 0:fs/nfft:(fs-fs/nfft);     % frequency axis after fft
freq = F([102 108 113 118]);    % carrier signal frequency
interval = 1;
ang_range = -89:interval:90;    % search from -89 to 90 with interval
phi = zeros(1,length(freq));
phase_all = [];     % phase difference of all the angles in ang_range
for i = 1:length(ang_range)
    angle0 = ang_range(i);
    td = -sin(angle0/180*pi)*d; % time difference of arrival
    phi =  2*pi*td/100./(v_speed./freq);
    phi = wrapping(phi);
    phase_all = [phase_all; phi];
end
size(phase_all)                 % size = searching angles by freqs
%% create error_table
error_table = [];
error = zeros(size(ang_range));
for num = 1:length(ang_range)
    phase = phase_all(num,:);
    for x = 1:length(ang_range)
        diff_phase = phase_all(x,:);
        m = length(diff_phase);
        error_temp = wrapping(diff_phase - phase);
        error(x) = sum(abs(error_temp));
    end
    error_table = [error_table; error];
end
figure;plot(ang_range,error_table(35+90,:));    % the error table at 35 [deg]
title('Search-based AOA Estimation [\theta = 35]')
xlabel('Searching Angle');ylabel('Error');
%% fern matrix creation
fern_num = 5;
fern_point = [47   -29     7   -69    -4] + 90;
rn = combnk(fern_point,2);  % create all possible combinations
r1 = rn(:,1);
r2 = rn(:,2);
bin = zeros(1,2^10);
fern_mat = [];
result = zeros(1,length(ang_range));
for i = 1:length(ang_range)
    value = sign(error_table(i + (180-searching_area)/2,r1)-error_table(i + (180-searching_area)/2,r2)+0.001);
    result(i) = (value+1)/2*fliplr(2.^(0:9))';
    bin(1+result(i)) = bin(1+result(i)) + 1;
    fern_mat = [fern_mat ;value];
end
% size(fern_mat)
figure;plot(bin);
xlabel('Cluster Index in Dec');ylabel('Num of Cluster Elements');

max(bin)
sum(bin)/length(nonzeros(bin))

figure;plot(ang_range,error_table(71,:))
hold on;plot(fern_point-90,error_table(71,fern_point),'r*')
% value = sign(error_table(91,r1)-error_table(91,r2)+0.001);
title('Search-based AOA Estimation')
xlabel('Horizontal Angle Theta');ylabel('Error');

% save fern_raw.mat r1 r2 fern_mat raw_phase
% save fern_vec.mat r1 r2 fern_mat vec_phase

%% error table test
test_num = -33 + 90;
result0 = error_table(test_num,r1)-error_table(test_num,r2);
value0 = sign(result0)

for i = 1:180
    value = sign(error_table(i,r1)-error_table(i,r2));
    result(i) = value0*value';
end
figure;plot(ang_range,(result));
title('Pattern Recognition for Angles');
xlabel('Candidate Angle Value [Deg]');ylabel('Cross Correlation Value');
temp = find(result==10)-90
% error_table(76,(find(result==10)-71+90))
% figure;plot(error_table(76,:))
%% test AOA algorithm
angle_real = -55.5;              
td = -sin(angle_real/180*pi)*d;
phi_obs =  2*pi*td/100./(v_speed./freq);
    
noise = 0.01;   % standard deviation of the noise
phi_obs = phi_obs + randn(size(phi_obs))*noise;
phase_test = wrapping(phi_obs);

angle_est = alg_aoa_fern(phase_test,r1,r2,fern_mat, phase_all)

disp(['The estimated angle is ', num2str(angle_est,'%2.2f'), '[deg] | error = ', ...
    num2str(abs(angle_real-angle_est),'%2.2f'), '[deg]']);
%% Extra part: best fern selection test
e = [];
r_all = [];
range_max = 140;
fern_points = 5;
error_table_n = error_table(90-floor(range_max/2)+1:90+floor(range_max/2),:);
for sel_i = 1:100
    fern_point = randperm(range_max,fern_points);
    rn = combnk(fern_point,2);
    fern_num = length(rn(:,1));
    r1 = rn(:,1);
    r2 = rn(:,2);
    bin = zeros(1,2^fern_num);
    fern_mat = [];
    value = [];
    for i = 1:range_max
        value = sign(error_table_n(i,r1)-error_table_n(i,r2)+0.001);
        result(i) = (value+1)/2*fliplr(2.^(0:fern_num-1))';
        bin(1+result(i)) = bin(1+result(i)) + 1;
        fern_mat = [fern_mat ;value];
end
% size(fern_mat)
% figure;plot(bin);
p = bin/sum(bin);
entropy = - sum(nonzeros(p).*log2(nonzeros(p))) / ceil(log2(length(p)));
e(sel_i) = entropy;
r_all = [r_all; fern_point];
% save fern.mat r1 r2 fern_mat
end
% r_all(37,:)
figure;plot(e);xlabel('Simulation Time');ylabel('Efficiency')
%% visualize the selected fern
[a1 a2] = max(e);
a1
fern_point = r_all(a2,:);
rn = combnk(fern_point,2);
r1 = rn(:,1);
r2 = rn(:,2);
bin = zeros(1,2^fern_num);
fern_mat = [];
value = [];
for i = 1:range_max
    value = sign(error_table_n(i,r1)-error_table_n(i,r2)+0.001);
    result(i) = (value+1)/2*fliplr(2.^(0:fern_num-1))';
    bin(1+result(i)) = bin(1+result(i)) + 1;
    fern_mat = [fern_mat ;value];
end
% size(fern_mat)
figure;plot(bin);
length(value)
a1
length(bin)
length(nonzeros(bin))
max(bin)
sum(bin)/length(nonzeros(bin))

