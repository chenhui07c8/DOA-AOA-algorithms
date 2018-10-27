close all; 
clear all;
clc;
%% 1 Set all the parameters
global tx_ref;
global nfft fs F freq v d;
global theta_max projections changes;
nfft = 1024;
fs = 96000;
F = 0:fs/nfft:(fs-fs/nfft);
P1_Initialization_phone
b_array = [183 199 215]; %96K 17 18.5 20; 
b_two = [1 3];  % use only two frequencies to fit the paper. (can use three actually..)
freq = F(b_array(b_two));
v = 343;
d = 0.025/sqrt(2);
theta_max = 89;
search_area = [-theta_max, theta_max];

error_s = zeros(1,7);
error_m = zeros(1,7);
anglexs = zeros(3,1000); % Dual Freq, Grid Search, PDP 
angleys = zeros(3,1000);
success = zeros(3,15);
angle_cell = cell(1,15);
%% 2 Load Data & Preprocess (find block peak)
for file_i = 1:15
    file_ind = file_i
    if(file_i<10)     % zero ahead 
        file_ind = ['data0' num2str(file_ind) '.mat'];
    else
        file_ind = ['data' num2str(file_ind) '.mat'];
    end
    load(file_ind);
    [cc,gesture_area,peak_pos] = Sub_Locating_Active_Area(y(1,:),tx_ref,0.0008);
    pulse_pos = peak_pos(5:1004);   % discard first and last several.
    g1 = zeros(size(freq));
    g2 = zeros(size(freq));
    tic
    phi_obs = [];
    phi_rf = [];    % phase for range refinement.
    phi_r = [];
    mangle1 = [];
    mangle2 = [];
    for i = 1:length(pulse_pos)
        part = y([1 2], pulse_pos(i):pulse_pos(i)+p_len-1);
        fs = 96000;
        signal_e = 0;
        noise_e = 0;
% check the SNR        
%         for snr_i = 1:999
%             signal_e = signal_e + sum(y(2,peak_pos(snr_i):peak_pos(snr_i)+600-1).^2);
%             noise_e = noise_e + sum(y(2,peak_pos(snr_i)+600:peak_pos(snr_i)+1200-1).^2);
%         end
%         snr(file_i) = 20*log10(sqrt((signal_e-noise_e) / noise_e));
        
        hamming_mat = [hamming(p_len/3) hamming(p_len/3)];
        Y1 = fft(part(:,1:200)'.*hamming_mat,nfft);     % fft for frequency 1
        Y2 = fft(part(:,201:400)'.*hamming_mat,nfft);
        Y3 = fft(part(:,401:600)'.*hamming_mat,nfft);
        Y = [Y1(b_array(1),:) ;Y2(b_array(2),:);Y3(b_array(3),:)];
        Y1 = Y(:,1); Y2 = Y(:,2); 
        for fre_i = 1:length(b_array)
            g_temp = Y(fre_i,2).*(transpose(Y(fre_i,1)'));
            g1(fre_i) = angle(g_temp);
        end
        phi_obs = [phi_obs; g1];  % phase for horizontal angle
    end
    toc

%% AOA estimation
    [projections,changes,pro_start, pro_stop] = get_projections(theta_max, d, v, freq);

    tic
    for angle_i = 1:length(phi_obs)
        anglexs(3,angle_i) = alg_projection_constrained(phi_obs(angle_i,b_two), search_area);
    end
    toc
    
    tic
    for angle_i = 1:length(phi_obs)
        ref = phi_obs(angle_i,b_two);
        anglexs(2,angle_i) = alg_aoa_search(ref,freq,search_area,d);
    end
    toc
    
    tic
    for angle_i = 1:length(phi_obs)
        ref = phi_obs(angle_i,b_two);
        anglexs(1,angle_i) = alg_disambiguity_freq(freq(1),freq(2),ref(1),ref(2),d);
    end
    toc
    results_all{file_i} = anglexs;
end
% save F8_1027.mat
%% Plot figure
% load F8_1027.mat
for i = 1:15
    anglexs = results_all{i};
    for s_i = 1:3
        success(s_i,i) = sum(abs(anglexs(s_i,:)) > 2.5);
    end
end
snr_set = -18:-4;
sr = (1000-success)/1000;
figure;
plot(snr_set,sr([1],:)','-.or','LineWidth',2);hold on;
plot(snr_set,sr([2],:)','--og','LineWidth',2);hold on;
plot(snr_set,sr([3],:)','-ob','LineWidth',2);

for i = 1:15
    anglexs = results_all{i};
    for s_i = 1:3
        success(s_i,i) = sum(abs(anglexs(s_i,:)) > 5);
    end
end
snr_set = -18:-4;
sr = (1000-success)/1000;
hold on;
plot(snr_set,sr([1],:)','-.^r','LineWidth',2);hold on;
plot(snr_set,sr([2],:)','--^g','LineWidth',2);hold on;
plot(snr_set,sr([3],:)','-^b','LineWidth',2);% legend('Search','Fine Search','PDP')
xlabel('SNR')
ylabel('Percentage')
set(gca,'FontSize',12)
legend('P(e<2.5^{o}) Dual Freq','P(e<2.5^{o}) Grid Search','P(e<2.5^{o}) PDP',...
'P(e<5^{o}) Dual Freq','P(e<5^{o}) Grid Search','P(e<5^{o}) PDP')
axis([-18 -4 0 1])
% grid on;
%% 3d plotting with synchronized signal
rmse = zeros(3,15);
for i = 1:15
anglex = results_all{i};
rmse(1,i) = sqrt(mean(anglex(1,:).^2));
rmse(2,i) = sqrt(mean(anglex(2,:).^2));
rmse(3,i) = sqrt(mean(anglex(3,:).^2));
end
figure;plot(rmse','LineWidth',2);
legend('Dual Freq','Grid Search','PDP')
xlabel('SNR')
ylabel('RMSE of the Error [Deg]')
set(gca,'FontSize',12)

