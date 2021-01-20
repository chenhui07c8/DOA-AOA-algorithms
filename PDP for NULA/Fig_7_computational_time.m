close all;
clear all;
clc;
%%
global v S_lambda;
nLoops  = 1;
freq    = 20000; % frequency (Hz)
Fs      = 4*freq;  % sampling frequency (Hz)
v = 343;
D = 1.2;
delta = 1.25;
lambda = v/freq;
S = [0 1 1+delta]*D*lambda; % sensor location, in [lambda]
S_lambda    = S/(lambda/2); % Array element positions

comb_3d = combnk(1:length(S),2);
d_3d = abs(S(comb_3d(:,1)) - S(comb_3d(:,2)));
ind_2d = [1 2];
comb_2d = comb_3d(ind_2d,:);
d_2d = abs(S(comb_2d(:,1)) - S(comb_2d(:,2)));

thetam = 70;
ang_range = (-thetam):(thetam);
Ns = floor(2*freq*sin(thetam/180*pi)*d_3d/v);

Nr      = length(S);  %length(arrp) % number of sensors
Np      = 1;  % number of snapshot

snrdB   = 0:5:60;
% snrdB = 50;
nSnr    = length(snrdB);
rvec1 = zeros(1,length(snrdB)); % pdp 3d
rvec2 = zeros(1,length(snrdB)); % pdp 2d
rvec3 = zeros(1,length(snrdB)); % music
rvec4 = zeros(1,length(snrdB)); % SD (spatial diversity)
rvec5 = zeros(1,length(snrdB)); % MLE
rvec6 = zeros(1,length(snrdB)); % 2q-order large
time_alg = zeros(1,6);
doa_center = 40;
trials = 1000;
for i_snr = 1
    snr = snrdB(i_snr)
    ang_range = doa_center + randn(1,1000)*1;
    results1 = zeros(size(ang_range));
    results2 = zeros(size(ang_range));
    results3 = zeros(size(ang_range));
    results4 = zeros(size(ang_range));
    results5 = zeros(size(ang_range));
    results6 = zeros(size(ang_range));
    
    nLoops = 1;
    for i_angle = 1
        temp_angle = ang_range(i_angle);
        for iLoop = 1:nLoops       
            theta   = temp_angle/180*pi;       
            a      = exp(1j*pi*S_lambda*sin(theta))';
            a      = sqrt(Nr)*a/norm(a);              
            t      = (0:Np-1)'/Fs;
            s      = exp(1j*pi*freq*t );
            s      = s/norm(s);
            y0     = a*transpose(s);  
            y      = awgn(y0, snr, 'measured');
            % covariance matrix
            estR   = zeros(Nr, Nr);
            for k=1:Np
                estR = estR + y(:, k)*y(:, k)';
            end
        end
        
        
        % pdp 3d
        
        test = y(comb_3d(:,1),:).*transpose(y(comb_3d(:,2),:)');
        p_test = atan2(imag(test),real(test))';
        % p_test    % -2.6278   -1.7140    0.9138
        [P, U] = alg_get_projections(thetam, d_3d, v, freq);
        
        tic
        for trials = 1:trials
            doa = alg_projection_constrain(p_test, P, U, freq, d_3d, thetam);
        end
        time_alg(1) = toc*1000;

        % 2q order
        tic
        for trials = 1:trials
            [r1, r2] = alg_2qorder_triplet(p_test, d_3d, lambda);
        end
        time_alg(6) = toc*1000

       
        
        % SD, spatial diversity
        tic
        for trials = 1:trials
            s12 = p_test(1)/freq/2/pi;
            s23 = p_test(3)/freq/2/pi;
            u1 = d_3d(1)/(d_3d(3)-d_3d(1));
            t12 = u1*(s23-s12+[-1 0 1]/freq);
            [~,k] = min(abs(t12));
            t12m = t12(k);
            k12 = round(freq*(t12m-s12));
            t12_hat = s12+k12/freq;
            output = real(asin(t12_hat*v/d_3d(1))/pi*180);
        end
        time_alg(4) = toc*1000


        % pdp 2d
        test = y(comb_2d(:,1),:).*transpose(y(comb_2d(:,2),:)');
        p_test = atan2(imag(test),real(test))';
        [P, U] = alg_get_projections(thetam, d_2d, v, freq);
        tic
        for trials = 1:trials
            doa = alg_projection_constrain(p_test, P, U, freq, d_2d, thetam);
        end
        time_alg(2) = toc*1000
        
        
        % MLE
        tic
        for trials = 1:trials
            output = alg_mle(estR, S_lambda, thetam);
        end
        time_alg(5) = toc*1000
        
        % musictic
        for trials = 1:trials
            output = alg_music_complex(estR, thetam);
        end
        time_alg(3) = toc*1000
        
    end
%     rvec1(i_snr) = mean((results1).^2);
%     rvec2(i_snr) = mean((results2).^2);
%     rvec3(i_snr) = mean((results3).^2);
%     rvec4(i_snr) = mean((results4).^2);
%     rvec5(i_snr) = mean((results5).^2);
%     rvec6(i_snr) = mean((results6).^2);


end


% figure;plot(results1)
% 
% figure;plot(results2)
%% crlb
snr = 10.^(snrdB/10);
% snr = 10.^(20/10);
d = S/(lambda/2);
M = length(d);
result = 0;
for i = 1:M
    result = result + (d(i)-mean(d))^2;
end
result = result/M;

CRB = 1/2./snr/result;
crlb = (CRB/pi*180*cos(30/180*pi))
%% crlb 2
snr = 10.^(snrdB/10);
% snr = 10.^(20/10);
d = S/(lambda/2);
M = length(d);
result = 0;
for i = 1:M
    result = result + (d(i)-mean(d))^2;
end
result = result/M;

snr2 = M*snr; 
CRB = 1/2/pi^2./snr2/result;
crlb = (CRB/pi*180*cos(30/180*pi))

% crlb = (M*snr+1)/2/pi^2/M^2./snr/result
%%
% save fig_5_2_data.mat;
% load fig_5_2_data.mat;
% save fig_5_2_data_40.mat;
% load fig_5_2_data_40.mat;
% save fig_5_2_data_40_1000.mat;
% load fig_5_2_data_40_1000.mat;

%%
linewidth = 2;
markersize = 7;
figure;semilogy(snrdB,(rvec2),'-+c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec3),'-^r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec4),'-*m','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec6),'-dg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec1),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec5),'-sk','LineWidth',linewidth,'MarkerSize',markersize);

% hold on;semilogy(snrdB,crlb,'-sk','LineWidth',linewidth,'MarkerSize',markersize);

axis([0 60 10^-5.5 10^4])
xlabel('SNR [dB]');
ylabel('MSE [Deg]');
% legend('PDP-3d','PDP-2d','MUSIC','SD')
legend('PDP-2D','MUSIC','SD','2Q-Order','PDP-3D','CRLB')
grid on;
set(gca,'FontSize',12)
% set(gcf,'position',[100,100,420*1.2,400*1.2])

