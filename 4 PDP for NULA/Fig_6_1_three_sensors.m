close all;
clear all;
clc;
%%
global v S_lambda;
nLoops  = 1;
freq    = 20000; % frequency (Hz)
Fs      = 4*freq;  % sampling frequency (Hz)
v = 343;
D = 2.4;
delta = 1.25;
lambda = v/freq;
S = [0 1 1+delta]*D*lambda/2; % sensor location, in [lambda]
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

snrdB   = 0:2.5:60;
% snrdB = 50;
nSnr    = length(snrdB);
rvec1 = zeros(1,length(snrdB)); % pdp 3d
rvec2 = zeros(1,length(snrdB)); % pdp 2d
rvec3 = zeros(1,length(snrdB)); % music
rvec4 = zeros(1,length(snrdB)); % SD (spatial diversity)
rvec5 = zeros(1,length(snrdB)); % MLE
rvec6 = zeros(1,length(snrdB)); % 2q-order
rvec7 = zeros(1,length(snrdB)); % MLE 1

time_all = zeros(1,7);
doa_center = 40;
for i_snr = 1:length(snrdB)
    snr = snrdB(i_snr)
%     ang_range = doa_center + randn(1,10000)*1;
    ang_range = doa_center + rand(1,1000)-0.5;

    results1 = zeros(size(ang_range));
    results2 = zeros(size(ang_range));
    results3 = zeros(size(ang_range));
    results4 = zeros(size(ang_range));
    results5 = zeros(size(ang_range));
    results6 = zeros(size(ang_range));
    results7 = zeros(size(ang_range));

    nLoops = 1;
    for i_angle = 1:length(ang_range)
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
        
        tic;
        doa = alg_projection_constrain(p_test, P, U, freq, d_3d, thetam);
        time_all(1) = time_all(1)+toc;
        r1 = doa - temp_angle;
        results1(i_angle) = r1;
        
        
        % 2q order
        tic;
        [r1, r2] = alg_2qorder_triplet(p_test, d_3d, lambda);
        time_all(6) = time_all(6)+toc;
        results6(i_angle) = r2 - temp_angle;
        
        
        % SD, spatial diversity
        tic;
        doa = alg_disambiguity_spatial(p_test, freq, d_3d);
        time_all(4) = time_all(4)+toc;

        r4 = doa-temp_angle;
        results4(i_angle) = r4;

        
        % pdp 2d
        test = y(comb_2d(:,1),:).*transpose(y(comb_2d(:,2),:)');
        p_test = atan2(imag(test),real(test))';
        [P, U] = alg_get_projections(thetam, d_2d, v, freq);
        tic;
        doa = alg_projection_constrain(p_test, P, U, freq, d_2d, thetam);
        time_all(2) = time_all(2) + toc;
        r2 = doa - temp_angle;
        results2(i_angle) = r2;
        
        % MLE
        tic;
        output = alg_mle(estR, S_lambda, thetam, .2);
        time_all(5) = time_all(5)+toc;
        r5 = output - temp_angle;
        results5(i_angle) = r5;
        
        % MLE
        tic;
        output = alg_mle(estR, S_lambda, thetam, 1);
        time_all(7) = time_all(7)+toc;
        r7 = output - temp_angle;
        results7(i_angle) = r7;
        
        
        % music
        tic;
        output = alg_music_complex(estR, thetam, .2);
        time_all(3) = time_all(3)+toc;
        r3 = output-temp_angle;
        results3(i_angle) = r3;


    end
    rvec1(i_snr) = mean((results1).^2);
    rvec2(i_snr) = mean((results2).^2);
    rvec3(i_snr) = mean((results3).^2);
    rvec4(i_snr) = mean((results4).^2);
    rvec5(i_snr) = mean((results5).^2);
    rvec6(i_snr) = mean((results6).^2);
    rvec7(i_snr) = mean((results7).^2);


end


% figure;plot(results1)
% 
% figure;plot(results2)

%%
% save fig_6_1_data.mat;
% load fig_6_1_data.mat;
% time_all/length(snrdB)
%% crlb
d = S/(lambda/2);
M = length(d);
snr = M*10.^(snrdB/10);
result = 0;
for i = 1:M
    result = result + (d(i)-mean(d))^2;
end
result = result/M;
CRB = 1/2/pi^2/result./snr;
% CRB = 1./snr2/result;

crlb = (CRB*(180/pi/cos(39.5/180*pi))^2);
% hold on; semilogy(snrdB, crlb, 'LineWidth',10);

plot_ind = 1:17;
linewidth = 1.5;
markersize = 6;
figure;
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(rvec2(plot_ind)),'-.+b','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'-sg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec6(plot_ind)),'-.+m','LineWidth',linewidth,'MarkerSize',markersize);
hold on; semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);

hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

% axis([0 60 10^-5.5 10^4])
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
legend('PDP','SD','MUSIC (0.2^\circ)','MLE (0.2^\circ)','MLE (1^\circ)','2Q-Order','CRLB')
% legend('PDP-3D','PDP-2D','MUSIC','SD','MLE','2Q-Order','CRLB')
grid on;
set(gca,'FontSize',12)
% set(gcf,'position',[100,100,420*1.2,400*1.2])
%%

linewidth = 1.5;
markersize = 7;
figure;semilogy(snrdB,(rvec2),'-c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec3),'-r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec4),'-m','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec6),'-.g','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec1),'-b','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB,(rvec5),'-.k','LineWidth',linewidth,'MarkerSize',markersize);

% hold on;semilogy(snrdB,crlb,'-sk','LineWidth',linewidth,'MarkerSize',markersize);
hold on; semilogy(snrdB, crlb,'--b','LineWidth',2);

% axis([0 60 10^-5.5 10^4])
xlabel('SNR [dB]');
ylabel('MSE [Deg]');
% legend('PDP-3d','PDP-2d','MUSIC','SD')
legend('PDP-2D','MUSIC','SD','2Q-Order','PDP-3D','MLE','CRLB')
grid on;
set(gca,'FontSize',12)

% print -dpng -r600 sim-3sensors.png