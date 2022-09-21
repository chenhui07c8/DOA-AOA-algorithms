close all;
clear all;
clc;
%%
global v S_lambda;
rng(1);
nLoops  = 1;
freq    = 20000; % frequency (Hz)
Fs      = 4*freq;  % sampling frequency (Hz)
v = 343;
lambda = v/freq;
% S = [0 1 3 6 13 20 27 31 35 36]*lambda/2;


% d0 = 0.2;   % antenna spacing in wavelength.
% S = [0 1 6]*lambda/2*d0;
% % S = [0 1 6 10 23 26]*lambda/2*d0;
% % S = [0 1 6 10 23 26 34 41]*lambda/2*d0;
% 
% S_lambda    = S/(lambda/2); % Array element positions
% S_integer = S_lambda/d0;

S_integer = [0 1 6];
% S_integer = [0 1 6 10 23];
% S_integer = [0 1 6 10 23 26 34 41];

d0 = 0.4;
S = S_integer*lambda/2*d0;
S_lambda = S/(lambda/2); % Array element positions


comb_3d = combnk(1:length(S),2);
% comb_3d = comb_3d(4:12,:);
d_all = abs(S(comb_3d(:,1)) - S(comb_3d(:,2)));

% comb_pdp = comb_3d([1 3 6 7 10],:);
% d_pdp = abs(S(comb_pdp(:,1)) - S(comb_pdp(:,2)));

% comb_pdp2 = comb_3d(7:10,:);
% d_pdp2 = abs(S(comb_pdp2(:,1)) - S(comb_pdp2(:,2)));

thetam = 70;
ang_range = (-thetam):(thetam);
Ns = floor(2*freq*sin(thetam/180*pi)*d_all/v);
Nr      = length(S);  %length(arrp) % number of sensors
Np      = 1;  % number of snapshot

snrdB   = (-10:2.5:30);
% snrdB = 15;
% snrdB = 60;
nSnr    = length(snrdB);
rvec1 = zeros(1,length(snrdB)); % pdp
rvec2 = zeros(1,length(snrdB)); % MUSIC 0.01
rvec3 = zeros(1,length(snrdB)); % MLE 0.01
rvec4 = zeros(1,length(snrdB)); % 2Q-Order or MLE 0.5
rvec5 = zeros(1,length(snrdB)); % Two-Step
rvec6 = zeros(1,length(snrdB)); % EM-ESPRIT 5

time_all = zeros(1,5);

doa_center = 40;
for i_snr = 1:length(snrdB)
    snr = snrdB(i_snr)
%     ang_range = -29:0.05:31;
%     ang_range = ones(1,50)*30;
    ang_range = doa_center + rand(1000,1)-0.5;
    results1 = zeros(size(ang_range));
    results2 = zeros(size(ang_range));
    results3 = zeros(size(ang_range));
    results4 = zeros(size(ang_range));
    results5 = zeros(size(ang_range));
	results6 = zeros(size(ang_range));

    nLoops = 1;
    for i_angle = 1:length(ang_range)
        temp_angle = ang_range(i_angle);
        for iLoop = 1:nLoops       
            theta   = temp_angle/180*pi;       
            a      = exp(1j*pi*S_lambda*sin(theta))';
            a      = sqrt(Nr)*a/norm(a);              
            t      = (0:Np-1)'/Fs;
            s      = exp(1j*pi*freq*t);
            s      = s/norm(s);
            y0     = a*transpose(s);  
            y      = awgn(y0, snr, 'measured');
%             y = y0;
            % covariance matrix
            estR   = zeros(Nr, Nr);
            for k=1:Np
                estR = estR + y(:, k)*y(:, k)';
            end
        end

        
        % two-step
        comb_ref = [ones(1,length(S)-1); 2:length(S)]';
        test = y(comb_ref(:,1),:).*transpose(y(comb_ref(:,2),:)');
        p_test_oc = atan2(imag(test),real(test))';
        doa = alg_offset_correction([0 p_test_oc], S_lambda);
        results5(i_angle) = doa - temp_angle;
        
        
        % pdp all
        test = y(comb_3d(:,1),:).*transpose(y(comb_3d(:,2),:)');
        p_test = atan2(imag(test),real(test))';
        [P, U] = alg_get_projections(thetam, d_all, v, freq);
%         size(P)
        tic;
        doa = alg_projection_constrain_mirror(p_test, P, U, freq, d_all, thetam);
%         time_all(1) = time_all(1)+toc;
        r1 = doa - temp_angle;
        results1(i_angle) = r1;
%     end 
        
        
        % pdp sub
%         test = y(comb_pdp(:,1),:).*transpose(y(comb_pdp(:,2),:)');
%         p_test = atan2(imag(test),real(test))';
%         [P, U] = alg_get_projections(thetam, d_pdp, v, freq);
% %         doa = alg_projection(p_test, P, U, freq, d_pdp);
%         doa = alg_projection_constrain_mirror(p_test, P, U, freq, d_pdp, thetam);
% 
%         r4 = doa - temp_angle;
%         results4(i_angle) = r4;
        
        % music
        tic;
        output = alg_music_complex(estR, thetam, 0.2);
        r2 = output-temp_angle;
        results2(i_angle) = r2;
        
        % mle fine
        % tic;
        output = alg_mle_fine(estR, S_lambda, thetam, 0.01);
        r3 = output - temp_angle;
        results3(i_angle) = r3;

        % mle coarse
        tic;
        output = alg_mle(estR, S_lambda, thetam, 0.5);
%         time_all(4) = time_all(4)+toc;
        r4 = output - temp_angle;
        results4(i_angle) = r4;
        

        % EM-esprit
        a = exp(1j*pi*(S_integer+1)'*d0*sin(theta));
        Y = awgn(a, snr, 'measured');
%         theta0 = output;
        theta0 = temp_angle+(rand(1,1)-0.5)*2*2;

        DOA = alg_em_esprit(theta0, Y, s, S_integer+1, d0, 1, 20);
        r6 = DOA(end);
        results6(i_angle) = r6 - temp_angle;
        

        

    end
    rvec1(i_snr) = mean((results1).^2);
    rvec2(i_snr) = mean((results2).^2);
    rvec3(i_snr) = mean((results3).^2);
    rvec4(i_snr) = mean((results4).^2);
    rvec5(i_snr) = mean((results5).^2);
    rvec6(i_snr) = mean((results6).^2);
%     rvec8(i_snr) = mean((results8).^2);


end

%%
% save fig_8_4_data_3.mat;
% load fig_8_4_data_3.mat;

% save fig_8_4_data_5.mat;
% load fig_8_4_data_5.mat;

% save fig_8_4_data_8.mat;
% load fig_8_4_data_8.mat;

% save fig_8_4_data_3_res2.mat;
% load fig_8_4_data_3_res2.mat;

% save fig_8_4_data_5_res2.mat;
% load fig_8_4_data_5_res2.mat;

% save fig_8_4_data_8_res2.mat;
% load fig_8_4_data_8_res2.mat;

d = S/(lambda/2);
M = length(d);
SNR = M*10.^(snrdB/10);
result = 0;
for i = 1:M
    result = result + (d(i)-mean(d))^2;
end
result = result/M;
CRB = 1/2/pi^2/result./SNR;
crlb = (CRB*(180/pi/cos(39.5/180*pi))^2);

plot_ind = 1:17;
linewidth = 1.5;
markersize = 6;

figure;
semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec2(plot_ind)),'--^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize);

hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-dg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec6(plot_ind)),'>m-.','LineWidth',linewidth,'MarkerSize',markersize*1.2);
hold on; semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);

hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
% semilogy(snrdB,(rvec4),'-+c','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;
% hold on;semilogy(snrdB,(rvec5),'-^m','LineWidth',linewidth,'MarkerSize',markersize);
% % semilogy(snrdB,(crlb),'-sk','LineWidth',linewidth,'MarkerSize',markersize);
% hold on; semilogy(snrdB, crlb,'--b','LineWidth',2);

% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
% legend('PDP-3d','PDP-2d','MUSIC','SD')
% legend('PDP-3d','PDP-2d','MUSIC','SD','2q Small Dis','2q Large Dis')
% legend('PDP-SUB','MUSIC','2Q-Order','PDP-ALL','MLE','CLRB')
legend('PDP','MUSIC (0.01^\circ)','MLE (0.01^\circ)','MLE (0.5^\circ)','Two-Step','EM-ESPRIT','CRLB');
grid on; set(gca,'FontSize',12)
% set(gcf,'position',[100,100,420*1.2,400*1.2])
% axis([5 30 10^-0.7 10^1.5])

% print -dpng -r600 sim-10sensors.png