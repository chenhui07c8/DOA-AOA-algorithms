close all;
clear all;
clc;
%%
global v S_lambda;
nLoops  = 1;
freq    = 20000; % frequency (Hz)
Fs      = 4*freq;  % sampling frequency (Hz)
v = 343;
lambda = v/freq;

% S_integer = [0 10 21];
S_integer = [0 10 21 33 46];
% S_integer = [0 10 21 33 46 60 75 91];

d0 = 0.5;
S = S_integer*lambda/2*d0;
S_lambda = S/(lambda/2); % Array element positions


comb_3d = combnk(1:length(S),2);
d_all = abs(S(comb_3d(:,1)) - S(comb_3d(:,2)));


thetam = 70;
ang_range = (-thetam):(thetam);
Ns = floor(2*freq*sin(thetam/180*pi)*d_all/v);
Nr      = length(S);  %length(arrp) % number of sensors
Np      = 1;  % number of snapshot

snrdB   = 0:2.5:40;
% snrdB = 35;
nSnr    = length(snrdB);
rvec1 = zeros(1,length(snrdB)); % pdp
rvec2 = zeros(1,length(snrdB)); % music
rvec3 = zeros(1,length(snrdB)); % MLE
rvec4 = zeros(1,length(snrdB)); % 2q -order
rvec5 = zeros(1,length(snrdB)); % Two-Stage
rvec6 = zeros(1,length(snrdB)); % ESPRIT

doa_center = 40;
for i_snr = 1:length(snrdB)
    rng(1);
    snr = snrdB(i_snr);
    disp(snr);
%     snr = 20
    ang_range = doa_center + rand(1,100)-0.5;
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

        % pdp all
        test_pdp = y(comb_3d(:,1),:).*transpose(y(comb_3d(:,2),:)');
        p_test = atan2(imag(test_pdp),real(test_pdp))';
        [P, U] = alg_get_projections(thetam, d_all, v, freq);
        tic;
%         doa = alg_projection(p_test, P, U, freq, d_all);  % original PDP algorithm
        doa = alg_projection_constrain_mirror(p_test, P, U, freq, d_all, thetam); % adding constraints, improves performance in low SNR
        r1 = doa - temp_angle;
        results1(i_angle) = r1;
        
        % mle fine
        tic;
        output = alg_mle_fine(estR, S_lambda, thetam, 0.01);
        r3 = output - temp_angle;
        results3(i_angle) = r3;
        
        % mle
%         tic;
%         output = alg_mle(estR, S_lambda, thetam, 0.5);
%         time_all(4) = time_all(4)+toc;
%         r4 = output - temp_angle;
%         results4(i_angle) = r4;

        % 2q order
        comb = combnk(1:length(S),2);
        coe = 1;
        psi_ind = zeros(1,length(S)-coe);
        for j = 1:(length(S)-coe)
            psi_ind(j) = find((comb(:,1) == j)&(comb(:,2)==(j+coe)));
        end
        tic;
        [~, r4] = alg_2qorder(p_test, psi_ind, d_all, lambda);
        results4(i_angle) = r4 - temp_angle;


        % two-step
        comb_ref = [ones(1,length(S)-1); 2:length(S)]';
        test = y(comb_ref(:,1),:).*transpose(y(comb_ref(:,2),:)');
        p_test_oc = atan2(imag(test),real(test))';
        doa = alg_offset_correction([0 p_test_oc], S_lambda);
        results5(i_angle) = doa - temp_angle;
        
        % EM-esprit
        a = exp(1j*pi*(S_integer+1)'*d0*sin(theta));
        Y = awgn(a, snr, 'measured');
        theta0 = temp_angle+(rand(1,1)-0.5)*2*2;
        DOA = alg_em_esprit(theta0, Y, s, S_integer+1, d0, 1, 20);
        r6 = DOA(end);
        results6(i_angle) = r6 - temp_angle;

        % music
        tic;
        output = alg_music_complex(estR, thetam, 0.01);
        r2 = output-temp_angle;
        results2(i_angle) = r2;
    end
    rvec1(i_snr) = mean((results1).^2);
    rvec2(i_snr) = mean((results2).^2);
    rvec3(i_snr) = mean((results3).^2);
    rvec4(i_snr) = mean((results4).^2);
    rvec5(i_snr) = mean((results5).^2);
    rvec6(i_snr) = mean((results6).^2);
end


% figure;plot(results1)
% 
% figure;plot(results2)

%%
% save fig_8_3_data_3.mat;
% load fig_8_3_data_3.mat;

% save fig_8_3_data_5.mat;
% load fig_8_3_data_5.mat;

% save fig_8_3_data_8.mat;
% load fig_8_3_data_8.mat;

% time_all / length(snr)
%%
% CRLB
result = 0;
d = S/(lambda/2);
M = length(d);
snr = M*10.^(snrdB/10);
for i = 1:M
    result = result + (d(i)-mean(d))^2;
end
result = result/M;

CRB = 1/2/pi^2/result./snr;
crlb = CRB*(180/pi/cos(39.5/180*pi))^2;
% hold on;semilogy(snrdB, crlb);
% plot_ind = 9:21;
% plot_ind = 9:17;
% plot_ind = 1:13;    
plot_ind = 1:length(snrdB);    

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
% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
legend('PDP','MUSIC (0.01^\circ)','MLE (0.01^\circ)','2Q-Order','Two-Step','EM-ESPRIT','CRLB');

grid on; set(gca,'FontSize',12)

% set(gcf,'position',[100,100,420*1.2,400*1.2])


% print -dpng -r600 sim-4sensors.png