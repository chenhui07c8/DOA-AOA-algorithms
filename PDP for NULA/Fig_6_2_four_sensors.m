close all;
clear all;
clc;
%%
global v S_lambda;
nLoops  = 1;
freq    = 20000; % frequency (Hz)
Fs      = 4*freq;  % sampling frequency (Hz)
v = 343;
% D = 1.2;
% delta = 1.25;
lambda = v/freq;
S = [0 10 21 33]*lambda/2*0.5;
% S = [0 10 21 33 45]*lambda/2*0.8;

S_lambda    = S/(lambda/2); % Array element positions

comb_3d = combnk(1:length(S),2);
% comb_3d = comb_3d(4:12,:);
d_all = abs(S(comb_3d(:,1)) - S(comb_3d(:,2)));
% 
% comb_pdp = comb_3d([1 3 6 7 10],:);
% d_pdp = abs(S(comb_pdp(:,1)) - S(comb_pdp(:,2)));

% comb_pdp2 = comb_3d(7:10,:);
% d_pdp2 = abs(S(comb_pdp2(:,1)) - S(comb_pdp2(:,2)));

thetam = 70;
ang_range = (-thetam):(thetam);
Ns = floor(2*freq*sin(thetam/180*pi)*d_all/v);
Nr      = length(S);  %length(arrp) % number of sensors
Np      = 1;  % number of snapshot

snrdB   = 0:2.5:60;
% snrdB = 15;
nSnr    = length(snrdB);
rvec1 = zeros(1,length(snrdB)); % pdp - all pairs
rvec2 = zeros(1,length(snrdB)); % pdp - sub
rvec3 = zeros(1,length(snrdB)); % music
rvec4 = zeros(1,length(snrdB)); % MLE
rvec6 = zeros(1,length(snrdB)); % music 1 deg
rvec7 = zeros(1,length(snrdB)); % MLE 1 deg
rvec5 = zeros(1,length(snrdB)); % 2q -order
time_all = zeros(1,7);
rng(7);
doa_center = 40;
for i_snr = 1:length(snrdB)
    snr = snrdB(i_snr)
%     ang_range = -29:0.05:31;
%     ang_range = ones(1,50)*30;
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

        % pdp all
        test = y(comb_3d(:,1),:).*transpose(y(comb_3d(:,2),:)');
        p_test = atan2(imag(test),real(test))';
        [P, U] = alg_get_projections(thetam, d_all, v, freq);
        tic;
        doa = alg_projection_constrain_mirror(p_test, P, U, freq, d_all, thetam);
        time_all(1) = time_all(1)+toc;
        r1 = doa - temp_angle;
        results1(i_angle) = r1;
        
        
        % 2q order
        comb = combnk(1:length(S),2);
        psi_ind = zeros(1,length(S)-1);
        for j = 1:(length(S)-1)
            psi_ind(j) = find((comb(:,1) == j)&(comb(:,2)==(j+1)));
        end
        tic;
        [~, r5] = alg_2qorder_more(p_test, psi_ind, d_all, lambda);
        time_all(5) = time_all(5)+toc;
        results5(i_angle) = r5 - temp_angle;

        % pdp sub
%         test = y(comb_pdp(:,1),:).*transpose(y(comb_pdp(:,2),:)');
%         p_test = atan2(imag(test),real(test))';
%         [P, U] = alg_get_projections(thetam, d_pdp, v, freq);
% %         doa = alg_projection(p_test, P, U, freq, d_pdp);
%         tic;
%         doa = alg_projection_constrain_mirror(p_test, P, U, freq, d_pdp, thetam);
%         time_all(2) = time_all(2)+toc;
%         r2 = doa - temp_angle;
%         results2(i_angle) = r2;
        
        % music
        tic;
        output = alg_music_complex(estR, thetam, 0.2);
        time_all(3) = time_all(3)+toc;
        r3 = output-temp_angle;
        results3(i_angle) = r3;
        
        % mle
        tic;
        output = alg_mle(estR, S_lambda, thetam, 0.2);
        time_all(4) = time_all(4)+toc;
        r4 = output - temp_angle;
        results4(i_angle) = r4;
        
        % music
%         tic;
%         output = alg_music_complex(estR, thetam, 1);
%         time_all(6) = time_all(6)+toc;
%         r6 = output-temp_angle;
%         results6(i_angle) = r6;
%         
        % mle
        tic;
        output = alg_mle(estR, S_lambda, thetam, 1);
        time_all(7) = time_all(7)+toc;
        r7 = output - temp_angle;
        results7(i_angle) = r7;

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
% save fig_6_2_data_four.mat;
% load fig_6_2_data_four.mat;

time_all / length(snr)
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
plot_ind = 1:21;
linewidth = 1.5;
markersize = 6;
figure;

semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB,(rvec2),'-dg','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec3(plot_ind)),'-^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec4(plot_ind)),'--dr','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec7(plot_ind)),'-*r','LineWidth',linewidth,'MarkerSize',markersize);
% hold on;semilogy(snrdB(plot_ind),(rvec6(plot_ind)),'--^c','LineWidth',linewidth,'MarkerSize',markersize);
hold on;semilogy(snrdB(plot_ind),sqrt(rvec5(plot_ind)),'-.+m','LineWidth',linewidth,'MarkerSize',markersize);

% semilogy(snrdB,(crlb),'-sk','LineWidth',linewidth,'MarkerSize',markersize);
hold on; semilogy(snrdB(plot_ind), sqrt(crlb(plot_ind)),'--k','LineWidth',2);

hold on;semilogy(snrdB(plot_ind),sqrt(rvec1(plot_ind)),'-ob','LineWidth',linewidth,'MarkerSize',markersize);

% axis([0 60 10^-7 10^4])
xlabel('SNR [dB]');
ylabel('RMSE [Deg]');
legend('PDP','MUSIC (0.2^\circ)','MLE (0.2^\circ)','MLE (1^\circ)','2Q-Order','CLRB');
% legend('PDP','MUSIC (0.2^\circ)','MLE (0.2^\circ)','MUSIC (1^\circ)','MLE (1^\circ)','2Q-Order','CLRB')
grid on; set(gca,'FontSize',12)
% set(gcf,'position',[100,100,420*1.2,400*1.2])

% print -dpng -r600 sim-4sensors.png