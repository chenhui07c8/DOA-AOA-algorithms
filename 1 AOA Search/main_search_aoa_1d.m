% Alg-1: Search based phase unwrapping algorithm
% Author: hui.chen@kaust.edu.sa
% Last revised on 2018-09-19
close all;
clear all;
clc;
%% 1d phase difference calculation.
s1 = [0.8839 0];    % sensor 1 location in a 2-d plane, in [cm]
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
%% search algorithm: estimate aoa from observed phase difference
angle_real = 55.3;              
td = -sin(angle_real/180*pi)*d;
phi_obs =  2*pi*td/100./(v_speed./freq);
noise = 0.1;   % standard deviation of the noise
phi_obs = phi_obs + randn(size(phi_obs))*noise;

angle_est = alg_aoa_search(phi_obs, phase_all);
disp(['The estimated angle is ', num2str(angle_est,'%2.2f'), '[deg] | error = ', ...
    num2str(abs(angle_real-angle_est),'%2.2f'), '[deg]']);





