%% 1 Initialization: fhss signal
% global dir fs b_len h_len num_hops p_len chan_min chan_bw chan_max ...
% num_bits num_pulses v carr_gain bw_carr f_carr tx_ref
global tx_ref;
global nfft fs F freq v d;
dir = 'C:\Users\chenh0c\Desktop\Demo_2018\signal\';
s1 = [1.25/sqrt(2) 0];     % sensor position definition, 1d aoa in [cm]
s2 = [-1.25/sqrt(2) 0];
d = norm(s1-s2)/100;    % in meter
r = 150;
nfft = 1024;
fs = 96000;
F = 0:fs/nfft:(fs-fs/nfft);
b_array = [183 199 215]; %96K 17 18.5 20; 
%     b_array = [193 209 225]; %96K 18 19.5 21;
freq = F(b_array);
fs = 96000;            % sampling frequency
b_len = 1920;           % length of one block ( hops--> pulse--> block )
h_len = 400/2;            % length of hop in samples
num_hops = 3;           % hop numbers in the pulse    
p_len = num_hops*h_len; % pulse length
chan_min = 17000;       % channel minimum freq (Hz)
chan_bw = 3000;         % channel bandwidth (Hz)
chan_max = chan_min+chan_bw;                      % channel maximum freq (Hz)
num_bits = 24;          % quantization bits
num_pulses = 2000;       
v = 343;                % speed of sound
carr_gain = ones(1,num_hops);                     % gain for channel (no attenuation case all 1)                                                               
bw_carr  = chan_bw/(num_hops-1);                  % bandwidth of each carrier (Hz)
f_carr   = chan_min+bw_carr*(0:num_hops-1);       % carrier frequencies (Hz)                                            
t = 0:1/fs:(h_len-1)/fs;
tx_ref = zeros(1,b_len);
% tx_time = 0:length(tx_ref)figure;plot()
win = hamming(h_len)';
for i = 1:num_hops
    tx_ref(1+(i-1)*h_len:i*h_len) = win.*sin(2*pi*f_carr(i)*t);
end
% figure;plot(tx_ref)
% figure;plot(abs(fft(tx_ref,2^17)))
tx = zeros(1,b_len*num_pulses);
for i = 1:num_pulses
    tx(1+(i-1)*b_len:i*b_len) = tx_ref;
end
% figure;plot(tx(1:10000))
%%
% filename = strcat(dir,'tx_phone_18K_21K_3hops_20s.wav'); 
% audiowrite(filename,tx,fs,'BitsPerSample',num_bits);