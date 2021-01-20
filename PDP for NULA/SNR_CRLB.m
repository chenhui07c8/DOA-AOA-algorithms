

snr = 100000;
10*log10(snr)
% linear_snr = log10(snr)
d = [0 1 2.25];
M = length(d);


result = 0;
for i = 1:M
    result = result + (d(i)-1/M*(sum(d(1:end))))^2;
end

CRB = 1/2/M/snr

log10(CRB/pi*180/cos(30/180*pi))




% sin(41/180*pi)-sin(40/180*pi)
% sin(1/180*pi)
% cos(40/180*pi)