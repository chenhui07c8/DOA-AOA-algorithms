close all;
clear all;
clc;
%%
% calculate Lmin and Ld for 4 setups L
% 3 sensor, 3 pair, 2 pair.
% 5 sensor, 10 pair, 5 pair.
%%
global v S_lambda;
Lset = zeros(4,2);
freq    = 20000; % frequency (Hz)
v = 343;
D = 1.2;
delta = 1.25;
lambda = v/freq;
thetam = 70;

S = [0 1 1+delta]*D*lambda; % sensor location, in [lambda]
% S = [0 1 3 6 13 20 27 31 35 36]*lambda;
S = [0 3 6 13 20 27 31 35 36]*lambda/2;


comb_3d = combnk(1:length(S),2);
d = abs(S(comb_3d(:,1)) - S(comb_3d(:,2)));
ang_range = (-thetam):(thetam);
[P, U] = alg_get_projections(thetam, d, v, freq);

Lmin = zeros(1,length(ang_range));
for i = 1:length(ang_range)
    doa = ang_range(i);
    psi = wrapping(2*pi*freq*sin(doa/180*pi)*d/v);
    rx_p = project_nd(d, psi);
    L = vecnorm(P-rx_p,2,2);
    L = sort(L);
    Lmin(i) = L(2);
end

Lset(1,1) = min(Lmin);
Lset(1,2) = norm(d);


Ns = floor(2*sin(thetam/180*pi)*d/lambda)
dr = round(d/lambda*100);
G = dr(1);
for i = 2:length(dr)
    G = gcd(G, dr(i));
end
di = dr/G
if(sum((dr/G)<Ns) == length(dr) && (d(1)/(dr(1)/G))>=(lambda/(2*sin(thetam))))
    disp('Can not be used...')
else
    disp('Can be used')
end
%%
comb_3d = combnk(1:length(S),2);
comb_2d = comb_3d([1 3],:);
d = abs(S(comb_2d(:,1)) - S(comb_2d(:,2)));
ang_range = (-thetam):(thetam);
[P, U] = alg_get_projections(thetam, d, v, freq);

Lmin = zeros(1,length(ang_range));
for i = 1:length(ang_range)
    doa = ang_range(i);
    psi = wrapping(2*pi*freq*sin(doa/180*pi)*d/v);
    rx_p = project_nd(d, psi);
    L = vecnorm(P-rx_p,2,2);
    L = sort(L);
    Lmin(i) = L(2);
end

Lset(2,1) = min(Lmin);
Lset(2,2) = norm(d);
% figure;plot(ang_range, Lmin)


Ns = floor(2*sin(thetam/180*pi)*d/lambda)
dr = round(d/lambda*100);
G = dr(1);
for i = 2:length(dr)
    G = gcd(G, dr(i));
end
di = dr/G
if(sum((dr/G)<Ns) == length(dr) && (d(1)/(dr(1)/G))>=(lambda/(2*sin(thetam))))
    disp('Can not be used...')
else
    disp('Can be used')
end
%%
S = [0 10 21 33 45]*lambda/2*0.8;
comb = combnk(1:length(S),2);
d = abs(S(comb(:,1)) - S(comb(:,2)));
thetam = 70;
ang_range = (-thetam):(thetam);
[P, U] = alg_get_projections(thetam, d, v, freq);

Lmin = zeros(1,length(ang_range));
for i = 1:length(ang_range)
    doa = ang_range(i);
    psi = wrapping(2*pi*freq*sin(doa/180*pi)*d/v);
    rx_p = project_nd(d, psi);
    L = vecnorm(P-rx_p,2,2);
    L = sort(L);
    Lmin(i) = L(2);
end

Lset(3,1) = min(Lmin);
Lset(3,2) = norm(d);


Ns = floor(2*sin(thetam/180*pi)*d/lambda)
dr = round(d/lambda*100);
G = dr(1);
for i = 2:length(dr)
    G = gcd(G, dr(i));
end
di = dr/G
if(sum((dr/G)<Ns) == length(dr) && (d(1)/(dr(1)/G))>=(lambda/(2*sin(thetam))))
    disp('Can not be used...')
else
    disp('Can be used')
end
%% 5 pairs.
S = [0 10 21 33 45]*lambda/2*0.8;
comb = combnk(1:length(S),2);
comb2 = comb([1 3 6 7 10],:);
d = abs(S(comb2(:,1)) - S(comb2(:,2)));
thetam = 70;
ang_range = (-thetam):(thetam);
[P, U] = alg_get_projections(thetam, d, v, freq);

Lmin = zeros(1,length(ang_range));
for i = 1:length(ang_range)
    doa = ang_range(i);
    psi = wrapping(2*pi*freq*sin(doa/180*pi)*d/v);
    rx_p = project_nd(d, psi);
    L = vecnorm(P-rx_p,2,2);
    L = sort(L);
    Lmin(i) = L(2);
end
Lset(4,1) = min(Lmin);
Lset(4,2) = norm(d);


Ns = floor(2*sin(thetam/180*pi)*d/lambda)
dr = round(d/lambda*100);
G = dr(1);
for i = 2:length(dr)
    G = gcd(G, dr(i));
end
di = dr/G
if(sum((dr/G)<Ns) == length(dr) && (d(1)/(dr(1)/G))>=(lambda/(2*sin(thetam))))
    disp('Can not be used...')
else
    disp('Can be used')
end
