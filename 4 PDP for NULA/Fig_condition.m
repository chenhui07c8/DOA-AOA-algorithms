% Search based phase unwrapping algorithm tests
close all;
clear all;
clc;
%%
% project_3d([a b c], p)
%% Simulation Model Building for Phase Unwrapping
global v;
v = 343;
freq = 20000;
D = 6;
delta = 0.2;
lambda = v/freq;
s = [0 1 1+delta]*D*lambda; % sensor location, in [lambda]
% s = [0 1 3 6 13 20 27 31 35 36]*lambda/2;
s = [0 10 21 33 45]*lambda/2*.4;
% s = [0 1 6 10 23 26 34 41 53 55]*lambda/2;
comb = combnk(1:length(s),2);
% comb = comb(1:8,:);
% comb = comb(4:12,:);
d = abs(s(comb(:,1)) - s(comb(:,2)));
thetam = 70;
ang_range = (-thetam):(thetam);
Ns = floor(2*sin(thetam/180*pi)*d/lambda);

lambda = v/freq;
angle_train = [];
phase_real = [];
phase_wrapped = [];
phi = zeros(1,length(freq)-1);
p = [];
for x_i = -70:1:70
    angle0 = x_i;
    % angle = -33.3
    phi = 2*pi*freq*sin(angle0/180*pi)*d/v;
    phase_real = [phase_real; (phi)];
    phase_wrapped = [phase_wrapped; wrapping(phi)];
    % phi = sin(phi);
end
%
[P, U] = alg_get_projections(thetam, d, v, freq);
b = sort(vecnorm(P,2,2));
Lmin = b(2)


phi_max = 2*pi*freq*sin(thetam/180*pi)*d/v;
        
Lmin_doa_all = zeros(size(ang_range));
for doa_i = 1:length(ang_range)
    angle0 = ang_range(doa_i);
    phi = 2*pi*freq*sin(angle0/180*pi)*d/v;
    psi = wrapping(phi);
    rx_p = project_nd(d, psi);
    [m, n] = size(P);
    if(m == 1)
        dis_p = norm(P - rx_p);
        Lmin_doa = sqrt(length(d))*pi;
    else
        dis_p = vecnorm(P - rx_p, 2, 2);
        [B,I] = sort(dis_p);
        for ind = 2:length(I)
            hat_target = psi + U(I(ind),:);
            if(abs(hat_target(end)) < phi_max(end))
                Lmin_doa = B(ind);
                break;
            end
        end
    end
    Lmin_doa_all(doa_i) = Lmin_doa;
end

L_all_constrain = min(Lmin_doa_all)
        
        
%% condition
max(phase_real);

Ns
d./Ns;
% d/lambda*10000;
dr = round(d/lambda*100);
G = dr(1);
for i = 2:length(dr)
    G = gcd(G, dr(i));
end

dr/G
Ns
if(sum((dr/G)<Ns) == length(dr) && (d(1)/(dr(1)/G))>=(lambda/(2*sin(thetam))))
    disp('Can not be used...')
else
    disp('Can be used')
end

