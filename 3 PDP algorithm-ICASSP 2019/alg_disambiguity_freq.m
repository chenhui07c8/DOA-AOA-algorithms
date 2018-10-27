% unwrap the phase using 2 freq components and corresponding phase
% difference
function data_out = alg_disambiguity_freq(f1,f2,phi1,phi2,d)
mu = f1/(f2-f1);
% step 1
k = -1:1;
phi_cand = mu*(-phi1 + phi2 + 2*pi*k);  % phi1-phi2 /// phi2-phi1
% 2*pi*f1*d/c
% step 2
[phi_u1,kuv] = min(abs(phi_cand));
phi_u1 = phi_cand(kuv);
% step 3
ku = (phi_u1-phi1)/2/pi;
ku = round(ku);
% step 4
data_out = phi1+2*pi*ku;

% 
% thres = 4.5;
%         if(data_out>thres) 
%             data_out = data_out-2*pi;
%         end
%         if(data_out<-thres)
%             data_out = data_out+2*pi;
%         end
%         
%         phase(i) = data_out;
%        
        data_out = data_out/2/pi*343./f1;
        data_out = real(asin(data_out*sqrt(2)/0.025)/pi*180); 
%         
        
end