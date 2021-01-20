function output = alg_music_complex(estR, theta_b, resolution)
    global S_lambda;
    
    Rxx = estR;
    [EV,D]=eig(Rxx);%%%% 
    EVA=diag(D)';
    [EVA,I]=sort(EVA);
    EVA=fliplr(EVA);
    EV=fliplr(EV(:,I));
    ang_range0 = -theta_b:resolution:theta_b;
    derad = pi/180;         % deg -> rad
    kelm = length(S_lambda);
    dint = S_lambda;
    iwave = 1;
    for iang = 1:length(ang_range0)
        phim=derad*ang_range0(iang);
        a=exp(-1j*pi*dint*sin(phim)).';
        L=iwave;    
        En=EV(:,L+1:kelm);
        SP(iang)=(a'*a)/(a'*En*En'*a);
    end
    % figure;plot(ang_range2, abs(SP))
    [~, b] = max(abs(SP));
    coarse = ang_range0(b);
    
    ang_range1 = (-resolution:0.02:resolution) + coarse;
    for iang = 1:length(ang_range1)
        phim=derad*ang_range1(iang);
        a=exp(-1j*pi*dint*sin(phim)).';
        L=iwave;    
        En=EV(:,L+1:kelm);
        SP(iang)=(a'*a)/(a'*En*En'*a);
    end
    [~, b] = max(abs(SP));
    fine = ang_range1(b);
    
    ang_range1 = (-0.02:0.001:0.02) + fine;
    for iang = 1:length(ang_range1)
        phim=derad*ang_range1(iang);
        a=exp(-1j*pi*dint*sin(phim)).';
        L=iwave;    
        En=EV(:,L+1:kelm);
        SP(iang)=(a'*a)/(a'*En*En'*a);
    end
    [~, b] = max(abs(SP));
    fine1 = ang_range1(b);
    
    output = fine1;
end