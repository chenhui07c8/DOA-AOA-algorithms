function output = alg_music_complex(estR, thetam, resolution)
    global S_lambda;
    
    Rxx = estR;
    [EV,D]=eig(Rxx);%%%% 
    EVA=diag(D)';
    [EVA,I]=sort(EVA);
    EVA=fliplr(EVA);
    EV=fliplr(EV(:,I));
    derad = pi/180;         % deg -> rad
    kelm = length(S_lambda);
    dint = S_lambda;
    iwave = 1;

%     ang_range0 = -thetam:resolution:thetam;
%     for iang = 1:length(ang_range0)
%         phim=derad*ang_range0(iang);
%         a=exp(-1j*pi*dint*sin(phim)).';
%         L=iwave;    
%         En=EV(:,L+1:kelm);
%         SP(iang)=(a'*a)/(a'*En*En'*a);
%     end
%     % figure;plot(ang_range0, abs(SP))
%     [~, b] = max(abs(SP));
%     coarse = ang_range0(b);
%     output = coarse;

    ang_range1 = (-thetam:0.2:thetam);
    for iang = 1:length(ang_range1)
        phim=derad*ang_range1(iang);
        a=exp(-1j*pi*dint*sin(phim)).';
        L=iwave;    
        En=EV(:,L+1:kelm);
        SP(iang)=(a'*a)/(a'*En*En'*a);
    end
    [~, b] = max(abs(SP));
    coarse = ang_range1(b);
    
    ang_range1 = (-0.2:resolution:0.2) + coarse;
    for iang = 1:length(ang_range1)
        phim=derad*ang_range1(iang);
        a=exp(-1j*pi*dint*sin(phim)).';
        L=iwave;    
        En=EV(:,L+1:kelm);
        SP(iang)=(a'*a)/(a'*En*En'*a);
    end
    [~, b] = max(abs(SP));
    fine1 = ang_range1(b);
% %     
    output = fine1;
end