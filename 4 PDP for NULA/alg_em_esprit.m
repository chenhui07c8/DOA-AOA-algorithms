% theta_initial: initial guess of the DOA
% Y: received symbol
% s: transmitted symbol
% S: array position
% d: antenna interval
% Nsource: number of sources
% K: max iterations.

% [ref]: El Kassis, Carine, José Picheral, and Chafic Mokbel. "EM-ESPRIT 
% algorithm for direction finding with nonuniform arrays." In 2007 IEEE/SP 
% 14th Workshop on Statistical Signal Processing, pp. 453-457. IEEE, 2007.

function output = alg_em_esprit(theta_initial, Y, s, S_int, d, Nsource, K)
    
    L = S_int;
    Lvec = zeros(1, L(end));
    Lvec(L) = 1;
    Lvecp = ones(1, L(end)) - Lvec;
    Lp = find(Lvecp==1);
    E =  eye(L(end));
    G = E(:, L);
    Gp = E(:, Lp);
    N = length(S_int);
    hat_theta = zeros(1, K);
    hat_theta(1) = theta_initial;

    for k = 2:K
        % E-step
        Ap = exp(1j*pi*(Lp)'*d*sind(hat_theta(k-1)));
        hatX = G*Y + Gp*Ap*s;
        eigen_values = sort(eig(Y*Y'), 'descend');
        sigma0 = min(eigen_values(1:Nsource));
%         sigma = 1000;
        % M-step
        R = (hatX*hatX')/length(hatX) + Gp*Gp'*sigma0;
        Nv = L(end);
        % Eigen decomposition
        [Vi,Li] = eig(R);
        [L0,I] = sort(diag(Li),'descend');
        V = Vi(:,I);
        Vs = V(:,1:Nsource);
        Vs1=Vs(1:Nv-1,:);
        Vs2=Vs(2:Nv,:);
        xsi = Vs1\Vs2;
%   (Vs1'*Vs1)^-1*Vs1'*Vs2
%         xsi=linsolve(Vs(1:Nv-1,:),Vs(2:Nv,:));
        % DOA estimation
        doa = asind((angle(eig(xsi))/(pi*d)));
        if abs(doa-hat_theta(k-1)) < 1e-4
            hat_theta(end) = doa;
            break;
        end
        hat_theta(k) = doa;
    end

    output = real(hat_theta);

end