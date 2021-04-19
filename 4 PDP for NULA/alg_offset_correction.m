function output = alg_offset_correction(phi, S_lambda)
% phi = [0 p_test_oc];
    Nr = length(S_lambda);
    Nbf = 90;
    bf_matrix = zeros(Nr, Nr);
    for i = 1:Nbf
        for j = 1:Nr
            angle_temp = ((i-1)/Nbf-0.5);
            bf_matrix(i,j) = 1/Nr*exp(-1j*pi*2*angle_temp*S_lambda(j));
        end
    end

    X = bf_matrix*exp(1j*phi'); % n*n muls
    [m,n] = max(abs(X));
    initial = asind(2*(n-1)/Nbf-1);
    psi = pi*S_lambda*sind(initial);
    eta = sum(wrapping(phi-psi))/sum(S_lambda);
    output = real(asind(sind(initial) + eta/pi));

end