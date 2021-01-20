function output = alg_disambiguity_spatial(p_test, freq, d_3d)
    global v;
    s12 = p_test(1)/freq/2/pi;
    s23 = p_test(3)/freq/2/pi;
    u1 = d_3d(1)/(d_3d(3)-d_3d(1));
    t12 = u1*(s23-s12+[-1 0 1]/freq);
    [~,k] = min(abs(t12));
    t12m = t12(k);
    k12 = round(freq*(t12m-s12));
    t12_hat = s12+k12/freq;
    output = real(asin(t12_hat*v/d_3d(1))/pi*180);
end