function [number factors time] = factor_test_finding(n_bits) 
    ratio = log(2)/log(10);

    n_digits = ceil(ratio*n_bits)

    number = int2str(randi(9));

    for i=2:n_digits,
        number = strcat(number,int2str(randi(10)-1));
    end
    
    number = sym(number)

    tic;
    factors = factor(number)
    time = toc;
 
