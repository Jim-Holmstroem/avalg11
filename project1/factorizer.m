function [] = factorizer(filename,start,goal)
    file = fopen(strcat('testdata/',filename),'w')

    for it=start:goal,
        number = sym(int2str(it))
        factors = factor(number)
        fprintf(file,'%s = %s\n',char(number),char(factors))
    end


    

