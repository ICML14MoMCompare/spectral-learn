function [ ] = interpret_pautomac(infile, outfile)

[wa,A,n] = load_model(char(infile));

outfilep = fopen(char(outfile),'w');

fprintf(outfilep, '%d %d\n', n,size(A,2));

for source=0:(n-1)
    prob = wa.a1(source+1);
    fprintf(outfilep,'%d %f\n', source, prob);
end

for source=0:(n-1)
    for symbol=0:(size(A,2)-1)
        for target=0:(n-1)
            operator = wa.(A(symbol+1));
            prob = operator(source+1,target+1);
            fprintf(outfilep,'%d %d %d %f\n',source,target,symbol,prob);
        end
    end
end

for source=0:(n-1)
    prob = wa.ainf(source+1);
    fprintf(outfilep,'%d %f\n', source, prob);
end
    


end

