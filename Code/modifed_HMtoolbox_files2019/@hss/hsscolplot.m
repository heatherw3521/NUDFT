function out = hsscolplot(H, level, ranks)
%plots the HSS columns for H at level l (l = 0 is root).
%ranks =1 : ranks of the cols are depicted. ranks = 0 turns this off. 

%%
% start by getting the matrix at the given level: 

For j = 1:level
    [Hll, Hlr] = split(Hl); 
    [Hrl, Hrr] = split(Hr);
    if (j==level)
        



    


