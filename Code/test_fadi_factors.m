%test if fadi is working correctly: 

Ct = buildcauchy(ridx, cidx, N, nodes, a, b); 
lra = H.U*Ct(J, :);

norm(Ct -lra)