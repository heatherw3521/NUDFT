function H = hss_build_hss_tree_nudft(m, n, nodes, rou, block_size)
%constructs a tree for compressing transformed nudft-type matrices. 
%the indices are split according to node distribution, rou = angle(roots of
%unity), nodes = angle(nodes), where angle ranges between pi/N to 2*pi+2*pi/N. 
% exp(1i j/N), j = 1...N, where N = # cols of nudft matrix. 

    gap = pi/n; %diff(rou(1:2))/2;  
    nu = ceil(m/n); 
    H = build_hss_tree_nudft_rec(m, n, nodes, rou, block_size, gap, nu);
    H.topnode = 1;
    
end

function H = build_hss_tree_nudft_rec(m, n, nodes, rou, block_size, gap, nu)
H = hss();

H.topnode  = 0;

if n > block_size && m > nu*block_size %
    %split rou down the middle: 
    n1 = ceil(n/2);
    n2 = n-n1; 

    %find split that separates rou(1:n1) from nodes(n1+1:end)
    % and rou(n1+1:end) from nodes(1:m1); 
    jj = (1:m)'; 
    jj = jj(nodes <= (rou(n1+1)-gap)); 
    m1 = length(jj); %what to do if there are no nodes? 
    m2 = m-m1; 

    H.ml = m1; H.mr = m2;
    H.nl = n1; H.nr = n2;
    
    H.A11 = build_hss_tree_nudft_rec(m1, n1, nodes(1:m1), rou(1:n1), block_size, gap, nu);
    H.A22 = build_hss_tree_nudft_rec(m2, n2, nodes(m1+1:end), rou(n1+1:end), block_size, gap, nu);
    
    H.leafnode = 0;
else
    H.leafnode = 1;
    H.D = zeros(m, n);
    H.U = zeros(m, 0);
    H.V = zeros(n, 0);
end
end

