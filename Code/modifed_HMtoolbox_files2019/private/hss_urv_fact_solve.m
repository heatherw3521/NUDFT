function x = hss_urv_fact_solve(F, b)
% HSS_URV_FACT_SOLVE computes the least-squares solution of the linear system 
%              A x = b where A is a HSS matrix and b a block column
%

F = apply_unitaries(F, b);
x = urv_solve_recursion(F, [], []);

end

function F = apply_unitaries(F, b)
    if ~F.leafnode
        F.left  = apply_unitaries(F.left,  b(1:F.left.m,:));
        F.right = apply_unitaries(F.right, b(1+F.left.m:end,:));
        c = [F.left.c2;F.right.c2];
    else
        c = b;
    end
    if ~F.topnode && ~isempty(F.Om)
        c = F.Om' * c;
        c = c(1:size(F.Q,1),:);
    end
    c = F.Q' * c;
    F.c1 = c(1:size(F.D11,1),:);
    F.c2 = c(size(F.D11,1)+1:end,:);
end

function x = urv_solve_recursion(F, z, y2)

if ~F.topnode
    y1 = F.D11 \ (F.c1 - F.D12 * y2 - F.U1 * z);
    x = F.P*[y1;y2];
else
    x = F.D11 \ F.c1;
end


if ~F.leafnode
    yl2 = x(1:size(F.left.D12,2), 1:size(x,2));
    yr2 = x(size(F.left.D12,2)+1:end, 1:size(x,2));
    zl = F.B12 * (F.right.V' * yr2);
    zr = F.B21 * (F.left.V' * yl2);
    if ~F.topnode
        zl = zl + F.Rl * z;
        zr = zr + F.Rr * z;
    end
    x = [urv_solve_recursion(F.left, zl, yl2);
         urv_solve_recursion(F.right, zr, yr2)];
end

end