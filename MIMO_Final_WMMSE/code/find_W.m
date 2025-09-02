function W = find_W(U,H,V,I,K,d)
    W = cell(I, K);
    for i = 1:I
        for k = 1:K
            U_ik = U{i, k};
            H_ik = H{i, k, k};
            V_ik = V{i, k};
            
            E = eye(d) - U_ik' * H_ik * V_ik;

            W{i, k} = inv(E);
        end
    end
end