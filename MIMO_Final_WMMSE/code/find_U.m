function U = find_U(H,V,sigma2,R,I,K,d)
    U = cell(I, K); % I * K * (R * d)
    for i = 1:I
        for k = 1:K
            J = sigma2 * eye(R);
            for j = 1:K
                for m = 1:I
                    H_j = H{i, k, j};
                    V_mj = V{m, j};
                    J = J + H_j * V_mj * V_mj' * H_j';
                end
            end
            H_ik = H{i, k, k};
            V_ik = V{i, k};
            U{i, k} = (J \ H_ik) * V_ik; % size: R × d
        end
    end
end