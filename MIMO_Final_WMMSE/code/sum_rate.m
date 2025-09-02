function system_rate = sum_rate(H,V,sigma2,R,I,K,alpha1)
    system_rate = 0;
    for i = 1:I
        for k = 1:K
            H_ik = H{i, k, k};
            V_ik = V{i, k};
            d = size(V_ik, 2);
            J = sigma2 * eye(R); % J = sigma (H_ikj * V_lj * V_lj' * H_ikj' + sigma2 * I)
            for j = 1:K
                for m = 1:I
                    if j ~= k || m ~= i
                        H_ij = H{i, k, j};
                        V_mj = V{m, j};
                        J = J + H_ij * V_mj * V_mj' * H_ij';
                    end
                end
            end
            signal = H_ik * V_ik;
            R_ik = log2(real(det(eye(R) + (signal * signal') / J)));
            system_rate = system_rate + alpha1(i, k) * R_ik;
        end
    end
end