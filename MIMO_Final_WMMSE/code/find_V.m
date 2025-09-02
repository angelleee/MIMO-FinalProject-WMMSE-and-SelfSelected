function V = find_V(alpha1,H,U,W,T,I,K,P)
    d = size(W{1,1},1);
    V = cell(I, K);
    
    tol = 1e-6;
    max_iter = 50;
    
    for k = 1:K
        mu_low = 0;
        mu_high = 1e3;
        V_tmp = cell(I,1);
        
        for iter = 1:max_iter % find mu_k*
            mu_k = (mu_low + mu_high) / 2;
            
            power_k = 0;
            for i = 1:I
                A = mu_k * eye(T);
                for j = 1:K
                    for l = 1:I
                        H_ljk = H{l, j, k};
                        U_lj = U{l, j};
                        W_lj = W{l, j};
                        alpha_lj = alpha1(l,j);
                        A = A + alpha_lj * (H_ljk') * U_lj * W_lj * (U_lj') * H_ljk;
                    end
                end
                H_ikk = H{i, k, k};
                U_ik = U{i, k};
                W_ik = W{i, k};
                alpha_ik = alpha1(i,k);
                B = alpha_ik * (H_ikk') * U_ik * W_ik;
                
                V_tmp{i} = A \ B;
                power_k = power_k + norm(V_tmp{i}, 'fro')^2;
            end
            
            if abs(power_k - P) < tol
                break;
            elseif power_k > P
                mu_low = mu_k;  % power too large-->higher mu
            else
                mu_high = mu_k; % power too small-->lower mu
            end
        end
        
        for i = 1:I
            V{i,k} = V_tmp{i};
        end
    end
end
