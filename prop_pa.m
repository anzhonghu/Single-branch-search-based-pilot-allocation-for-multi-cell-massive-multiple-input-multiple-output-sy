%%%output:Psi_pro
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_pro = zeros(K, K*(L-1));
b = zeros(L-1, K);
j0 = 1;%BS
for l0 = 2 : L
    Psi_pro(:, K*(l0-2)+1:K*(l0-1)) = eye(K);
end
Psi_pro_temp = Psi_pro;
Psi_pro_temp_s = Psi_pro;
for l0 = 2 : L
    branch_sel_t;
    Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0;
end
s_count = 0;
sum_s_stor = zeros(2^6, 1);
for s2 = 1 : 2
    l0 = 2;
    if 1 ~= s2
        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1));
    else
        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
    end
    for s3 = 1 : 2
        l0 = 3;
        if 1 ~= s3
            Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1));
        else
            Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
        end
        for s4 = 1 : 2
            l0 = 4;
            if 1 ~= s4
                Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1));
            else
                Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
            end
            for s5 = 1 : 2
                l0 = 5;
                if 1 ~= s5
                    Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1));
                else
                    Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
                end
                for s6 = 1 : 2
                    l0 = 6;
                    if 1 ~= s6
                        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1));
                    else
                        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
                    end
                    for s7 = 1 : 2
                        l0 = 7;
                        if 1 ~= s7
                            Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro_temp_s(:, K*(l0-2)+1:K*(l0-1));
                        else
                            Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
                        end
                        s_count = s_count + 1;
                        for jj_opt = 1 : L%j
                            for kk_opt = 1 : K%k
                                %%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%uplink%%%%%%%%%%%%
                                xxx_opt = 0;
                                for jjj_opt = 1 : L%l
                                    if jj_opt == jjj_opt
                                    else
                                        if 1 == jj_opt
                                            X_opt = Psi_pro_temp(:, K*(jjj_opt-2)+1 : K*(jjj_opt-2) + K);
                                            for kkk_opt = 1 : K
                                                xxx_opt = xxx_opt + X_opt(kk_opt, kkk_opt) * Betaq(kkk_opt, (jj_opt-1)*L+jjj_opt);
                                            end
                                        else
                                            if 1 == jjj_opt
                                                X_opt = Psi_pro_temp(:, K*(jj_opt-2)+1 : K*(jj_opt-2) + K)';
                                                for kkk_opt = 1 : K
                                                    xxx_opt = xxx_opt + X_opt(kk_opt, kkk_opt) * Betaq(kkk_opt, (jj_opt-1)*L+jjj_opt);
                                                end
                                            else
                                                X_opt = Psi_pro_temp(:, K*(jj_opt-2)+1 : K*(jj_opt-2) + K)' * Psi_pro_temp(:, K*(jjj_opt-2)+1 : K*(jjj_opt-2) + K);
                                                for kkk_opt = 1 : K%k'
                                                    xxx_opt = xxx_opt + X_opt(kk_opt, kkk_opt) * Betaq(kkk_opt, (jj_opt-1)*L+jjj_opt);
                                                end
                                            end
                                        end
                                    end
                                end
                                xxx_opt = Betaq(kk_opt, (jj_opt-1)*L+jj_opt) / xxx_opt;
                                sum_s_stor(s_count, 1) = sum_s_stor(s_count, 1) + log2(1+xxx_opt);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%downlink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                xxx_opt = 0;
                                for jjj_opt = 1 : L%l
                                    if jj_opt == jjj_opt
                                    else
                                        if 1 == jj_opt
                                            X_opt = Psi_pro_temp(:, K*(jjj_opt-2)+1 : K*(jjj_opt-2) + K);
                                            for kkk_opt = 1 : K
                                                xxx_opt = xxx_opt + X_opt(kk_opt, kkk_opt) * Betaq(kk_opt, (jjj_opt-1)*L+jj_opt) / Betaq(kkk_opt, (jjj_opt-1)*L+jjj_opt);
                                            end
                                        else
                                            if 1 == jjj_opt
                                                X_opt = Psi_pro_temp(:, K*(jj_opt-2)+1 : K*(jj_opt-2) + K)';
                                                for kkk_opt = 1 : K
                                                    xxx_opt = xxx_opt + X_opt(kk_opt, kkk_opt) * Betaq(kk_opt, (jjj_opt-1)*L+jj_opt) / Betaq(kkk_opt, (jjj_opt-1)*L+jjj_opt);
                                                end
                                            else
                                                X_opt = Psi_pro_temp(:, K*(jj_opt-2)+1 : K*(jj_opt-2) + K)' * Psi_pro_temp(:, K*(jjj_opt-2)+1 : K*(jjj_opt-2) + K);
                                                for kkk_opt = 1 : K%k'
                                                    xxx_opt = xxx_opt + X_opt(kk_opt, kkk_opt) * Betaq(kk_opt, (jjj_opt-1)*L+jj_opt) / Betaq(kkk_opt, (jjj_opt-1)*L+jjj_opt);
                                                end
                                            end
                                        end
                                    end
                                end
                                xxx_opt = 1 / xxx_opt;
                                sum_s_stor(s_count, 1) = sum_s_stor(s_count, 1) + log2(1+xxx_opt);
                            end
                        end
                    end
                end
            end
        end
    end
end
[~, indopt] = max(sum_s_stor);
n_sopt = zeros(L-2, 1);
indopt = indopt-1;
for n_opt = L-1 : -1 : 1
    if indopt / 2^(n_opt-1) >= 1
        n_sopt(n_opt, 1) = 1;
        indopt = indopt - 2^(n_opt-1);
    else
        n_sopt(n_opt, 1) = 0;
    end
end
for l0 = 2 : L
    if 0 ~= n_sopt(-l0+8, 1)
        branch_sel_t;
        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0;
    else
        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
    end
end
Psi_pro = Psi_pro_temp;