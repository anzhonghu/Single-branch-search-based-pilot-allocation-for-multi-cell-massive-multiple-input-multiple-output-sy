%%%%%%%%%%%%%uplink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j_fz = 1 : L
    if l0 == j_fz%%%j
        for k_zl0 = 1 : K%%%%%k
            finalsumkzl0 = 0;
            for l_final = 1 : L%%%l
                if 1 == l_final
                    X_final = Zl0_temp';
                else
                    X_final = Zl0_temp' * Psi_pro(:, K*(l_final-2)+1:K*(l_final-1));
                end
                if l0 == l_final
                else
                    for k_final = 1 : K%%%k'
                        finalsumkzl0 = finalsumkzl0 + X_final(k_zl0, k_final) * Betaq(k_final, (l0-1)*L+l_final) / Betaq(k_zl0, (l0-1)*L+l0);
                    end
                end
            end
            f_zl0(f_count, 1) = f_zl0(f_count, 1) + gamma(l0, k_zl0)^2 / (1+gamma(l0, k_zl0)) * finalsumkzl0;
        end
    else
        if 1 == j_fz
            X_final = Zl0_temp;
        else
            X_final = Psi_pro(:, K*(j_fz-2)+1:K*(j_fz-1))' * Zl0_temp;
        end
        for k_zl0 = 1 : K%%%%%%%%%%%%%%5k
            finalsumkzl0 = 0;
            for k_final = 1 : K%%%k'
                finalsumkzl0 = finalsumkzl0 + X_final(k_zl0, k_final) * Betaq(k_final, (j_fz-1)*L+l0) / Betaq(k_zl0, (j_fz-1)*L+j_fz);
            end
            f_zl0(f_count, 1) = f_zl0(f_count, 1) + gamma(j_fz, k_zl0)^2 / (1+gamma(j_fz, k_zl0)) * finalsumkzl0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%downlink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j_fz = 1 : L
    if l0 == j_fz%%%j
        for k_zl0 = 1 : K%%%%%k
            finalsumkzl0 = 0;
            for l_final = 1 : L%%%l
                if 1 == l_final
                    X_final = Zl0_temp';
                else
                    X_final = Zl0_temp' * Psi_pro(:, K*(l_final-2)+1:K*(l_final-1));
                end
                if l0 == l_final
                else
                    for k_final = 1 : K%%%k'
                        finalsumkzl0 = finalsumkzl0 + X_final(k_zl0, k_final) * Betadown(k_zl0, (l_final-1)*L+l0) / Betadown(k_final, (l_final-1)*L+l_final);
                    end
                end
            end
            f_zl0(f_count, 1) = f_zl0(f_count, 1) + gamma_down(l0, k_zl0)^2 / (1+gamma_down(l0, k_zl0))  * finalsumkzl0;
        end
    else
        if 1 == j_fz
            X_final = Zl0_temp;
        else
            X_final = Psi_pro(:, K*(j_fz-2)+1:K*(j_fz-1))' * Zl0_temp;
        end
        for k_zl0 = 1 : K%%%%%%%%%%%%%%5k
            finalsumkzl0 = 0;
            for k_final = 1 : K%%%k'
                finalsumkzl0 = finalsumkzl0 + X_final(k_zl0, k_final) * Betadown(k_zl0, (l0-1)*L+j_fz) / Betadown(k_final, (l0-1)*L+l0);
            end
            f_zl0(f_count, 1) = f_zl0(f_count, 1) + gamma_down(j_fz, k_zl0)^2 / (1+gamma_down(j_fz, k_zl0)) * finalsumkzl0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_zl0(f_count, 1) = -f_zl0(f_count, 1);