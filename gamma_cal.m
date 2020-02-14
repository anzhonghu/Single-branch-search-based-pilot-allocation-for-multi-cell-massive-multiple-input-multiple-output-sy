%%%%%%%%%%%%%%%%%%%uplink%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = zeros(L, K);
for l_gamma = 1 : L%%%j
    for k_gamma = 1 : K%%%k
        for lgama = 1 : L%%%%l
            if lgama == l_gamma
            else
                if 1 == lgama
                    X_forgammma = Psi_pro_temp(:, K*(l_gamma-2)+1:K*(l_gamma-1))';
                    for kgama = 1 : K%%%%k'
                        gamma(l_gamma, k_gamma) = gamma(l_gamma, k_gamma) + X_forgammma(k_gamma, kgama) * Betaq(kgama, (l_gamma-1)*L+lgama) / Betaq(k_gamma, (l_gamma-1)*L+l_gamma);
                    end
                else
                    if 1 == l_gamma
                        X_forgammma = Psi_pro_temp(:, K*(lgama-2)+1:K*(lgama-1));
                        for kgama = 1 : K%%%%k'
                            gamma(l_gamma, k_gamma) = gamma(l_gamma, k_gamma) + X_forgammma(k_gamma, kgama) * Betaq(kgama, (l_gamma-1)*L+lgama) / Betaq(k_gamma, (l_gamma-1)*L+l_gamma);
                        end
                    else
                        X_forgammma = Psi_pro_temp(:, K*(l_gamma-2)+1:K*(l_gamma-1))' * Psi_pro_temp(:, K*(lgama-2)+1:K*(lgama-1));
                        for kgama = 1 : K%%%%k'
                            gamma(l_gamma, k_gamma) = gamma(l_gamma, k_gamma) + X_forgammma(k_gamma, kgama) * Betaq(kgama, (l_gamma-1)*L+lgama) / Betaq(k_gamma, (l_gamma-1)*L+l_gamma);
                        end
                    end
                end
            end
        end
        if gamma(l_gamma, k_gamma) < 1e-10
            gamma(l_gamma, k_gamma) = 1e10;
        else
            gamma(l_gamma, k_gamma) = 1 / gamma(l_gamma, k_gamma);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%downlink%%%%%%%%%%%%%%
gamma_down = zeros(L, K);
for l_gamma = 1 : L%%%j
    for k_gamma = 1 : K%%%k
        for lgama = 1 : L%%%%l
            if lgama == l_gamma
            else
                if 1 == lgama
                    X_forgammma = Psi_pro_temp(:, K*(l_gamma-2)+1:K*(l_gamma-1))';
                    for kgama = 1 : K%%%%k'
                        gamma_down(l_gamma, k_gamma) = gamma_down(l_gamma, k_gamma) + X_forgammma(k_gamma, kgama) * Betadown(k_gamma, (lgama-1)*L+l_gamma) / Betadown(kgama, (lgama-1)*L+lgama);
                    end
                else
                    if 1 == l_gamma
                        X_forgammma = Psi_pro_temp(:, K*(lgama-2)+1:K*(lgama-1));
                        for kgama = 1 : K%%%%k'
                            gamma_down(l_gamma, k_gamma) = gamma_down(l_gamma, k_gamma) + X_forgammma(k_gamma, kgama) * Betadown(k_gamma, (lgama-1)*L+l_gamma) / Betadown(kgama, (lgama-1)*L+lgama);
                        end
                    else
                        X_forgammma = Psi_pro_temp(:, K*(l_gamma-2)+1:K*(l_gamma-1))' * Psi_pro_temp(:, K*(lgama-2)+1:K*(lgama-1));
                        for kgama = 1 : K%%%%k'
                            gamma_down(l_gamma, k_gamma) = gamma_down(l_gamma, k_gamma) + X_forgammma(k_gamma, kgama) * Betadown(k_gamma, (lgama-1)*L+l_gamma) / Betadown(kgama, (lgama-1)*L+lgama);
                        end
                    end
                end
            end
        end
        if gamma_down(l_gamma, k_gamma) < 1e-10
            gamma_down(l_gamma, k_gamma) = 1e10;
        else
            gamma_down(l_gamma, k_gamma) = 1 / gamma_down(l_gamma, k_gamma);
        end
    end
end