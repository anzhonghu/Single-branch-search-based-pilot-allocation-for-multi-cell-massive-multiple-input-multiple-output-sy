flag_k_sel = zeros(K, 1);
flag_k_sel_unbeli = zeros(K, 1);
Zl0 = zeros(K, K);
Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1));
gamma_cal;
for kq = 1 : K
    if 1 == kq
        f_zl0 = zeros(K, 1);
        Zl0_temp = zeros(K, K);
        f_zl0(1, 1) = 0;%%%%%output: f_zl0(1, 1)
        f_count = 1;
        for k_t = 2 : K
            f_count = f_count + 1;
            Zkq = zeros(K, K);
            Zkq(kq, kq) = -1;
            Zkq(kq, k_t) = 1;
            Zl0_temp = Zkq;
            Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0_temp;
            fl0f_count1cal_temp;%%%%%output: f_zl0(f_count, 1)
        end
        [~, f_ind] = sort(f_zl0,'descend');
        if 1 == f_ind(1, 1)
            %Zl0 = Zl0;
        else
            Zkq = zeros(K, K);
            Zkq(kq, kq) = -1;
            if f_ind(1, 1) - 1 < kq
                k_t = f_ind(1, 1) - 1;
            else
                k_t = f_ind(1, 1);
            end
            Zkq(kq, k_t) = 1;
            Zl0 = Zl0 + Zkq;
            flag_k_sel(k_t, 1) = 1;
        end
    else%%%%%%%%%%%%%%1~=kq
        if kq < K
            f_zl0 = zeros(K, 1);
            Zl0_temp = Zl0;
            x_zl0 = sum(Zl0_temp);
            x_count = 0;
            for k_x = 1 : K
                if x_zl0(1, k_x) < 0
                    x_count = x_count + 1;
                else
                end
            end
            if x_count > K - kq
            else
                Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0_temp;
                fl011cal_temp;%%%%%output: f_zl0(1, 1)
            end
            f_count = 1;
            for k_t = 1 : K
                if k_t == kq
                else
                    f_count = f_count + 1;
                    if 1 == flag_k_sel(k_t, 1)
                    else
                        Zkq = zeros(K, K);
                        Zkq(kq, kq) = -1;
                        Zkq(kq, k_t) = 1;
                        Zl0_temp = Zl0;
                        Zl0_temp = Zl0_temp + Zkq;
                        x_zl0 = sum(Zl0_temp);
                        x_count = 0;
                        for k_x = 1 : K
                            if x_zl0(1, k_x) < 0
                                x_count = x_count + 1;
                            else
                            end
                        end
                        if x_count > K - kq
                        else
                            Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0_temp;
                            fl0f_count1cal_temp;%%%%%output: f_zl0(f_count, 1)
                        end
                    end
                end
            end
            [~, f_ind] = sort(f_zl0,'descend');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for kk_unu = 1 : K
                if 1 == f_ind(kk_unu,1)
                    if 1 == flag_k_sel(kq, 1)
                    else
                        Zl0_temp = Zl0;
                        flag_k_sel_unbeli(kq, 1) = 1;
                        break;
                    end
                else
                    if f_ind(kk_unu,1) - 1 < kq
                        k_t = f_ind(kk_unu,1) - 1;
                    else
                        k_t = f_ind(kk_unu,1);
                    end
                    if 1 == flag_k_sel(k_t, 1) || 1 == flag_k_sel_unbeli(k_t, 1)
                    else
                        Zkq = zeros(K, K);
                        Zkq(kq, kq) = -1;
                        Zkq(kq, k_t) = 1;
                        Zl0_temp = Zl0 + Zkq;
                        flag_k_sel(k_t, 1) = 1;
                        break;
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if 1 == f_ind(1,1)
%                 if 1 == flag_k_sel(kq, 1)
%                     for kk_unu = 2 : K
%                         if f_ind(kk_unu,1) - 1 < kq
%                             k_t = f_ind(kk_unu,1) - 1;
%                         else
%                             k_t = f_ind(kk_unu,1);
%                         end
%                         if 1 == flag_k_sel(k_t, 1) || 1 == flag_k_sel_unbeli(k_t, 1)
%                         else
%                             Zkq = zeros(K, K);
%                             Zkq(kq, kq) = -1;
%                             Zkq(kq, k_t) = 1;
%                             Zl0_temp = Zl0 + Zkq;
%                             flag_k_sel(k_t, 1) = 1;
%                             break;
%                         end
%                     end
%                 else
%                     Zl0_temp = Zl0;
%                     flag_k_sel_unbeli(kq, 1) = 1;
%                 end
%             else
%                 if f_ind(1,1) - 1 < kq
%                     k_t = f_ind(1,1) - 1;
%                 else
%                     k_t = f_ind(1,1);
%                 end
%                 if 1 == flag_k_sel(k_t, 1) || 1 == flag_k_sel_unbeli(k_t, 1)
%                     for kk_unu = 2 : K
%                         if 1 == f_ind(kk_unu,1)
%                             if 1 == flag_k_sel(kq, 1)
%                             else
%                                 Zl0_temp = Zl0;
%                                 flag_k_sel_unbeli(kq, 1) = 1;
%                                 break;
%                             end
%                         else
%                             if f_ind(kk_unu,1) - 1 < kq
%                                 k_t = f_ind(kk_unu,1) - 1;
%                             else
%                                 k_t = f_ind(kk_unu,1);
%                             end
%                             if 1 == flag_k_sel(k_t, 1) || 1 == flag_k_sel_unbeli(k_t, 1)
%                             else
%                                 Zkq = zeros(K, K);
%                                 Zkq(kq, kq) = -1;
%                                 Zkq(kq, k_t) = 1;
%                                 Zl0_temp = Zl0 + Zkq;
%                                 flag_k_sel(k_t, 1) = 1;
%                                 break;
%                             end
%                         end
%                     end
%                 else
%                     Zkq = zeros(K, K);
%                     Zkq(kq, kq) = -1;
%                     Zkq(kq, k_t) = 1;
%                     Zl0_temp = Zl0 + Zkq;
%                     flag_k_sel(k_t, 1) = 1;
%                 end
%             end
            Zl0 = Zl0_temp;
        else%%%%%kq = K;
            f_zl0 = zeros(K, 1);
            Zl0_temp = Zl0;
            x_zl0 = sum(Zl0_temp);
            x_count = 0;
            for k_x = 1 : K
                if x_zl0(1, k_x) < 0
                    x_count = x_count + 1;
                else
                end
            end
            if 0 == x_count
                Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0_temp;
                gamma_cal;
                fl011cal_temp;%%%%%output: f_zl0(1, 1)
            else
            end
            f_count = 1;
            for k_t = 1 : K
                if k_t == kq
                else
                    f_count = f_count + 1;
                    Zkq = zeros(K, K);
                    Zkq(kq, kq) = -1;
                    Zkq(kq, k_t) = 1;
                    Zl0_temp = Zl0;
                    Zl0_temp = Zl0_temp + Zkq;
                    x_zl0 = sum(Zl0_temp);
                    x_count = 0;
                    for k_x = 1 : K
                        if x_zl0(1, k_x) < 0
                            x_count = x_count + 1;
                        else
                        end
                    end
                    if 0 == x_count
                        Psi_pro_temp(:, K*(l0-2)+1:K*(l0-1)) = Psi_pro(:, K*(l0-2)+1:K*(l0-1)) + Zl0_temp;
                        gamma_cal;
                        fl0f_count1cal_temp;%%%%%output: f_zl0(f_count, 1)
                    else
                    end
                end
            end
            if max(f_zl0) > 0
                [~, f_ind] = sort(f_zl0,'descend');
                if 1 == f_ind(1,1)
                else
                    Zkq = zeros(K, K);
                    Zkq(kq, kq) = -1;
                    if f_ind(1,1) - 1 < kq
                        k_t = f_ind(1,1) - 1;
                    else
                        k_t = f_ind(1,1);
                    end
                    Zkq(kq, k_t) = 1;
                    Zl0 = Zl0 + Zkq;
                end
            else
                Zl0 = zeros(K, K);
            end
        end
    end
end