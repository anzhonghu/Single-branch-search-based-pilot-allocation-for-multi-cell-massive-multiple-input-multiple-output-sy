%Semiblind channel estimation for rapid fading multicell MIMO systems
%This file is for the simulation calculation
clear;
close all;
%%constants
L = 7;%cell number
K = 10;%user number2,3
tao = K;%pilot length
eta_updown = 3 / 7;
r = 1600;%center to edge distance(m)
rc = r * 0.8;
rh = 100;%minimum terminal radius of the cell(m)
ra = rc / rh - 1;
gamma_decay = 3.8;%decay exponent
height = 32;
i_ant = 196;
sigma = 8;%in dB
Num = 1e3;%iteration 2e2
ant_s = [64; 100; 144; 196; 400];%rho_up, rho_down, in dB
SNR = 20;
amp = 10 ^ (SNR*0.05);
%%position of every base
base(1:7,1) = [0;(1i * 2 * rc);(sqrt(3) * rc + 1i * rc);(sqrt(3) * rc - 1i * rc);(-1i * 2 * rc);(-sqrt(3) * rc - 1i * rc);(-sqrt(3) * rc + 1i * rc);];
D = zeros(K,K*L*L);
Dq = zeros(K,K*L*L);
Pl_down = zeros(K,K*L);
pos = zeros(K, L);
phi_UT = zeros(K, 1);
x_temp=zeros(K,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate_store = zeros(length(ant_s), 4);
rate_storett = zeros(length(ant_s), 4);
Psi_ran = zeros(K, K*(L-1));
Psi_ran_s = zeros(K, K*(L-1));
l1 = 1;%BS
for l2 = 2 : L%user
    p_x_store = zeros(K, 1);
    for k = 1 : K
        flag = 0;
        if 1 == k
            p_x = ceil(rand * K);
            Psi_ran_s(k, (l2-2)*K+p_x) = 1;
        else
            while flag < k-1
                p_x = ceil(rand * K);
                for kk = 1 : k-1
                    if p_x == p_x_store(kk, 1)
                        break;
                    else
                        flag = flag + 1;
                    end
                end
                if flag < k-1
                    flag = 0;
                else
                end
            end
            Psi_ran_s(k, (l2-2)*K+p_x) = 1;
        end
        p_x_store(k, 1) = p_x;
    end
end
%Here the iteration begins
for nn = 1 : length(ant_s)
    sum_rate_t = zeros(1, 4);
    sum_rate_td = zeros(1, 4);
    i_ant = ant_s(nn, 1);
    H = zeros(i_ant,K*L*L);
    G = zeros(i_ant,K*L*L);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1 : Num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%position of every terminal, unifrom distribute
        dis(1:K,1:L) = (rem(rand(K,L) * ra, ra) + 1) * rh;
        ang(1:K,1:L) = rand(K,L) * 2 * pi;
        pos(1:K,1:L) = dis .* (exp(1i * ang));
        for ll = 1 : L - 1
            pos(:,ll+1) = pos(:,ll+1) + base(ll+1,1);
        end
        shadow_amp = sqrt(10.^(randn(1,K*L) * sigma * 0.1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %channel matrix, large scale fading and small scale fading
        dl_sum = zeros(L, 1);
        for l1 = 1 : L%BS
            for l2 = 1 : L%user
                H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
                x = ((abs(pos(:,l2)-base(l1,1))).^2 + height^2).^(0.5);
                D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  diag(((x*0.01).^(-0.5*gamma_decay))) * diag(sqrt(shadow_amp(:,(l2-1)*K+1:l2*K)));
            end
        end
        for l1 = 1 : L%BS
            for l2 = 1 : L%user
                Dq(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) =  D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) / D(:,(l2-1)*L*K+(l2-1)*K+1:(l2-1)*L*K+l2*K);
            end
        end
        amp_use = amp;%sqrt(rho),or sqrt(rho_d)
        %%%%%%%%%%%%%%%%%%%%
        Betaq = zeros(K, L*L);
        Betadown = zeros(K, L*L);
        for lll1 = 1 : L
            for lll2 = 1 : L
                for kkk1 = 1 : K
                    Betaq(kkk1, (lll1-1)*L+lll2) = Dq(kkk1,(lll1-1)*L*K+(lll2-1)*K+kkk1)^4;
                    Betadown(kkk1, (lll1-1)*L+lll2) = D(kkk1,(lll1-1)*L*K+(lll2-1)*K+kkk1)^4;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%
        pa_23;%%%output:Psi_23
        prop_pa;%%%output:Psi_pro
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate error rate
        for n_count = 1 : 3
            Dja_t_s = zeros(K, K*L);
            Hjjhat_t_s = zeros(i_ant,K*L);
            switch n_count
                case 1
                    Psi_ran = Psi_pro;
                case 2
                    Psi_ran = Psi_23;
                case 3
                    Psi_ran = Psi_ran_s;
                otherwise
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%uplink%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j = 1 : L
                Yqj_t = H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
                for ll = 1 : L
                    if ll == j
                    else
                        D_temp = zeros(K, K);
                        if 1 == j
                            D_n = Dq(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K) * Psi_ran(:, (ll-2)*K+1 : (ll-1)*K)';
                        else
                            if 1 == ll
                                D_n = Dq(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K) * Psi_ran(:, (j-2)*K+1 : (j-1)*K);
                            else
                                D_n = Dq(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K) * Psi_ran(:, (ll-2)*K+1 : (ll-1)*K)' * Psi_ran(:, (j-2)*K+1 : (j-1)*K);
                            end
                        end
                        for kk = 1 : K
                            D_temp(:, kk) = D_n(:,kk) / Dq(kk,(j-1)*L*K+(j-1)*K+kk);
                        end
                        Yqj_t = Yqj_t + H(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K) * D_temp;
                    end
                end
                if 1 == j
                    for kk = 1 : K
                        Yqj_t(:, kk) = Yqj_t(:, kk) + (randn(i_ant,1)+1i*randn(i_ant,1)) / (amp_use * Dq(kk,(j-1)*L*K+(j-1)*K+kk) *  sqrt(2));
                    end
                else
                    for kk = 1 : K
                        Yqj_t(:, kk) = Yqj_t(:, kk) + (randn(i_ant,K)+1i*randn(i_ant,K)) * Psi_ran(:, (j-2)*K+kk) / (amp_use * Dq(kk,(j-1)*L*K+(j-1)*K+kk) *  sqrt(2));
                    end
                end
                Dja_t = zeros(K, K);
                for kk = 1 : K
                    ajk_t = 0;
                    for ll = 1 : L
                        if ll == j
                        else
                            if 1 == j
                                XD_t = Dq(:,(j-1)*L*K+(ll-1)*K+1 : (j-1)*L*K+ll*K)* Psi_ran(:, (ll-2)*K+1 : (ll-1)*K)';
                            else
                                if 1 == ll
                                    XD_t = Dq(:,(j-1)*L*K+(ll-1)*K+1 : (j-1)*L*K+ll*K) * Psi_ran(:, (j-2)*K+1 : (j-1)*K);
                                else
                                    XD_t = Dq(:,(j-1)*L*K+(ll-1)*K+1 : (j-1)*L*K+ll*K)* Psi_ran(:, (ll-2)*K+1 : (ll-1)*K)' * Psi_ran(:, (j-2)*K+1 : (j-1)*K);
                                end
                            end
                            ajk_t = ajk_t + norm(XD_t(:, kk))^2;
                        end
                    end
                    ajk_t = ajk_t + 1 / amp_use^2;
                    ajk_t = ajk_t / Dq(kk,(j-1)*L*K+(j-1)*K+kk)^2 + 1;
                    ajk_t = 1 / ajk_t;
                    Dja_t(kk, kk) = ajk_t;
                end
                Dja_t_s(:, (j-1)*K+1 : j*K) = Dja_t;
                Hjjhat_t = Yqj_t * Dja_t;
                Hjjhat_t_s(:, (j-1)*K+1 : j*K) = Hjjhat_t;
                alphaj_t = 0;
                for ll = 1 : L
                    if ll == j
                        for  kk = 1 : K
                            alphaj_t = alphaj_t + (1 - Dja_t(kk, kk)) * Dq(kk,(j-1)*L*K+(j-1)*K+kk)^2;
                        end
                    else
                        alphaj_t = alphaj_t +  trace(Dq(:,(j-1)*L*K+(ll-1)*K+1 : (j-1)*L*K+ll*K)^2);
                    end
                end
                alphaj_t = alphaj_t + 1;
                for k = 1 : K
                    Ajk_t = Hjjhat_t * Dq(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)^2 * Hjjhat_t' + alphaj_t * eye(i_ant);
                    fjk_t = 1 / amp_use * Dq(k,(j-1)*L*K+(j-1)*K+k) * (Hjjhat_t(:, k))' / Ajk_t;
                    xxx_t = (Hjjhat_t' * Hjjhat_t * Dq(:,(j-1)*L*K+(j-1)*K+1 : (j-1)*L*K+j*K)^2 + alphaj_t * eye(K))^(-1) * Hjjhat_t' * H(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
                    Sjk_t =  Dq(k,(j-1)*L*K+(j-1)*K+k)^4 * (abs(xxx_t(k,k)))^2;
                    Injk_t = 0;
                    Injk1_t = 0;
                    for ll = 1 : L
                        if ll == j
                            for  kk = 1 : K
                                if kk == k
                                else
                                    Injk_t = Injk_t + amp_use^2 * (abs(fjk_t * H(:,(j-1)*L*K+(j-1)*K+kk)))^2 * Dq(kk,(j-1)*L*K+(j-1)*K+kk)^2;
                                end
                            end
                        else
                            xxx_t = (Hjjhat_t' * Hjjhat_t * Dq(:,(j-1)*L*K+(j-1)*K+1 : (j-1)*L*K+j*K)^2 + alphaj_t * eye(K))^(-1) * Hjjhat_t' * H(:,(j-1)*L*K+(ll-1)*K+1:(j-1)*L*K+ll*K)...
                                * Dq(:,(j-1)*L*K+(ll-1)*K+1 : (j-1)*L*K+ll*K);
                            Injk_t = Injk_t + Dq(k,(j-1)*L*K+(j-1)*K+k)^2 * (norm(xxx_t(k,:)))^2;
                        end
                    end
                    Injk_t = Injk_t + (norm(fjk_t))^2;
                    sinr_t = Sjk_t / Injk_t;
                    sum_rate_t(1, n_count) = sum_rate_t(1, n_count) + log2(1 + sinr_t);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%downlink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            G = zeros(i_ant, K*L);
            for j = 1 : L
                delta_pow = 0;
                for k = 1 : K
                    delta_pow = delta_pow + D(k,(j-1)*L*K+(j-1)*K+k)^(-2);
                end
                amp_use_d = amp_use * sqrt(delta_pow) / sqrt(K);
                Dja_t = Dja_t_s(:, (j-1)*K+1 : j*K);
                deltaj = K / amp_use_d^2 + trace(Dja_t^2) / amp_use^2 + trace((Dja_t-eye(K))^2 * D(:,(j-1)*L*K+(j-1)*K+1 : (j-1)*L*K+j*K)^2);
                for l = 1 : L
                    if l == j
                    else
                        if 1 == j
                            Psi_jl_use =  Psi_ran(:, (l-2)*K+1 : (l-1)*K);
                        else
                            if 1 == l
                                Psi_jl_use =  Psi_ran(:, (j-2)*K+1 : (j-1)*K)';
                            else
                                Psi_jl_use = Psi_ran(:, (j-2)*K+1 : (j-1)*K)' * Psi_ran(:, (l-2)*K+1 : (l-1)*K);
                            end
                        end
                        deltaj = deltaj + trace(Dja_t^2 * Psi_jl_use * D(:,(j-1)*L*K+(l-1)*K+1 : (j-1)*L*K+l*K)^2 * Psi_jl_use');
                    end
                end
                Hjjhat_t = Hjjhat_t_s(:, (j-1)*K+1 : j*K);
                Fj = (Hjjhat_t * D(:,(j-1)*L*K+(j-1)*K+1 : (j-1)*L*K+j*K)^2 * Hjjhat_t' + deltaj * eye(i_ant)) \ Hjjhat_t * D(:,(j-1)*L*K+(j-1)*K+1 : (j-1)*L*K+j*K);
                ksij = amp_use_d / sqrt(trace(Fj * Fj'));
                G(:, (j-1)*K+1 : j*K) = ksij * Fj;
            end
            for j = 1 : L
                for k = 1 : K
                    Sjk_t = D(k,(j-1)*L*K+(j-1)*K+k)^2 * (abs((H(:,(j-1)*L*K+(j-1)*K+k))' * G(:, (j-1)*K+k)))^2;
                    Injk_t = 1;
                    for kk = 1 : K
                        if kk == k
                        else
                            Injk_t = Injk_t + D(k,(j-1)*L*K+(j-1)*K+k)^2 * (abs((H(:,(j-1)*L*K+(j-1)*K+k))' * G(:, (j-1)*K+kk)))^2;
                        end
                    end
                    for ll = 1 : L
                        if ll == j
                        else
                            Injk_t = Injk_t + (norm(D(k,(ll-1)*L*K+(j-1)*K+k) * (H(:,(ll-1)*L*K+(j-1)*K+k))' * G(:, (ll-1)*K+1:ll*K)))^2;
                        end
                    end
                    sinr_t = Sjk_t / Injk_t;
                    sum_rate_td(1, n_count) = sum_rate_td(1, n_count) + log2(1 + sinr_t);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sprintf('%d,%d',nn,jj)
    end
    sum_rate_t = sum_rate_t / Num;
    sum_rate_td = sum_rate_td / Num;
    rate_store(nn, :) = (sum_rate_t + sum_rate_td) * eta_updown;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
set(h1,'PaperType','A4');
xx = axes('FontSize',16);
plot(ant_s, rate_store(:,3),'b:o','LineWidth',2,'MarkerSize',10)
hold on
plot(ant_s, rate_store(:,2),'m:x','LineWidth',2,'MarkerSize',10)
plot(ant_s, rate_store(:,1),'k-^','LineWidth',2,'MarkerSize',10)
xlim([min(ant_s), max(ant_s)])
grid on
le = legend('Random pilot allocation', 'Pilot allocation [19]','Proposed pilot allocation', 'Location','Southeast');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',ant_s)
xlabel('Number of the BS antennas','Fontsize',20,'Fontname','Times')
ylabel('Spectral efficiency (bps/Hz)','Fontsize',20,'Fontname','Times')
print(h1,'-dpdf','BSantenna')


