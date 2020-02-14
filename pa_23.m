%%%output:Psi_23
Psi_23 = zeros(K, K*(L-1));
Interf_23 = zeros(K, 1);
Sign_23 = zeros(K, 1);
for k_23 = 1 : K
    for l_23 = 2 : L
        Interf_23(k_23, 1) = Interf_23(k_23, 1) + D(k_23,(l_23-1)*K+k_23)^4;
    end
    Sign_23(k_23, 1) = D(k_23,k_23)^4;
end
[~, indint23] = sort(Interf_23);%increasing
[~, indsig23] = sort(Sign_23);%increasing
for k_23 = 1 : K
    Psi_23(indsig23(k_23,1), indint23(k_23,1)) = 1;
end
for l0 = 3 : L
    Psi_23(:, K*(l0-2)+1:K*(l0-1)) = Psi_23(:, 1:K);
end