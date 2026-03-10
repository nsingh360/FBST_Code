function [X1,X2,X3,X4,b_11,b_12,b_21,b_22,l1,r1,l2,r2] = eval_canonical_vecs(T)
%% Evaluate the canonical fundamental system required to perform fast inversion of Toeplitz system
sz = size(T);
n = sz(1);
c = T(:,1);
r = T(1,2:end);
ak = [flip(c);r.']; % coefficient vector for the Toeplitz matrix

power = 0;
smallest_power_of_2 = 2^power;  
% Loop until the smallest power of 2 greater than the input is found
while smallest_power_of_2 < n
    power = power + 1;
    smallest_power_of_2 = 2^power;
end
N = smallest_power_of_2;

k = linspace(0,2*N-1,2*N);
Sk = exp(1i*2*pi*k/(2*N)); % 2N interpolation nodes
Sk_bad = []; % list for bad points

B_prev = eye(2);
leftpadz = @(p1,p2) [zeros(1,max(0,numel(p2) - numel(p1)));p1];

%% coefficients of the initial polynomial matrix
b_11 = [1];
b_12 = [0];
b_21 = [0];
b_22 = [1];

l1 = [];
l2 = [];
r1 = [];
r2 = [];

for i=1:N

    % store old_coeff
    b_11_old = b_11;
    b_12_old = b_12;
    b_21_old = b_21;
    b_22_old = b_22;

    % evaluate l*,s* at 2*(i-1)+1 node
    % s_1 =  Sk(2*(i-1)+1);
    % lr_1 = [s_1^n -eval_fk(ak,s_1,n)]*[polyval(b_11_old,s_1) polyval(b_12_old,s_1); polyval(b_21_old,s_1) polyval(b_22_old,s_1)];
    % l1(i) = lr_1(1);
    % r1(i) = lr_1(2);
    % 
    % % evaluate ;**,s** at 2*(i-1)+2 node
    % s_2 = Sk(2*(i-1)+2);
    % lr_2 = [s_2^n -eval_fk(ak,s_2,n)]*[polyval(b_11_old,s_2) polyval(b_12_old,s_2); polyval(b_21_old,s_2) polyval(b_22_old,s_2)];
    % l2(i) = lr_2(1);
    % r2(i) = lr_2(2);
    [lr_1, lr_2, lambda_1, lambda_2, s_1, s_2] = select_interpolation_points(b_11_old, b_12_old, b_21_old, b_22_old, Sk, ak, n);
    Sk = Sk(Sk~=s_1);
    Sk = Sk(Sk~=s_2);
    l1(i) = lr_1(1);
    r1(i) = lr_1(2);
    l2(i) = lr_2(1);
    r2(i) = lr_2(2);

    temp_mat = [lr_1(1) 0 lr_1(2) 0;0 lr_1(1) 0 lr_1(2); lr_2(1) 0 lr_2(2) 0;0 lr_2(1) 0 lr_2(2)];
    temp_b = -[s_1*lr_1.';s_2*lr_2.'];

    % coeff_Bij = temp_mat \ temp_b;

    % lambda_1 = [lr_1(1) lr_1(1)*s_1 lr_1(2) lr_1(2)*s_1]';
    % lambda_2 = [lr_2(1) lr_2(1)*s_2 lr_2(2) lr_2(2)*s_2]';

    [Q,R] = qr([lambda_1 lambda_2]);

    % check for bad points
    qq = [Q(2,3) Q(2,4); Q(4,3) Q(4,4)];
    if rcond(qq) < 1e-10
        Sk_bad(end+1) = s_1;
        Sk_bad(end+1) = s_2;
    else
        Coeff_Bij_unitary =  [Q(1,3) Q(1,4); Q(3,3) Q(3,4)]/qq;
        Coeff_Bij_unitary = Coeff_Bij_unitary.';
        coeff_Bij = Coeff_Bij_unitary(:);

        % multiply coefficients
        b_11 = leftpadz(conv(b_11_old,[1;coeff_Bij(1)]), conv(b_12_old,[coeff_Bij(3)])) + leftpadz(conv(b_12_old,[coeff_Bij(3)]),conv(b_11_old,[1;coeff_Bij(1)]));
        b_12 = leftpadz(conv(b_11_old,[coeff_Bij(2)]), conv(b_12_old,[1;coeff_Bij(4)])) + leftpadz(conv(b_12_old,[1;coeff_Bij(4)]),conv(b_11_old,[coeff_Bij(2)]));
        b_21 = leftpadz(conv(b_21_old,[1;coeff_Bij(1)]), conv(b_22_old,[coeff_Bij(3)])) + leftpadz(conv(b_22_old,[coeff_Bij(3)]),conv(b_21_old,[1;coeff_Bij(1)]));
        b_22 = leftpadz(conv(b_21_old,[coeff_Bij(2)]), conv(b_22_old,[1;coeff_Bij(4)])) + leftpadz(conv(b_22_old,[1;coeff_Bij(4)]),conv(b_21_old,[coeff_Bij(2)]));
    end

end

if ~isempty(Sk_bad)
    for i=1:(length(Sk_bad)/2)
    
        % store old_coeff
        b_11_old = b_11;
        b_12_old = b_12;
        b_21_old = b_21;
        b_22_old = b_22;
    
        % evaluate l*,s* at 2*(i-1)+1 node
        % s_1 =  Sk(2*(i-1)+1);
        % lr_1 = [s_1^n -eval_fk(ak,s_1,n)]*[polyval(b_11_old,s_1) polyval(b_12_old,s_1); polyval(b_21_old,s_1) polyval(b_22_old,s_1)];
        % l1(i) = lr_1(1);
        % r1(i) = lr_1(2);
        % 
        % % evaluate ;**,s** at 2*(i-1)+2 node
        % s_2 = Sk(2*(i-1)+2);
        % lr_2 = [s_2^n -eval_fk(ak,s_2,n)]*[polyval(b_11_old,s_2) polyval(b_12_old,s_2); polyval(b_21_old,s_2) polyval(b_22_old,s_2)];
        % l2(i) = lr_2(1);
        % r2(i) = lr_2(2);
        [lr_1, lr_2, lambda_1, lambda_2, s_1, s_2] = select_interpolation_points(b_11_old, b_12_old, b_21_old, b_22_old, Sk_bad, ak, n);
        Sk_bad = Sk_bad(Sk_bad~=s_1);
        Sk_bad = Sk_bad(Sk_bad~=s_2);
        l1(i) = lr_1(1);
        r1(i) = lr_1(2);
        l2(i) = lr_2(1);
        r2(i) = lr_2(2);
    
        temp_mat = [lr_1(1) 0 lr_1(2) 0;0 lr_1(1) 0 lr_1(2); lr_2(1) 0 lr_2(2) 0;0 lr_2(1) 0 lr_2(2)];
        temp_b = -[s_1*lr_1.';s_2*lr_2.'];
    
        % coeff_Bij = temp_mat \ temp_b;
    
        % lambda_1 = [lr_1(1) lr_1(1)*s_1 lr_1(2) lr_1(2)*s_1]';
        % lambda_2 = [lr_2(1) lr_2(1)*s_2 lr_2(2) lr_2(2)*s_2]';
    
        [Q,R] = qr([lambda_1 lambda_2]);
    
        Coeff_Bij_unitary =  [Q(1,3) Q(1,4); Q(3,3) Q(3,4)]*inv( [Q(2,3) Q(2,4); Q(4,3) Q(4,4)]);
        Coeff_Bij_unitary = Coeff_Bij_unitary.';
        coeff_Bij = Coeff_Bij_unitary(:);
    
        % multiply coefficients
        b_11 = leftpadz(conv(b_11_old,[1;coeff_Bij(1)]), conv(b_12_old,[coeff_Bij(3)])) + leftpadz(conv(b_12_old,[coeff_Bij(3)]),conv(b_11_old,[1;coeff_Bij(1)]));
        b_12 = leftpadz(conv(b_11_old,[coeff_Bij(2)]), conv(b_12_old,[1;coeff_Bij(4)])) + leftpadz(conv(b_12_old,[1;coeff_Bij(4)]),conv(b_11_old,[coeff_Bij(2)]));
        b_21 = leftpadz(conv(b_21_old,[1;coeff_Bij(1)]), conv(b_22_old,[coeff_Bij(3)])) + leftpadz(conv(b_22_old,[coeff_Bij(3)]),conv(b_21_old,[1;coeff_Bij(1)]));
        b_22 = leftpadz(conv(b_21_old,[coeff_Bij(2)]), conv(b_22_old,[1;coeff_Bij(4)])) + leftpadz(conv(b_22_old,[1;coeff_Bij(4)]),conv(b_21_old,[coeff_Bij(2)]));
    end
end

% evaluate the canonical vecs from coefficients
u = flip(b_21(2:end));
v = flip(b_22(2:end));

u_c1 = u;
u_r1 = zeros(N,1);
u_r1(1) = u(1);
X1 = toeplitz(u_c1, u_r1);

v_r1 = ones(N,1);
v_r1(2:end) = flip(v(2:N));
v_c1 = zeros(N,1);
v_c1(1) = v_r1(1);
X2 = toeplitz(v_c1, v_r1);

v_c1 = v;
v_r1 = zeros(N,1);
v_r1(1) = v(1);
X3 = toeplitz(v_c1, v_r1);

u_r1 = zeros(N,1);
u_r1(2:end) = flip(u(2:N));
u_c1 = zeros(N,1);
u_c1(1) = u_r1(1);
X4 = toeplitz(u_c1, u_r1);

end

function [lr_1, lr_2, lambda_1, lambda_2, s1, s2] = select_interpolation_points(b_11, b_12, b_21, b_22, S_k, ak, n)

len = length(S_k);

l1_norm_max = -999;
l12_ortho = -999;

for i=1:len
    s_i = S_k(i);
    lr_i = [s_i^n -eval_fk(ak,s_i,n)]*[polyval(b_11,s_i) polyval(b_12,s_i); polyval(b_21,s_i) polyval(b_22,s_i)];

    lambda_i = [lr_i(1) lr_i(1)*s_i lr_i(2) lr_i(2)*s_i]';
    if l1_norm_max < norm(lambda_i)
        l1_norm_max = norm(lambda_i);
        lambda_1 = lambda_i;
        lr_1 = lr_i;
        s1 = s_i;
    end
end

for i=1:len
    s_i = S_k(i);
    lr_i = [s_i^n -eval_fk(ak,s_i,n)]*[polyval(b_11,s_i) polyval(b_12,s_i); polyval(b_21,s_i) polyval(b_22,s_i)];

    lambda_i = [lr_i(1) lr_i(1)*s_i lr_i(2) lr_i(2)*s_i]';
    comp =  norm(lambda_i - (lambda_1'*lambda_i/(norm(lambda_1)^2))*lambda_1) ;
    
    if l12_ortho < comp && s_i~=s1
        l12_ortho = comp;
        lambda_2 = lambda_i;
        lr_2 = lr_i;
        s2 = s_i;
    end
end

end


function [val] = eval_fk(a,sk,n)
poly_val = polyval(a,sk);
val = poly_val/sk^(n-1);
end