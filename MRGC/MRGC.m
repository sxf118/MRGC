function U = MRGC(X, num_bases, alpha, beta, num_iters)

normal = 1;
k = 9;

%% Set initial value
% myeps = 10e-11;
threshold = 10e-6;
v = size(X, 2);
for i = 1:v
	X{i} = X{i}';
end
n = size(X{1}, 1);
distx = cell(1, v);
sorted_distx = cell(1, v);
idx = cell(1, v);
B = cell(1,v);
W = cell(1, v);
D_W = cell(1, v);
L_W = cell(1, v);
rr = zeros(1, n);
gamma = zeros(1, v);
w = ones(1, v) / v;

%% Normalization
if normal
    for i = 1:v
        for j = 1:n
            X{i}(j,:) = ( X{i}(j,:) - mean(X{i}(j,:)) ) / std(X{i}(j,:));
        end
    end
end

%% Initialize W, B, U, F
for i = 1:v
    distx{i} = L2_distance_1( X{i}', X{i}' );
    [sorted_distx{i}, idx{i}] = sort(distx{i}, 2);
    W{i} = zeros(n);
    for j = 1:n
        di = sorted_distx{i}(j, 2:k+2);
        id = idx{i}(j, 2:k+2);
        W{i}(i,id) = ( di(k+1) - di ) / ( k * di(k+1) - sum(di(1:k)) + eps );
    end
%     options = [];
%     options.NeighborMode = 'KNN';
%     options.k = k;
%     W{i} = constructW(X{i}, options);
    B{i} = rand(size(X{i}, 2), num_bases(i))-0.5;
	B{i} = B{i} - repmat(mean(B{i},1), size(B{i},1),1);
    B{i} = B{i} * diag(1./sqrt(sum(B{i} .* B{i})));
end
U = zeros(n);
for i = 1:v
    U = U + W{i};
end
U = U / v;
for i = 1:v
    U(i,:) = U(i,:) / sum(U(i,:));
end
% U = (U + U') / 2;
% D_U = diag(sum(U));
% L_U = D_U - U;

%% Update Itertions
disp('Iteration begin!');
t = 1;
while t <= num_iters
    for i = 1:v
        W{i} = (W{i} + W{i}') / 2;
        D_W{i} = diag(sum(W{i}));
        L_W{i} = D_W{i} - W{i};
        
        % update S
        if t ==1
            S{i} = learn_coefficients(B{i}, X{i}', alpha(i), beta(i), L_W{i});
        else
            S{i} = learn_coefficients(B{i}, X{i}', alpha(i), beta(i), L_W{i}, S{i});
        end
        S{i}(isnan(S{i}))=0;
		
		% update B
        B{i} = learn_basis(X{i}', S{i}, 1);
        
        % update W
        W{i} = zeros(n);
        distx{i} = L2_distance_1(S{i}, S{i});
        [sorted_distx{i}, idx{i}] = sort(distx{i}, 2);
        for j = 1:n
            id = idx{i}(j, 2:k+2);
            di = distx{i}(j, id);
			numerator = alpha(i) * (di(k+1) - di) + 2 * w(i) * U(j,id(:)) - 2 * w(i) * U(j,id(k+1));
            denominator1 = k * di(k+1) - sum(di(1:k));
			denominator1 = alpha(i) * denominator1;
			denominator2 = 2 * w(i) * sum(U(j,id(1:k))) - 2 * k * w(i) * U(j, id(k+1));
            W{i}(j,id) = max(numerator/(denominator1+denominator2+eps), 0);
            rr(j) = k * di(k+1) - sum(di(1:k)) - 2 * k * w(i) * U(j,id(k+1)) - 2 * w(i);
        end
        gamma(v) = mean(rr);
        
        % update w
        US = U - W{i};
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end
        w(i) = 0.5 / sqrt(distUS);
    end
    
    % update U
    U = zeros(n);
    for i=1:n
        idxw = zeros();
        for j = 1:v
            W_row = W{j}(i,:);
            idxw = [idxw, find(W_row>0)];
        end
        idxw = unique(idxw(2:end));
        
        for j = 1:v
            W_row = W{j}(i,:);
            wi = W_row(idxw);
            q(j,:) = wi;
        end
        
        U(i,idxw) = SloutionToP19(q, v);
        clear q;
    end
    
    obj(t) = mgetobj(X, B, S, W, U, L_W, alpha, beta, gamma, w);
    fprintf('Iteration %d done\n', t);
    
    if(t > 1 && abs((obj(t-1)-obj(t))/obj(t-1)) <= threshold)
        break;
    end
    t = t + 1;
end
