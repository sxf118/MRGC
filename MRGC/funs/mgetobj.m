function re = mgetobj(X, B, S, W, U, L, alpha, beta, gamma, omega)
re = 0;
v = size(X, 2);
%     size(X{1})
%     size(B{1})
%     size(S{1})
%     size(L)
for i = 1:v
    re = re + sum( sum( (X{i}'-B{i}*S{i}).^2 ) ) + alpha(i)*trace(S{i}*L{i}*S{i}') + beta(i)*sum(sum(abs(S{i}))) + gamma(i)*norm(W{i}, 'fro')^2 + omega(i)*norm(U-W{i}, 'fro')^2;
end
