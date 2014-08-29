% Parameter M adjacency matrix where M_i,j represents the link 
%     from 'j' to 'i', such that for all 'j' sum(i, M_i,j) = 1
% Parameter d damping factor
% Parameter v_quadratic_error quadratic error for v
% Return v, a vector of ranks such that v_i is the i-th rank from [0, 1]
 
function [v] = rank(M, d, v_quadratic_error)
 
N = size(M, 2); % N is equal to half the size of M
v = ones(N,1)*0.5; %rand(N, 1);
v = v ./ norm(v, 2);
last_v = ones(N, 1) * inf;
M_hat = (d .* M) + (((1 - d) / N) .* ones(N, N));

i =0;
while(norm(v - last_v, 2) > v_quadratic_error)
        last_v = v;
        v = M_hat * v;
        v = v ./ norm(v, 2);
        i = i + 1;
        if i == 30
            %warning('Page rank did not converge');
            break
        end
end
 
end

