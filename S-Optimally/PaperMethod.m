function [Error_y_F, Error_cf_F, Time_F, Criteria_F] = PaperMethod(Y, A, number_selection, Subset_Size, n_simulation, A_test, Y_test, C_real)
% number of coeffient
n_coef = size(A, 2);
n_selection = length(number_selection);

Error_cf_F = zeros(n_simulation, length(number_selection));
Error_y_F = zeros(n_simulation, length(number_selection));
Time_F = zeros(n_simulation, length(number_selection));
Criteria_F = zeros(n_simulation, length(number_selection));
S_S = Subset_Size;
parfor s = 1:n_simulation
fprintf('simulation %i started\n', s);
%Initial Point set
A_Tread = A;
Y_Thread = Y;
S_max = -1;
tic
indicies = 1:size(A, 1);
initial_index = randsample(indicies, 1);
A_I = A_Tread(initial_index, 1);
k = 1;
selected_Indicies = zeros(max(number_selection), 1);
selected_Indicies(k) = initial_index;
indicies(indicies == initial_index) = [];
while k < n_coef
   Subset_Indicies = randsample(indicies, min(length(indicies),S_S));           
    S_max = -1;
    selected = 0;
    X_prime = pinv(A_I' * A_I);
    for i = 1:length(Subset_Indicies)
        r = (A_Tread(Subset_Indicies(i), 1:k))';
        gamma = A_Tread(Subset_Indicies(i), k + 1);
        c = A_Tread(selected_Indicies(1:k), k + 1);    
        b = X_prime * r;
        g = X_prime * A_I' * c;
        alpha = (c' * A_I + gamma * r') * (eye(k) - (b * r') / (1 + r' * b)) * (g + gamma * b);
        demn = 1;
        modifier = 1e-7;
        for n = 1:k
            demn = modifier * demn * (norm(A_I(:, n)) ^ 2 + r(n)^2);
        end
        S = ((1 + r' * b) / (demn)) * (c' * c + gamma ^ 2 - alpha) / (c' * c + gamma ^ 2);
        if i == 1
            S_max = S;
            selected = Subset_Indicies(i); 
        else
            if S > S_max
                selected = Subset_Indicies(i);   
                S_max = S;
            end
        end
        
    end
    k = k + 1;
    selected_Indicies(k) = selected;
    indicies(indicies == selected) = [];
    A_I = A_Tread(selected_Indicies(1:k), 1:k);
end
A_S = A(selected_Indicies(1:n_coef),1:n_coef);
X_prime = pinv(A_S' * A_S);
for m = 1:n_selection
    %Greedy Algorithm
    while k < number_selection(m)
        Subset_Indicies = randsample(indicies, min(length(indicies),S_S));
        S_max = -1;
        selected = 0;
        for j = 1:length(Subset_Indicies)
            r = (A_Tread(Subset_Indicies(j), 1:n_coef))';
            demn = 1;
            modifier = 1e-7;
            for n = 1:n_coef
                demn = modifier * demn * (norm(A_S(:, n)) ^ 2 + r(n)^2);
            end
            S = (1 + r' * X_prime * r) / demn;
            if j == 1
                S_max = S;
                selected = Subset_Indicies(j); 
            else
                if S > S_max
                    selected = Subset_Indicies(j);   
                    S_max = S;
                end
            end
            
        end
        k = k + 1;
        selected_Indicies(k) = selected;
        indicies(indicies == selected) = [];
        phi = A_Tread(selected_Indicies(k), 1:n_coef);
        X_prime = X_prime - (X_prime * (phi' * phi) * X_prime) / (1 + phi * X_prime * phi');
        A_S = A_Tread(selected_Indicies(1:k), 1:n_coef);
    end
    
    % Compute the Coefficients
    Y_temp =  Y_Thread(selected_Indicies(1:number_selection(m)));
    A_S = A_Tread(selected_Indicies(1:number_selection(m)), :);
    C = lsqr(A_S, Y_temp, 1e-6, 1e6);
%       C = pinv(A_S' * A_S) * A_S' * Y;
    %compute the l2 norm of error
    Y_hat = A_test * C;
    Error1 = abs(Y_test - Y_hat);
    Error_y_F(s, m) = norm(Error1) / norm(Y_test);
    Criteria_F(s, m) = S_max
    Error2 = abs(C_real - C);
    Error_cf_F(s, m) = norm(Error2) / norm(C_real);

    Time_F(s, m) = toc;
end
fprintf('simulation %i finished\n', s);
end
end