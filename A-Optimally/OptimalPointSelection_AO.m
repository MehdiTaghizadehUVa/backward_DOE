function [Error_y_B, Error_cf_B, Time_B, Criteria_B]  = OptimalPointSelection_AO(Y, A, NumberofSelection, Subset_Size, n_simulation, A_test, Y_test, C_real)
% number of coeffient

n_selection = length(NumberofSelection);

Criteria_B = zeros(n_simulation, length(NumberofSelection));
Error_cf_B = zeros(n_simulation, length(NumberofSelection));
Error_y_B = zeros(n_simulation, length(NumberofSelection));
Time_B = zeros(n_simulation, length(NumberofSelection));
S_S = Subset_Size;
for j = 1:n_simulation
    fprintf('simulation %i started\n', j);
    tic
    A_S = A;
    AV_min = 2;
    Y_Thread = Y;
    NS_Thread = NumberofSelection;
    Opt_Indicies = 1:size(A_S, 1);
    X_prime = pinv(A_S' * A_S);
    for m = 1:n_selection
        while size(A_S,1) > NS_Thread(n_selection + 1 - m)
            indicies = 1:size(A_S, 1);
            Subset_Indicies = randsample(indicies, min(length(indicies),S_S));       
            AV_min = 2;
            selected = 0;
            for ind = 1:length(Subset_Indicies)
                r = (A_S(Subset_Indicies(ind), :))';
                AV = trace((X_prime * (r * r') * X_prime) /  (1 - r' * X_prime * r));
                if ind == 1
                    AV_min = AV;
                    selected = Subset_Indicies(ind); 
                else
                    if AV_min > AV
                        selected = Subset_Indicies(ind);   
                        AV_min = AV;
                    end
                end
            end
            r = (A_S(selected, :))';
            X_prime = X_prime + (X_prime * (r * r') * X_prime) / (1 - r' * X_prime * r);
            A_S(selected, :) = [];
            Opt_Indicies(selected) = [];
        end
        % Compute the Coefficients
        Y_temp = Y_Thread(Opt_Indicies);
        C = lsqr(A_S, Y_temp, 1e-6, 1e6);
%       C = pinv(A_S' * A_S) * A_S' * Y;
        %compute the l2 norm of error
        Y_hat = A_test * C;
        Error1 = abs(Y_test - Y_hat);
        Error_y_B(j, m) = norm(Error1) / norm(Y_test);

        Criteria_B(j, m) = AV_min;

        Error2 = abs(C_real - C);
        Error_cf_B(j, m) = norm(Error2) / norm(C_real);

        Time_B(j, m) = toc;
    end
    fprintf('simulation %i finished\n', j);
end
end