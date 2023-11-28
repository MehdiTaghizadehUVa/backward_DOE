function [Error_cf, Error_y] = MonteCarlo_Fast(A, Y, number_selection, n_simulation, A_test, Y_test, C_real)

E = zeros(n_simulation,length(number_selection));
for k = 1:length(number_selection)
    n_s = number_selection(k);
    parfor s = 1:n_simulation
        A_thread = A;
        Indicies = randsample(1:length(Y), n_s);
         %Compute Experiment Matrix
        A_temp = A_thread(Indicies, :);
        % Compute the Coefficients
        Y_temp = Y(Indicies,:);
        C = lsqr(A_temp, Y_temp, 1e-6, 1e6);
        %C = pinv(A' * A) * A' * Y;
        %compute the l2 norm of error
        Error1 = abs(C_real - C) ;
        Error_cf(s, k) = norm(Error1) / norm(C_real);

        Y_hat = A_test * C;
        Error2 = abs(Y_test - Y_hat) ;
        Error_y(s, k) = norm(Error2) / norm(Y_test);
    end
end
end