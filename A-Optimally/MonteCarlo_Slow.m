function Error_y = MonteCarlo_Slow(dim, order, number_selection, n_simulation)
X_test = -1 + 2 * rand(4000, dim);
Y_test = IshigamiFunction(X_test);
A_test = legendrePoly(order, X_test);
X = -1 + 2 * rand(2000, dim);
Y = IshigamiFunction(X);
A = legendrePoly(order, X);
Error_y = zeros(n_simulation,length(number_selection));
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
        Y_hat = A_test * C;
        Error = abs(Y_test - Y_hat) ;
        Error_y(s, k) = norm(Error) / norm(Y_test);
    end
end
end