function Output = BareholeFunction(R)
x(:, 1) = ((R(:, 1) + 1) * (0.15 - 0.05) / 2) + 0.05;
x(:, 2) = ((R(:, 2) + 1) * (50000 - 100) / 2) + 100;
x(:, 3) = ((R(:, 3) + 1) * (115600 - 63700) / 2) + 63700;
x(:, 4) = ((R(:, 4) + 1) * (1100 - 990) / 2) + 990;
x(:, 5) = ((R(:, 5)+ 1) * (116 - 63.1) / 2) + 63.1;
x(:, 6) = ((R(:, 6) + 1) * (820 - 700) / 2) + 700;
x(:, 7) = ((R(:, 7) + 1) * (1680 - 1120) / 2) + 1120;
x(:, 8) = ((R(:, 8) + 1) * (12045 - 9855) / 2) + 9855;

num = (2 * pi * x(:, 3) .* (x(:, 4) - x(:, 6)));
dem1 = log(x(:, 2) ./ x(:, 1));
dem2 = (1 + (2 * x(:, 7) .* x(:, 3)) ./ (dem1 + x(:, 1) .^2 .* x(:, 8)) + x(:, 3) ./ x(:, 5));
Output = num ./ dem1 ./ dem2;
end

