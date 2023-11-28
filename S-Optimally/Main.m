clear all
clc
addpath ./Negin_codes/Poly-functions
addpath ./Functions
rng(1,"twister");

number_Sample = 2000;
number_Test = 4000;
%Dimension
dim = 3;
% Polynomial Order
order = 12;
%--------------------------------------------------------------------------
load('3D-12O-Ishigami-Data.mat')

% Number of coeffient
n_coef = factorial(dim + order) / (factorial(dim) * factorial(order));
%Training Data
X = sample.rv(1:number_Sample, :);
Y = bsxfun(@times, sample.w(1:number_Sample, :), IshigamiFunction(X));
% Training model matrix
A = sample.wlhs(1:number_Sample, :);
%%Test data
X_test = sample.rv(number_Sample+1:number_Sample + number_Test, :);
Y_test = bsxfun(@times, sample.w(number_Sample+1:number_Sample + number_Test, :), IshigamiFunction(X_test));
%Test mide matrix
A_test = sample.wlhs(number_Sample+1:number_Sample + number_Test, :);
%--------------------------------------------------------------------------
n_simulation = 300;
NumberofSelection = [470, 500:50:900];
%--------------------------------------------------------------------------
C_real = lsqr(A_test, Y_test, 1e-6, 1e6);
%Monte Carlo Method
% [Error_cf, Error_y] = MonteCarlo_Fast(A, Y, NumberofSelection, n_simulation, A_test, Y_test, C_real);
%Paper Method
Subset_Size = number_Sample;
[Error_y_F, Error_cf_F, Time_F, Criteria_F] = PaperMethod(Y, A, NumberofSelection, Subset_Size, n_simulation, A_test, Y_test, C_real);
%Our Method
Subset_Size = 30;
[Error_y_B, Error_cf_B, Time_B, Criteria_B] = OptimalPointSelection(Y, A, NumberofSelection, Subset_Size, n_simulation, A_test, Y_test, C_real);

save('3D-12O-Ishigami-Result-S.mat')
exit