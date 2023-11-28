% % This is an example to show improvement achived by preconditioning
addpath ./Poly-functions
% %Setting 
%2-Pol_dimension-20-order polynomial
Pol_dim=3;
Pol_order=8;
Num_bas=func_PCE_NumPols(Pol_dim,Pol_order);
Pol_Type='Legendre';
% NumofCoe=50; %Number of nonzero coefficients for the coefficient vector
sam_size_low=1000;
% %%%%%%%%%%%%%%%%%%%%%%

       rng('shuffle')

       Samp_rand =unifrnd(-1,1,[sam_size_low,Pol_dim]); %Randomly selects samples 
       
        Psi=ones(sam_size_low,Num_bas);%construct the measurement matrix
        for i=1:Num_bas
        multi_index = func_PCE_MultiIndex(i-1, Pol_dim); 
            for d=1:Pol_dim
            c=func_PCE_1DPoly(multi_index(d), Pol_Type);
                for j=1:sam_size_low
                Psi(j,i) = polyval(c,Samp_rand(j,d))*Psi(j,i);
                end      
            end
        end
 
        %Construct the mesh
x = -5:.08:5; 
y = -5:.08:5; 
[X,Y] = meshgrid(x,y);

x1 = -5:.01:5; 
y1 = -5:.01:5; 
[X1,Y1] = meshgrid(x1,y1);

%Evaluate the integral 

for i=1:sam_size_low
  a=Samp_rand(i,1);
  b=Samp_rand(i,2);
  c=Samp_rand(i,3);
F=-20*(1+0.2*a)*exp(-0.2*(1+0.2*b)*0.25*(X.^2+Y.^2).^0.5)-exp(0.5*(cos(2*pi*(1+0.3*c)*X)+cos(2*pi*(1+0.2*c)*Y)));
F1=-20*(1+0.2*a)*exp(-0.2*(1+0.2*b)*0.25*(X1.^2+Y1.^2).^0.5)-exp(0.5*(cos(2*pi*(1+0.3*c)*X1)+cos(2*pi*(1+0.2*c)*Y1)));
b_low(i,1)=trapz(y,trapz(x,F,2));
b_high(i,1)=trapz(y1,trapz(x1,F1,2));
end

