function Output = AckleyFunction(Samples)
%Construct the mesh
x1 = -5:.01:5; 
y1 = -5:.01:5; 
[X1,Y1] = meshgrid(x1,y1);

number_Sample = size(Samples, 1);
Output = zeros(number_Sample, 1);

parfor i=1:number_Sample
    Thread_samples = Samples;
    a=Thread_samples(i,1);
    b=Thread_samples(i,2);
    c=Thread_samples(i,3);
    
    F1=-20*(1+0.2*a)*exp(-0.2*(1+0.2*b)*0.25*(X1.^2+Y1.^2).^0.5)-exp(0.5*(cos(2*pi*(1+0.3*c)*X1)+cos(2*pi*(1+0.2*c)*Y1)));
    Output(i,1)=trapz(y1,trapz(x1,F1,2));
end
end