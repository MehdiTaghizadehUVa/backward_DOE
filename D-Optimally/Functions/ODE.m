function fsol = ODE(R)
solinit=bvpinit(linspace(0,1,1000),[0 1]);
    sol = bvp4c(@odefun,@mybvpbc,solinit);
    m=max(size(sol.y));
    if mod(m,2)==1
    fsol=sol.y(1,(m+1)/2);
    else
    fsol=0.5*(sol.y(1,m/2)+sol.y(1,m/2+1));
    end   
        function f= odefun(t,Y)
            a=1;
            b=0;
            for i=1:size(R,2)
            a= a+1*cos(2*pi*i*t)*R(i)/(i*pi);  
            b=b-2*sin(2*pi*i*t)*R(i);
            end
        f=[Y(2); (-2-b*Y(2))/a];
        end
    
        function res= mybvpbc(yl,yr)
          res=[yl(1);yr(1)];
        end
    
end