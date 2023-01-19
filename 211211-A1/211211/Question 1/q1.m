
prompt= 'Enter the non-linear equation in x: ';
str = input(prompt,'s');
str1='@(x)';
str2=append(str1,str);
f=str2func(str2);
n=input('Enter the number corresponding to the method you wish to use.\n1.Bisection\n2.False-position\n3.Fixed Point Iteration\n4.Newton-Raphson\n5.Secant\nchoice: ');

maxe_s= input('Enter the maximum approximate error in percentage: ');

maxite =input('Enter the maximum number of iteration you want: ');
 subplot(2,1,1);
fplot(f);
title('Plot of f(x) vs x');
xlabel('x');ylabel('f(x)');


if n==1
    x_l = input('Enter first guess value: ');
    x_u = input('Enter second guess value: ');
    %Bisection (f,x_l,x_u,maxe_s,maxite);
    fprintf('Root= %f\n',Bisection (f,x_l,x_u,maxe_s,maxite));
    
  elseif n==2
        x_l = input('Enter first guess value: ');
        x_u = input('Enter second guess value: ');
        %FalsePosition(f,x_l,x_u,maxe_s,maxite);
        fprintf('Root= %f\n',FalsePosition(f,x_l,x_u,maxe_s,maxite));
   elseif n==3
          x_l = input('Enter first guess value: ');
              prompt= 'Enter the second  non-linear equation in x: ';
               stR = input(prompt,'s');
                     stR1='@(x)';
                 stR2=append(stR1,stR);
                    g=str2func(stR2);
             % FixedPointIteration(f,g,x_l,maxe_s,maxite);
             fprintf('Root= %f\n',fixedPoint(f,g,x_l,maxe_s,maxite));
     elseif n==4
              x_l = input('Enter first guess value: ');
              prompt= 'Enter the derivative of  non-linear equation in x: ';
               stR = input(prompt,'s');
                     stR1='@(x)';
                 stR2=append(stR1,stR);
                    g=str2func(stR2);
             % NewtonRaphson(f,g,x_l,maxe_s,maxite);
             fprintf('Root= %f\n',NewtonRaphson(f,g,x_l,maxe_s,maxite));
       elseif n==5
                  x_l = input('Enter first guess value: ');
                  x_u = input('Enter second guess value: ');
                  %Secant(f,x_l,x_u,maxe_s,maxite);
                  fprintf('Root= %f\n',Secant(f,x_l,x_u,maxe_s,maxite));
end

% fplot(f);
% title('Plot of f(x) vs x');
% xlabel('x');ylabel('f(x)');
        
     %Bisection method for finding roots
   function root =Bisection (f,x_l,x_u,maxe_s,maxite)
%    fplot(f);
%     title('Plot of f(x) vs x');
%     xlabel('x');ylabel('f(x)');
      if(x_l>x_u)
        x_temp=x_l;
        x_l=x_u;
        x_u=x_temp;
      end
        if f(x_l)==0
          'x_l is one of the root';
         elseif f(x_u)==0
             'x_u is one of the root';
         end
             test=f(x_l)*f(x_u);
        if(test>0)
             error('No root is present');
         else 
              ite=1;e_a=100;
              x_r=x_l;
              x_r_prev=x_r;
              x_r=(x_l+x_u)/2;
              temp=f(x_l)*f(x_r);
              if temp<0
                      x_u=x_r;
                  else
                      x_l=x_r;
                  end
              for ite=2:maxite
                  x_r_prev=x_r;
                  x_r=(x_l+x_u)/2;
                  ite=ite+1;
                  %temp= eval(subs(y,x,x_l))*eval(subs(y,x,x_r));
                  temp=f(x_l)*f(x_r);
                  if x_r~=0,e_a=abs((x_r-x_r_prev)/x_r)*100,end
                  if temp<0
                      x_u=x_r;
                  else
                      x_l=x_r;
                  end
                if e_a<= maxe_s,break,end  
                p(ite)=e_a;
              end
              root=x_r;

              fprintf('Number of iteration = %f\n',(ite-1));

              
%               fplot(f);
%                title('Plot of f(x) vs x');
%                xlabel('x');ylabel('f(x)');
           % fplot(e_a,ite);
           subplot(2,1,2);
           plot(p);
           title('Plot of relative approximate error vs iteration number');
           xlabel('Iteration number');ylabel('Relative approximate error');
         end
   end
    
             %False- Position method
   function root= FalsePosition (f,x_l,x_u,maxe_s,maxite)
   fplot(f);
    title('Plot of f(x) vs x');
    xlabel('x');ylabel('f(x)');
       if(x_l>x_u)
        x_temp=x_l;
        x_l=x_u;
        x_u=x_temp;
       end
        if f(x_l)==0
          'x_l is one of the root';
         elseif f(x_u)==0
             'x_u is one of the root';
         end
             test=f(x_l)*f(x_u);
        if(test>0)
             error('No root is present');
         else 
              
              for ite=1:maxite 
                  
                  x_r=x_l-(f(x_l)*(x_u-x_l)/(f(x_u)-f(x_l)));
                  ite=ite+1;
                  e_a=abs((x_r - x_l)/x_r)*100 

                   if e_a<= maxe_s,break,end  
                 p(ite)=e_a;
              

                  
                  temp=f(x_l)*f(x_r);
                 
                  if temp<0
                      x_u=x_r;
                  else
                      x_l=x_r;
                  end
              end
               
              root=x_r;
                fprintf('Number of iteration = %f\n',(ite-1));

              subplot(2,1,2);
           plot(p);
           title('Plot of relative approximate error vs iteration number');
           xlabel('Iteration number');ylabel('Relative approximate error');
         end
   end
   
               %fixed point iteration
     function root= fixedPoint(f,g,x_l,maxe_s,maxite)
   fplot(f);
    title('Plot of f(x) vs x');
    xlabel('x');ylabel('f(x)');
   ite=1;e_a=100;
   for ite=1:maxite
       x_l_prev=x_l;
       
       x_l=g(x_l_prev);
       ite=ite+1;
       if x_l~=0,e_a=abs((x_l - x_l_prev)/x_l)*100,end
       if e_a<= maxe_s,break,end
        p(ite)=e_a;
   end
    root=x_l;
      fprintf('Number of iteration = %f\n',(ite-1));

    subplot(2,1,2);
    plot(p);
           title('Plot of relative approximate error vs iteration number');
           xlabel('Iteration number');ylabel('Relative approximate error');
   end 
             % Newton Raphson method
   function root= NewtonRaphson(f,g,x_l,maxe_s,maxite)
   fplot(f);
    title('Plot of f(x) vs x');
    xlabel('x');ylabel('f(x)');
   ite=1;e_a=100;
   for ite=1:maxite
       x_l_prev=x_l;
       
       x_l=x_l_prev- (f(x_l_prev)/g(x_l_prev));
       ite=ite+1;
       if x_l~=0,e_a=abs((x_l - x_l_prev)/x_l)*100,end
       if e_a<= maxe_s,break,end
        p(ite)=e_a;
   end
    root=x_l;
      fprintf('Number of iteration = %f\n',(ite-1));

    subplot(2,1,2);
    plot(p);
           title('Plot of relative approximate error vs iteration number');
           xlabel('Iteration number');ylabel('Relative approximate error');
   end  
   
               % Secant method
   function root= Secant(f,x_l,x_u,maxe_s,maxite)
  
   fplot(f);
    title('Plot of f(x) vs x');
    xlabel('x');ylabel('f(x)');
   ite=1;e_a=100;
      g=f(x_u)-f(x_l);
       x_r=x_u-(f(x_u)*(x_u-x_l)/g);
       h=f(x_r);
       x_l=x_u;
       x_u=x_r;
   for ite=2:maxite
       g=f(x_u)-f(x_l);
       x_r=x_u-(f(x_u)*(x_u-x_l)/g);
       h=f(x_r);
       ite=ite+1;
       x_l=x_u;
       x_u=x_r;
       if x_u~=0,e_a=abs((x_u - x_l)/x_u)*100,end
       if e_a<=maxe_s, break ,end
        p(ite)=e_a;
   end
   root=x_u;
     fprintf('Number of iteration = %f\n',(ite-1));

    subplot(2,1,2);
   plot(p);
           title('Plot of relative approximate error vs iteration number');
           xlabel('Iteration number');ylabel('Relative approximate error');
   end