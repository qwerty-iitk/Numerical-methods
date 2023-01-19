f=fopen('Input.txt');
F=fopen('Output.txt','w');
m= input('Enter the number corresponding to the method you wish to use.\n1. Direct power method \n2. Inverse power method \n3. Shifted-power method \n4. QR method: ');
n=fscanf(f,'%f\n',1);
C=readmatrix('Input.txt'); %the augmented matrix i.e C=[A|b]
A= C(:,1:n);  %the equation being A X= b
f=fopen('Input2.txt');
i_max=fscanf(f,'%f\n',1);
e_max=fscanf(f,'%f\n',1);
X=ones(n,1); %column matrix with all entries 1

                        %%% Direct power method
 lambda=1;
 lambda_new=1;
if (m==1)
    for i=1:i_max
        lambda=lambda_new;
        X_new = A*X;
        lambda_new= max(abs(X_new));
        e_a=abs((lambda_new-lambda)/lambda_new)*100;
        X=X_new/lambda_new;
        if (e_a<=e_max) 
            break;
        end
    end
    X= X/norm(X);
    fprintf(F,'\t\t\tDirect power method\n\n');
    fprintf(F,'Iteration\n');
    fprintf(F,' %d\n\n',i);
    fprintf(F,'Eigenvalue \n');
    fprintf(F,' %f\n\n',lambda_new);
    fprintf(F,'Eigenvector \n');
    for j=1:n
            fprintf(F,' %f\n',X(j,1));
    end 
end


                    %%% Inverse power method
if (m==2)
    B=inv(A);
   for i=1:i_max
        lambda=lambda_new;
        X_new = B*X;
        lambda_new= max(abs(X_new));
        e_a=abs((lambda_new-lambda)/lambda_new)*100;
        X=X_new/lambda_new;
        if (e_a<=e_max) 
            break;
        end
    end
    X= X/norm(X);
    fprintf(F,'\t\t\tInverse power method\n\n');
    fprintf(F,'Iteration\n');
    fprintf(F,' %d\n\n',i);
    fprintf(F,'Eigenvalue \n');
    fprintf(F,' %f\n\n',1/lambda_new); % to find minimum eigen value 
    fprintf(F,'Eigenvector \n');
    for j=1:n
            fprintf(F,' %f\n',X(j,1));
    end  
end

                    %%% Shifted-power method
 if (m==3)
      
             S=fscanf(f,'%f\n',1);
            I=eye(n); %return identity matrix of n by n 
            B=(A-S*I);
            %inverse power method to estimate the eigen value
            B=inv(B);
             for i=1:i_max
                lambda=lambda_new;
                X_new = B*X;
                lambda_new= max(abs(X_new));
                e_a=abs((lambda_new-lambda)/lambda_new)*100;
                X=X_new/lambda_new;
                if (e_a<=e_max) 
                    break;
                end
             end
             X= X/norm(X);
             fprintf(F,'\t\t\tShifted power method\n\t\t\t\n\n');
            fprintf(F,'Iteration\n');
            fprintf(F,' %d\n\n',i);
            fprintf(F,'Eigenvalue \n');
            lambda_new=((1/lambda_new)+S);
            fprintf(F,' %f\n\n',lambda_new); % to find minimum eigen value 
            fprintf(F,'Eigenvector \n');
            for j=1:n
                    fprintf(F,' %f\n',X(j,1));
            end 
     
 end
 
  
                 %%% QR method
 if (m==4)
     lambda=0;
       lambda_new=0;
       
   for k=1:i_max
       Q=zeros(n);
       lambda=lambda_new;
             %Gram Schmidt Process for finding Q and R
% %        Q(:,1)=A(:,1)/norm(A(:,1)); %finding Q
% %        fprintf('%d\n',Q);    
       for i=1:n
           Q(i,1)=A(i,1)/norm(A(:,1)); %finding first column of Q
       end  
       for i=2:n
           B=zeros(n,1);  
           B= A(:,i);
        for j=1:i-1
           B= B-(Q(:,j)'*B)*Q(:,j);     
        end
         Q(:,i)=B/norm(B);
       end
       R=zeros(n);
       for i=1:n     %finding R
           for j=i:n
               R(i,j)= Q(:,i)'*A(:,j);
           end
       end
      A=R*Q;
      lambda_new =max(abs(A));
      lambda_new=max(abs(lambda_new));
%       fprintf('%f \n',lambda_new);
       e_a=abs((lambda_new-lambda)/lambda_new)*100;
       if e_a<=e_max
           break;
       end
   end
    fprintf(F,'\t\t\tQR method\n\n');
    fprintf(F,'Iteration\n');
    fprintf(F,' %d\n\n',k);
    fprintf(F,'Eigenvalues \n');
        for i=1:n
            fprintf(F,' %f\n',A(i,i));
        end 
 end
fclose(F);