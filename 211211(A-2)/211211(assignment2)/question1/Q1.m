f=fopen('Input.txt');
F2=fopen('Output.txt','w');
n=fscanf(f,'%f\n',1);  %number of equation
C=readmatrix('Input.txt'); %the augmented matrix i.e C=[A|b]
A= C(:,1:n);  %the equation being A X= b
B= C(:,n+1);
D=repmat(C,1); %making duplicate of matrix C
 for i=1:n         % identity matrix I(n*n)
    for j=1:n
        I(i,j)=0;
        if i==j 
            I(i,j)=1;
        end
    end
 end
m= input('Enter the number corresponding to the method you wish to use.\n1. Gauss elimination (GE; without pivoting)\n2. Gauss elimination (GE; with pivoting)\n3. LU decomposition by using Doolittle (without pivoting)\n4. LU decomposition by using Crout method (without pivoting)\n5. Cholesky decomposition (for symmetric positive definite matrix) : ');
               
if(m==1)           
                                %%%  Gauss elimination (GE; without pivoting) %%%
       for j=1:(n-1) %coverting augmented matrix C into upper diagonal matrix
           for i=j+1:n
               mult=C(i,j)/C(j,j);
               for k=j:n+1
               C(i,k)=C(i,k)- mult*C(j,k);
%                C
               end
           end
       end
       fprintf(F2,'\t\t\tGauss elimination without pivoting\n\n');
       fprintf(F2,'\nUpper Triangular Matrix U\n\n');
             for i=1:n
                for j=1:n
                fprintf(F2,'%f\t\t',C(i,j));
                end 
                fprintf(F2,'\n');
              end
%            fprintf('%f ',C);
            for i=n:-1:1        % back substitution
                for j=1:n-1
                    if i==n+1-j
                        break;
                    end
                    C(i,n+1)=C(i,n+1)-C(i,n+1-j)*X(n+1-j);
                end
                X(i)=C(i,n+1)/C(i,i);
            end
              fprintf(F2,'Solution Matrix X\n\n'); %Printing the output
              for i=1:n
                 fprintf(F2,'%f\n',X(i)); 
              end
end

if m==2                   %%%%   Gauss elimination (GE; with pivoting)  %%%
        for j=1:(n-1) %converting augmented matrix C into upper diagonal matrix
          r=j; %row no of the current pivot element
          Big=abs(C(j,j));%current pivot element
              for k=r+1:n   %pivotting
               temp=abs(C(k,r));
               if temp>Big
                   Big=temp;   %     largest element in the pivot column
                   r=k;
               end
              end
              if(r~=j)
                C([j r],:)=C([r j],:); % replacing current row with pivotted row
                I([j r],:)=I([r j],:);
              end   
%               fprintf('%f ',C);
%               fprintf('\n%f ',I);
            for i=j+1:n
               mult=C(i,j)/C(j,j);
               C(i,:)=C(i,:)- mult*C(j,:);
            end
        end
      fprintf(F2,'\t\t\tGauss elimination with pivoting\n\n');
      fprintf(F2,'\nUpper Triangular Matrix U\n\n');
         for i=1:n
            for j=1:n
            fprintf(F2,'%f\t\t',C(i,j));
            end 
            fprintf(F2,'\n');
          end
      fprintf(F2,'Permutation Matrix P\n\n');
          for i=1:n
            for j=1:n
            fprintf(F2,'%f\t\t',I(i,j));
            end 
            fprintf(F2,'\n');
          end
        for i=n:-1:1        % back substitution
            for j=1:n-1
                if i==n+1-j
                    break;
                end
                C(i,n+1)=C(i,n+1)-C(i,n+1-j)*X(n+1-j);
            end
            X(i)=C(i,n+1)/C(i,i);
        end
        %   fprintf('\n%f ',X); 
          fprintf(F2,'Solution Matrix X\n\n'); %Printing the output
          for i=1:n
             fprintf(F2,'%f\n',X(i)); 
          end
end


if (m==3)
            %%%  LU decomposition by using Doolittle (without pivoting)
for i=1:n-1
    for j=i+1:n
        L(i,j)=0;
        U(j,i)=0;
    end
end
for i=1:n
 for k=1:i-1   % Finding L
      L(i,k)=A(i,k);
      for j=1:k-1
          L(i,k)= L(i,k)-L(i,j)*U(j,k);
      end
      L(i,k) = L(i,k)/U(k,k);
  end
 for k=i:n % Finding U
      U(i,k) = A(i,k);
      for j=1:i-1
          U(i,k)= U(i,k)-L(i,j)*U(j,k);
      end
  end
end
  for i=1:n
      L(i,i)=1;
  end
%     fprintf('\n%f ',L);
    fprintf(F2,'\t\t\tLU decomposition by using Doolittle (without pivoting)\n\n');
    fprintf(F2,'Lower Diagonal Matrix L\n\n'); 
          for i=1:n
            for j=1:n
            fprintf(F2,' %f\t\t ',L(i,j));
            end 
            fprintf(F2,'\n');
          end
     fprintf(F2,'Upper Diagonal Matrix U\n\n'); 
          for i=1:n
            for j=1:n
            fprintf(F2,' %f\t\t ',U(i,j));
            end 
            fprintf(F2,'\n');
          end
          %forward substitution to find Y= [L|B]
          for i=1:n
              for j=1:n
                  if  i==j
                      break;
                  end
                  B(i,1)=B(i,1)-L(i,j)*Y(j,1);
              end
              Y(i,1)=B(i,1)/L(i,i);
          end
          %Backward substitution to find X= [U|Y] 
          for i=n:-1:1       
                for j=n:-1:1
                    if i==j
                        break;
                    end
                    Y(i,1)=Y(i,1)-U(i,j)*X(j);
                end
                X(i)=Y(i,1)/U(i,i);
          end
        fprintf(F2,'Solution Matrix X\n\n'); %Printing the output
          for i=1:n
             fprintf(F2,'%f\n',X(i)); 
          end    
end

                            %%% LU decomposition by using Crout method (without pivoting)
if m==4
L=zeros(size(A));
U=zeros(size(A));
L(:,1)=A(:,1);
U(1,:)=A(1,:)/L(1,1);
U(1,1)=1;
for k=2:n
for j=2:n
    for i=j:n
        L(i,j)=A(i,j)-(L(i,1:j-1)*U(1:j-1,j));
%  L
    end
    U(k,j)=(A(k,j)-(L(k,1:k-1)*U(1:k-1,j)))/L(k,k);
end
end
    fprintf(F2,'\t\t\tLU decomposition by using Crout method (without pivoting)\n\n');
    fprintf(F2,'Lower Diagonal Matrix L\n\n'); 
          for i=1:n
            for j=1:n
            fprintf(F2,' %f\t\t ',L(i,j));
            end 
            fprintf(F2,'\n');
          end
     fprintf(F2,'Upper Diagonal Matrix U\n\n'); 
          for i=1:n
            for j=1:n
            fprintf(F2,' %f\t\t ',U(i,j));
            end 
            fprintf(F2,'\n');
          end

           %forward substitution to find Y= [L|B]
          for i=1:n
              for j=1:n
                  if  i==j
                      break;
                  end
                  B(i,1)=B(i,1)-L(i,j)*Y(j,1);
              end
              Y(i,1)=B(i,1)/L(i,i);
          end        
         %Backward substitution to find X= [U|Y]
          for i=n:-1:1       
                for j=n:-1:1
                    if i==j
                        break;
                    end
                    Y(i,1)=Y(i,1)-U(i,j)*X(j);
                end
                X(i)=Y(i,1)/U(i,i);
          end
        fprintf(F2,'Solution Matrix X\n\n'); %Printing the output
          for i=1:n
             fprintf(F2,'%f\n',X(i)); 
          end  
          
end
if(m==5)
                %%% Cholesky decomposition (for symmetric positive definite matrix)
      L= zeros(n,n);
     if(A==A.')  %Only if A is a symmetric matrix
       for k=1:n
          for i=1: k-1
              sum=0;
              for j=1:i-1
                  sum=sum+L(i,j)*L(k,j);
              end
              L(k,i)=(A(k,i)-sum)/L(i,i);
          end
          sum1=0;
          for j=1:k-1
              sum1=sum1+L(k,j)*L(k,j);
          end 
          L(k,k)=sqrt(A(k,k)-sum1);
       end
      end
      fprintf(F2,'\t\t\tCholesky decomposition\n\n');
      fprintf(F2,'Cholesky factor Matrix L\n\n'); 
          for i=1:n
            for j=1:n
            fprintf(F2,' %f\t\t ',L(i,j));
            end 
            fprintf(F2,'\n');
          end
          for i=1:n % forward substitution to find Y=[L|B]
              for j=1:n
                  if  i==j
                      break;
                  end
                  B(i,1)=B(i,1)-L(i,j)*Y(j,1);
              end
              Y(i,1)=B(i,1)/L(i,i);
          end
          U=L.';
          % backward substitution to find X=[U|Y]
          for i=n:-1:1       
                for j=n:-1:1
                    if i==j
                        break;
                    end
                    Y(i,1)=Y(i,1)-U(i,j)*X(j);
                end
                X(i)=Y(i,1)/U(i,i);
          end
          fprintf(F2,'Solution Matrix X\n\n'); %Printing the output
          for i=1:n
             fprintf(F2,'%f\n',X(i)); 
          end
end

fclose(F2);
fclose(f);