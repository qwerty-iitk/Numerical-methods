fileID = fopen("input.txt",'r') ;
formatSpec = '%c' ;
func = fscanf(fileID,formatSpec,6);
func = char(func);
f = inline(func,'x','y');
formatSpec='%f' ;
x0 = fscanf(fileID,formatSpec,1) ;
y0 = fscanf(fileID,formatSpec,1) ;
xf = fscanf(fileID,formatSpec,1) ;
h = fscanf(fileID,formatSpec,1) ;
n = (xf-x0)/h ;
method = fscanf(fileID,formatSpec,1) ;
x(1)=x0 ;
y(1)=y0 ;
fclose(fileID) ;
if method == 1
 for i = 2:n+1
        x(i) = x(i-1)+h;
        y(i) = y(i-1)+h*f(x(i-1),y(i-1));
 end 
        
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','Euler Forward Method');
    title('y vs t');
    xlabel('t');
    ylabel('y');
    legend;
    hold on;
    fID = fopen('euler_forward_output.txt','wt');
    fprintf(fID,'t                 y\n');
    for i = 1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
end
if method == 2
     for i = 2:n+1
        x(i) = x(i-1)+h;
        z(i) = y(i-1)+h*f(x(i-1),y(i-1)) ;
        y(i) = y(i-1)+h*f(x(i-1),z(i));
     end 
      plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','2nd order RK method ');
    title('y vs t');
    xlabel('t');
    ylabel('y');
    legend;
    hold on;
    fID = fopen('RK_method_2nd_order.txt','wt');
    fprintf(fID,'t                 y\n');
    for i = 1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
end
if method == 3
    for i = 2:n+1
        x(i) = x(i-1)+h;
        phi0 = f(x(i-1),y(i-1));
        phi1 = f(x(i-1)+0.5*h,y(i-1)+0.5*h*phi0);
        phi2 = f(x(i-1)+0.5*h,y(i-1)+0.5*h*phi1);
        phi3 = f(x(i-1)+h,y(i-1)+h*phi2);
        y(i) = y(i-1)+h*(phi0/6+(phi1+phi2)/3+phi3/6);
    end
    plot(x,y,'-o','MarkerEdgeColor','r','DisplayName','4th order Runge Kutta');
    title('y vs x');
    xlabel('x');
    ylabel('y');
    legend;
    hold on;
    fID = fopen('runga_kutta_output.txt','wt');
    fprintf(fID,'x                 y\n');
    for i = 1:n+1
        fprintf(fID,'%f     %f\n',x(i),y(i));
    end
end