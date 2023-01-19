function bairstow(p,f, r, s, threshold, iter)
%{
p is ploymial in string format.
r and s are inital guesses.
threshold is the maximum error in perentage.
iter in the maximum iterations to be performed.
%} 

%x = linspace(r-30, r+30, 1000);
%y = zeros(1, 1000);
%fplot(f);
%title('Plot of f(x) vs x');
%xlabel('x');ylabel('f(x)');

x=[-10:0.1:10];
y=subs(f,x);
plot(x,y), grid on, xlabel('x'), ylabel('f(x)'), title('f(x) vs x')


syms x;
n = polynomialDegree(str2sym(p), x);
a = coeffs(str2sym(p), x, 'All');
roots = zeros(n,1);

while(n >= 3)
    count = 1;
    b = zeros(1,n+1);
    c = zeros(1,n+1);
    err1 = 100;
    err2 = 100;
    
    while((err1 > threshold || err2 > threshold) && (iter - count >= 0))
        b(n+1) = a(n+1);
        b(n) = a(n) + r*b(n+1);
        c(n+1) = b(n+1);
        c(n) = b(n) + r*c(n+1);
        for i = n-1:-1:1
            b(i) = a(i) + r*b(i+1) + s*b(i+2);
            c(i) = b(i) + r*c(i+1) + s*c(i+2);
        end
        
        det = c(3)*c(3) - c(2)*c(4);
        dr = (b(1)*c(4) - b(2)*c(3)) / det;
        ds = (b(2)*c(2) - b(1)*c(3)) / det;
        r = r + dr;
        s = s + ds;
        err1 = abs(dr/r)*100; 
        err2 = abs(ds/s)*100; 
        count = count + 1;
    end
    
    disp(b);
    disp(c);
    disp(['values of r and s: ', num2str(r), ', ', num2str(s)]);
    disc = r*r + 4*s;
    roots(n) = (r + sqrt(disc)) / 2;
    roots(n-1) = (r - sqrt(disc)) / 2;
    
    a = b(3:(n+1));
    n = n-2;
end

if(n == 1)
    roots(1) = -a(1)/a(2);
elseif(n == 2)
    disc = a(2)*a(2) - 4*a(1)*a(3);
    roots(1) = (-a(2) + sqrt(disc)) / (2*a(3));
    roots(2) = (-a(2) - sqrt(disc)) / (2*a(3));
end

disp('Root of equation are: ');
disp(roots);
end