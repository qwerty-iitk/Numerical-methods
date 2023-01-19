function muller(f, x0, x1, x2, threshold, iter)
x = linspace(x0-30, x0+30, 1000);
y = zeros(1, 1000);
c=1;
for v = x
    y(1,c) = feval(f,v);
    c = c+1;
end
plot(x, y);
xlabel('x')
ylabel('f(x)')
title('f(x) vs x');

while(1)
    x_new = x2;
    h0 = x1 - x0;
    h1 = x2 - x1;
    d0 = (f(x1) - f(x0)) / h0;
    d1 = (f(x2) - f(x1)) / h1;
    a = (d1 - d0) / (h1 + h0);
    b = a*h1 + d1;
    c = f(x2);
    disc = sqrt(b*b - 4*a*c);
   
    if(abs(b+disc) > abs(b-disc))
        den = b + disc;
    else
        den = b - disc;
    end
    diff = -2*c/den;
    x_new = x_new + diff;
    err = (abs(diff)/x_new) * 100;
    iter = iter - 1;
    if ~(err > threshold && iter > 0 && f(x_new) ~= 0)
        break;
    end
    x0 = x1;
    x1 = x2;
    x2 = x_new;
end
disp('\n');
disp(['Root of equation is ',num2str(x2), ' and final error is ', num2str(err)]);
end