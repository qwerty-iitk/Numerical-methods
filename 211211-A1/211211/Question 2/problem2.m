str = input('Give an equation in x: ', 's');
f = str2func(strcat('@(x)',str));
method = input('\n Please select the method that you want to use: \n 1. Muller \n 2. Brainstow \n Enter the number of method that you want to use: ');

switch method
    case 1
        x0 = input('Enter value of x0: ');
        x1 = input('Enter value of x1: ');
        x2 = input('Enter value of x2: ');
        threshold = input('Enter the value of maximum error in %: ');
        max_iter = input('Maximum Iteration to be performed: ');
        muller(f, x0, x1, x2, threshold, max_iter)
    case 2
        r = input('Enter value of r: ');
        s = input('Enter the value of s: ');
        threshold = input('Enter the value of maximum error in %: ');
        max_iter = input('Maximum Iteration to be performed: ');
        bairstow(str,f, r, s, threshold, max_iter);
    otherwise
        disp('Wrong value selected')
end