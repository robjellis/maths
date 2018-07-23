function res = fibonacci(n)

% RJ version

if n < 2
    res = 1;
else
    res = fibonacci(n-1) + fibonacci(n-2); 
end

end