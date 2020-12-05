clear
clc

% Number of nodes
n = 5;
A = zeros(n);

% Create the adjacent matrix
for i=1:n
    for j=i:n
        if i==j
            continue;
        else
            A(i,j) = randi([0,1]);
            A(j,i) = A(i,j);
        end
    end
end

% e vector all ones
e = ones(1,n);

A
B = A*A;

% The last vector
c3 = (A.*B) * e' ./2