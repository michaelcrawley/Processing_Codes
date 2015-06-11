function soln = TDMS(diag,RHS,sub,sup)

%TDMS solves a tridiagonal matrix system and returns the solution in the
%form of a vector. The inputs are a vector of the diagonal elements,
%right-hand-side elements, sub-diagonal elements, and super-diagonal
%elements. An example usage is shown below:
%
%Solution Vector = TDMS(diagonal, right-hand side, sub-diagonal, super-diagonal);
%

%reformatting the input vectors if necessary
if size(diag,1) == 1
    diag = diag';
end
n = length(diag);

if size(RHS,1) == 1
    RHS = RHS';
end

if size(sup,1) == 1
    sup = sup';
end
if length(sup) < n
    sup = [sup; 0];
end

if size(sub,1) == 1
    sub = sub';
end
if length(sub) < n
    sub = [0; sub];
end

if (length(RHS) + length(sup) + length(sub))/3 ~= length(diag)
    fprintf('Error: Dimension of the inputs does not match');
    soln = zeros(size(diag));
else
    %Solving the tridiagonal system
    for i = 2:n
        diag(i) = diag(i) - sub(i)/diag(i-1)*sup(i-1);
        RHS(i) = RHS(i) - sub(i)/diag(i-1)*RHS(i-1);
    end
    
    for i = n-1:-1:1
        RHS(i) = RHS(i) - sup(i)/diag(i+1)*RHS(i+1);
    end
    
    soln = RHS./diag;
end