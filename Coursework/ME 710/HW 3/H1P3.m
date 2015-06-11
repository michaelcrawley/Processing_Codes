function [Tnumeric x] = H1P3(Tl,Tr, N)
    x = 0:1/(N-1):1; %create N nodes along the x-direction
    h = mean(diff(x)); %determine step size along the x-direction

    %create coefficients for discretization of general equation
    coefs = repmat([1 -2 1]/h^2,N-2,1); 
    %create sparse matrix for solver
    A = spdiags(coefs,-1:1,N-2,N-2); 
    B = zeros(1,N-2);
    B(1) = -Tl/h^2; %set left boundary condition
    B(end) = -Tr/h^2; %set right boundary condition
    
    %solve for temperature at internal nodes
    Tnumeric = [Tl B/A Tr]; 
    Tanalytic = Tl + (Tr-Tl)*x;
    
    fid=fopen('H1P3.txt','w');
    fprintf(fid,['ME 710 Homework 1, Problem 3\nCompleted by Michael Crawley on ',date,'\n\nx/L\t\tNumeric Soln\tAnalytic Soln\n']);
    fprintf(fid,'%1.2f\t\t %3.5f\t\t %3.5f\n',[x; Tnumeric; Tanalytic]);
    fclose(fid);
end