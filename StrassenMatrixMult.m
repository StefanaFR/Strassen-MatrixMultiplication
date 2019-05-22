% CSCI445
% Stefana Rusu
% Assignment:

% A function that computes the product of two matrices A & B 
% using Strassen's Matrix Multiplication algorithm.
% Input: 2 matrices (2^n x 2^n)
% Output: 1 matrix (2^n x 2^n)

% Test Data: 4x4 matrix
% In the run command line, type: StrassenMatrixMult([21 96 75 74; 44 58 1 30; 32 87 40 85; 100 7 7 21], [22 29 97 12; 21 92 54 13; 86 25 82 12; 38 47 26 77])

function matrixC = StrassenMatrixMult(matrixA, matrixB)
    n = length(matrixA); % find the length of the largest array dimension
    
    % If the matrices are of size 1x1
    if (n == 1)
        % Perform normal multiplication on the matrices 
        matrixC = matrixA * matrixB;
        
    % Otherwise
    else
        n = n/2;
       
        %Divide matrix A into an n/2 x n/2 sub-matrix
        A11 = matrixA(1:n, 1:n); 
        A12 = matrixA(1:n, (n+1):end);
        A21 = matrixA((n+1):end, 1:n);
        A22 = matrixA((n+1):end, (n+1):end);
        
        %Divide matrix B into an n/2 x n/2 sub-matrix
        B11 = matrixB(1:n, 1:n);
        B12 = matrixB(1:n, (n+1):end);
        B21 = matrixB((n+1):end, 1:n);
        B22 = matrixB((n+1):end, (n+1):end);
        
        % Compute the 7 multiplications - I used a variation of this
        % computation different from the book
        P1 = StrassenMatrixMult(A11 + A22, B11 + B22);
        P2 = StrassenMatrixMult(A21 + A22, B11);
        P3 = StrassenMatrixMult(A11, B12 - B22);
        P4 = StrassenMatrixMult(A22, B21 - B11);
        P5 = StrassenMatrixMult(A11 + A12, B22);
        P6 = StrassenMatrixMult(A21 - A11, B11 + B12);
        P7 = StrassenMatrixMult(A12 - A22, B21 + B22);

        % Compute the four submatrices that will go back into C
        C11 = P1 + P4 - P5 + P7;
        C12 = P3 + P5;
        C21 = P2 + P4;
        C22 = P1 - P2 + P3 + P6;
        matrixC = [ C11, C12; C21, C22 ]; % Matrix C is returned
        
    end
    
end