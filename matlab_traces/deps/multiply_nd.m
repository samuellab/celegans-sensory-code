function Y = multiply_nd(A, X, d)
% Y = MULTIPLY_ND(A, X, d)
%
%   Performs matrix multiplication between (M x K) 2D array A and
%   ND array X. Dimension d of X must have size K, and the result Y will
%   have the same size of X, with the exception of dimension d which is now
%   size M instead of K.
%
%   A = randn(6,4);
%   X = randn(2,3,4,5);
%   Y = multiply_ND(A, X, 3);
%
%   size(Y) % [2, 3, 8, 4]
%
% Y = MULTIPLY_ND(A, X)
%
%   If A is MxN, the transformation will be applied to the first dimension
%   of X that has size N.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

SA = size(A);
SX = size(X);

N = ndims(X);
M = SA(1);
K = SA(2);

type = element_class(X);

if nargin < 3
    d = find(SX == K, 1);
end

X_float = double(X);
A_float = double(A);

% Reshape X to prepare for standard matrix multiplication.
fwd_permutation = [d, 1:d-1, d+1:N];
X_permuted = permute(X_float, fwd_permutation);
X_2D = reshape(X_permuted, K, []);

% Do the matrix multiplication.
Y_2D = A_float * X_2D;

% Now reshape Y to match the original dimensions.
[~, rev_permutation] = sort(fwd_permutation);
Y_permuted = reshape(Y_2D, [M, SX([1:d-1, d+1:N])]);
Y_float = permute(Y_permuted, rev_permutation);

Y = cast(Y_float, type);