function A = toeplitz_matrix ( n, x )

%*****************************************************************************80
%
%% toeplitz_matrix() returns a TOEPLITZ matrix.
%
%  Formula:
%
%    A(I,J) = X(N+J-I)
%
%  Example:
%
%    N = 5, X = ( 1, 2, 3, 4, 5, 6, 7, 8, 9 )
%
%    5  6  7  8  9
%    4  5  6  7  8
%    3  4  5  6  7
%    2  3  4  5  6
%    1  2  3  4  5
%
%  Properties:
%
%    A is generally not symmetric: A' /= A.
%
%    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
%
%    A is Toeplitz: constant along diagonals.
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    21 October 2007
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer N, the order of A.
%
%    real X(2*N-1), the diagonals of A, with X(1) being
%    the A(N,1) entry, X(N) being the main diagonal value of A,
%    and X(2*N-1) being the A(1,N) entry.
%
%  Output:
%
%    real A(N,N), the matrix.
%
  A = zeros ( n, n );

  for i = 1 : n
    A(i,1:n) = x(n-i+1:n-i+n);
  end

  return
end
