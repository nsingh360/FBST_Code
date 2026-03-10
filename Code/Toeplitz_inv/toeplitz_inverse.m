function B = toeplitz_inverse ( n, x )

%*****************************************************************************80
%
%% toeplitz_inverse() computes the inverse of a Toeplitz matrix.
%
%  Discussion:
%
%    This function uses the fact that if T is a Toeplitz matrix,
%    and J is the "exchange" matrix, then H=J*T is a Hankel matrix.
%
%    By lucky chance, the vector x defining T is the same vector x
%    used to define J*T.
%
%    Hence, we can use a known algorithm for the inverse of a Hankel matrix.
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    16 May 2020
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Miroslav Fiedler,
%    Hankel and Loewner Matrices
%
%  Input:
%
%    integer N, the order of the matrix.
%
%    real X(2*N-1), the vector that defines the matrix.
%
%  Output:
%
%    real B(N,N), the inverse matrix.
%
  x = x(:);
%
%  Get the Exchange matrix J.
%
  J = exchange_matrix ( n );
%
%  Define the Toeplitz matrix A.
%
  A = toeplitz_matrix ( n, x );
%
%  Define the corresponding Hankel matrix H = J * A.
%
  H = J * A;
%
%  Solve two linear systems.
%
  p = [ x(n+1:2*n-1); 0.0 ];
  u = H \ p;

  q = zeros ( n, 1 );
  q(n) = 1.0;
  v = H \ q;
%
%  Construct four matrices.
%
  z1 = zeros ( n, 1 );
  w1 = [ v(2:n); z1 ];
  M1 = hankel_matrix ( n, w1 );

  z2 = zeros ( n-1, 1 );
  w2 = [ z2; u ];
  M2 = toeplitz_matrix ( n, w2 );

  z3 = zeros ( n, 1 );
  z3(1) = -1.0;
  w3 = [ u(2:n); z3 ];
  M3 = hankel_matrix ( n, w3 );

  z4 = zeros ( n-1, 1 );
  w4 = [ z4; v  ];
  M4 = toeplitz_matrix ( n, w4 );
%
%  Construct K, the inverse of the Hankel matrix.
%
  K = ( M1 * M2 - M3 * M4 );
%
%  Compute B, the inverse of the Toeplitz matrix.
%
  B = K * J;

  return
end

