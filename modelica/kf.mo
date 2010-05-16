

block KF

  parameter Real A[:,size(A,1)];
  parameter Real B[size(A,1),:];
  parameter Real C[1,size(A,1)];

  // Noise
  parameter Real V_2[1,1];
  Real V_2_inv[1,1];
  parameter Real V_1[size(A,1),size(A,1)];

  Real H[size(A,1),1];
  Real P[size(A,1),size(A,1)](start=P_o);

  input Real X[size(A,1),1];
  input Real U[size(B,2),1];

  // Real output
  input Real Y[1,1];

  parameter Real P_o[size(A,1),size(A,1)] = zeros(size(A,1),size(A,1));

  output Real X_e[size(A,1),1];
  output Real Y_e[1,1];

  replaceable NonLinearBase
    core(inputs = size(B,2),
	 outputs = size(C,1), 
	 states = size(A,1)) 
    extends NonLinearBase;

equation

  der(X) = A * X + B * (U - core.p.U_o) + H * (Y - core.p.Y_o - C * X);

  H = P * transpose(C) * V_2_inv;
  der(P) = P * transpose(A) + A * P - P * transpose(C) * V_2_inv * C * P  + V_1;
  identity(1) = V_2 * V_2_inv;

  core.U = U;
  core.X = X_e;

  X_e = X + core.p.X_o;
  Y_e = C * X + core.p.Y_o;

end KF;
