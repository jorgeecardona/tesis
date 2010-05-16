block Luenberger

  parameter Real A[:,size(A,1)];
  parameter Real C[1,size(A,1)];
  parameter Real B[size(A,1),:];

  parameter Real K[size(A,1),1];

  input Real X[size(A,1),1];
  input Real U[size(B,2),1];

  // Real output
  input Real Y[1,1];

  output Real X_e[size(A,1),1];
  output Real Y_e[1,1];

  replaceable NonLinearBase
    core(inputs = size(B,2), 
	 outputs = size(C,1), 
	 states = size(A,1)) 
    extends NonLinearBase;

equation

  der(X) = A * X + B * (U - core.p.U_o) + K * (Y -  Y_e);

  X_e = X + core.p.X_o;
  Y_e = C * X + core.p.Y_o;

  // No es parte del observador pero sirve para comparar
  core.U = U;
  core.X = X_e;

end Luenberger;

block LuenbergerD

  parameter Real period = 0.01;

  parameter Real A[:,size(A,1)];
  parameter Real C[1,size(A,1)];
  parameter Real B[size(A,1),:];

  parameter Real K[size(A,1),1];

  input Real X[size(A,1),1];
  input Real Xd[size(A,1),1];
  input Real U[size(B,2),1];

  // Real output
  input Real Y[1,1];

  output Real X_e[size(A,1),1];
  output Real Y_e[1,1];

  replaceable NonLinearBase
    core(inputs = size(B,2), 
	 outputs = size(C,1), 
	 states = size(A,1)) 
    extends NonLinearBase;

  Real NextSample(start = 0);

equation

  when (NextSample > 1) then
    reinit(X, Xd);
    reinit(NextSample, 0);
  end when;

  Xd = pre(X) + period * (A * X + B * (U - core.p.U_o) + K * (Y -  Y_e));
  der(X) = zeros(size(A,1),1);
  der(NextSample) = 1/period;

  X_e = X + core.p.X_o;
  Y_e = C * X + core.p.Y_o;

  // No es parte del observador pero sirve para comparar
  core.U = U;
  core.X = X_e;

end LuenbergerD;
