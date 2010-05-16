
partial model HG

  parameter Integer inputs, outputs, states;

  // Input
  Real U[inputs,1];

  // Real output
  Real Y[outputs,1];

  // Estimations
  Real X_e[states,1];
  Real Y_e[outputs,1];

  // Core Block
  replaceable NonLinearBase
    core(inputs  = inputs,
	 outputs = outputs,
	 states  = states)
    extends NonLinearBase;

  // Jacobian
  Real Jac[states, states];
  Real Jac_inv[states, states];
  
  // Delta and K parameter
  parameter Real lambda;
  Real Delta[states, states];
  parameter Real K[states,1];

equation
  
  der(X_e) = core.F  +  Jac_inv * Delta * K * (Y - Y_e);

  for i in 1:size(K,1) loop
    for j in 1:size(K,1) loop
      Delta[i,j] = if i==j then lambda ^ i else 0;
    end for;
  end for;

  identity(size(K,1)) = Jac * Jac_inv;

  core.U = U;
  core.X = X_e;
  core.H = Y_e;

end HG;

partial model HGD

  parameter Integer inputs, outputs, states;

  // Input
  Real U[inputs,1];

  // Real output
  Real Y[outputs,1];

  // Estimations
  Real X_e[states,1];
  Real Xd[states,1];
  Real Y_e[outputs,1];

  // Core Block
  replaceable NonLinearBase
    core(inputs  = inputs,
	 outputs = outputs,
	 states  = states)
    extends NonLinearBase;

  // Jacobian
  Real Jac[states, states];
  Real Jac_inv[states, states];
  
  // Delta and K parameter
  parameter Real lambda;
  Real Delta[states, states];
  parameter Real K[states,1];

  Real NextSample(start = 0);
  parameter Real period = 0.01;

equation
  
  when (NextSample > 1) then
    reinit(X_e, Xd);
    reinit(NextSample, 0);
  end when;
  Xd = pre(X_e) + period * (core.F  +  Jac_inv * Delta * K * (Y - Y_e));
  der(X_e) = zeros(size(K,1),1);
  der(NextSample) = 1/period;

  for i in 1:size(K,1) loop
    for j in 1:size(K,1) loop
      Delta[i,j] = if i==j then lambda ^ i else 0;
    end for;
  end for;

  identity(size(K,1)) = Jac * Jac_inv;

  core.U = U;
  core.X = X_e;
  core.H = Y_e;

end HGD;

