


model SandBox
  import Modelica.Blocks.Sources;
  import Modelica.Blocks.Nonlinear;
  
  Table1D vm(fileName="../20100125_vm.csv", offset=0);
  Table1D ve(fileName="../20100125_ve.csv", offset=0);
  Table1D theta(fileName="../20100125_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100125_ref.csv", offset=0);
  Table1D T(fileName="../20100125_T.csv", offset=0.36);

  inner Params p;

  EKF
    K(
      A=
      [                   0,              1.0,                 0;
       -52.6532512081627, -4.93532096464933, -2768.15571542308;
       -0.00121325046945302,                 0, -0.06106225795865],
      
      B=
      [         0,                    0;
       0,                    0;
       0.00161355, -0.00181174183759376],

      C=
      [1.0, 0, 0],

      V_1 = 
      [3.48945380789953, 0, 0;
       0, 189.496757341876, 0;
       0, 0, 0.00175590859584379],

      V_2 = 
      [0.0100000000000000],

      X_e(start=[0;0;0.1]),

      redeclare PlantaBase core
      );

  Planta P;

  Nonlinear.Limiter u(uMax={p.W_max,p.W_max}, uMin ={0,0});
equation
  
  // Control PID para el motor esclavo
  vm.u = time;
  ve.u = time;
  theta.u = time;
  ref.u = time;
  T.u = time;

  u.u[1] = 5;
  u.u[2] = 5 + 19.1 * (ref.y - P.Y[1,1]);

  P.U[1,1] = u.y[1];
  P.U[2,1] = u.y[2];

  // Conexiones al Observador
  K.U = P.U;
  K.Y = P.Y;

end SandBox;


model SandBoxI
  import Modelica.Blocks.Sources;
  import Modelica.Blocks.Nonlinear;
  
  Table1D vm(fileName="../20100206b_vm.csv", offset=0);
  Table1D ve(fileName="../20100206b_ve.csv", offset=0);
  Table1D theta(fileName="../20100206b_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100206b_ref.csv", offset=0);
  Table1D T(fileName="../20100206b_T.csv", offset=0);

  inner Params p;

  EKF
    K(
      A=
      [                   0,              1.0,                 0;
       -52.6532512081627, -4.93532096464933, -2768.15571542308;
       -0.00121325046945302,                 0, -0.06106225795865],
      
      B=
      [         0,                    0;
       0,                    0;
       0.00161355, -0.00181174183759376],

      C=
      [1.0, 0, 0],

      V_1 = 
      [3.48945380789953, 0, 0;
       0, 189.496757341876, 0;
       0, 0, 0.00175590859584379],

      V_2 = 
      [0.0100000000000000],

      X_e(start=[0;0;0.1]),

      redeclare PlantaBase core
      );

  IndiceError indice(startTime = 10);

equation
  
  // Control PID para el motor esclavo
  vm.u = time;
  ve.u = time;
  theta.u = time;
  ref.u = time;
  T.u = time;

  // Conexiones al Observador
  K.U = [vm.y;ve.y];
  K.Y = [theta.y];

  indice.u[1] = T.y;
  indice.u[2] = K.core.T_aux;

end SandBoxI;