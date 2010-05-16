


model SandBox
  import Modelica.Blocks.Sources;
  import Modelica.Blocks.Nonlinear;
  
  Table1D vm(fileName="../20100125_vm.csv", offset=0);
  Table1D ve(fileName="../20100125_ve.csv", offset=0);
  Table1D theta(fileName="../20100125_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100125_ref.csv", offset=0);
  Table1D T(fileName="../20100125_T.csv", offset=0.36);

  inner Params p;

  Luenberger
    L(

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

      K = 
      [25.0036167773897 ; 158.117735578911 ; -0.485143653721667],
            
      X(start=[1;1;1]),

      redeclare PlantaBase core
      );

  Planta P;

  Nonlinear.Limiter u(uMax={p.W_max,p.W_max}, uMin ={0,0});
  Sources.Pulse pulse(amplitude = {1}, period = {2}, offset = {-0.5});
  IndiceError indice(startTime = 50);

equation
  
  // Control PID para el motor esclavo
  vm.u = time;
  ve.u = time;
  theta.u = time;
  ref.u = time;
  T.u = time;

  u.u[1] = 5;
  u.u[2] = 5 + 19.1 * (ref.y - P.Y[1,1]);
//  u.u[2] = 5 + 19.1 * ((if time < 10 then 0 else pulse.y[1]) - P.Y[1,1]);

  P.U[1,1] = u.y[1];
  P.U[2,1] = u.y[2];

  // Conexiones al Observador
  L.U = [vm.y;ve.y];//P.U;
  L.Y = [theta.y];//P.Y;

  indice.u[1] = T.y;//P.core.T_aux;
  indice.u[2] = L.core.T_aux;

end SandBox;


model SandBoxI
  import Modelica.Blocks.Sources;
  import Modelica.Blocks.Nonlinear;
  
  Table1D vm(fileName="../20100206b_vm.csv", offset=0);
  Table1D ve(fileName="../20100206b_ve.csv", offset=0);
  Table1D theta(fileName="../20100206b_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100206b_ref.csv", offset=0);
  Table1D T(fileName="../20100206b_T.csv", offset=0.0);

  inner Params p;

  Luenberger
    L(
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

      K = 
      [25.0036167773897 ; 158.117735578911 ; -0.485143653721667],
            
      X(start=[1;1;1]),

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
  L.U[1,1] = vm.y;
  L.U[2,1] = ve.y;
  L.Y[1,1] = theta.y;

  indice.u[1] = T.y;
  indice.u[2] = L.core.T_aux;

end SandBoxI;



model SandBoxD
  import Modelica.Blocks.Sources;
  import Modelica.Blocks.Nonlinear;
  
  Table1D ref(fileName="../20100125_ref.csv", offset=0);
  Table1D T(fileName="../20100125_T.csv", offset=0.36);

  inner Params p(k_f = 20, I=30);

  PlantaD P;
equation
  
  // Control PID para el motor esclavo
  ref.u = time;
  T.u = time;


  P.U[1,1] = 5;
  P.U[2,1] = 5 + 19.1 * (ref.y - P.Y[1,1]);

end SandBoxD;

model SandBoxI
  import Modelica.Blocks.Sources;
  import Modelica.Blocks.Nonlinear;
  
  Table1D vm(fileName="../20100125_vm.csv", offset=0);
  Table1D ve(fileName="../20100125_ve.csv", offset=0);
  Table1D theta(fileName="../20100125_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100125_ref.csv", offset=0);
  Table1D T(fileName="../20100125_T.csv", offset=0.36);

  inner Params p;

  Luenberger 
    L(
      
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

      K = 
      [25.0036167773897 ; 158.117735578911 ; -0.485143653721667],
            
      X(start=[1;1;1]),

      redeclare PlantaBase core
      );

  Planta P;

  parameter Real startTime = 30;

  Nonlinear.Limiter u(uMax={p.W_max,p.W_max}, uMin ={0,0});
  Sources.Pulse pulse(amplitude = {0.7}, period = {8}, offset = {-0.35}, startTime={startTime-10});
  IndiceError indice(startTime = startTime);

equation
  
  // Control PID para el motor esclavo
  vm.u = time;
  ve.u = time;
  theta.u = time;
  ref.u = time;
  T.u = time;

  u.u[1] = 5;
  u.u[2] = 5 + 19.1 * ((if time < startTime-10 then 0 else pulse.y[1]) - P.Y[1,1]);

  P.U[1,1] = u.y[1];
  P.U[2,1] = u.y[2];

  // Conexiones al Observador
  L.U = P.U;
  L.Y = P.Y;

  indice.u[1] = P.core.T_aux;
  indice.u[2] = L.core.T_aux;

end SandBoxI;
