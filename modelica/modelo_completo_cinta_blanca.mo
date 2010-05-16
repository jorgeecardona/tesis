

// Parameters
class Params
  import Modelica.Math.*;

  parameter Real X_o[3,1] = [0.02824; 0; 0.13233];
  parameter Real U_o[2,1] = [5; 4.46];
  parameter Real Y_o[1,1] = [0.02824];

  parameter Real Y 	  = 3.744;
  parameter Real A 	  = 31e-6;
  parameter Real rho_e 	  = 1.041;
  parameter Real rho_o 	  = 1.2157;
  parameter Real L_k   	  = 2.50;
  parameter Real r   	  = 0.05;

  parameter Real l_1   	  = 0.27;
  parameter Real l_2   	  = 0.384;
  parameter Real l_3   	  = 0.41;
  parameter Real l_4   	  = 0.50;
  parameter Real l_5   	  = 0.092;

  parameter Real l_1_  	  = sqrt(l_1^2 + l_5^2);
  parameter Real theta_   = atan(l_5 / l_1);

  parameter Real g   	  = 9.81;
  parameter Real m        = 1.013;
  parameter Real m_1 	  = 0.5;

  parameter Real k_f 	  = 8.9182;
  parameter Real I	  = 13.9720;

  parameter Real W_max 	  = 10;
end Params;

// Functions 
function rho
  input Real L, M;
  input Real A;
  
  output Real rho;
algorithm
  rho := M / (A * L);
//  rho := if(rho < 10) then 10 else rho;
end rho;


function L
  
  input Real beta_1, beta_2, theta;
  input Real L_k, r, l_1_, l_2, l_3, l_4, theta_;
  
  output Real L;
algorithm
  L := L_k + 2 * r * (beta_1 + 1/tan(beta_1) + beta_2 + 1/tan(beta_2)) + (1/sin(beta_1) + 1/sin(beta_2)) * (l_4 - l_1_*sin(theta + theta_));
//  L := if (L < 0) then 0 else L; 
end L;

function T
  
  input Real rho;
  input Real A, Y, rho_o; 
  
  output Real T;
algorithm
  T := 1000000 * Y * A * (1000 * rho_o / rho - 1);
//  T := if (T < 0) then 0 else T;
end T;


function beta_1
  
  input Real theta;
  input Real r, l_1_, l_2, l_3, l_4, theta_;
  
  output Real beta_1;

  Real a, b, c, F_a, F_b, F_c;

algorithm

  a := 0;

  b := 3.14159265;

  F_a := 
    (l_3 - l_2) * sin(a) - l_4 * cos(a) + l_1_ * sin(a + theta + theta_) - 2*r;
    
  F_b := 
    (l_3 - l_2) * sin(b) - l_4 * cos(b) + l_1_ * sin(b + theta + theta_) - 2*r;

  while abs(a - b) > 1e-10 loop
    
    c := (a + b)/2;
    
    F_c :=
      (l_3 - l_2) * sin(c) - l_4 * cos(c) + l_1_ * sin(c + theta + theta_) - 2*r;
    
    if (F_a * F_c >=0) then
      a := c;
      F_a := F_c;
    elseif (F_b * F_c >=0) then
	b := c;
      F_b := F_c;
    else
      terminate("Error in algorithm.")
      break;
    end if;

  end while;
  
  beta_1 := (a + b)/2;

end beta_1;


function beta_2

  input Real theta;
  input Real r, l_1_, l_2, l_3, l_4, theta_;

  output Real beta_2;

  Real a, b, c, F_a, F_b, F_c;

algorithm

  a := 0;

  b := 3.14159265;

  F_a := 
    l_2 * sin(a) - l_4 * cos(a) - l_1_ * sin(a - theta - theta_) - 2*r;
    
  F_b := 
    l_2 * sin(b) - l_4 * cos(b) - l_1_ * sin(b - theta - theta_) - 2*r;

  while abs(a - b) > 1e-10 loop
    
    c := (a + b)/2;

    F_c :=
      l_2 * sin(c) - l_4 * cos(c) - l_1_ * sin(c - theta - theta_) - 2*r;
    
    if (F_a * F_c >=0) then
      a := c;
      F_a := F_c;
    elseif (F_b * F_c >=0) then
      b := c;
      F_b := F_c;
    else
      terminate("Error in algorithm.")
      break;
    end if;

  end while;
  
  beta_2 := (a + b)/2;

end beta_2;



model PlantaBase
  extends NonLinearBase(X(start = p.X_o));

  Real beta_1_aux, beta_2_aux, rho_aux, T_aux, L_aux;
  
  outer Params p;

equation

  // State equations
  F[1,1] = X[2,1];

  F[2,1] = (p.l_1_ * T_aux * (sin(X[1,1] + p.theta_ + beta_1_aux) - sin(X[1,1] + p.theta_- beta_2_aux)) - p.l_1_ * p.m * p.g * cos(X[1,1] + p.theta_) - abs(p.k_f)/10 * X[2,1]) / (abs(p.I)/100 + p.l_1_^2  * p.r);

  F[3,1] = p.A * p.r * (1000 * p.rho_e * U[1,1] - rho_aux * U[2,1]);

  // Output equation
  H[1,1] = X[1,1];

  // Auxiliary equations
  beta_1_aux = beta_1(X[1,1], p.r, p.l_1_, p.l_2, p.l_3, p.l_4, p.theta_);
  beta_2_aux = beta_2(X[1,1], p.r, p.l_1_, p.l_2, p.l_3, p.l_4, p.theta_);
  T_aux    = T(rho_aux, p.A, p.Y, p.rho_o);
  rho_aux  = rho(L_aux, X[3,1], p.A);
  L_aux    = L(beta_1_aux, beta_2_aux, X[1,1], p.L_k, p.r, p.l_1_, p.l_2, p.l_3, p.l_4, p.theta_);

end PlantaBase;

model PlantaRuidoBase
  extends NonLinearBase(X(start = p.X_o));

  Real beta_1_aux, beta_2_aux, rho_aux, T_aux, L_aux;
  
  outer Params p;

  Normal N_1(seed=1);
  Normal N_2(seed=2);
  Normal N_3(seed=3);
  Normal N_4(seed=4);

equation

  // State equations
  F[1,1] = N_1.y + X[2,1];

  F[2,1] = N_2.y + (p.l_1 * T_aux * (sin(X[1,1] + beta_1_aux) - sin(X[1,1] - beta_2_aux)) - p.l_1 * p.M * p.g * cos(X[1,1]) - p.k_f * X[2,1]) / (p.I + p.l_1^2 * p.m_1);

  F[3,1] = N_3.y + p.A * p.r * (1000 * p.rho_e * U[1,1] - rho_aux * U[2,1]);


  // Output equation
  H[1,1] = N_4.y + X[1,1];


  // Auxiliary equations
  beta_1_aux = beta_1(pre(beta_1_aux), X[1,1], p.r, p.l_1, p.l_2, p.l_3, p.l_4);
  beta_2_aux = beta_2(pre(beta_2_aux), X[1,1], p.r, p.l_1, p.l_2, p.l_3, p.l_4);
  T_aux    = T(rho_aux, p.A, p.Y, p.rho_o);
  rho_aux  = rho(L_aux, X[3,1], p.A);
  L_aux    = L(beta_1_aux, beta_2_aux, X[1,1], p.L_k, p.r, p.l_1, p.l_2, p.l_3, p.l_4);

end PlantaRuidoBase;


model Planta
  extends NonLinear(inputs = 2,
		    states = 3,
		    outputs = 1,
		    X(start = p.X_o),
		    redeclare PlantaBase core);

  outer Params p;
end Planta;

model PlantaD
  extends NonLinearDiscrete(inputs = 2,
			    states = 3,
			    outputs = 1,
			    periodTime = 0.09,
			    X(start = p.X_o),
			    redeclare PlantaBase core);

  outer Params p;
end PlantaD;


model PlantaRuido
  extends NonLinear(inputs = 2,
		    states = 3,
		    outputs = 1,
		    X(start = p.X_o),
		    redeclare PlantaRuidoBase 
		    core(
			 N_1(d=dev[1,1]),
			 N_2(d=dev[2,1]),
			 N_3(d=dev[3,1]),
			 N_4(d=dev[4,1])));

  parameter Real dev[states + outputs,1] = ones(states + outputs,1);
  outer Params p;
end PlantaRuido;

model SandBox
  import Modelica.Blocks.Sources;

  inner Params p;
  Planta P;

  Sources.Step 
    s1(startTime = 0), 
    s2(startTime = 100), 
    s3(startTime = 200), 
    s4(startTime = 300);

  PID pid(K_p = -15, K_i = -40);

  Limiter 
    Lim_1(umin = 0, umax = p.W_max), 
    Lim_2(umin = 0, umax = p.W_max);

equation

  // Entrada del motor maestro
  P.U[1,1] = Lim_1.y;
  
  Lim_1.u = 
    (p.W_max/2) * s1.y + 
    (p.W_max/3) * s2.y - 
    (2 * p.W_max/3) * s3.y + 
    (p.W_max/3) * s4.y;

  // Control PID para el motor esclavo
  P.U[2,1] = Lim_2.y;
  Lim_2.u = pid.y;
  pid.u = P.Y[1,1] - p.Y_o[1,1];


end SandBox;