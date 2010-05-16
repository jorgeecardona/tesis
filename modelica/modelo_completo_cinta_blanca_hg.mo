
model PlantaHG
  //extends HGD(
  //	      period = period,

   extends HG(
 	      inputs = 2, outputs = 1, states =3,
	      
	      lambda = 1,
	      
	      K = [26; 240; 800],
	      
	      redeclare PlantaBase core,
	      
	      X_e(start = [0.1;0;0.1])
	      
	      );

  outer Params p;

  // Real dbeta1_dtheta;
  // Real dbeta2_dtheta;
  // Real dL_dtheta;

  // Real drho_dtheta;
  // Real drho_dM;

  // Real dT_dtheta;
  // Real dT_dM;

  // Real df_dtheta;
  // Real df_dM;

  // Real df2_dtheta;
  // Real df2_domega;
  // Real df2_dM;

equation


  // dbeta1_dtheta =  - p.l_1_ * cos(core.beta_1_aux + X_e[1,1] + p.theta_) / ((p.l_3 - p.l_2) * cos(core.beta_1_aux) + p.l_4 * sin(core.beta_1_aux) + p.l_1_ * cos(core.beta_1_aux + X_e[1,1] + p.theta_));

  // dbeta2_dtheta =   p.l_1_ * cos(core.beta_2_aux - X_e[1,1] - p.theta_) / (p.l_2 * cos(core.beta_2_aux) + p.l_4 * sin(core.beta_2_aux) - p.l_1_ * cos(core.beta_2_aux - X_e[1,1] - p.theta_));

  // dL_dtheta = 2 * p.r * ( -dbeta1_dtheta / tan(core.beta_1_aux)^2 - dbeta2_dtheta /tan(core.beta_2_aux)^2) + ( -dbeta1_dtheta * cos(core.beta_1_aux) / sin(core.beta_1_aux)^2 - dbeta2_dtheta * cos(core.beta_2_aux) / sin(core.beta_2_aux)^2)  * (p.l_4 - p.l_1_ * sin(X_e[1,1] + p.theta_)) + (1/sin(core.beta_1_aux) + 1/sin(core.beta_2_aux))*(-p.l_1_ * cos(X_e[1,1] + p.theta_));

  // drho_dtheta = - (X_e[3,1] / (p.A * core.L_aux^2)) * dL_dtheta;

  // drho_dM = 1 / (p.A * core.L_aux);

  // dT_dtheta = - (drho_dtheta * p.Y * 1e6 * p.A * p.rho_o *1000) / core.rho_aux^2;

  // dT_dM = - (drho_dM * p.Y * 1e6 * p.A * p.rho_o * 1000) / core.rho_aux^2;

  // df_dtheta = dT_dtheta * (sin(X_e[1,1] + p.theta_ + core.beta_1_aux) - sin(X_e[1,1] + p.theta_ - core.beta_2_aux)) + core.T_aux * (cos(X_e[1,1] + p.theta_+ core.beta_1_aux) * (1 + dbeta1_dtheta) - cos(X_e[1,1] + p.theta_ - core.beta_2_aux) * (1 - dbeta2_dtheta));

  // df_dM = dT_dM * (sin(X_e[1,1] + p.theta_ + core.beta_1_aux) - sin(X_e[1,1] + p.theta_ - core.beta_2_aux));

  // df2_dtheta = (1/(p.I/100 + p.l_1_^2 * p.m_1)) * (p.l_1_ * df_dtheta + p.l_1_ * p.m * p.g * sin(X_e[1,1] + p.theta_));

  // df2_domega = - (p.k_f/10) / (p.I/100 + p.l_1_^2 * p.m_1);

  // df2_dM = (p.l_1_ * df_dM) / (p.I/100 + p.l_1_^2 * p.m_1);

  // Jac = [1, 0, 0; 0, 1, 0; df2_dtheta , df2_domega, df2_dM];

  Jac[1,1] =  1 ;
  Jac[1,2] =  0 ;
  Jac[1,3] =  0 ;
  Jac[2,1] =  0 ;
  Jac[2,2] =  1 ;
  Jac[2,3] =  0 ;
  Jac[3,1] =  (p.g)*(p.l_1_)*(p.m)*sin((p.theta_) + X_e[1,1])/((p.I/100) + (p.m_1)*(p.l_1_)^2) + (p.A)*(p.Y*1e6)*(p.l_1_)*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.I/100) + (p.m_1)*(p.l_1_)^2) - (p.A)*(p.Y*1e6)*(p.l_1_)*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.I/100) + (p.m_1)*(p.l_1_)^2) + (p.L_k)*(p.Y*1e6)*(p.l_1_)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - (p.L_k)*(p.Y*1e6)*(p.l_1_)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) + (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) + (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*sin((p.theta_) + X_e[1,1])/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*sin((p.theta_) + X_e[1,1])/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) - (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*sin((p.theta_) + X_e[1,1])/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*sin((p.theta_) + X_e[1,1])/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_1_aux)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_2_aux)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*core.beta_1_aux*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 + (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)))*core.beta_2_aux*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)) - 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))) - 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_1_aux)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_2_aux)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*core.beta_1_aux*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*(1 - (p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux)))*core.beta_2_aux*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)) + 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))) + 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))) + (p.Y*1e6)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*cos(core.beta_2_aux)*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))*sin(core.beta_2_aux)^2) + (p.Y*1e6)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*cos(core.beta_1_aux)*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))*sin(core.beta_1_aux)^2) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^3*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*cos(core.beta_2_aux)*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))*sin(core.beta_2_aux)^2) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^3*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*cos(core.beta_1_aux)*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))*sin(core.beta_1_aux)^2) - (p.Y*1e6)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*cos(core.beta_2_aux)*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))*sin(core.beta_2_aux)^2) - (p.Y*1e6)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*cos(core.beta_1_aux)*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))*sin(core.beta_1_aux)^2) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^3*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*cos(core.beta_2_aux)*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))*sin(core.beta_2_aux)^2) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^3*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*cos(core.beta_1_aux)*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))*sin(core.beta_1_aux)^2) - 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 + tan(core.beta_2_aux)^2)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))*tan(core.beta_2_aux)^2) - 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 + tan(core.beta_1_aux)^2)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))*tan(core.beta_1_aux)^2) + 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 + tan(core.beta_2_aux)^2)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux)*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_2)*cos(core.beta_2_aux) + (p.l_4)*sin(core.beta_2_aux) - (p.l_1_)*cos((p.theta_) + X_e[1,1] - core.beta_2_aux))*tan(core.beta_2_aux)^2) + 2*(p.Y*1e6)*(p.r)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*(1 + tan(core.beta_1_aux)^2)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux)*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]*((p.I/100) + (p.m_1)*(p.l_1_)^2)*((p.l_1_)*cos((p.theta_) + X_e[1,1] + core.beta_1_aux) + (p.l_4)*sin(core.beta_1_aux) + ((p.l_3) - (p.l_2))*cos(core.beta_1_aux))*tan(core.beta_1_aux)^2) ;
    Jac[3,2] =  -(p.k_f/10)/((p.I/100) + (p.m_1)*(p.l_1_)^2) ;
    Jac[3,3] =  (p.L_k)*(p.Y*1e6)*(p.l_1_)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)) - (p.L_k)*(p.Y*1e6)*(p.l_1_)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)) + (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) + (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) + (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) - (p.Y*1e6)*(p.l_1_)*(p.l_4)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_1_aux)) - (p.Y*1e6)*(p.rho_o*1000)*(p.A)^2*(p.l_1_)^2*sin((p.theta_) + X_e[1,1])*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*sin(core.beta_2_aux)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_1_aux)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_2_aux)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*core.beta_1_aux*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)) - 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*core.beta_2_aux*sin((p.theta_) + X_e[1,1] + core.beta_1_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_1_aux)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)*tan(core.beta_2_aux)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*core.beta_1_aux*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)) + 2*(p.Y*1e6)*(p.l_1_)*(p.r)*(p.rho_o*1000)*(p.A)^2*core.beta_2_aux*sin((p.theta_) + X_e[1,1] - core.beta_2_aux)/(X_e[3,1]^2*((p.I/100) + (p.m_1)*(p.l_1_)^2)) ;
end PlantaHG;

model SandBox
  import Modelica.Blocks.Nonlinear;
  import Modelica.Blocks.Sources;
  
  Table1D vm(fileName="../20100125_vm.csv", offset=0);
  Table1D ve(fileName="../20100125_ve.csv", offset=0);
  Table1D theta(fileName="../20100125_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100125_ref.csv", offset=0);
  Table1D T(fileName="../20100125_T.csv", offset=0.36);

  inner Params p;

  PlantaHG K;

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
  K.U = [vm.y;ve.y];
  K.Y = [theta.y];

end SandBox;


model SandBoxI
  import Modelica.Blocks.Nonlinear;
  import Modelica.Blocks.Sources;

  Table1D vm(fileName="../20100206b_vm.csv", offset=0);
  Table1D ve(fileName="../20100206b_ve.csv", offset=0);
  Table1D theta(fileName="../20100206b_theta.csv", offset=0.36);
  Table1D ref(fileName="../20100206b_ref.csv", offset=0);
  Table1D T(fileName="../20100206b_T.csv", offset=0);
  
  inner Params p;

  PlantaHG K;

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
