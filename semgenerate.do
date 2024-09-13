//----------------------------------------------------------------------------//
// Title  : SEM simulation program
// Purpose: Program for simulating data using structural equations with latent variables 
// Author : Steffen Erickson
// Date   : 6/10/24
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
// Generates Structural Equation Models using LISERAL notation and 
// Equations from Bollen, 1989
/*
	/* Inputs */ 
	Endogenous 
	matrix Lambda_y...................(p x m) measurement loading matrix 
	matrix Theta_y....................(p x p) measurement error matrix
	matrix Beta.......................(m x m) latent structural coefficients
	matrix Psi........................(m x m) latent variance/covariance
	vector alpha .....................(m x 1) latent mean vector 
	
	Exogenous
	matrix Lambda_x...................(q x n) measurement loading matrix 
	matrix Theta_x....................(q x q) measurement error matrix 
	matrix Gamma......................(m x n) latent structural coefficients
	matrix Phi........................(n x n) latent variance/covariance
	vector kappa......................(n x 1) latent mean vector 
	
	/* Outputs */ 
	Matrix Sigma......................(m x n) observed covariance matrix 
	vector upsilon_y..................(p x 1) observed endogenous mean vector 
	vector upsilon_x..................(q x 1) observed exogenous mean vector 


*/
//----------------------------------------------------------------------------//

mata 
// ---------------- Define Structure -----------------------------------------//
struct myproblem  {
	real matrix 	Lambda_y 
	real matrix 	Theta_y  
	real matrix 	Beta     
	real matrix 	Psi      
	real matrix 	Lambda_x 
	real matrix 	Theta_x  
	real matrix 	Gamma    
    real matrix 	Phi  
	real colvector 	kappa
	real colvector 	alpha 
	struct derived 	scalar d
}
struct derived { 
	real matrix 	Sigma 
	real matrix 	YY 
	real matrix 	XX 
	real matrix 	YX
	real matrix 	XY 
	real colvector  upsilon_y
	real colvector  upsilon_x
	real matrix     IminusBeta 
}
// ---------------- Population Covariance Matrix -----------------------------//			
struct derived generate_sem_cov(real matrix Lambda_x,
								real matrix Theta_x, 
								real matrix Phi, 
								| real matrix Lambda_y,
								real matrix Theta_y,
								real matrix Beta,
								real matrix Gamma, 
								real matrix Psi) 
{
	struct myproblem scalar pr
	initialize_objects_cov(pr.d)
	
	pr.Lambda_x = Lambda_x 
	pr.Theta_x  = Theta_x  
	pr.Phi      = Phi 
	pr.Lambda_y = Lambda_y 
	pr.Theta_y  = Theta_y  
	pr.Beta     = Beta      
	pr.Gamma    = Gamma
	pr.Psi      = Psi     
	pr.d.IminusBeta = I(rows(pr.Beta),rows(pr.Beta))- pr.Beta
	
	if (args() <= 3) {
		get_xx(pr)
		pr.d.Sigma = pr.d.XX
		_makesymmetric(pr.d.Sigma)
	}
	else {
		get_yy(pr)
		get_xx(pr)
		get_yx(pr)
		get_xy(pr)
		pr.d.Sigma = (pr.d.YY,pr.d.YX\pr.d.XY,pr.d.XX)
		_makesymmetric(pr.d.Sigma)
	}
	return(pr.d)
}

void initialize_objects_cov(struct derived scalar d)
{
	d.Sigma = J(0,0,.) 
	d.YY = J(0,0,.) 
	d.XX = J(0,0,.)
	d.YX = J(0,0,.) 
	d.XY = J(0,0,.) 
}

void get_xx(struct myproblem scalar pr)
{
  pr.d.XX = pr.Lambda_x*pr.Phi*pr.Lambda_x' + pr.Theta_x
}
void get_yy(struct myproblem scalar pr)
{
	pr.d.YY = pr.Lambda_y*luinv(pr.d.IminusBeta)*(pr.Gamma*pr.Phi*pr.Gamma' + pr.Psi)*luinv(pr.d.IminusBeta)'*pr.Lambda_y' + pr.Theta_y
}
void get_yx(struct myproblem scalar pr)
{
	pr.d.YX = pr.Lambda_y*luinv(pr.d.IminusBeta)*pr.Gamma*pr.Phi*pr.Lambda_x'
}
void get_xy(struct myproblem scalar pr) // the transpose of yx
{
	pr.d.XY = pr.Lambda_x*pr.Phi*pr.Gamma'*luinv(pr.d.IminusBeta)'pr.Lambda_y'
}

// ---------------- Population Mean Vector -----------------------------------//
struct derived generate_sem_mean(real matrix Lambda_x,
								 real vector kappa,
								| real matrix Lambda_y,
								real matrix Beta,
								real matrix Gamma, 
								real vector alpha) 
{
	struct myproblem scalar pr
	initialize_objects_mean(pr.d)
	row_to_col(kappa)
	row_to_col(alpha)
	
	pr.Lambda_x 	= Lambda_x 
	pr.kappa 		= kappa 
	pr.Lambda_y 	= Lambda_y 
	pr.Beta 		= Beta     
	pr.Gamma 		= Gamma   
	pr.alpha 		= alpha 
	pr.d.IminusBeta = I(rows(pr.Beta),rows(pr.Beta))- pr.Beta
	
	if (args() <= 2) {
		get_upsilon_x(pr)
	}
	else {
		get_upsilon_x(pr)
		get_upsilon_y(pr)
	}
	return(pr.d)
}

void initialize_objects_mean(struct derived scalar d)
{
	d.upsilon_y = J(0,1,.)
	d.upsilon_x = J(0,1,.)
}
void get_upsilon_x(struct myproblem scalar pr)
{
	pr.d.upsilon_x = pr.Lambda_x*pr.kappa
}
void get_upsilon_y(struct myproblem scalar pr)
{
	pr.d.upsilon_y = pr.Lambda_y*luinv(pr.d.IminusBeta)*(pr.alpha + pr.Gamma*pr.kappa)
}

real scalar is_colvec(z) return(anyof(("colvector","scalar"),orgtype(z)))
void row_to_col(real vector v) 
{
	if (is_colvec(v) == 0) v = v'
}
end 
