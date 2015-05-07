#' Modified Wilcoxon rank test for zero-inflated data
#'
#' @param X the covariates for zero proportion
#' @param Z the covariates for non-zero compoment
#' @param alpha the regression coeffcient
#' @param beta the regression coeffcient
#' @return test.table, eta, b, phi
#' @export
#' @examples
#' likelihood_ratio_test_mada(Data,Covariates,total.counts,marker.length)



Wilcoxon.zeros <- function(x, y, alpha = 0.05, perm = FALSE) {
	
	N = c(length(x), length(y));

	###########Our Test ###################
	n = c(sum(x != 0), sum(y != 0));
	prop = n/N; pmax = max(prop); pmean = mean(prop);
	Ntrun = round(pmax*N);
	Xtrun = c(x[x != 0], rep(0, Ntrun[1] - n[1]));
	Ytrun = c(y[y != 0], rep(0, Ntrun[2]- n[2]));
	rankdata = sum(Ntrun) + 1 - rank(c(Xtrun, Ytrun));
	r = sum(rankdata[1:Ntrun[1]]);
	s = r - Ntrun[1]*(sum(Ntrun) + 1)/2;

	v1 = pmean*(1 - pmean);
	var1 = N[1]*N[2]*pmean*v1*(sum(N)*pmean + 3*(1 - pmean)/2) + sum(N^2)*v1^2*5/4 + 2*sqrt(2/pi)*pmean*(sum(N)*v1)^(1.5)*sqrt(N[1]*N[2]);
	var1 = var1/4;
	var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
	vars = var1 + var2;
	w = s/sqrt(vars);

	if(perm == FALSE) {p = 2*(1 - pnorm(abs(w)));}
	else{
		numrep = 10000; permu.w = rep(0, numrep)
		Z = c(x, y);
		for(i in 1:numrep){
			ind = sample(sum(N), N[1]);
			x = Z[ind]; y = Z[-ind];
			n = c(sum(x != 0), sum(y != 0));
			prop = n/N; pmax = max(prop); pmean = mean(prop);
			Ntrun = round(pmax*N);
			Xtrun = c(x[x != 0], rep(0, Ntrun[1] - n[1]));
			Ytrun = c(y[y != 0], rep(0, Ntrun[2]- n[2]));
			rankdata = sum(Ntrun) + 1 - rank(c(Xtrun, Ytrun));
			r = sum(rankdata[1:Ntrun[1]]);
			s = r - Ntrun[1]*(sum(Ntrun) + 1)/2;

			v1 = pmean*(1 - pmean);
			var1 = N[1]*N[2]*pmean*v1*(sum(N)*pmean + 3*(1 - pmean)/2) + sum(N^2)*v1^2*5/4 + 2*sqrt(2/pi)*pmean*(sum(N)*v1)^(1.5)*sqrt(N[1]*N[2]);
			var1 = var1/4;
			var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
			vars = var1 + var2;
			permu.w[i] = s/sqrt(vars);
		}
		p = sum(abs(w) < abs(permu.w))/numrep;	
	}
	H = (p < alpha);

	result <- list(H = H, p.value = p, statistics=w)
	return(result)
}
