#' Modified Kruskal Wallis test for zero-inflated data
#'
#' @param X the covariates for zero proportion
#' @param Z the covariates for non-zero compoment
#' @param alpha the regression coeffcient
#' @param beta the regression coeffcient
#' @param sim.seed the seed for the random number generator
#' @return test.table, eta, b, phi
#' @export
#' @examples
#' likelihood_ratio_test_mada(Data,Covariates,total.counts,marker.length)


#x: a numeric vector of data values, or a list of numeric data vectors.
#g: a vector or factor object giving the group for the corresponding elements of x. Ignored if x is a list.

ANOVA.zeros <- function(x, g, alpha = 0.05, perm = FALSE){
	if(class(x) == "numeric" || class(x) == "data.frame"){
		s = levels(g); newx = as.vector(NULL, mode = "list")
		for(i in 1:s){
			newx = list(newx, x[g==leves(g)[i]])
		}
		names(newx) = names(g)
	}

	K = length(x); N = n = rep(0, K); xvec = numeric(0);
	for(i in 1:K){ 
		N[i] = length(x[[i]]); n[i] = sum(x[[i]] !=0)
		xvec = c(xvec, x[[i]]);
		}
	prop = n/N; pmax = max(prop);
	Ntrun = round(pmax*N);
	
	Xtrun.vec = numeric(0);
	for(i in 1:K){ 
		data = x[[i]];
		Xtrun.vec = c(Xtrun.vec, data[data != 0], rep(0, Ntrun[i] - n[i]));
		}
	rankdata = sum(Ntrun) + 1 - rank(Xtrun.vec); 
	
	r = sum(rankdata[1:Ntrun[i]]);
	for(i in 2:K){ 
		r = c(r, sum(rankdata[1:Ntrun[i] + sum(Ntrun[1:(i-1)])]))
		}
	s = r - Ntrun*(sum(Ntrun) + 1)/2; u = numeric(0);
	for(i in 1:(K - 1))
		u = c(u, N[i+1]*sum(s[1:i]) - sum(N[1:i])*s[i+1]);
	u = u/sum(N)^2;

	thetam = mean(prop);
	simun = matrix(0, nrow = 5000, ncol = K); simup = simun;
	for(ss in 1:K){
		simun[,ss] = rbinom(5000, N[ss], thetam);		
		simup[,ss] = simun[,ss]/N[ss];
	}
	simupmax = apply(simup, 1, max);
	varsimu = numeric(K - 1);
	varsimu[1] = N[2]^2*mean(simupmax^2*(simup[,1] - simup[,2])^2)*N[1]^2;
	varu2 = N[2]*N[1]*(N[1] + N[2]);
	for(ss in 2:(K-1)){
		varsimu[ss] = N[ss+1]^2*mean(simupmax^2*(apply(simun[,1:ss], 1, sum) - simup[,ss+1]*sum(N[1:ss]))^2);
		varu2 = c(varu2, N[ss+1]*sum(N[1:ss])*sum(N[1:(ss+1)]));
	}
	varsimu = varsimu/(sum(N))^2/4;

	varu2 = varu2*thetam^2*(thetam + 1/sum(N))/12/(sum(N))^2;
	varu = varsimu + varu2;

	w = sum(u^2/varu);

	if(perm){
		numrep = 10000; permu.w = rep(0, numrep)
		for(i in 1:numrep){
			ind = sample(xvec, sum(N));
			for(i in 1:K) n[i] = sum(xvec[sum(N[0:(i - 1)]) + 1:N[i]] !=0);
			prop = n/N; pmax = max(prop);
			Ntrun = round(pmax*N);
	
			Xtrun.vec = numeric(0);
			for(i in 1:K){ 
				data = xvec[sum(N[0:(i - 1)]) + 1:N[i]];
				Xtrun.vec = c(Xtrun.vec, data[data != 0], rep(0, Ntrun[i] - n[i]));
			}
			rankdata = sum(Ntrun) + 1 - rank(Xtrun.vec); 
	
			r = sum(rankdata[1:Ntrun[i]]);
			for(i in 1:length(x)){ 
				r = c(r, sum(rankdata[1:Ntrun[i] + sum(Ntrun[1:(i-1)])]))
			}
			s = r - Ntrun*(sum(Ntrun) + 1)/2;
			for(i in 1:(K - 1))
				u = c(u, N[i+1]*sum(s[1:i]) - sum(N[1:i])*s[i+1]);
			u = u/sum(N)^2;

			thetam = mean(prop);
			simun = matrix(0, nrow = 5000, ncol = K); simup = simun;
			for(ss in 1:K){
				simun[,ss] = rbinom(5000, N[ss], thetam);		
				simup[,ss] = simun[,ss]/N[ss];
			}
			simupmax = apply(simup, 1, max);
			varsimu = numeric(K - 1);
			varsimu[1] = N[2]^2*mean(simupmax^2*(simup[,1] - simup[,2])^2)*N[1]^2;
			for(ss in 2:(K-1)){
				varsimu[ss] = N[ss+1]^2*mean(simupmax^2*(apply(simun[,1:ss], 1, sum) - simup[,ss+1]*sum(N[1:ss]))^2);
			}
			varsimu = varsimu/(sum(N))^2/4;

			varu2 = c(N[1]*N[2], N[3]*sum(N))*(N[1] + N[2])*sum(N)*thetam^2*(sum(N)*thetam + 1)/12/(sum(N))^4;
			varu = varsimu + varu2;

			permu.w[i] = sum(u^2/varu);
			}
		pw = sum(abs(w) < abs(permu.w))/numrep;
		}

	else{
		pw = pchisq(w, K-1, lower.tail = F);}

	H = (pw < alpha);
	result <- list(H = H, p.value = pw, statistics=w)
	return(result)
}
