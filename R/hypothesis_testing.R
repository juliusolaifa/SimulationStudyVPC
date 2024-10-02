testVpc <- function( null, alt, beta, ns, Sigma = NULL,alpha = 0.05, num = 1,
                     X = NULL, covariate = NULL,family = "gaussian",
                     link = "identity", num_cores=1, type,
                     seed = NULL,...) {

    if(!is.null(seed)) {
      set.seed(seed)
    }

    data <- glmmVpc::batchGLMMData(beta=beta,ns=ns,Sigma=Sigma,num=num,
                                 X=X,family=family,link=link,...)

    fits <- glmmVpc::batchGLMMFit(formula=alt,dataMat=data, family=family)
    vpc.est <- glmmVpc::vpc(fits,x=covariate)
    hyp.test <- glmmVpc::vpc.test(vpc.est, null_formula=null, type=type)

}



