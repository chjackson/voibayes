
data {
    int GMSM;
    int NGMSM;
    int PMSM;
    int MSM;
    int I_gumcad;
    int<lower=0> obs_ons;
    int<lower=0> obs_natsal[4];
    int<lower=0> obs_gumcad[I_gumcad];
    int<lower=0> obs_uagum[2];
    int<lower=0> obs_gmshs[2,2];
    int<lower=0> denom_gmshs;
    int<lower=0> obs_sophid[2];
    int<lower=0> obs_handd;
    real pioptout_max;
}

transformed data {
//// could declare and define constants here
    vector[MSM] alpha;
    for (i in 1:MSM) alpha[i] = 1;
}

parameters {
    real<lower=0> mean_obs_ons;
    simplex[MSM] rho;
    real<lower=0,upper=1> aux_delta[PMSM];
    ////// These have different notation from JAGS, since can't define both stochastic and deterministic nodes in same array
    real<lower=0,upper=1> gamma[I_gumcad-1];
    real<lower=0,upper=1> aux_gumcad[3];
    real<lower=0,upper=1> prob_gmgum;
    real<lower=0,upper=1> prob_gmshs[2];
    real lor_pinodelta;
    real<lower=0> aux_sophid;
    real<lower=0,upper=1> aux_handd;
}

transformed parameters {
    vector<lower=0,upper=1>[MSM] pi;
    real<lower=0,upper=1> delta[PMSM];
    vector<lower=0,upper=1>[MSM] pidelta;
    real<lower=0,upper=1> pinodelta[MSM];
    real<lower=0,upper=1> delta_msm;
    real<lower=0,upper=1> rho_msm;
    vector<lower=0>[MSM] r;
    vector<lower=0>[MSM] mu;
    vector<lower=0>[MSM] mudelta;
    real<lower=0> munodelta[MSM];

    real<lower=0,upper=1> pigd; // prev of newly-diag infection 
    real<lower=0,upper=1> piga;
    real<lower=0,upper=1> piun;
    real<lower=0,upper=1> piop;
    real<lower=0,upper=1> aex;

    real<lower=0,upper=1> prev_unoffered;
    real<lower=0> or_gmshs;
    real<lower=0> or_pinodelta;
    real<lower=0> mean_obs_sophid;
    real<lower=0,upper=1> prob_handd;

    rho_msm = sum(rho[GMSM:PMSM]); // Relative MSM size.

    pigd = gamma[1]*gamma[2]*gamma[3]*gamma[4]; // Newly diagnosed HIV prevalence among GMSM. Clearer written as product of four probs?  p1*p2*p3*p4
    prev_unoffered = inv_logit(logit(gamma[I_gumcad-1]) + log(0.5) + aux_gumcad[1]*(log(1.5) - log(0.5))); 
    piun = gamma[1] * (1 - gamma[2]) * prev_unoffered;  // Undiagnosed HIV prevalence (due to unoffered HIV test) among GMSM

    aex = (pioptout_max - gamma[I_gumcad-1]) * aux_gumcad[2];
    piop = gamma[1] * gamma[2] * (1 - gamma[3]) * (gamma[4] + aex); 
    
    pinodelta[GMSM] = piun + piop;      //// True undiagnosed prevalence in people attending GUM in last year
    piga = (pinodelta[GMSM] + pigd) / gamma[1]; // HIV prevalence among previously undiagnosed GMSM (directly estimated from GUM Anon)

    ////// Assumption. OR for prevalence between NG/G same in GMSHS as in wider population
    or_gmshs = (prob_gmshs[NGMSM]/(1-prob_gmshs[NGMSM])) / (prob_gmshs[GMSM] / (1 - prob_gmshs[GMSM]));
    pinodelta[NGMSM] = pinodelta[GMSM] * or_gmshs  / (pinodelta[GMSM] * or_gmshs  + 1 - pinodelta[GMSM]);

    or_pinodelta = exp(lor_pinodelta); // PMSM:GMSM odds ratio of undiagnosed HIV prevalences
    pinodelta[PMSM] = pinodelta[GMSM] * or_pinodelta / (1 - pinodelta[GMSM] * (1 - or_pinodelta));

    for(grp in GMSM:PMSM)// Subgroup counter ("GMSM", "NGMSM", "PMSM")
	{
	    delta[grp] = aux_delta[grp]*(1 - pinodelta[grp]); // so that pinodelta < 1 - delta  in definition of pi
	    r[grp] = rho[grp] * mean_obs_ons;  // Absolute subgroup size
	    pi[grp] = pinodelta[grp] / (1 - delta[grp]); // HIV prevalence among NGMSM
	    pidelta[grp]   = pi[grp]  * delta[grp];    // Diagnosed HIV prevalence
	    mu[grp]    = pi[grp]  * r[grp];  // Number of PLWH
	    mudelta[grp] = delta[grp] * mu[grp]; // Number of diagnosed PLWH
	    munodelta[grp] = mu[grp] - mudelta[grp]; // Number of undiagnosed PLWH
	}
	pi[MSM] = dot_product(rho[GMSM:PMSM], pi[GMSM:PMSM])  / rho_msm;

	delta_msm = dot_product(rho[GMSM:PMSM], pidelta[GMSM:PMSM]) / rho_msm / pi[MSM];  // Proportion of PLWH diagnosed among MSM
	pidelta[MSM] = pi[MSM] * delta_msm;// Diagnosed HIV prevalence among MSM
 	pinodelta[MSM] = pi[MSM] * (1 - delta_msm);// Undiagnosed HIV prevalence among MSM.
	r[MSM] =   sum(r[GMSM:PMSM]); // Absolute MSM size.
	mu[MSM] =    sum(mu[GMSM:PMSM]);   // Number of PLWH among MSM
	mudelta[MSM] = sum(mudelta[GMSM:PMSM]);// Number of diagnosed PLWH among MSM
	munodelta[MSM] = mu[MSM] - mudelta[MSM];// Number of undiagnosed PLWH among MSM

        mean_obs_sophid = aux_sophid*mudelta[MSM]; // log expected number of MSM diagnosed 

	//// HANDD informs mudelta[GMSM].
	//// prob_handd is p(GUM | MSM diag), mudelta ratio is P(GMSM diag | MSM diag)
	//// aux_handd is proportion of GMSM diagnoses registered in GUM clinics: P(GUM|GMSM)
	prob_handd = aux_handd * mudelta[GMSM] / mudelta[MSM]; // Probability of a new MSM HIV diagnosis being a GMSM.

}



model {

    for(grp in GMSM:PMSM)// Subgroup counter ("GMSM", "NGMSM", "PMSM")
	{
	    aux_delta[grp] ~ uniform(0, 1);
	}
    rho ~ dirichlet(alpha);

    obs_ons ~ poisson(mean_obs_ons); // Likelihood of male population size (source: ONS)
    mean_obs_ons ~ lognormal(0, 1000); // Non-informative prior on male population log-size

    for(grp in GMSM:PMSM)// Subgroup counter ("GMSM", "NGMSM", "PMSM")
	{
	    obs_natsal[grp+1] ~ binomial(obs_natsal[1], rho[grp]);// Likelihood of relative subgroup size (source: NATSAL-2)
	    //// CJ: shouldn't this be multinomial?  Doesn't matter for model fitting, as long as rho are dirichlet
	}

    for(i in 2:I_gumcad)// GUM clinic event counter ("Previously undiagnosed", "HIV test offered", "HIV test accepted", "New HIV diagnosis")
	{
	    obs_gumcad[i] ~ binomial(obs_gumcad[i-1], gamma[i-1]);// Likelihood of GUM clinic HIV event (source: GUMCAD 2012)
	}
    for(i in 1:(I_gumcad-2)) {
	gamma[i] ~ uniform(0, 1);
    }
    gamma[I_gumcad-1] ~ uniform(0, pioptout_max);

    aux_gumcad[1] ~ uniform(0, 1); // Constraint on undiagnosed HIV prevalence (due to unoffered HIV test) among GMSM
    aux_gumcad[2] ~ uniform(0, 1); // Constraint on undiagnosed HIV prevalence (due to declined HIV test) among GMSM
    aux_gumcad[3] ~ uniform(0, 1); // Non-informative prior on normalised undiagnosed HIV prevalence among GMSM

    obs_uagum[2] ~ binomial(obs_uagum[1], piga); // Likelihood of HIV infection among previously undiagnosed GMSM (source: GUM Anon 2012-13)

    obs_gmshs[1, 1] ~ binomial(denom_gmshs, prob_gmgum);
    prob_gmgum ~ uniform(0, 1); 
    for(grp in GMSM:NGMSM)// Subgroup counter ("GMSM", "NGMSM")
	{
	    obs_gmshs[grp, 2] ~ binomial(obs_gmshs[grp, 1], prob_gmshs[grp]); // Likelihood of undiagnosed HIV infection (source: GMSHS 2008)
	    prob_gmshs[grp] ~ uniform(0, 1);
	}

    lor_pinodelta ~ normal(-2.3, sqrt(1/20.77)); // Informative prior on PMSM:GMSM log-odds ratio of undiagnosed HIV prevalences

    obs_sophid[2] ~ poisson(mean_obs_sophid); // Likelihood of a male HIV diagnosis being an MSM (source: SOPHID 2012)
    aux_sophid ~ normal(1, 0.018)T[0,]; // from Anne

    obs_handd ~ binomial(obs_sophid[2], prob_handd); // Likelihood of a new GMSM HIV diagnosis among MSM (source: HANDD 2012)
    //	prob_handd ~ uniform(0, 1); // if cutting inference
    aux_handd ~ uniform(0, 1); // Constraint on the probability of an MSM HIV diagnosis being a GMSM
}



generated quantities {
}
