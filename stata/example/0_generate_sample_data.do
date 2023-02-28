*You might need to install uirt_sim first:
* ssc install uirt_sim

mata
	real matrix multinormal(real rowvector mu, real matrix sigma, real scalar obs){
		return(mu:+ (cholesky(sigma)*rnormal(rows(sigma),obs,0,1))')
	}	
end

clear

set seed 31415

* generating multilevel data
mata:
	R_12=0.9
	ICC=0.2
	EWD=0.05
	cov_ICC_EWD=.04
	s_s=1-ICC
	s_g=1-ICC-EWD-2*cov_ICC_EWD
	s_gs=R_12-ICC-cov_ICC_EWD
	sigma_w=(s_s , s_gs \ s_gs , s_g)
	mu_w=(0,0)
	N_w=10000
	N_b=50
	W=multinormal(mu_w,sigma_w,N_w)
	mu_b=(0,0,0)
	sigma_b=( ICC,(ICC-0.0000000001), cov_ICC_EWD \ (ICC-0.0000000001), ICC, cov_ICC_EWD \ cov_ICC_EWD, cov_ICC_EWD, EWD )
	B=sort( ( J(N_b,1,(1::N_w/N_b)) , J(N_b,1,multinormal(mu_b,sigma_b,N_w/N_b)) ) ,1 )
	X=W:+(B[.,2],rowsum(B[.,3..4]))
	st_addobs(N_w)
	ind=st_addvar("double",("w_1","w_2","id_schl","schl_mean_in","schl_mean_out","eva","theta_in","theta_out"))
	st_store(.,ind,(W,B,X))
end

keep id_schl theta* schl_mean* eva
gen id_stud=_n
order id_stud id_schl theta* schl_mean* eva

* generating item responses
mata:
	a_pars0=J(5,1,1.2)\J(5,1,2.2)
	b_pars0=J(2,1,(-2,-1,0,1,2)')
	test_length = 20
	a_pars=J(test_length/10,1,a_pars0)
	b_pars=J(test_length/10,1,b_pars0)
	
	st_matrix("i_par_in",(a_pars,b_pars))
	st_matrixrowstripe("i_par_in",( ("i_in_":+strofreal((1::test_length))) , J(test_length,1,"2plm") ) )
	st_matrixcolstripe("i_par_in",(J(2,1,""),("a"\"b")))
	
	st_matrix("i_par_out",(a_pars,b_pars))
	st_matrixrowstripe("i_par_out",( ("i_out_":+strofreal((1::test_length))) , J(test_length,1,"2plm") ) )
	st_matrixcolstripe("i_par_out",(J(2,1,""),("a"\"b")))
end

uirt_sim, ipar(i_par_in) theta(theta_in)
uirt_sim, ipar(i_par_out) theta(theta_out)

*Changing one item to be our exogenous variable
rename i_in_3 student_gender

*Some grouping variable for the second measurement, schools nested within grouping_var_out
gen grouping_var_out=0
gen temp=schl_mean_out+eva
sum temp
replace temp=(temp-r(min))/(r(max)-r(min))
foreach schl of numlist 1/200{
	local uni=runiform()
	qui sum temp if id_schl==`schl'
	if(`uni'<r(mean)){
		qui replace grouping_var_out=1 if id_schl==`schl'
	}
}
drop temp

*Some grouping variable for the first measurement, crossed with school membership
gen grouping_var_in=0
gen rnorm=rnormal()
replace grouping_var_in=1 if theta_in>rnorm
drop rnorm

order id_stud id_schl student_gender grouping_var_in grouping_var_out theta_in theta_out schl_mean_in eva
compress
drop schl_mean_out

*saving data
save example_data_all,replace

*responses from the first exam
preserve
keep id_stud grouping_var_in i_in*
save responses_in,replace
restore

*responses from the second exam
preserve
keep id_stud grouping_var_out i_out*
save responses_out,replace
restore

*sneak peak at the grouping parameters on raw thetas:
tabstat theta_out,stats(mean sd N) by(grouping_var_out)
tabstat theta_in,stats(mean sd N) by(grouping_var_in)
