*export_exog
*ver 1.1
*2022.12.02
*everythingthatcounts@gmail.com

cap mata: mata drop export_exog()
cap mata: mata drop verify_isnumvar()
cap prog drop export_exog

program define export_exog
version 10
syntax [varlist], STudentid(str) SCHoolid(str) [Outpath(str)]

	
	
	m: st_local("st_isnumvar",verify_isnumvar("`studentid'"))
	if(strlen("`st_isnumvar'")){
		di as err "{p 0 2}Student variable must be numeric, `st_isnumvar' is not{p_end}"
		exit 109
	}
	m: st_local("sch_isnumvar",verify_isnumvar("`schoolid'"))
	if(strlen("`sch_isnumvar'")){
		di as err "{p 0 2}School variable must be numeric, `sch_isnumvar' is not{p_end}"
		exit 109
	}
	
	unab allvars: *
	if("`allvars'"=="`varlist'"){
		local varlist=""
	}
	if(strlen("`varlist'")){	
		m: st_local("exog_isnumvar",verify_isnumvar("`varlist'"))
		if(strlen("`exog_isnumvar'")){
			di as err "{p 0 2}All exog variables must be numeric, {p_end}"
			di as err "{p 0 2}the following item variables are strings: `exog_isnumvar'{p_end}"
			exit 109
		}
		unab exog: `varlist'
	}
	else{
		local exog=""
	}
	qui preserve
	
	keep `exog' `studentid' `schoolid'
	
	m: export_exog("`exog'" , "`studentid'" , "`schoolid'" , "`outpath'")
	
	qui restore

end

mata:

	function verify_isnumvar(string scalar items){
		itemlist	= tokens(items)'
		notoklist		=""
		for(i=1;i<=rows(itemlist);i++){
			if(st_isnumvar(itemlist[i])==0){
				notoklist	= notoklist+" "+itemlist[i]
			}
		}
		return(notoklist)
	}

	void export_exog(string scalar exog_vars , string scalar student_id, string scalar school_id, string scalar out_path){
	
		if(fileexists(pwd()+"pvreg_config.json")){
			display("pvreg_config.json found, working in append mode")
			X = fopen(pwd()+"pvreg_config.json", "rw")
			fseek(X,0,1)
		}
		else{
			_error("pvreg_config.json NOT found, please run export_egog after uirt estimates are exported")
		}
	
		wd = subinstr(pwd(),"\","/")
				
		stata("qui compress")
		stata("export delimited using "+char(34)+"student_exog.csv"+char(34)+", replace nolabel")
		
		fput(X,char(9)+char(34)+"path"+char(34)+": " +char(34)+wd+char(34)+",")
		fput(X,char(9)+char(34)+"student_exog"+char(34)+": " +char(34)+"student_exog.csv"+char(34)+",")
		
		fput(X,char(9)+char(34)+"student_id"+char(34)+": " +char(34)+student_id+char(34)+",")
		fput(X,char(9)+char(34)+"school_id"+char(34)+": " +char(34)+school_id+char(34)+",")
		fput(X,char(9)+char(34)+"fixed_effects"+char(34)+": " +char(34)+exog_vars+char(34)+",")
		
		fput(X,char(9)+char(34)+"keep_pv"+char(34)+": " +char(34)+"True"+char(34)+",")
		fput(X,char(9)+char(34)+"npv"+char(34)+": 15,")
		fput(X,char(9)+char(34)+"njobs"+char(34)+": 1,")
		
		
		if(strlen(out_path)){
			out_path = subinstr(out_path,"\","/")
			if(substr(out_path,-1)!="/"){
				out_path=out_path+"/"
			}
			fput(X,char(9)+char(34)+"out_path"+char(34)+": " +char(34)+out_path+char(34)+",")
		}
		
		fput(X,char(9)+char(34)+"out_files"+char(34)+": " +char(34)+"pvreg_results_"+char(34))
		
		
		fput(X,"}")
		
		fclose(X)
	}

end

