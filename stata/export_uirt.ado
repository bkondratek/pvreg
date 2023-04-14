*export_uirt.ado 
*ver 1.2
*2022.12.02
*everythingthatcounts@gmail.com

cap mata: mata drop export_uirt()
cap prog drop export_uirt

program define export_uirt
version 10
syntax [anything], [Filetype(numlist integer max=2 >=0) Studentid(str) PREFix(str)]

	if( strlen("`prefix'")>0 & (strlen("`filetype'")>0 | strlen("`studentid'")>0) ){
		di as err "Proper usage is either: export_uirt ,  Filetype(numlist integer max=2 >=0) Studentid(str)"
		di as err "                    or: export_uirt ,  PREFix(str))"
		exit 109
	}
	if( strlen("`prefix'")==0 & (strlen("`filetype'")==0 | strlen("`studentid'")==0) ){
		di as err "Proper usage is either: export_uirt ,  Filetype(numlist integer max=2 >=0) Studentid(str)"
		di as err "                    or: export_uirt ,  PREFix(str))"
		exit 109
	}
	if(strlen("`filetype'")==0){
	    local filetype=.
	}
	
	qui preserve
	m: export_uirt(`filetype' , "`studentid'", "`prefix'")
	qui restore

end

mata:
	void export_uirt(real scalar filetype , string scalar student_id, string scalar prefix){
	
		if(strlen(prefix)){
			file_prefix=prefix
		}
		else{	
			if(filetype){
				file_prefix="out"
			}
			else{
				file_prefix="in"
			}
		}
		
		if(strlen(prefix)==0){
			if(fileexists(pwd()+"pvreg_config.json")){
				display("pvreg_config.json found, working in append mode")
				X = fopen(pwd()+"pvreg_config.json", "rw")
				fseek(X,0,1)
			}
			else{
				display("pvreg_config.json not found, will be created")
				X = fopen(pwd()+"pvreg_config.json", "rw")
				fput(X,"{")
			}
		}
		
		wd = subinstr(pwd(),"\","/")
		wd = substr(wd,1,strlen(wd)-1)
		wd = substr(wd,strrpos(wd, "/")+1,strlen(wd))+"/"
	
		cmdstrip = st_global("e(cmdstrip)")
		if(strpos(cmdstrip,"gr(")){
			gr_var = substr(cmdstrip, strpos(cmdstrip,"gr(")+3,strlen(cmdstrip))
			gr_var = tokens(gr_var)[1]
			gr_var = subinstr(gr_var,",","")
			gr_var = subinstr(gr_var,")","")
		}
		else{
			gr_var = "group"
			ind = st_addvar("byte",gr_var)
			st_store(.,ind,J(st_nobs(),1,1))
		}

		group_vals=uniqrows(st_data(.,gr_var))
		group_vals = select(group_vals,group_vals :< .)
		group_mat_labels="group_":+strofreal(group_vals)
		group_N=J(rows(group_vals),1,0)
		for(i=1;i<=rows(group_vals);i++){
			group_N[i] = sum(st_data(.,gr_var):==group_vals[i])
		}
		
		string_vars = st_matrixcolstripe("e(b)")
		e_b = st_matrix("e(b)")
		e_v = st_matrix("e(V)")
		e_g_par = st_matrix("e(group_par)")
		e_g_N  = st_matrix("e(group_N)")
		e_g_labs = st_matrixcolstripe("e(group_N)")[.,2]
		
		e_c=st_matrix("e(item_cats)")
		e_c_items=st_matrixrowstripe("e(item_cats)")[.,2]
		
		groupN_itemCATS=J(rows(string_vars),1,.)
		for(g=1;g<=rows(group_vals);g++){
			ind = select( (1::rows(string_vars)) , (string_vars[.,2] :== "mean_theta" ) :* (string_vars[.,1] :== ("group_"+strofreal(group_vals[g])) ) ) 
			groupN_itemCATS[ind] = group_N[g]
		}
		for(i=1;i<=rows(e_c_items);i++){
			cats = select(e_c[i,.]',e_c[i,.]':!=.)
			ind = select( (1::rows(string_vars)) , (string_vars[.,2] :!= "3plm_c" ) :* (string_vars[.,1] :== e_c_items[i] ) ) 
			groupN_itemCATS[ind] = cats
		}
		
		
		if(strlen(prefix)==0){
			i_start=min( select( (1::rows(string_vars)) ,strpos(string_vars[.,2],"_theta"):==0  ) )
			items=""
			for(i=i_start;i<=rows(string_vars);i++){
				if(string_vars[i-1,1]!=string_vars[i,1]){
					items=items+" "+string_vars[i,1]
				}
			}
		
			stata("qui drop if "+gr_var+"==.")
			stata("keep " + student_id + " " + gr_var+ " " + items)
			stata("order " + student_id + " " + gr_var+ " " + items)
	
			stata("qui compress")
			stata("export delimited using "+char(34)+file_prefix+"_responses.csv"+char(34)+", replace nolabel")
			fput(X,char(9)+char(34)+"responses_"+file_prefix+char(34)+": " +char(34)+file_prefix+"_responses.csv"+char(34)+",")
		}

		//I shall always rescale
		//rescale_groups = 0
		
		//if( ( group_N' != e_g_N ) | ( group_mat_labels!=e_g_labs )  ){
		
			rescale_groups = 1
			
			sel_g=J(cols(e_g_N),1,0)
			sel_g_long=J(2*cols(e_g_N),1,0)
			for(g=1;g<=cols(e_g_N);g++){
				sel_group = select((1::rows(group_N)), group_mat_labels:==e_g_labs[g])
				if( rows(sel_group)==1 ){
					sel_g[g] = 1
					sel_g_long[2*g-1::2*g] = 1\1
					e_g_N[g] = group_N[sel_group]
					
				}
				else{
					e_g_N[g] = 0
				}
			}
			sel_g= select( (1::rows(sel_g)), sel_g)
			sel_g_long = select( (1::cols(e_b)), sel_g_long \ J(cols(e_b)-rows(sel_g_long),1,1) )
			
			string_vars = string_vars[sel_g_long,.]
			e_b = e_b[.,sel_g_long']
			e_v = e_v[sel_g_long,sel_g_long']
			e_g_par = e_g_par[.,sel_g']
			e_g_N = e_g_N[.,sel_g']
			groupN_itemCATS = groupN_itemCATS[sel_g_long]
		//}
		
					
		ind = select( (1::rows(string_vars)) , ( strpos(string_vars[.,2],"_theta"):>0  )  ) 
		string_vars[ind,1] = subinstr(string_vars[ind,1],"group",gr_var)
		
		stata("clear")
		
		st_addobs(cols(e_b))
		
		
		string_names = ("var","par")
		for(c=1;c<=2;c++){
			ind = st_addvar("str"+strofreal(max(strlen(string_vars[.,c]))),string_names[c])
			st_sstore(.,ind,string_vars[.,c])
		}
		
		ind = st_addvar("double","est")
		st_store(.,ind,e_b')
		
		ind = st_addvar("double","groupN_itemCATS")
		st_store(.,ind,groupN_itemCATS)

		ind = st_addvar("double", (string_vars[.,1]:+"_":+string_vars[.,2])' )
		st_store(.,ind,e_v)
		
		if(rescale_groups){
			n_g = e_g_N
			mean_g = e_g_par[1,.]
			var_g = e_g_par[2,.] :* e_g_par[2,.]
			
			w_g = n_g :/ sum(n_g)
			global_mean = w_g * mean_g'
			global_var = w_g * (var_g :+ (mean_g :* mean_g))' - global_mean^2
			global_sd = sqrt(global_var) 
			
			_sel = st_sdata(.,"par"):=="mean_theta"
			_sel = select((1::rows(_sel)),_sel)
			st_store(_sel, "est" , (st_data(_sel,"est") :- global_mean) :/ global_sd )
			
			_sel = st_sdata(.,"par"):=="sd_theta"
			_sel = select((1::rows(_sel)),_sel)
			st_store(_sel,"est",  st_data(_sel,"est") :/ global_sd)
			
			_sel = strpos(st_sdata(.,"par"),"_a"):>0
			_sel = select((1::rows(_sel)),_sel)
			st_store(_sel, "est",  st_data(_sel,"est") :* global_sd)
			
			_sel = strpos(st_sdata(.,"par"),"_b"):>0
			_sel = select((1::rows(_sel)),_sel)
			st_store(_sel,"est", (st_data(_sel,"est") :- global_mean) :/ global_sd )
			
			st_store(.,ind , st_data(.,ind) :/ global_var )
			
		}

		stata("qui compress")
		stata("export delimited using "+char(34)+file_prefix+"_estimates.csv"+char(34)+", replace nolabel")
		if(strlen(prefix)==0){
			fput(X,char(9)+char(34)+"estimates_"+file_prefix+char(34)+": " +char(34)+file_prefix+"_estimates.csv"+char(34)+",")
			fclose(X)
		}
			
	}

end

