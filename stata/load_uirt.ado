*load_uirt.ado
*ver 1.1
*2023.04.13
*everythingthatcounts@gmail.com

cap mata: mata drop load_uirt()
cap prog drop load_uirt

program define load_uirt, eclass
version 10
syntax [anything], Estimates(str)
	preserve
	qui import delimited "`estimates'",case(preserve) clear
	m: load_uirt()
	restore
end

mata
	void load_uirt(){
	
		var = st_sdata(.,"var") 
		par = st_sdata(.,"par")
		est = st_data(.,"est")
		groupN_itemCATS = st_data(.,"groupN_itemCATS")
		
		eret_b_colnames=(var,par)
		eret_b=est'
		
		gr_var=select(var,strpos(par,"theta"):>0)[1]
		gr_var=substr(gr_var,1, strrpos(gr_var,"_")-1)		
		ind = select( (1::rows(par)) , ( strpos(par,"_theta"):>0  )  ) 
		eret_b_colnames[ind,1] = subinstr(eret_b_colnames[ind,1],gr_var,"group")
		
		st_matrix("b",eret_b)
		st_matrixcolstripe("b",eret_b_colnames)
		st_matrixrowstripe("b", ("","y1"))

		
		V_vars=(var:+"_":+par)'
		st_matrix("V",st_data(.,V_vars))
		st_matrixcolstripe("V", eret_b_colnames)
		st_matrixrowstripe("V", eret_b_colnames)
		
		
		group_N=select(groupN_itemCATS,strpos(par,"mean_theta"):>0)'
		group_N_colnames=select(eret_b_colnames[.,1],strpos(par,"mean_theta"):>0)		
		group_N_colnames=J(rows(group_N_colnames),1,""),group_N_colnames
		
		st_matrix("group_N",group_N)
		st_matrixcolstripe("group_N",group_N_colnames)
		st_matrixrowstripe("group_N", ("","N"))
		
		group_par=J(2,cols(group_N),.)
		group_par_colnames=group_N_colnames
		for(g=1;g<=cols(group_par);g++){
			m=select(est, rowsum(eret_b_colnames :== (group_par_colnames[g,2] ,  "mean_theta" ) ) :==2  )
			sd=select(est, rowsum(eret_b_colnames :== (group_par_colnames[g,2] ,  "sd_theta" ) ) :==2  )
			group_par[1,g] = m
			group_par[2,g] = sd
		}
		
		st_matrix("group_par",group_par)
		st_matrixcolstripe("group_par",group_par_colnames)
		st_matrixrowstripe("group_par", ("","mean"\"","sd"))
		
		
		items=select(var,strpos(par,"_a"):>0)
		max_cat=0
		for(i=1; i<=rows(items); i++){
			i_cats=rows(select(par, (var:==items[i]) :* (strpos(par,"c")!=0) ))
			max_cat=max( (max_cat,i_cats) )
		}
		item_cats=J(rows(items),max_cat,.)
		for(i=1; i<=rows(items); i++){
			cats=select(groupN_itemCATS, (var:==items[i]) :* (strpos(par,"c")!=0) )
			item_cats[i,1..rows(cats)]=cats'
		}
		
		item_cats_colnames=J(max_cat,1,""),("cat_":+strofreal(1::max_cat))
		item_cats_rownames=J(rows(items),1,""),items
		
		st_matrix("item_cats",item_cats)
		st_matrixcolstripe("item_cats",item_cats_colnames)
		st_matrixrowstripe("item_cats",item_cats_rownames)
		
		
		models=J(rows(items),1,"")
		for(i=1;i<=rows(items);i++){
			model = select(par, var:==items[i])
			model = uniqrows(substr(model, 1, strpos(model,"_"):-1))
			models[i] = model
		}
		
		par_labels=substr(par,strpos(par,"_"):+1,strlen(par))
		par_labels=select(par_labels,par_labels:!="theta")
		par_labels=uniqrows(par_labels)
		
		if_c=(sum(par_labels:=="c")>0)
		if_b=(sum(par_labels:=="b")>0)
		max_bn=max(strtoreal(subinstr(select(par_labels,strpos(par_labels,"b")),"b","")))
		
		par_labels="a"
		if(if_b){
			par_labels=par_labels\"b"
		}
		if(if_c){
			par_labels=par_labels\"c"
		}
		if(max_bn!=.){
			par_labels=par_labels\(J(max_bn,1,"b"):+strofreal(1::max_bn))
		}
			
		item_pars = J(rows(items),rows(par_labels),.)
		for(i=1;i<=rows(items);i++){
			for(j=1;j<=rows(par_labels);j++){
				p = select(est, ( var:==items[i] ) :* ( par :== (models[i]+"_"+par_labels[j]) ) )
				if(rows(p)){
					item_pars[i,j]=p
				}
			}
		}
	
		st_matrix("item_par",item_pars)
		st_matrixrowstripe("item_par",(items,models))
		st_matrixcolstripe("item_par",(J(rows(par_labels),1,""),par_labels))
		
		items_list=""
		for(i=1;i<=rows(items);i++){
			items_list=items_list + " " + items[i]
		}
		
		eret_cmdstrip = "uirt "+items_list+",gr("+gr_var+")"
		
		stata("ereturn clear")
		stata("ereturn post b V")
		
		stata("ereturn matrix item_cats item_cats")
		stata("ereturn matrix group_N group_N")
		stata("ereturn matrix group_par group_par")
		stata("ereturn matrix item_par item_par")
		
		stata("ereturn local cmdstrip "+char(34)+eret_cmdstrip+char(34))
		stata("ereturn local cmd "+char(34)+"uirt"+char(34))
		
		display("uirt estimates loaded into memory")
		
	}
	

end