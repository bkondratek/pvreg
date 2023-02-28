*You might need to install uirt first:
* ssc install uirt

* fitting irt model to the first dataset
use responses_in,clear
uirt i_in*, not gr(grouping_var_in,slow) nit(200)
mat l e(group_par)
export_uirt, f(0) s(id_stud)

* fitting irt model to the second dataset
use responses_out,clear
uirt i_out*, not gr(grouping_var_out,slow) nit(200)
mat l e(group_par)
export_uirt, f(1) s(id_stud)

*exporting the exog data
use example_data_all,clear
export_exog student_gender , st(id_stud) sch(id_schl)  o("./")

