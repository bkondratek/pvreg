## Description
Three Stata programs:
1) ```export_uirt``` - after fitting an IRT model
with [uirt](https://github.com/bkondratek/uirt) you can use it to
export (i) IRT model parameter estimates and (ii) item data to CSV
files that are understood by ```pvreg```.  
2) ```load_uirt``` - allows you to load IRT model estimates 
saved by ```export_uirt``` back into Stata creating mockup.
[uirt](https://github.com/bkondratek/uirt) estimates 
3) ```export_exog``` - used to export the dataset with exogenous 
variables to CSV files that are understood by ```pvreg```.

Consecutively running ```export_uirt``` on two measurement occasions 
and then using ```export_exog``` will also produce a proper
JSON config file for ```pvreg```.

## Usage examples
Usage examples are provided in the
[1_fit_uirt_and_export_data.do](./example/1_fit_uirt_and_export_data.do) file
that you can find in the [example](./example) folder. 
After putting all three programs in your Stata personal ADO folder you can
run them like this:

```commandline
*installing uirt and uirt_sim in case not already installed
ssc install uirt
ssc install uirt_sim
*simulating some datasets in Stata for our example
*(make sure you cd to the 'example' folder on your machine before)
do 0_generate_sample_data.do
*running the example:
do 1_fit_uirt_and_export_data
```

