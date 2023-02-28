## pvreg
 Python package to run 3-level latent regression model for educational value added


## Table of Contents
[Description](#description)

[Installation](#installation)

[Usage examples](#usage-examples)

## Description
TBD


## Installation 

```commandline
pip install git+https://github.com/bkondratek/pvreg
```

## Usage examples
```pvreg``` should be executed as Python module. 

Typing:
```commandline
python -m pvreg -h
```
...will give you a hint about its arguments: 
```commandline
usage: pvreg [-h] [-c FILE] [-p PATH] [--out_path OUT_PATH] [--njobs NJOBS] [--npv NPV] [--keep_pv {True,False}] [--student_id STUDENT_ID] [--school_id SCHOOL_ID] [--estimates_in ESTIMATES_IN] [--responses_in RESPONSES_IN]
             [--estimates_out ESTIMATES_OUT] [--responses_out RESPONSES_OUT] [--student_exog STUDENT_EXOG] [--fixed_effects FIXED_EFFECTS] [-o OUT_FILES]

Generate PVs conditioned on multilevel regression

optional arguments:
  -h, --help            show this help message and exit
  -c FILE, --conf_file FILE
                        JSON config file
  -p PATH, --path PATH  Path where input data is stored
  --out_path OUT_PATH   Path where output is stored (default = args.path)
  --njobs NJOBS         Number of parallel jobs for multiprocessing
  --npv NPV             Number of plausible values drawn
  --keep_pv {True,False}
                        Whether to keep plausible values
  --student_id STUDENT_ID
                        Name of student id column (common in responses_# and student_exog df)
  --school_id SCHOOL_ID
                        name of school id column (in student_exog df)
  --estimates_in ESTIMATES_IN
                        CSV file with irt estimates for exam 0
  --responses_in RESPONSES_IN
                        CSV file with item responses for exam 0
  --estimates_out ESTIMATES_OUT
                        CSV file with irt estimates for exam 1
  --responses_out RESPONSES_OUT
                        CSV file with item responses for exam 1
  --student_exog STUDENT_EXOG
  --fixed_effects FIXED_EFFECTS
                        Comma-separated list of names of exog variables (in student_exog file) that are to be included in the fixed effects part of latent regression
  -o OUT_FILES, --out_files OUT_FILES
                        Prefix used for CSV files that are saved by pvreg
```


