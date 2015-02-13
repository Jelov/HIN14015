# HIN14015 fit macros
* _ctauErrorRange_step8: Contains ctau error ranges for analysis bins
* fit2DData.h, fit2DData_pbpb.cpp: Fit macros, need to be complied (Tested uner ROOTv5.28.00d)
* runBatch_***.sh: Make batch jobs and run fits for all analysis bins with options
* run.sh: Feed RooDataSet files to runBatch_***.sh, determine name of results
* extract.py: After all fitting jobs are done, all numbers are sorted into excel files by this script
* rfcp.sh: use extract.py and find if there is any missing fitting jobs
