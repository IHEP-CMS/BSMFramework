1. multicrab_RD.py
   python file that can handle multiple crab configurations per run, 
   where the different configurations change depending on different inputs
2. MuonEG_2016B_cfg.py
   example of crab configuration file where a parameter of the running configuration python file is set (see config.JobType.pyCfgParams =)
3. runBTagAnalyzer_cfg.py
   exmaple of running configuration python file

To test the implementation, you can try like this:
1. Modify the multicrab_RD.py file as mentioned before
2. python multicrab_RD.py
   this command usually create and submit a task subsequently, but you can wait it create the information for the task
   and then ctr+c the task.
   Or let it go and the kill the task later
3. Open the file datasetnames/crab_datasetnames/inputs/PSetDump.py
   and check whether in PSetDump.py you see the proper jec files and also the correct json files
