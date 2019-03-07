cp DATA2018/miniAOD_Run*.py .
cp DATA2018/crab_Run*_cfg.py .
crab submit -c crab_RunA_cfg.py
crab submit -c crab_RunB_cfg.py
crab submit -c crab_RunC_cfg.py
crab submit -c crab_RunD_cfg.py
rm -rf crab_Run*_cfg.py*
rm -rf miniAOD_Run*
