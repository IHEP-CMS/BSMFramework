cp DATA2016/miniAOD_Run*.py .
cp DATA2016/crab_Run*_cfg.py .
crab submit -c crab_RunB_cfg.py
crab submit -c crab_RunC_cfg.py
crab submit -c crab_RunD_cfg.py
crab submit -c crab_RunE_cfg.py
crab submit -c crab_RunF_cfg.py
crab submit -c crab_RunG_cfg.py
crab submit -c crab_RunH_cfg.py
rm -rf crab_Run*_cfg.py*
rm -rf miniAOD_Run*
