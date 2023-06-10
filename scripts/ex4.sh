python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/input/diffmodel/791loci \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffmodel/791loci \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params_diffModel.json \
-f structure 
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/input/diffmodel/13loci \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffmodel/13loci \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params_diffModel.json \
-f structure 
python -m clumppling.diffModel \
-i G:/My\ Drive/Projects/ImprClsAlign/output/diffmodel \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffmodel \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/default_params_diff.json