python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/input/diffind/ind200 \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffind/ind200 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/default_params.json \
-f structure \
--cd_mod_thre 0.16
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/input/diffind/ind100 \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffind/ind100 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/default_params.json \
-f structure \
--cd_mod_thre 0.16
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/input/diffind/ind50 \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffind/ind50 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/default_params.json \
-f structure \
--cd_mod_thre 0.16
python -m clumppling.diffInd \
-i G:/My\ Drive/Projects/ImprClsAlign/output/diffind \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffind \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/default_params_diff.json