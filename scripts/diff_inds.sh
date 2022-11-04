python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/diversity200/diversity200_Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002/HGDP_200 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.16
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/diversity100/diversity100_Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002/HGDP_100 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.16
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/diversity50/diversity50_Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002/HGDP_50 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.16
python -m clumppling.diffInd \
-i G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002 \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diversity \
--plot_separate N