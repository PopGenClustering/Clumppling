python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/Fortier2020_791loci/b10000/Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/HGDP_diffmodel/791loci \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.2
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/Fortier2020_13loci/b10000/Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/HGDP_diffmodel/13loci \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.2
python -m clumppling.diffModel \
-i G:/My\ Drive/Projects/ImprClsAlign/output/HGDP_diffmodel \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diffmodel \
--plot_separate N