python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/diversity200/diversity200_Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002/HGDP_200 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.18 \
--custom_cmap Y \
--cmap "#FF0000 #FFFF00 #00EAFF #AA00FF #FF7F00 #BFFF00 #0095FF #FF00AA #FFD400 #6AFF00 #0040FF #EDB9B9 #B9D7ED #E7E9B9 #DCB9ED #B9EDE0 #8F2323 #23628F #8F6A23 #6B238F #4F8F23 #000000 #737373 #CCCCCC"
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/diversity100/diversity100_Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002/HGDP_100 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.18 \
--custom_cmap Y \
--cmap "#FF0000 #FFFF00 #00EAFF #AA00FF #FF7F00 #BFFF00 #0095FF #FF00AA #FFD400 #6AFF00 #0040FF #EDB9B9 #B9D7ED #E7E9B9 #DCB9ED #B9EDE0 #8F2323 #23628F #8F6A23 #6B238F #4F8F23 #000000 #737373 #CCCCCC"
python -m clumppling.main \
-i G:/My\ Drive/Projects/ImprClsAlign/StructureHuman/diversity50/diversity50_Results \
-o G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002/HGDP_50 \
-p G:/My\ Drive/Projects/ImprClsAlign/Clumppling/scripts/params.xml \
-f structure \
--cd_mod_thre 0.18 \
--custom_cmap Y \
--cmap "#FF0000 #FFFF00 #00EAFF #AA00FF #FF7F00 #BFFF00 #0095FF #FF00AA #FFD400 #6AFF00 #0040FF #EDB9B9 #B9D7ED #E7E9B9 #DCB9ED #B9EDE0 #8F2323 #23628F #8F6A23 #6B238F #4F8F23 #000000 #737373 #CCCCCC"
python -m clumppling.diffInd \
-i G:/My\ Drive/Projects/ImprClsAlign/output/Rosenberg2002 \
-o G:/My\ Drive/Projects/ImprClsAlign/output/diversity \
--plot_separate N \
--custom_cmap Y \
--cmap "#FF0000 #FFFF00 #00EAFF #AA00FF #FF7F00 #BFFF00 #0095FF #FF00AA #FFD400 #6AFF00 #0040FF #EDB9B9 #B9D7ED #E7E9B9 #DCB9ED #B9EDE0 #8F2323 #23628F #8F6A23 #6B238F #4F8F23 #000000 #737373 #CCCCCC"