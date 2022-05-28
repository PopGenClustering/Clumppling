# input_path = "G:/My Drive/Projects/ImprClsAlign/StructureChicken/chicken/Rosenberg2001/Results"
# output_path = "G:/My Drive/Projects/ImprClsAlign/output/Rosenberg2001"
# prj_type="structure" 

# input_path = "G:/My Drive/Projects/ImprClsAlign/StructureHuman/FortierEtAl2020/para/Results/"
# output_path = "G:/My Drive/Projects/ImprClsAlign/output/Fortier2020/"
# prj_type="structure" 

# input_path = "G:/My Drive/Projects/ImprClsAlign/StructureHuman/PembertonEtAl2008/"
# output_path = "G:/My Drive/Projects/ImprClsAlign/output/Pemberton2008/"
# prj_type="fastStructure" 

##  community detection parameter (used if not leader-clustering)

class Params:
    def __init__(self, input_path,output_path,prj_type="structure"):
        
        self.input_path = input_path
        self.output_path = output_path
        self.prj_type = prj_type
        
        
        # mode detection parameters
        self.default_cd = True
        self.cd_mod_thre = None
        self.lc_flag = False     
        self.adaptive_thre_flag = True
        self.lc_cost_thre = 0
        
        # plotting parameters
        self.custom_cmap = False
        self.cmap = ""
        self.plot_flag_community_detection = False
        self.plot_flag_all_modes = True
        self.plot_flag_all_within_mode = False
        self.plot_flag_mode_within_K = False
        self.plot_flag_mode_across_K_multipartite = True
        self.plot_flag_mode_across_K_chains = False
        
        # other parameters
        # self.Qbar_flag = False
        self.reorder_inds = False
        
        
      
# default_cd = True
# cd_mod_thre = 0.0
# # method=("louvain",1.0,0.0)

# # method=("mcl",2)
# # method=("belief",None)
# # method = ("agdl",3)
# # method = ("paris",None)
# # method = ("eigenvector",None)
# # method=("belief",None)
# # method = ("mcode",None)
# Qbar_flag = False

# lc_flag = False     
# adaptive_thre_flag = True
# lc_cost_thre = 0

    
# # reorder_ind_flag = False
# reorder_cls_flag = False


# plot_flag_community_detection = False
# plot_flag_all_modes = True
# plot_flag_all_within_mode = False
# plot_flag_mode_within_K = False
# plot_flag_mode_across_K_multipartite = True
# plot_flag_mode_across_K_chains = False