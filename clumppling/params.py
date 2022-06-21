
class Params:
    def __init__(self, input_path,output_path,input_type="structure"):
        
        self.input_path = input_path
        self.output_path = output_path
        self.input_type = input_type
        
        
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
        
        