
class Params:
    def __init__(self, input_path,output_path,input_type="structure"):
        
        self.input_path = input_path
        self.output_path = output_path
        self.input_type = input_type
        
        # mode detection parameters
        self.default_cd = True
        self.cd_mod_thre = -1

        # colormap parameters
        self.custom_cmap = False
        self.cmap = ""

        # Advanced: ILP when number of clusters differ by 1
        self.enum_combK = True

        # Advanced: leader clustering parameters
        self.lc_flag = False     
        self.adaptive_thre_flag = False
        self.lc_cost_thre = 0.2
        
        # Advanced: plotting parameters
        self.plot_flag_all_modes = True
        self.plot_flag_all_within_mode = False
        self.plot_flag_mode_within_K = False
        self.plot_flag_mode_across_K_multipartite = True
        self.plot_flag_mode_across_K_chains = False
        
        
    def display(self):
        disp = []
        disp.append("========== [Parameters] ========== ")
        disp.append("Input path: {}".format(self.input_path))
        disp.append("Input data type: {}".format(self.input_type))
        disp.append("Output path: {}".format(self.output_path))
        disp.append("Mode detection method: {}".format("community detection" if not self.lc_flag else "leader clustering"))
        if not self.lc_flag:     
            disp.append("Using default community detection method: {}".format(self.default_cd))
            if self.default_cd:
                disp.append("Community detection modularity threshold: {}".format(self.cd_mod_thre))
        else:
            disp.append("Using adaptive threshold in leader clustering: {}".format(self.adaptive_thre_flag))
            if not self.adaptive_thre_flag:
                disp.append("Fixed threshold in leader clustering: {}".format(self.lc_cost_thre))
        disp.append("Enumerating all combinations of two clusters when aligning two replicates with K differing by one: {}".format(self.enum_combK))
        disp.append("Using customized colormap: {}".format(self.custom_cmap))
        if self.custom_cmap:
            disp.append("Customized colormap: {}".format(", ".join(self.cmap)))
        disp.append("Plotting all modes: {}".format(self.plot_flag_all_modes))
        disp.append("Plotting all replicates within each mode: {}".format(self.plot_flag_all_within_mode))
        disp.append("Plotting all modes of the same K: {}".format(self.plot_flag_mode_within_K))
        disp.append("Plotting alignment relationships between modes across K: {}".format(self.plot_flag_mode_across_K_multipartite))
        disp.append("Plotting the optimally aligned modes across K: {}".format(self.plot_flag_mode_across_K_chains))
        disp.append("")

        return "\n".join(disp)
