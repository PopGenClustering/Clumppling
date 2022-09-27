
class Params:
    def __init__(self, input_path,output_path,params_path,input_format="structure"):
        
        self.input_path = input_path
        self.output_path = output_path
        self.params_path = params_path
        self.input_format = input_format
        
        self.load_from_xml(self.params_path)

        # # mode detection parameters
        # self.default_cd = True
        # self.cd_mod_thre = -1

        # # colormap parameters
        # self.custom_cmap = False
        # self.cmap = ""

        # # Advanced: ILP when number of clusters differ by 1
        # self.enum_combK = True

        # # Advanced: leader clustering parameters
        # self.lc_flag = False     
        # self.adaptive_thre_flag = False
        # self.lc_cost_thre = 0.2
        
        # # Advanced: plotting parameters
        # self.plot_flag_all_modes = True
        # self.plot_flag_all_within_mode = False
        # self.plot_flag_mode_within_K = False
        # self.plot_flag_mode_across_K_multipartite = True
        # self.plot_flag_mode_across_K_chains = False
    

    def load_from_xml(self,path):

        # xml parser
        xml = __import__('xml.etree.ElementTree')
        params_xml = xml.etree.ElementTree.parse(path).getroot() 

        def find_element(param_name):
            elem = params_xml.find(param_name)
            return elem.text 
        
        def parse_bool(var):
            return var.lower() in ['true', 't', 'y', 'yes', '1']

        # mode detection parameters
        self.default_cd = parse_bool(find_element('default_cd'))
        self.cd_mod_thre = float(find_element('cd_mod_thre'))

        # colormap parameters
        self.custom_cmap = parse_bool(find_element('custom_cmap'))
        cmap = find_element('cmap')
        self.cmap = str(cmap) if cmap else ""

        # Advanced: ILP when number of clusters differ by 1
        self.enum_combK = parse_bool(find_element('enum_combK'))

        # Advanced: leader clustering parameters
        self.lc_flag = parse_bool(find_element('lc_flag'))     
        self.adaptive_thre_flag = parse_bool(find_element('adaptive_thre_flag'))    
        self.lc_cost_thre = float(find_element('lc_cost_thre'))
        
        # Advanced: plotting parameters
        self.plot_flag_all_modes = parse_bool(find_element('plot_flag_all_modes'))    
        self.plot_flag_all_within_mode = parse_bool(find_element('plot_flag_all_within_mode'))    
        self.plot_flag_mode_within_K = parse_bool(find_element('plot_flag_mode_within_K'))    
        self.plot_flag_mode_across_K_multipartite = parse_bool(find_element('plot_flag_mode_across_K_multipartite'))    
        self.plot_flag_mode_across_K_chains = parse_bool(find_element('plot_flag_mode_across_K_chains'))    
 
        
    def display(self):
        disp = []
        disp.append("========== [Parameters] ========== ")
        disp.append("Input path: {}".format(self.input_path))
        disp.append("Input data format: {}".format(self.input_format))
        disp.append("Output path: {}".format(self.output_path))
        disp.append("Parameter file path: {}".format(self.params_path))
        disp.append("")
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
