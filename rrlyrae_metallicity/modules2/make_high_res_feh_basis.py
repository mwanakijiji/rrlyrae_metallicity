import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

class LitMetallicities():
    '''
    Read in metallicities from high-res studies and basis data sets, and
    define matching of stars between data sets.
    '''
    
    def __init__(self):
    
        stem = "./rrlyrae_metallicity/src/high_res_feh/"

        
        ####################################
        #### calibration program stars #####

        # stand-in that consists of our program star names
        self.our_program_stars = pd.read_csv(stem + "our_program_stars_names_only.csv")

        
        #######################################################
        #### basis set and high-res studies of RRab stars #####
        # (N.b. some RRcs are among them too, but basis set should only be RRabs)
        
        # Fe/H from Layden+ 1994; this serves as the common basis
        self.layden_feh = pd.read_csv(stem + "layden_1994_abundances.dat")
        # RES: "rather low"
        
        # Fe/H Clementini+ 1995
        self.clementini_feh = pd.read_csv(stem + "clementini_1995_abundances.dat")

        # Fe/H Fernley+ 1996
        self.fernley_feh = pd.read_csv(stem + "fernley_1996_abundances.dat")
        # RES: 60,000, FeI & FeII, 5900-8100 A
        
        # log(eps) from Lambert+ 1996
        self.lambert_logeps = pd.read_csv(stem + "lambert_1996_abundances.dat")
        # RES: ~23,000, FeII + photometric models, 3600-9000 A
        
        # Fe/H from Wallerstein and Huang 2010, arXiv 1004.2017
        self.wallerstein_feh = pd.read_csv(stem + "wallerstein_huang_2010_abundances.dat")
        # RES: ~30,000, FeII
        
        # Fe/H from Chadid+ 2017 ApJ 835.2:187 (FeI and II lines)
        self.chadid_feh = pd.read_csv(stem + "chadid_2017_abundances.dat")
        # RES: 38000, FeI & FeII, 3400-9900 A

        # Fe/H from Liu+ 2013 Res Ast Astroph 13:1307
        self.liu_feh = pd.read_csv(stem + "liu_2013_abundances.dat")
        # RES: ~60,000, FeI (& FeII?), 5100-6400 A

        # Fe/H from Nemec+ 2013
        self.nemec_feh = pd.read_csv(stem + "nemec_2013_abundances.dat")
        # RES: ~65,000 or 36,000, FeI & FeII, 5150-5200 A

        # Fe/H from Fernley+ 1997
        self.fernley97_feh = pd.read_csv(stem + "fernley_1997_abundances.dat")
        # RES: 60,000, two FeII lines, 5900-8100 A

        # Fe/H from Solano+ 1997
        self.solano_feh = pd.read_csv(stem + "solano_1997_abundances.dat")
        # RES: 22,000 & 19,000, strong FeI lines, 4160-4390 & 4070-4490 A
        
        # Fe/H from Pancino+ 2015 MNRAS 447:2404
        self.pacino_feh = pd.read_csv(stem + "pacino_2015_abundances.dat") 
        # RES: >30,000, FeI (weighted average), 4000-8500 A

        # Fe/H from Sneden+ 2017
        self.sneden_feh = pd.read_csv(stem + "sneden_2017_abundances.dat")
        # RES: ~27,000 (at 5000 A), FeI & FeII, 3400-9000 A
        
        # convert Lambert's values, which are in terms of log(eps)
        # FeH = log(epsFe) - log(epsFe,sol)
        #     = log(epsFe) - log(NFe,sol/NH,sol)
        #     = log(epsFe) - 7.51 # value of 7.51 from Anstee+ 1997, MNRAS
        self.lambert_logeps['feh'] = np.subtract(self.lambert_logeps['log_eps_fe_spec'], 7.51) 
        
        # average the values in Chadid from FeI and FeII lines
        self.chadid_feh['feh'] = np.mean([self.chadid_feh['fehI'].values,self.chadid_feh['fehII'].values],axis=0)

        # note also that Sneden+ 1997, Liu+ 2013, and 
        ## ## INCLUDE SINGLE DATA PT FROM KOLENBERG+ 2010? (SEE CHADID+ 2017, FIG. 7)
        
        # FYI: average Fe/H values in Liu+ 2013 which were taken at different phases
        # liu_feh.groupby(liu_feh['name'], axis=0, as_index=False).mean()
        
        # FYI: average Fe/H values in Sneden+ 1997 which were taken at different epochs
        # sneden_feh.groupby(sneden_feh['name'], axis=0, as_index=False).mean()

        
        #######################################################
        #### basis set and high-res studies of RRc stars ######
        
        # Fe/H from Kemper+ 1982; this serves as the common basis
        self.kemper_feh = pd.read_csv(stem + "kemper_1982_abundances.dat")

        # Fe/H from Govea+ 2014
        ## ## note: Govea+ has abundances for each phase value, and this includes NLTE phases; how to get single Fe/H?
        self.govea_feh = pd.read_csv(stem + "govea_2014_abundances.dat")
        

        #####################
        
        # initialize arrays: essential info
        empir_spec_name_array = []
        star_name_array = []
        H_data_array = []
        K_data_array = []
        err_H_data_array = [] 
        err_K_data_array = []

        # initialize arrays: other info
        Hbet_data_array = []
        err_Hbet_data_array = []
        Hgam_data_array = []
        err_Hgam_data_array = []
        rHgam_data_array = [] # rescaled Hgamma
        err_rHgam_data_array = []
        Hdel_data_array = []
        err_Hdel_data_array = []
        Heps_data_array = []
        err_Heps_data_array = []

        
    def __call__(self):
        
        # make a list of all UNIQUE, EMPIRICAL spectrum names
        uniqueSpecNames = line_data.drop_duplicates(subset='empir_spec_name')['empir_spec_name']

        
    def matchmaker(self, input_table, basis_table):
        '''
        Find what stars are common to the two input tables, and return array of FeHs from the first and basis tables

        INPUTS:
        input_table: table I'm interested in checking for overlapping stars
        basis_table: table with the names for which I am looking for repeats in the other table (most likely our own program stars)

        OUTPUTS:
        dictionary with
        1. overlapping star names
        2. FeHs from the input_table
        3. FeHs from the basis_table
        4. residuals in FeH: FeH_input - FeH_basis
        '''

        self.input_table = input_table
        self.basis_table = basis_table
            
        input_FeH = [] # Fe/H of high-res study
        basis_FeH = [] # Fe/H of basis (ex. Layden 1994)
        name_array = [] # name of star

        for row in range(0,len(input_table)): # scan over each row in input table
            if (basis_table['name'] == input_table['name'][row]).any():
                input_FeH = np.append(input_FeH,input_table['feh'][row])
                basis_FeH = np.append(basis_FeH,basis_table.loc[basis_table['name'] == input_table['name'][row]]['feh'])
                name_array = np.append(name_array,input_table['name'][row])

        d = dict()
        d['name'] = name_array
        d['input_FeH'] = input_FeH
        d['basis_FeH'] = basis_FeH
        d['residuals_FeH'] = np.subtract(d['input_FeH'],d['basis_FeH'])

        ## ##
        '''
        if self.__offset: # if we want to apply Chadid+ 2017-style offsets to form a common Fe/H basis
            print('PLACEHOLDER')
            # if there needs to be an offset (like in Fig. 6 of Chadid+ 2017)
            chadid_y_125 = -0.10583621694962 # from Chadid line at Fe/H=-1.25
            this_y_125 = np.multiply(coeff[0],-1.25)+coeff[1] # y-value of this line at Fe/H=-1.25
            net_offset = chadid_y_125 - this_y_125 # offset needed to move line
            print('Y_offset to add to residuals in order to overlap with Chadid+ 2017 at Fe/H=-1.25:')
            print(net_offset)
            print('Number of overlapping stars:')
            print(len(residuals))
            line_offset = np.add(line,net_offset)
        '''    
            
        return d

    
class MetalBasisTypeSpecific(LitMetallicities):
    '''
    For a given RR Lyrae subtype, generate a metallicity basis and calculate
    the abundances of stars in the calibration program dataset.
    '''

    def __init__(self, plot_name, star_type="RRab", offset=False):
        super().__init__()
        self.__plot_name = plot_name
        self.__star_type = star_type
        self.__offset = offset

    def make_basis(self):
        '''
        Find what stars overlap with basis data set, and return star name, FeH values, residuals
                
        The functionality of LitMetallicities is inherited, and we just add Chadid+ 17-style offsets, a best-fit line, and plotting functionality

        INPUTS:
        input_table: table of likely high-res-derived Fe/H values of RRabs, which I want to cross-ref with Layden 94
        layden_table: the layden table, which serves as the basis set
        plot_name: file name for saving a plot of the results

        OUTPUTS:
        dictionary with
        1. overlapping star names
        2. Fe/Hs from the input_table
        3. Fe/Hs from Layden
        '''

        # define the basis data set (like Layden+ 1994 for RRabs, or Kemper+ 1982 for RRcs)
        if self.__star_type == "RRab":
            type_string = "ab"
            basis_set = self.layden_feh
            basis_string = "Layden RRab basis set" # string for plots
        elif self.__star_type == "RRc":
            type_string = "c"
            basis_set = self.kemper_feh
            basis_string = "Kemper RRc basis set"
        else:
            sys.exit("Error! No RR Lyrae subtype chosen.")
            
        # match high-res studies with the basis set
        dict_Lambert_1996 = self.matchmaker(self.lambert_logeps, basis_set) # Lambert+ 1996 (logeps has already been converted to Fe/H)
        dict_Nemec_2013 = self.matchmaker(self.nemec_feh, basis_set) # Nemec+ 2013
        dict_Chadid_2017 = self.matchmaker(self.chadid_feh, basis_set) # Chadid+ 2017
        dict_Fernley_1997 = self.matchmaker(self.fernley97_feh, basis_set) # Fernley+ 1997
        dict_Solano_1997 = self.matchmaker(self.solano_feh, basis_set) # Solano+ 1997
        dict_Wallerstein_2010 = self.matchmaker(self.wallerstein_feh, basis_set) # Wallerstein 2010
        
        # for Liu+ 2013, we need to group multiple Fe/H values by star name
        # (the grouping is done here rather than further up because a bug causes the grouped column to disappear)
        self.liu_feh_grouped = self.liu_feh.groupby(self.liu_feh['name'], axis=0, as_index=False).mean()
        dict_Liu_2013 = self.matchmaker(self.liu_feh_grouped, basis_set) # Liu+ 2013

        # merge dictionaries of Fe/H values
        dict_collect = [dict_Lambert_1996, dict_Nemec_2013, dict_Liu_2013, dict_Chadid_2017,
                        dict_Fernley_1997, dict_Solano_1997, dict_Wallerstein_2010]
        dict_merged = {}
        for key in dict_Lambert_1996:
            dict_merged[key] = tuple(dict_merged[key] for dict_merged in dict_collect)

        # rename some things for neatness
        basis_data_merged = np.hstack(dict_merged['basis_FeH'])
        highres_data_merged = np.hstack(dict_merged['input_FeH'])
        residuals_data_merged = np.hstack(dict_merged['residuals_FeH'])
        names_merged = np.hstack(dict_merged['name'])

        # find best-fit line (note that user may have used a flag to make Fe/H values be offset)
        limits = [-3.0,0.5]
        m_merged_highres, b_merged_highres = np.polyfit(basis_data_merged, highres_data_merged, 1)
        line_highres = np.multiply(m_merged_highres,limits)+b_merged_highres # make best-fit line for high-res Fe/H
        m_merged_resid, b_merged_resid = np.polyfit(basis_data_merged, residuals_data_merged, 1)
        line_resid = np.multiply(m_merged_resid,limits)+b_merged_resid # make best-fit line for residuals
            
        # save a plot (high_res vs. basis on top; residuals vs. basis on bottom)
        plt.clf()
        fig, axs = plt.subplots(2, 1, figsize=(10,10), sharex=True)
        axs[0].plot([limits[0],limits[1]],[limits[0],limits[1]], linestyle='--') # make 1-to-1 line
        axs[0].plot([limits[0],limits[1]],np.add(np.multiply(m_merged_highres,[limits[0],limits[1]]),b_merged_highres), linestyle='--') # best-fit line
        axs[0].scatter(basis_data_merged, highres_data_merged) # input vs. basis
        axs[0].set_xlim(limits[0], limits[1])
        axs[0].set_ylabel("Fe/H, high-res")
        axs[0].set_title("m = "+str(m_merged_highres)+", b = "+str(b_merged_highres)+"; (blue line: 1-to-1; orange line: best fit)")
        axs[1].axhline(y=0, linestyle='--') # dashed line at y=0
        axs[1].scatter(basis_data_merged, residuals_data_merged) # input vs. basis
        axs[1].set_xlabel("Fe/H, "+basis_string)
        axs[1].set_ylabel('Fe/H Residuals: high-res minus basis set')
        axs[1].set_title("m = "+str(m_merged_resid)+", b = "+str(b_merged_resid)+"; (blue line: zero)")
        fig.suptitle('Finding remapping relation between\nhigh-res studies and basis dataset\n('+type_string+' subtype)')
        #fig.tight_layout()
        plt.savefig('remapping_'+self.__plot_name, overwrite=True)
        plt.clf()
        
            
        # return 
        # 1. overlapping Layden94 values
        # 2. FeH values from lit source
        # 3. Residuals between 1. and 2. (see Chadid+ 2017 ApJ 835:187, Figs. 5, 6, 7)
        # 4. coefficients of best-fit line
        # 5. offset in y to bring lit FeH values to match Chadid+ 2017 at FeH=-1.25 (see Chadid+ 2017 Figs. 5, 6)
        # 6. Residuals (from 3.) minus the offset (from 5.)  (see Chadid+ 2017 Fig. 7)
        # 7. The names of the stars (in same order as arrays for 1., 2., 3., 4.)
        
        d = dict()
        d['basis_FeH'] = basis_data_merged
        d['input_FeH'] = highres_data_merged
        d['residuals'] = residuals_data_merged
        d['coeff_merged_highres'] = [m_merged_highres, b_merged_highres] # best-fit line coeffs for high-res vs. basis 
        d['coeff_merged_resid'] = [m_merged_resid, b_merged_resid] # best-fit line coeffs for (residuals: high-res minus basis) vs. basis 
        #d['net_offset'] = net_offset # this needs to be ADDED to high-res study data to make it match Chadid
        #d['residuals_shifted'] = np.add(residuals,net_offset)
        d['name'] = names_merged
        
        return d
    
    
    def calc_FeH_program_stars(self):
        '''
        Calculate metallicities for the program stars which form the basis of the
        metallicity calibration, by using the remapping relationships

        INPUTS:
        basis_set: basis set used for either RRab (such as Layden 1994) or RRc (such as Kemper+ 1982)
        '''

        # retrieve our own program stars and remap those of the right type
        if self.__star_type == "RRab":
            type_string = "ab"
            basis_set = self.layden_feh
            basis_string = "Layden RRab basis set" # string for plots
        elif self.__star_type == "RRc":
            type_string = "c"
            basis_set = self.kemper_feh
            basis_string = "Kemper RRc basis set"

        # retrieve our stars here, and extract only those which conform to the right type
        program_stars_subset = self.our_program_stars.loc[self.our_program_stars['type'] == type_string].reset_index()
        
        # find matches with the basis set
        program_stars_subset_matched = self.matchmaker(program_stars_subset, basis_set)

        # make the Fe/H basis and return the coefficients of the linear remapping
        map_info = self.make_basis()
        
        # remap metallicities via
        # [Fe/H]_highres = m*[Fe/H]_basis_set + b   
        program_stars_subset_matched['mapped_FeH'] = np.add(np.multiply(program_stars_subset_matched['basis_FeH'],
                                                                        map_info['coeff_merged_highres'][0]),
                                                            map_info['coeff_merged_highres'][1])
        
        print(program_stars_subset_matched.keys())
        # save a plot of calibration program stars Fe/H
        # post-mapped Fe/H vs. pre-mapped (i.e., basis set) Fe/H
        limits = [-3.0,0.5]
        plt.clf()
        fig, axs = plt.subplots(1, 1, figsize=(10,10))
        axs.plot([limits[0],limits[1]],[limits[0],limits[1]], linestyle='--') # make 1-to-1 line
        axs.scatter(program_stars_subset_matched['basis_FeH'], program_stars_subset_matched['mapped_FeH']) # input vs. basis
        axs.set_xlim(limits[0], limits[1])
        axs.set_ylabel("Fe/H, high-res")
        axs.set_xlabel("Fe/H, "+basis_string)
        #axs.set_title("m = "+str(m_merged_highres)+", b = "+str(b_merged_highres))

        fig.suptitle('Calculated Fe/H of calibration program stars\n('+type_string+' subtype)')
        #fig.tight_layout()
        plt.savefig('calculated_FeH_'+self.__plot_name, overwrite=True)
        plt.clf()
    
        '''
        # write out
        convert_to_df = pd.DataFrame.from_dict(dict_our_program_stars['name']) # initialize
        convert_to_df.columns = ['name'] # rename the column
        convert_to_df['mapped_feh'] = pd.DataFrame.from_dict(dict_our_program_stars['mapped_feh']) # add the remapped Fe/H
        no_return = convert_to_df.to_csv(write_loc + "mapped_feh.csv") # write out ## ## note 2 things: 1., this should be appeneded to our .csv with EWs; 2. there is no phase info here yet
        '''

    '''
    ## ##
    EXAMPLE COMMANDS IN THE HIGHER-LEVEL SCRIPT:
    test_rrab = MetalBasisTypeSpecific(plot_name='name_here',offset=True).calc_FeH_program_stars()
    test_rrc = MetalBasisTypeSpecific(plot_name='name_here',star_type="RRc").calc_FeH_program_stars()
    '''
