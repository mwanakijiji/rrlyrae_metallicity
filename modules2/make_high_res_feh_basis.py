class lit_metallicities():
    
    ##############################################################################
    # STEP X: READ IN LITERATURE METALLICITY VALUES AND RESCALE
    ##############################################################################
    
    def __init__(self):
    
        stem = "~/Documents/PythonPrograms/all_Python_code/2018_03_31_rrlyrae_rescale_a_la_chadid/"
        
        # Fe/H from Layden+ 1994
        layden_feh = pd.read_csv(stem + "layden_1994_abundances.dat",delimiter=';')
        # RES: "rather low"
        
        # Fe/H Clementini+ 1995
        clementini_feh = pd.read_csv(stem + "clementini_1995_abundances.dat")

        # Fe/H Fernley+ 1996
        fernley_feh = pd.read_csv(stem + "fernley_1996_abundances.dat")
        # RES: 60,000, FeI & FeII, 5900-8100 A
        
        # log(eps) from Lambert+ 1996
        lambert_logeps = pd.read_csv(stem + "lambert_1996_abundances.dat")
        # RES: ~23,000, FeII + photometric models, 3600-9000 A
        
        # Fe/H from Wallerstein and Huang 2010, arXiv 1004.2017
        wallerstein_feh = pd.read_csv(stem + "wallerstein_huang_2010_abundances.dat")
        # RES: ~30,000, FeII
        
        # Fe/H from Chadid+ 2017 (FeI and II lines)
        chadid_feh = pd.read_csv(stem + "chadid_2017_abundances.dat")
        # RES: 38000, FeI & FeII, 3400-9900 A

        # Fe/H from Liu+ 2013
        liu_feh = pd.read_csv(stem + "liu_2013_abundances.dat")
        # RES: ~60,000, FeI (& FeII?), 5100-6400 A

        # Fe/H from Nemec+ 2013
        nemec_feh = pd.read_csv(stem + "nemec_2013_abundances.dat")
        # RES: ~65,000 or 36,000, FeI & FeII, 5150-5200 A

        # Fe/H from Fernley+ 1997
        fernley97_feh = pd.read_csv(stem + "fernley_1997_abundances.dat",delimiter=';')
        # RES: 60,000, two FeII lines, 5900-8100 A

        # Fe/H from Solano+ 1997
        solano_feh = pd.read_csv(stem + "solano_1997_abundances.dat",delimiter=';')
        # RES: 22,000 & 19,000, strong FeI lines, 4160-4390 & 4070-4490 A
        
        # Fe/H from Pacino+ 2015
        pacino_feh = pd.read_csv(stem + "pacino_2015_abundances.dat") 
        # RES: >30,000, FeI (weighted average), 4000-8500 A

        # Fe/H from Sneden+ 2017
        sneden_feh = pd.read_csv(stem + "sneden_2017_abundances.dat")
        # RES: ~27,000 (at 5000 A), FeI & FeII, 3400-9000 A
        
        # convert Lambert's values, which are in terms of log(eps)
        # FeH = log(epsFe) - log(epsFe,sol)
        #     = log(epsFe) - log(NFe,sol/NH,sol)
        #     = log(epsFe) - 7.51 # value of 7.51 from Anstee+ 1997, MNRAS
        lambert_logeps['feh'] = np.subtract(lambert_logeps['log_eps_fe_spec'], 7.51) 
        
        # average the values in Chadid from FeI and FeII lines
        chadid_feh['feh'] = np.mean([chadid_feh[' fehI'].values,chadid_feh[' fehII'].values],axis=0)
        
        ## ## INCLUDE SINGLE DATA PT FROM KOLENBERG+ 2010? (SEE CHADID+ 2017, FIG. 7)
        
        # FYI: average Fe/H values in Liu+ 2013 which were taken at different phases
        # liu_feh.groupby(liu_feh['name'], axis=0, as_index=False).mean()
        
        # FYI: average Fe/H values in Sneden+ 1997 which were taken at different epochs
        # sneden_feh.groupby(sneden_feh['name'], axis=0, as_index=False).mean()
        
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
        
        
    # fcn: find stars that overlap with Layden 1994, and return (x,y,z)=(FeH_Lay94,FeH_input-FeH_Lay94,starname)
    def find_match_Layden(input_table, layden_table, plot_name, offset=False):
    
        inputFeH = []
        laydenFeH = []
        nameArray = []
    
        for row in range(0,len(input_table)): # scan over each row in input table
            if (layden_table['name'] == input_table['name'][row]).any():
            #if (input_table['name'][row] isin layden_table['name']): # if there's a star name that matches
                inputFeH = np.append(inputFeH,input_table['feh'][row])
                laydenFeH = np.append(laydenFeH,layden_table.loc[layden_table['name'] == input_table['name'][row]]['feh'])
                nameArray = np.append(nameArray,input_table['name'][row])
    
        residuals = np.subtract(inputFeH,laydenFeH)
    
        # best-fit line
        coeff = np.polyfit(laydenFeH, residuals, 1)
        limits = [-3.0,0.5]
        line = np.multiply(coeff[0],limits)+coeff[1]
    
        # if there needs to be an offset (like in Fig. 6 of Chadid+ 2017)
        chadid_y_125 = -0.10583621694962 # from Chadid line at Fe/H=-1.25
        this_y_125 = np.multiply(coeff[0],-1.25)+coeff[1] # y-value of this line at Fe/H=-1.25
        net_offset = chadid_y_125 - this_y_125 # offset needed to move line
        print('Y_offset to overlap with Chadid+ 2017 at Fe/H=-1.25:')
        print(net_offset)
        print('Number of overlapping stars:')
        print(len(residuals))
        line_offset = np.add(line,net_offset)
    
        # save a plot
        plt.scatter(laydenFeH, np.subtract(inputFeH,laydenFeH))
        plt.plot([-3.0,0.5], [0., 0.], linestyle='--')
        plt.plot(limits, line)
        plt.plot(limits, line_offset)
        #plt.xlim([-3.0,0.5])
        #plt.ylim([-0.6,0.6])
        plt.xlabel('[Fe/H]_Lay94')
        plt.ylabel('[Fe/H]_input - [Fe/H]_Lay94')
        plt.title('residuals between '+str(plot_name)+' and Lay94\ny=mx+b, m='+str(coeff[0])+', b='+str(coeff[1])+'\n offset '+str(net_offset))
        plt.savefig(plot_name+'_test_180708.png')
        #plt.show()
        plt.clf()
            
        # return 
        # 1. overlapping Layden94 values
        # 2. FeH values from lit source
        # 3. Residuals between 1. and 2.(see Chadid+ 2017 Figs. 5, 6, 7)
        # 4. coefficients of best-fit line
        # 5. offset in y to bring lit FeH values to match Chadid+ 2017 at FeH=-1.25 (see Chadid+ 2017 Figs. 5, 6)
        # 6. Residuals (from 3.) minus the offset (from 5.)  (see Chadid+ 2017 Fig. 7)
        # 7. The names of the stars (in same order as arrays for 1., 2., 3., 4.)
        
        d = dict()
        d['laydenFeH'] = laydenFeH
        d['inputFeH'] = inputFeH
        d['residuals'] = residuals
        d['coeff'] = coeff
        d['net_offset'] = net_offset
        d['residuals_shifted'] = np.subtract(inputFeH,residuals)
        d['name'] = nameArray
        
        return d


# now actually find the matches between datasets and apply the offsets

# find matches: Fernley 1996 ## ## WAIT-- FERNLEY 97 INCLUDES THESE
#dict_Fernley_96 = lit_metallicities.find_match_Layden(fernley_feh,layden_feh,'Fernley_96', offset=True)

# find matches: Lambert 1996
dict_Lambert_96 = lit_metallicities.find_match_Layden(lambert_logeps,layden_feh,'Lambert_96', offset=True)

# find matches: Nemec 2013
dict_Nemec_2013  = lit_metallicities.find_match_Layden(nemec_feh,layden_feh,'Nemec_2013', offset=True)

# find matches: Liu 2013
liu_feh2 = liu_feh.groupby(liu_feh['name'], axis=0, as_index=False).mean()
dict_Liu_2013  = lit_metallicities.find_match_Layden(liu_feh2,layden_feh,'Liu_2013', offset=True)

# find matches: Chadid 2017
dict_Chadid_2017  = lit_metallicities.find_match_Layden(chadid_feh,layden_feh,'Chadid_2017', offset=True)

# find matches: Fernley 1997
dict_Fernley_1997  = lit_metallicities.find_match_Layden(fernley97_feh,layden_feh,'Fernley_1997', offset=True)

# find matches: Solano 1997
dict_Solano_1997  = lit_metallicities.find_match_Layden(solano_feh,layden_feh,'Solano_1997', offset=True)

## ## SNEDEN_17 DOES NOT OVERLAP WITH LAYDEN! FIX THIS

# find matches: Wallerstein+ 2010
dict_Wallerstein_2010  = lit_metallicities.find_match_Layden(wallerstein_feh,layden_feh,'Wallerstein_2010', offset=True)

## ## IS THE BELOW NEEDED?
# find matches between Wallerstein and Chadid
# Chadid stars that appear in Wallerstein
#chadid_winnow = lit_metallicities.chadid_feh[chadid_feh['star'].isin(wallerstein_feh['star'])]
# Wallerstein stars that appear in Chadid
#wallerstein_winnow = lit_metallicities.wallerstein_feh[wallerstein_feh['star'].isin(chadid_feh['star'])]


# merge the metallicity dictionaries

dict_collect = [dict_Lambert_96, dict_Nemec_2013, dict_Liu_2013, dict_Chadid_2017, 
            dict_Fernley_1997, dict_Solano_1997, dict_Wallerstein_2010]
dict_merged = {}
for k in dict_Lambert_96.iterkeys():
    dict_merged[k] = tuple(dict_merged[k] for dict_merged in dict_collect)
    
## ## CAUTION: TEST TO SEE IF THE CONTENT IN THE KEYS IS IN ORDER (I.E., MAKE A PLOT AND SEE IF ITS THE SAME IF DATASETS ARE OVERLAID INDIVIDUALLY)

# rescale_lit_metallicities to find high-res Fe/H
