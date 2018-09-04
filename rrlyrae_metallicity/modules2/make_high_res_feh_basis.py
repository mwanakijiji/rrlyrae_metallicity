import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class lit_metallicities():
    
    ##############################################################################
    # STEP X: READ IN LITERATURE METALLICITY VALUES AND RESCALE
    ##############################################################################
    
    def __init__(self):
    
        stem = "./rrlyrae_metallicity/src/high_res_feh/"
        
        # Fe/H from Layden+ 1994
        self.layden_feh = pd.read_csv(stem + "layden_1994_abundances.dat",delimiter=';')
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
        self.fernley97_feh = pd.read_csv(stem + "fernley_1997_abundances.dat",delimiter=';')
        # RES: 60,000, two FeII lines, 5900-8100 A

        # Fe/H from Solano+ 1997
        self.solano_feh = pd.read_csv(stem + "solano_1997_abundances.dat",delimiter=';')
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
        self.chadid_feh['feh'] = np.mean([self.chadid_feh[' fehI'].values,self.chadid_feh[' fehII'].values],axis=0)
        
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

        ## ##
        #self.stem = "~/Documents/PythonPrograms/all_Python_code/2018_03_31_rrlyrae_rescale_a_la_chadid/"
        # log(eps) from Lambert+ 1996
        #lambert_logeps = pd.read_csv(self.stem + "lambert_1996_abundances.dat")
        # RES: ~23,000, FeII + photometric models, 3600-9000 A
        
    # fcn: find stars that overlap with Layden 1994, and return (x,y,z)=(FeH_Lay94,FeH_input-FeH_Lay94,starname)
    def find_match_Layden(self, input_table, layden_table, plot_name, offset=False):
    
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
        print('Y_offset to add to residuals in order to overlap with Chadid+ 2017 at Fe/H=-1.25:')
        print(net_offset)
        print('Number of overlapping stars:')
        print(len(residuals))
        line_offset = np.add(line,net_offset)
    
        # save a plot
        '''
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
        '''
            
        # return 
        # 1. overlapping Layden94 values
        # 2. FeH values from lit source
        # 3. Residuals between 1. and 2. (see Chadid+ 2017 ApJ 835:187, Figs. 5, 6, 7)
        # 4. coefficients of best-fit line
        # 5. offset in y to bring lit FeH values to match Chadid+ 2017 at FeH=-1.25 (see Chadid+ 2017 Figs. 5, 6)
        # 6. Residuals (from 3.) minus the offset (from 5.)  (see Chadid+ 2017 Fig. 7)
        # 7. The names of the stars (in same order as arrays for 1., 2., 3., 4.)
        
        d = dict()
        d['laydenFeH'] = laydenFeH
        d['inputFeH'] = inputFeH
        d['residuals'] = residuals
        d['coeff'] = coeff
        d['net_offset'] = net_offset # this needs to be ADDED to high-res study data to make it match Chadid
        d['residuals_shifted'] = np.add(residuals,net_offset)
        d['name'] = nameArray
        
        return d


def make_basis():
    
    # find the matches between datasets and apply the offsets

    # find matches: Fernley 1996 ## ## WAIT-- FERNLEY 97 INCLUDES THESE
    #dict_Fernley_96 = lit_metallicities.find_match_Layden(fernley_feh,layden_feh,'Fernley_96', offset=True)

    x = lit_metallicities() # instantiate the class

    # find matches: Lambert 1996
    dict_Lambert_96 = x.find_match_Layden(x.lambert_logeps,x.layden_feh,'Lambert_96', offset=True)

    # find matches: Nemec 2013
    dict_Nemec_2013  = x.find_match_Layden(x.nemec_feh,x.layden_feh,'Nemec_2013', offset=True)
    
    # find matches: Liu 2013
    x.liu_feh2 = x.liu_feh.groupby(x.liu_feh['name'], axis=0, as_index=False).mean()
    dict_Liu_2013  = x.find_match_Layden(x.liu_feh2,x.layden_feh,'Liu_2013', offset=True)
    
    # find matches: Chadid 2017
    dict_Chadid_2017  = x.find_match_Layden(x.chadid_feh,x.layden_feh,'Chadid_2017', offset=True)

    # find matches: Fernley 1997
    dict_Fernley_1997  = x.find_match_Layden(x.fernley97_feh,x.layden_feh,'Fernley_1997', offset=True)

    # find matches: Solano 1997
    dict_Solano_1997  = x.find_match_Layden(x.solano_feh,x.layden_feh,'Solano_1997', offset=True)

    ## ## SNEDEN_17 DOES NOT OVERLAP WITH LAYDEN! FIX THIS

    # find matches: Wallerstein+ 2010
    dict_Wallerstein_2010  = x.find_match_Layden(x.wallerstein_feh,x.layden_feh,'Wallerstein_2010', offset=True)

    ## ## IS THE BELOW NEEDED?
    # find matches between Wallerstein and Chadid
    # Chadid stars that appear in Wallerstein
    #chadid_winnow = x.chadid_feh[chadid_feh['star'].isin(wallerstein_feh['star'])]
    # Wallerstein stars that appear in Chadid
    #wallerstein_winnow = x.wallerstein_feh[wallerstein_feh['star'].isin(chadid_feh['star'])]
    
    # merge the metallicity dictionaries
    dict_collect = [dict_Lambert_96, dict_Nemec_2013, dict_Liu_2013, dict_Chadid_2017, 
            dict_Fernley_1997, dict_Solano_1997, dict_Wallerstein_2010]
    dict_merged = {}
    for key in dict_Lambert_96:
        dict_merged[key] = tuple(dict_merged[key] for dict_merged in dict_collect)
        
    # plot merged data and fit linreg line
    m_merged, b_merged = np.polyfit(np.hstack(dict_merged['laydenFeH']), np.hstack(dict_merged['residuals_shifted']), 1)

    
    ### make a plot like Chadid+ 2017 Fig. 5 (residuals between high res study FeHs and Layden94 FeH vs. Layden94 FeH)
    
    plt.clf() # clear plot space
    f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, figsize=(10,25), sharex=True)
    #
    quantx = np.array(dict_Fernley_1997['laydenFeH'],dtype=float)
    quanty = np.array(dict_Fernley_1997['residuals'], dtype=float)
    ax1.scatter(quantx, quanty)
    ax1.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax1_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax1.plot(quantx, np.add(np.multiply(quantx,ax1_coeff[0]),ax1_coeff[1]), linestyle='--') # regression line
    ax1.set_title('Fernley 97')
    #
    quantx = np.array(dict_Lambert_96['laydenFeH'],dtype=float)
    quanty = np.array(dict_Lambert_96['residuals'], dtype=float)
    ax2.scatter(quantx, quanty)
    ax2.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax2_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax2.plot(quantx, np.add(np.multiply(quantx,ax2_coeff[0]),ax2_coeff[1]), linestyle='--') # regression line
    ax2.set_title('Lambert 96')
    #
    quantx = np.array(dict_Nemec_2013['laydenFeH'],dtype=float)
    quanty = np.array(dict_Nemec_2013['residuals'], dtype=float)
    ax3.scatter(quantx, quanty)
    ax3.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax3_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax3.plot(quantx, np.add(np.multiply(quantx,ax3_coeff[0]),ax3_coeff[1]), linestyle='--') # regression line
    ax3.set_title('Nemec 13')
    #
    quantx = np.array(dict_Liu_2013['laydenFeH'],dtype=float)
    quanty = np.array(dict_Liu_2013['residuals'], dtype=float)
    ax4.scatter(quantx, quanty)
    ax4.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax4_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax4.plot(quantx, np.add(np.multiply(quantx,ax4_coeff[0]),ax4_coeff[1]), linestyle='--') # regression line
    ax4.set_title('Liu 13')
    #
    quantx = np.array(dict_Chadid_2017['laydenFeH'],dtype=float)
    quanty = np.array(dict_Chadid_2017['residuals'], dtype=float)
    ax5.scatter(quantx, quanty)
    ax5.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax5_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax5.plot(quantx, np.add(np.multiply(quantx,ax5_coeff[0]),ax5_coeff[1]), linestyle='--') # regression line
    ax5.set_title('Chadid 17')
    #
    quantx = np.array(dict_Solano_1997['laydenFeH'],dtype=float)
    quanty = np.array(dict_Solano_1997['residuals'], dtype=float)
    ax6.scatter(quantx, quanty)
    ax6.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax6_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax6.plot(quantx, np.add(np.multiply(quantx,ax6_coeff[0]),ax6_coeff[1]), linestyle='--') # regression line
    ax6.set_title('Solano 97')
    #
    quantx = np.array(dict_Wallerstein_2010['laydenFeH'],dtype=float)
    quanty = np.array(dict_Wallerstein_2010['residuals'], dtype=float)
    ax7.scatter(quantx, quanty)
    ax7.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax7_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax7.plot(quantx, np.add(np.multiply(quantx,ax7_coeff[0]),ax7_coeff[1]), linestyle='--') # regression line
    ax7.set_title('Wallerstein 10')
    ax7.set_xlabel('FeH_Layden94')
    #
    ax1.set_xlim([-2.9,0.4])
    ax1.set_ylim([-0.6,0.6])
    ax2.set_ylim([-0.6,0.6])
    ax3.set_ylim([-0.6,0.6])
    ax4.set_ylim([-0.6,0.6])
    ax5.set_ylim([-0.6,0.6])
    ax6.set_ylim([-0.6,0.6])
    ax7.set_ylim([-0.6,0.6])
    
    plt.ylabel('FeH_highres')
    #plt.title('chadid_fig5')
    plt.tight_layout()
    plt.savefig('chadid_fig5_imitation_test.pdf')
    
    plt.clf()

    ### make a plot like Chadid+ 2017 Fig. 6 (same as Fig. 5, but by shifting in y to match Chadid at FeH_Layden94=-1.25; except that Chadid+ 17 also
    ### reprocesses Clementini and Pancino FeH, which we dont do)
    plt.clf() # clear plot space
    f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, figsize=(10,25), sharex=True)
    ax1.scatter(dict_Fernley_1997['laydenFeH'], dict_Fernley_1997['residuals_shifted'])
    ax1.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax1.set_title('Fernley 97')
    ax2.scatter(dict_Lambert_96['laydenFeH'], dict_Lambert_96['residuals_shifted'])
    ax2.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax2.set_title('Lambert 96')
    ax3.scatter(dict_Nemec_2013['laydenFeH'], dict_Nemec_2013['residuals_shifted'])
    ax3.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax3.set_title('Nemec 13')
    ax4.scatter(dict_Liu_2013['laydenFeH'], dict_Liu_2013['residuals_shifted'])
    ax4.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax4.set_title('Liu 13')
    ax5.scatter(dict_Chadid_2017['laydenFeH'], dict_Chadid_2017['residuals_shifted'])
    ax5.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax5.set_title('Chadid 17')
    ax6.scatter(dict_Solano_1997['laydenFeH'], dict_Solano_1997['residuals_shifted'])
    ax6.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax6.set_title('Solano 97')
    ax7.scatter(dict_Wallerstein_2010['laydenFeH'], dict_Wallerstein_2010['residuals_shifted'])
    ax7.plot(dict_Fernley_1997['laydenFeH'], np.zeros(np.shape(dict_Fernley_1997['laydenFeH'])), linestyle=':') # zero line
    ax7.set_title('Wallerstein 10')
    ax7.set_xlabel('FeH_Layden94')
    ax1.set_xlim([-2.9,0.4])
    ax1.set_ylim([-0.6,0.6])
    ax2.set_ylim([-0.6,0.6])
    ax3.set_ylim([-0.6,0.6])
    ax4.set_ylim([-0.6,0.6])
    ax5.set_ylim([-0.6,0.6])
    ax6.set_ylim([-0.6,0.6])
    ax7.set_ylim([-0.6,0.6])
    plt.ylabel('FeH_highres_residuals_shifted')
    #plt.title('chadid_fig6')
    plt.tight_layout()
    plt.savefig('chadid_fig6_imitation_test.pdf')
    plt.clf()

    ### make a plot like Chadid+ 2017 Fig. 7
    plt.clf() # clear plot space
    f, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,9), sharex=True)
    # subplot 1
    to_plot_x_all = dict_merged['laydenFeH'] # consolidate stuff to find regression coefficients
    to_plot_y_all = np.add(dict_merged['residuals_shifted'],dict_merged['laydenFeH'])
    [ax1.scatter(dict_merged['laydenFeH'][p], np.add(dict_merged['residuals_shifted'][p],dict_merged['laydenFeH'][p])) for p in range(0,len(dict_merged['laydenFeH']))]
    coeff2 = np.polyfit(np.concatenate(to_plot_x_all).ravel(), np.concatenate(to_plot_y_all).ravel(), 1)
    ax1.plot(np.arange(-3,1,step=0.1),np.arange(-3,1,step=0.1),linestyle='--') # one-to-one
    ax1.set_xlim([-2.9,0.4])
    ax1.set_ylim([-2.9,0.4])
    ax1.set_ylabel('FeH_highres_shifted')
    ax1.set_title('m = '+str(coeff2[0])+', b = '+str(coeff2[1]))
    # subplot 2
    [ax2.scatter(dict_merged['laydenFeH'][p], dict_merged['residuals_shifted'][p]) for p in range(0,len(dict_merged['laydenFeH']))]
    reg_x, reg_y = zip(*sorted(zip(np.hstack(dict_merged['laydenFeH']), np.add(np.multiply(np.hstack(dict_merged['laydenFeH']),m_merged),b_merged)))) # sort values in x, so that plot line doesn't zigzag
    ax2.plot(reg_x, reg_y, linestyle='--') # regression line
    ax2.plot(reg_x, np.zeros(np.shape(reg_y)), linestyle=':') # zero line
    ax2.set_ylim([-0.6,0.6])
    ax2.set_xlabel('FeH_Layden94')
    ax2.set_ylabel('FeH_highres_residuals_shifted')
    plt.suptitle('chadid_fig7')
    plt.savefig('chadid_fig7_imitation_test.pdf')
    plt.clf()

    
    ## ## CAUTION: TEST TO SEE IF THE CONTENT IN THE KEYS IS IN ORDER (I.E., MAKE A PLOT AND SEE IF ITS THE SAME IF DATASETS ARE OVERLAID INDIVIDUALLY)

    # rescale_lit_metallicities to find high-res Fe/H
