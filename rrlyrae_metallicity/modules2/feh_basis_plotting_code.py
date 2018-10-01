    ########################################
    # BEGIN BUNCH OF PLOTS FOR RRABS
    ########################################
    '''
    plt.clf()
    remapped_Lambert = remap_metal(dict_Lambert_96['laydenFeH'], dict_Lambert_96['inputFeH'])
    plt.scatter(dict_Lambert_96['laydenFeH'], dict_Lambert_96['inputFeH'], color='orange')
    plt.scatter(dict_Lambert_96['laydenFeH'], remapped_Lambert['vals_mapped'])
    plt.plot(dict_Lambert_96['laydenFeH'], dict_Lambert_96['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Lambert['coeffs']))
    labels = dict_Lambert_96['name']
    labels_x = dict_Lambert_96['laydenFeH']
    labels_y = dict_Lambert_96['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_lambert.pdf')
    
    plt.clf()
    remapped_Nemec = remap_metal(dict_Nemec_2013['laydenFeH'], dict_Nemec_2013['inputFeH'])
    plt.scatter(dict_Nemec_2013['laydenFeH'], dict_Nemec_2013['inputFeH'], color='orange')
    plt.scatter(dict_Nemec_2013['laydenFeH'], remapped_Nemec['vals_mapped'])
    plt.plot(dict_Nemec_2013['laydenFeH'], dict_Nemec_2013['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Nemec['coeffs']))
    labels = dict_Nemec_2013['name']
    labels_x = dict_Nemec_2013['laydenFeH']
    labels_y = dict_Nemec_2013['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_nemec.pdf')

    plt.clf()
    remapped_Liu = remap_metal(dict_Liu_2013['laydenFeH'], dict_Liu_2013['inputFeH'])
    plt.scatter(dict_Liu_2013['laydenFeH'], dict_Liu_2013['inputFeH'], color='orange')
    plt.scatter(dict_Liu_2013['laydenFeH'], remapped_Liu['vals_mapped'])
    plt.plot(dict_Liu_2013['laydenFeH'], dict_Liu_2013['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Liu['coeffs']))
    labels = dict_Liu_2013['name']
    labels_x = dict_Liu_2013['laydenFeH']
    labels_y = dict_Liu_2013['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_liu.pdf')

    plt.clf()
    remapped_Chadid = remap_metal(dict_Chadid_2017['laydenFeH'], dict_Chadid_2017['inputFeH'])
    plt.scatter(dict_Chadid_2017['laydenFeH'], dict_Chadid_2017['inputFeH'], color='orange')
    plt.scatter(dict_Chadid_2017['laydenFeH'], remapped_Chadid['vals_mapped'])
    plt.plot(dict_Chadid_2017['laydenFeH'], dict_Chadid_2017['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Chadid['coeffs']))
    labels = dict_Chadid_2017['name']
    labels_x = dict_Chadid_2017['laydenFeH']
    labels_y = dict_Chadid_2017['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_chadid.pdf')

    plt.clf()
    remapped_Fernley = remap_metal(dict_Fernley_1997['laydenFeH'], dict_Fernley_1997['inputFeH'])
    plt.scatter(dict_Fernley_1997['laydenFeH'], dict_Fernley_1997['inputFeH'], color='orange')
    plt.scatter(dict_Fernley_1997['laydenFeH'], remapped_Fernley['vals_mapped'])
    plt.plot(dict_Fernley_1997['laydenFeH'], dict_Fernley_1997['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Fernley['coeffs']))
    labels = dict_Fernley_1997['name']
    labels_x = dict_Fernley_1997['laydenFeH']
    labels_y = dict_Fernley_1997['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_fernley.pdf')

    plt.clf()
    remapped_Solano = remap_metal(dict_Solano_1997['laydenFeH'], dict_Solano_1997['inputFeH'])
    plt.scatter(dict_Solano_1997['laydenFeH'], dict_Solano_1997['inputFeH'], color='orange')
    plt.scatter(dict_Solano_1997['laydenFeH'], remapped_Solano['vals_mapped'])
    plt.plot(dict_Solano_1997['laydenFeH'], dict_Solano_1997['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Solano['coeffs']))
    labels = dict_Solano_1997['name']
    labels_x = dict_Solano_1997['laydenFeH']
    labels_y = dict_Solano_1997['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_solano.pdf')

    plt.clf()
    remapped_Wallerstein = remap_metal(dict_Wallerstein_2010['laydenFeH'], dict_Wallerstein_2010['inputFeH'])
    plt.scatter(dict_Wallerstein_2010['laydenFeH'], dict_Wallerstein_2010['inputFeH'], color='orange')
    plt.scatter(dict_Wallerstein_2010['laydenFeH'], remapped_Wallerstein['vals_mapped'])
    plt.plot(dict_Wallerstein_2010['laydenFeH'], dict_Wallerstein_2010['laydenFeH'], linestyle = ':')
    plt.xlabel('Layden 94 Fe/H')
    plt.ylabel('High-res Fe/H')
    plt.title('coeffs '+str(remapped_Wallerstein['coeffs']))
    labels = dict_Wallerstein_2010['name']
    labels_x = dict_Wallerstein_2010['laydenFeH']
    labels_y = dict_Wallerstein_2010['inputFeH']
    for point in range(0,len(labels)): 
        plt.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    plt.savefig('test_wallerstein.pdf')
    
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
    labels = dict_Fernley_1997['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax1.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    ax1.set_title('Fernley 97')
    #
    quantx = np.array(dict_Lambert_96['laydenFeH'],dtype=float)
    quanty = np.array(dict_Lambert_96['residuals'], dtype=float)
    ax2.scatter(quantx, quanty)
    ax2.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax2_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax2.plot(quantx, np.add(np.multiply(quantx,ax2_coeff[0]),ax2_coeff[1]), linestyle='--') # regression line
    labels = dict_Lambert_96['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax2.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    ax2.set_title('Lambert 96')
    #
    quantx = np.array(dict_Nemec_2013['laydenFeH'],dtype=float)
    quanty = np.array(dict_Nemec_2013['residuals'], dtype=float)
    ax3.scatter(quantx, quanty)
    ax3.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax3_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax3.plot(quantx, np.add(np.multiply(quantx,ax3_coeff[0]),ax3_coeff[1]), linestyle='--') # regression line
    labels = dict_Nemec_2013['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax3.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    ax3.set_title('Nemec 13')
    #
    quantx = np.array(dict_Liu_2013['laydenFeH'],dtype=float)
    quanty = np.array(dict_Liu_2013['residuals'], dtype=float)
    ax4.scatter(quantx, quanty)
    ax4.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax4_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax4.plot(quantx, np.add(np.multiply(quantx,ax4_coeff[0]),ax4_coeff[1]), linestyle='--') # regression line
    labels = dict_Liu_2013['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax4.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    ax4.set_title('Liu 13')
    #
    quantx = np.array(dict_Chadid_2017['laydenFeH'],dtype=float)
    quanty = np.array(dict_Chadid_2017['residuals'], dtype=float)
    ax5.scatter(quantx, quanty)
    ax5.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax5_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax5.plot(quantx, np.add(np.multiply(quantx,ax5_coeff[0]),ax5_coeff[1]), linestyle='--') # regression line
    labels = dict_Chadid_2017['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax5.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    ax5.set_title('Chadid 17')
    #
    quantx = np.array(dict_Solano_1997['laydenFeH'],dtype=float)
    quanty = np.array(dict_Solano_1997['residuals'], dtype=float)
    ax6.scatter(quantx, quanty)
    ax6.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax6_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax6.plot(quantx, np.add(np.multiply(quantx,ax6_coeff[0]),ax6_coeff[1]), linestyle='--') # regression line
    labels = dict_Solano_1997['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax6.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    ax6.set_title('Solano 97')
    #
    quantx = np.array(dict_Wallerstein_2010['laydenFeH'],dtype=float)
    quanty = np.array(dict_Wallerstein_2010['residuals'], dtype=float)
    ax7.scatter(quantx, quanty)
    ax7.plot(quantx, np.zeros(np.shape(quantx)), color='k', linestyle=':') # zero line
    ax7_coeff = np.polyfit(quantx.ravel(), quanty.ravel(), 1)
    ax7.plot(quantx, np.add(np.multiply(quantx,ax7_coeff[0]),ax7_coeff[1]), linestyle='--') # regression line
    ax7.set_title('Wallerstein 10')
    labels = dict_Wallerstein_2010['name']
    labels_x = quantx
    labels_y = quanty
    for point in range(0,len(labels)): 
        ax7.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
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
    coeff2 = np.polyfit(np.concatenate(to_plot_x_all).ravel(), np.concatenate(to_plot_y_all).ravel(), 1) # regression coefficients to make high-res Fe/H values
    ax1.plot(np.arange(-3,1,step=0.1),np.arange(-3,1,step=0.1),linestyle='--') # one-to-one
    ax1.set_xlim([-2.9,0.4])
    ax1.set_ylim([-2.9,0.4])
    ax1.set_ylabel('FeH_highres_shifted')

    for p in range(0,len(dict_merged['laydenFeH'])): # add labels
        for label, x, y in zip(dict_merged["name"][p], dict_merged["laydenFeH"][p], np.add(dict_merged['residuals_shifted'][p],dict_merged['laydenFeH'][p])):
            ax1.annotate(label, xy = (x, y))
    '''
    labels = dict_Nemec_2013['name']
    labels_x = dict_Nemec_2013['laydenFeH']
    labels_y = dict_Nemec_2013['inputFeH']
    for point in range(0,len(labels)): 
        ax1.annotate(
            labels[point],
            xy=(labels_x[point], labels_y[point]), xytext=(labels_x[point]+0.1, labels_y[point]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    '''
    ax1.set_title('m = '+str(coeff2[0])+', b = '+str(coeff2[1]))
    # subplot 2
    [ax2.scatter(dict_merged['laydenFeH'][p], dict_merged['residuals_shifted'][p]) for p in range(0,len(dict_merged['laydenFeH']))]

    '''
    ax2.annotate(
            np.ravel(dict_merged['name'][0].values()),
            xy=(np.ravel(dict_merged['laydenFeH'][0].values()), np.ravel(dict_merged['residuals_shifted'][0].values())),
            xytext=(np.ravel(dict_merged['laydenFeH'][0].values()), np.ravel(dict_merged['residuals_shifted'][0].values())),
            textcoords='data')
    '''

    for p in range(0,len(dict_merged['laydenFeH'])): # add labels
        for label, x, y in zip(dict_merged["name"][p], dict_merged["laydenFeH"][p], dict_merged["residuals_shifted"][p]):
            ax2.annotate(label, xy = (x, y))
    
    '''
    for p in range(0,len(dict_merged['laydenFeH'])): # add labels
        ax2.annotate(
            dict_merged['name'][p],
            xy=(dict_merged['laydenFeH'][p], dict_merged['residuals_shifted'][p]),
            xytext=(dict_merged['laydenFeH'][p]+0.1, dict_merged['residuals_shifted'][p]+0.06),
            textcoords='data', ha='right', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=1),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    '''
    reg_x, reg_y = zip(*sorted(zip(np.hstack(dict_merged['laydenFeH']), np.add(np.multiply(np.hstack(dict_merged['laydenFeH']),m_merged_resid_shifted),b_merged_resid_shifted)))) # sort values in x, so that plot line doesn't zigzag
    ax2.plot(reg_x, reg_y, linestyle='--') # regression line
    ax2.plot(reg_x, np.zeros(np.shape(reg_y)), linestyle=':') # zero line
    ax2.set_ylim([-0.6,0.6])
    ax2.set_xlabel('FeH_Layden94')
    ax2.set_ylabel('FeH_highres_residuals_shifted')

    '''
    labels = dict_Nemec_2013['name']
    labels_x = dict_Nemec_2013['laydenFeH']
    labels_y = dict_Nemec_2013['inputFeH']
    for point in range(0,len(labels)): 
    '''

    plt.suptitle('chadid_fig7')
    plt.savefig('chadid_fig7_imitation_test.pdf')
    plt.clf()

    ## ## CAUTION: TEST TO SEE IF THE CONTENT IN THE KEYS IS IN ORDER (I.E., MAKE A PLOT AND SEE IF ITS THE SAME IF DATASETS ARE OVERLAID INDIVIDUALLY)


    
    # print high-res metallicities of our program stars: are the distributions bimodal?
    '''
    match_our_lambert = lit_metal.find_match_gen(lit_metal.lambert_logeps, lit_metal.our_program_stars) # find overlaps between high-res and our program stars
    match_our_nemec = lit_metal.find_match_gen(lit_metal.nemec_feh, lit_metal.our_program_stars)
    match_our_liu = lit_metal.find_match_gen(lit_metal.liu_feh2, lit_metal.our_program_stars)
    match_our_chadid = lit_metal.find_match_gen(lit_metal.chadid_feh, lit_metal.our_program_stars)
    match_our_fernley = lit_metal.find_match_gen(lit_metal.fernley97_feh, lit_metal.our_program_stars)
    match_our_solano = lit_metal.find_match_gen(lit_metal.solano_feh, lit_metal.our_program_stars)
    match_our_wallerstein = lit_metal.find_match_gen(lit_metal.wallerstein_feh, lit_metal.our_program_stars)
    all_high_res = [match_our_lambert,
                    match_our_nemec,
                    match_our_liu,
                    match_our_chadid,
                    match_our_fernley,
                    match_our_solano,
                    match_our_wallerstein]
    all_high_res_names = ['lambert',
                    'nemec',
                    'liu',
                    'chadid',
                    'fernley',
                    'solano',
                    'wallerstein']        
    # now plot the FeH values, star by star
    plt.clf() # clear plot space

    # for a given star in our own dataset, find all FeH values for that star in all the high-res studies
    for row in range(0,len(lit_metal.our_program_stars)):
        print('----START NEW STAR---')
        this_star_name = []
        this_star_feh = []
        this_dataset_name = []
        # check each dataset for a match
        for dataset in range(0,len(all_high_res)):
            ix = np.isin(all_high_res[dataset]['name'], lit_metal.our_program_stars['name'][row]) # index of first argument which appears in second (i.e., the name of the star)
            if (len(np.where(ix==True)[0]) == 1): # if there is one name match (protects against return of an empty array)

                this_star_name = np.concatenate((this_star_name,all_high_res[dataset]['name'][ix]))
                this_star_feh = np.concatenate((this_star_feh,all_high_res[dataset]['feh'][ix]))
                this_dataset_name = np.concatenate((this_dataset_name,[all_high_res_names[dataset]]))

        print(this_star_name)
        print(this_star_feh)
        print(this_dataset_name)
            
            

            
            #print('Matches with ')
            #print(lit_metal.our_program_stars['name'][row])
            #print('in dataset ')
            #print(dataset)
            #if (lit_metal.our_program_stars['name'][row] == all_high_res[dataset]['name']).any():
            #    print(all_high_res[dataset]['name'])
            #    print(all_high_res[dataset]['feh'])
            #    #all_high_res[row]['name']
            #    #inputFeH = np.append(inputFeH,input_table['feh'][row])
            #    #nameArray = np.append(nameArray,input_table['name'][row])
            
        #print('----')
        #print(this_star_name)
        #print(this_star_feh)
        #print(this_dataset_name)
    '''
    
    ########################################
    # END BUNCH OF PLOTS
    ########################################
    
