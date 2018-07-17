from astropy.io import fits

def ca_corrxn(corrxn_map):
# READ IN CA ABSORPTION AND MAKE CORRECTIONS

    image, header = fits.getdata(corrxn_map,0,header=True)

    # plot merged data and fit linreg line
    m_merged,b_merged = polyfit(dict_merged['laydenFeH'], dict_merged['residuals_shifted'], 1)
    plt.scatter(dict_merged['laydenFeH'], dict_merged['residuals_shifted'])
    plt.plot(dict_merged['laydenFeH'], np.add(np.multiply(dict_merged['laydenFeH'],m_merged),b_merged))
    plt.show()


    # CALCULATE FINAL FEH VALUES FOR OUR OWN STARS, AND WRITE OUT 


    ## APPLY OFFSETS to all datasets and overlap

    plt.plot(dict_Lambert_96['laydenFeH'], dict_Lambert_96['residuals_shifted'])
    plt.plot(dict_Nemec_2013['laydenFeH'], dict_Nemec_2013['residuals_shifted'])
    plt.plot(dict_Liu_2013['laydenFeH'], dict_Liu_2013['residuals_shifted'])
    plt.plot(dict_Chadid_2017['laydenFeH'], dict_Chadid_2017['residuals_shifted'])
    plt.plot(dict_Fernley_1997['laydenFeH'], dict_Fernley_1997['residuals_shifted'])
    plt.plot(dict_Solano_1997['laydenFeH'], dict_Solano_1997['residuals_shifted'])
    plt.plot(dict_Wallerstein_2010['laydenFeH'], dict_Wallerstein_2010['residuals_shifted'])
    plt.show()
