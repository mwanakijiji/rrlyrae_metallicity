import glob
import pandas as pd
import os
import numpy as np
import sys
import configparser
import matplotlib.pyplot as plt

# configuration data
config = configparser.ConfigParser() # for parsing values in .ini file
config.read("../rrlyrae_metallicity/modules2/config.ini")

# functions and classes in the pipeline
pipeline_modules_dir = str(config["data_dirs"]["DIR_HOME"] + 'modules2')
#print(pipeline_modules_dir)

#import rrlyrae_metallicity

# append the project directory to the paths, so
# that we can import pipeline functionality and test it
sys.path.append(str(config["data_dirs"]["DIR_HOME"]))

from modules2 import scrape_ew_and_errew


class TestScraper(scrape_ew_and_errew.Scraper):

    
    def setup_class(self):

        # make a directory for testing, which will be removed later
        self.test_dir = os.path.dirname(__file__) + "/"

        
    def teardown_class(self):
        
        ## ## POPULATE THIS LATER
        print('yadayada')

        
    def test_scrape_data(self):

        # check some rows in the table for consistency with the original file

        # check for any bad flags

        # check that the wavel for each of the lines makes sense
        
        pass



class TestfindHK(scrape_ew_and_errew.findHK):

    
    def setup_class(self):

        # make a directory for testing, which will be removed later
        self.test_dir = os.path.dirname(__file__) + "/"

        
    def teardown_class(self):
        ## ## POPULATE THIS LATER
        print('yadayada')


    def test_find_HK(self):

        # read in a couple fake tables and check that the resulting data
        # of all passed-on quantities is right

    
    


t = TestScraper() # instantiate
t.test_get_list()
