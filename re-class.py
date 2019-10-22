# author
# @ Andy Tzanidakis
# attempt to re-classify all Ia SN

from datetime import datetime
from datetime import *
import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import os
import requests
import json
from datetime import datetime
from astropy.table import Table
from astropy.time import Time
import matplotlib
import SNID_GROWTH_Flex as snid
import glob
from os import listdir
import os

# Define user login credentials:
user = "AndyTzanidakis"#input("Username: ")
passw = "Kepler314**"#getpass.getpass("Enter your password: ")

# Where will the data be stored from this python run? format: month/day/year (i.e 07022019)
date_directory= "09172019"#input("Create a new date-directory (i.e 07202019): ")
mkdir_date_dir = os.system("mkdir data/%s"%date_directory) # make a new directory where we will store the datapath

data = ascii.read("temp_data/ZTF_II_list.ascii")

src_number = 0

for name in data:
    name = name[0]
    src_number = src_number + 1
    print ("Source number: %s"%src_number)
    print ("Now downloading spectra for %s"%name)
    spectra_path = snid.fetch_ZTF_spectrum(name, user, passw, specfilter=True, program_idx=0) # fetch the ZTF spectrum path for only the best spectra available!
    print ("Attempted to download spectra!")
    print ("Here is the path I found: %s"%spectra_path)
    print ("_______")

    if spectra_path!=None:
        print ("Now in the downloading phase...")
        spec_download = snid.download_spectra(spectra_path, name, date_directory)
        print ("Downloading was completed!")
        print ("_____________")

        try:
            snid_outcome = snid.SNID_spectra(spec_download, name, date_directory, user, passw)
        except:
            print ("Sorry we had a problem running SNID on %s"%name)

        if snid_outcome==True:
            # Maybe do this interacively... show the TOP 5 fits... for each one if yes it will store or else it will reject!
            out_files_path = glob.glob("data/%s/%s/spectra/*.ascii"%(date_directory, name)) # print all the .output files!
            print ("Output files: %s"%out_files_path)

            ############ SNID OUTPUT SUMMARY #############
            print ("############## SNID OUTPUT SUMMARY ##############")

            # I'm feeding in the total path to the output file
            check_fit = snid.SNID_fit_check(out_files_path[0], name, date_directory, user, passw) # dislay the summary statistics of this SNID fit

            print ("________")
            print (check_fit)
            print ("________")

            matc = np.unique(check_fit)
            print ("Here are our unique pairs: %s"%matc)
            print ("Here is the shape of the unique pairs: %s"%np.shape(matc))

            print ("Shape of the unique pairs: %s" %(len(matc[0])))

            sh = np.where((check_fit!="IIn") & (len(matc[0])==2) )
            print ("LEN %s"%len(sh[0]))

            if len(sh[0])==0:
                print ("Looks like all the sources are IIn... Now writing to the normal list")
                snid.show_fits(out_files_path[0], name, date_directory, user, passw)
                normal = open('normals.txt', 'a') # this will be the list of sources that are typically Ia-norm...
                normal.write("%s, "%name)
                normal.close()

            elif len(sh[0])>0:
                print ("No there's at least one match that's not a SN II.. Now adding to the re-classified list")
                passe = open('re_class_list.txt', 'a') # this will be the list of source that we need to re-classify
                passe.write("%s, " %name)
                passe.close()


        elif snid_outcome==False or snid_outcome==None:
            print ("SNID couldn't run on this spectrum... Now moving to the next source")
            pass
    else:
        print ("No spectra were found in this source! moving on now!")
        pass
