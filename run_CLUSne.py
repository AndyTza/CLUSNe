import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy.io import ascii
import os
import requests
import json
import getpass
from astropy.io import ascii
import numpy as np
from astropy.table import Table
import warnings
from datetime import datetime
import pandas
import time
import subprocess
import logging
import webbrowser
from datetime import *
import glob
from os import listdir
import os
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
import SNID_GROWTH_Flex as snid

# Define user login credentials:
user = "AndyTzanidakis"#input("Username: ")
passw = "Kepler314**"#getpass.getpass("Enter your password: ")

# Where will the data be stored from this python run? format: month/day/year (i.e 07022019)
date_directory= "09102019"#input("Create a new date-directory (i.e 07202019): ")
mkdir_date_dir = os.system("mkdir data/%s"%date_directory) # make a new directory where we will store the datapath

# Chose what program you want to run
print ("###### Welcome to CLUSne ######")
print ("1. Single Source SNID & Posting")
print ("2. N-Source Automated Source SNID & Posting")

program_index = input("Select program [1,2]: ")

# Single source posting and SNID program
if float(program_index)==1:
    print ("-------------")
    print ("You have selected single source SNID & Commenting!")

    try:
        ZTF_target = input("Input ZTF name: ") # define ZTF name
        # Fetch & Download all spectra for the queried source
        spectra = snid.fetch_ZTF_spectrum(ZTF_target, user, passw, specfilter=False, program_idx=0) # NOTE: incorporate other programs such as RCF...
    except:
        logging.error("Invalid ZTF ID. Please try again.")
        ZTF_target = input("Input ZTF name: ") # define ZTF name
        # Fetch & Download all spectra for the queried source
        spectra = snid.fetch_ZTF_spectrum(ZTF_target, user, passw, specfilter=False, program_idx=0) # NOTE: incorporate other programs such as RCF...

    for spec in spectra:
        spec_download = snid.download_spectra(spec, ZTF_target, date_directory) # download all spectra to data path
        print (spec_download)
        snid.SNID_spectra(spec_download, ZTF_target, date_directory, user, passw) # Run SNID

    out_files_path = glob.glob("data/%s/%s/spectra/*.ascii"%(date_directory, ZTF_target))

    display_menu=True
    while display_menu==True:
        print ("Display the output SNID fits!")

        for i in enumerate(out_files_path):
            print ("%s) %s"%(i[0], i[1].split("/")[4])) # show the available outputs

        temp_1, temp_2 = [], []
        for j in enumerate(out_files_path):
            temp_1.append(j[0]) # which number on the menu
            temp_2.append(j[1]) # which output file

        menu_choice = input("Please select the spectrum number you would like to display [spec # or q]: ")

        if menu_choice=='q':
            display_menu=False
        else:
            try:
                print ("Great you have selected : %s"%temp_2[int(mo)].split("/")[4])
            except :
                pass

            file_plot = temp_2[int(menu_choice)] #data/date/ZTF_id/spectrum/.ascii
            snid.show_fits(file_plot, ZTF_target, date_directory, user, passw)


            continue

        if display_menu==False: # quit the forloop
            break



elif float(program_index)==2:
    # Define start-end dates for scanning in CLU
    start_date = input("Start Date <yearr-month-day>: ") # 06-11-2018 - 12-01-2018 OK
    end_date = input("End Date <yearr-month-date>: ")
