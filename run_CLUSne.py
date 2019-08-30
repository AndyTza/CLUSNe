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
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

# Define user login credentials:
user = input("Username: ")
passw = getpass.getpass("Enter your password: ")

# Where will the data be stored from this python run? format: month/day/year (i.e 07022019)
date_directory= input("Create a new date-directory (MDYr): ")

# choose what option to run
print ("###### Welcome to CLUSne ######")
print ("1. Single Source SNID & Posting")
print ("2. N-Source Automated Source SNID & Posting")
program_index = input("Select program [1,2]: ")

if float(program_index)==1:
    print ("test")

elif float(program_index)==2:
    # Define start-end dates for scanning in CLU
    start_date = input("Start Date <yearr-month-day>: ") # 06-11-2018 - 12-01-2018 OK
    end_date = input("End Date <yearr-month-date>: ")
