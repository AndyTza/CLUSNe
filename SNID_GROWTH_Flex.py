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
import glob
from datetime import *
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

# ABOUT: SNID_GROWTH_Flex.py will be used to impliment a new feature where for N>1 spectra, we query all -- and find for which one does med_rlap have the highest score
# Author: Anastasios Tzanidakis (atzanida@caltech.edu)

def download_marshall_spectra(user, passw, start_date, end_date, date_dir_name):
    """Fetch list of all ZTF targets with spectra within a given time period.

        Input
        ------
        user (str): username of GROWTH marshall account
        passw (str): password of GROWTH marshall account
        start & end date (str): start and end date of query. NOTE: format <yearr-month-day>

        Output
        -------
        Return the list of ZTF_id within the given time-stamp

        """
    # Contact GROWTH Marshall....
    programidx = 0 # CLU
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', auth=(user, passw), data = {'programidx': str(programidx)})
    sources_clu = json.loads(r.text) # preliminary information on CLU objects

    start_time = [datetime.strptime('%s'%start_date, '%Y-%m-%d')]
    end_time = [datetime.strptime('%s'%end_date, '%Y-%m-%d')]
    program_targets = []
    for source in sources_clu:
        in_time = False
        date_created = datetime.strptime(source['creationdate'], '%Y-%m-%d')

        for i in range(len(start_time)):
            # Check if source within given time-stamp
            if date_created>= start_time[i] and date_created <= end_time[i]: # creation dates WITHIN time stamp
                in_time = True

                s = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(user, passw), data = {'sourceid': str(source['id'])})
                sourceDict = json.loads(s.text)

                ZTF_names = sourceDict['name'] # write names to local variable

                # Avoid duplicate and bogus classified soures!
                if source['classification']=='bogus' or source['classification']=='duplicate' or source['classification']=='Q' or source['classification']=='None':
                    print ("%s does not have a correct classification... Now adding to flagged!"%ZTF_names)
                    flag = open('Flagged_spectra_%s.txt'%date_dir_name, 'a')
                    flag.write("%s " %ZTF_names)
                    flag.close()
                    break

                program_targets.append(ZTF_names) # append names

    return (program_targets)

def fetch_ZTF_spectrum(target_id, user, passw, specfilter=False, program_idx=0):
    """Returns list of all available spectra(data url paths) for a given ZTF target.

    Input
    ------
    target_id: ZTF name of source (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)
    specfilter: True; Choose a single spectrum based of of the peak magnitude (PM) (+7 days after PM or -4 before PM) (bool)
    specfilter: False; Choose multiple spectra uploaded in the Marshall (bool)

    Output
    ------
    list_spec: List of download paths (from marshall) in format: spectra/data/ZTFid_date_inst_vn.ascii (list)
    """

    programidx = program_idx # CLU:0
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', auth=(user, passw), data = {'programidx': str(programidx)})
    sources_clu = json.loads(r.text) # preliminary information on CLU objects
    list_spec, sd, dat_url, dat_inst = [], [], [], []

    for i in enumerate(sources_clu):
        name = i[1]['name'] # name generated from marshall

        if name == target_id: # if you find the ZTF_id you have queried
            marsh_info = i[1] # fetch marshall information!
            s = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(user, passw), data = {'sourceid': str(marsh_info['id'])})
            sourceDict = json.loads(s.text)

            # Spectral and Photometric information --
            spec = sourceDict['uploaded_spectra'] # information on spectroscopy
            photometry_info = sourceDict['uploaded_photometry'] # information on photometry

            N_spec = len(spec) # number of available spectra! -- check for uniqueness? [NOTE!]
            print ("Number of available spectra for %s : %s"%(name, N_spec))

            if specfilter==False: # Download all available spectra!
                for spectra_path in spec:
                    list_spec.append(str(spectra_path['datapath'])) # append datapaths...
                return (list_spec)

            elif specfilter==True:
                if N_spec==0:
                    # empty spectrum!
                    print ("Found empty spectrum!")
                    return (None)

            # Fetch the spectrum that will be the easisest to classify by +/- 7 days of max light
                if N_spec==1: # only one spectrum available
                    for spectra_path in spec:
                        print ("Only one spectrum found...")
                        datapath = str(spectra_path['datapath'])
                        ext = datapath.split("/")[2].split(".")[1]
                        if ext=="fits":
                            print ("It was .fits and now making it into .ascii")
                            datapath = datapath.replace(".fits", ".ascii")
                            return (datapath)
                        elif ext=="ascii":
                            return (datapath)

                if N_spec>1: # More than one spectrum available
                    # Choose best spectrum based via photometric maximum
                    print ("Found multiple spectra for %s. Now looking for the best likely spectrum available..."%target_id)
                    for path in enumerate(spec):
                            acc = path[1] # where the actual path is
                            indx_path = path[0] # index of the enumerated info

                            # Assign last modified and magnitude values to local variables!
                            last_mod = [photometry_info[i]['obsdate'] for i in range (0, len(photometry_info))] # append last motified dates of each target
                            magnitudes = [photometry_info[i]['magpsf'] for i in range (0, len(photometry_info))] # append associated magnitudes
                            mag_filter = [photometry_info[i]['filter'] for i in range (0, len(photometry_info))] # append assocaited magnitude filters

                            mag_filter = np.asanyarray(mag_filter) # convert magnitude filters to np.array
                            mag_filter = mag_filter.astype(str) # str
                            red_filter = np.where(mag_filter=='r') # Select only magnitudes with thte "red:r" filter!

                            magnitudes = np.asanyarray(magnitudes) # convert magnitudes to np.array
                            magnitudes = magnitudes.astype(float) # convert magnitudes to float
                            magnitudes = magnitudes[red_filter] # select only the red magnitudes

                            last_mod = np.asanyarray(last_mod) # convert last_motified dates to np.array
                            last_mod = last_mod.astype(str) # convert last_modidied dates to strings
                            last_mod = last_mod[red_filter] # sort last mod dates by magnitude sorting

                            minima = np.where(magnitudes==min(magnitudes)) # for not it will ignore outliers

                            magnitudes = magnitudes[minima] # brightest magnitude
                            last_mod = last_mod[minima] # associated date with brightest magnitude

                            # Find peak magnitude
                            con_date = last_mod[0].split('-') # splited list of dates by peak magnitude!
                            con_date = "".join(con_date) # collapse date: YrMonDay
                            photometric_max_date = datetime.strptime("%s"%con_date, '%Y%m%d') # date of peak photometry
                            #print (photometric_max_date)
                            # split the spectrum data path to find the date assocaited with the spectra
                            comp_spectra = acc['datapath'] # format: spectra/ZTF(id)_date_inst.ascii
                            comp_spectra = comp_spectra.split('/') # format ['spectra','data', 'ZTFid_date...']
                            comp_spectra = comp_spectra[2] # format: ZTFid_date_inst.ascii
                            comp_spectra = (comp_spectra.split('_'))[1] # format: date!
                            spectrum_date = datetime.strptime("%s"%comp_spectra, '%Y%m%d') # date of spectrum
                            sd.append(spectrum_date) # sd is the list of all spectra dates
                            dat_url.append(str(acc['datapath'])) # data url, is the url of data paths!
                            p = str(acc['datapath'])
                            p = p.split("/")[2] # ZTFid_date_inst_v.ascii
                            p = p.split("_")[2]
                            dat_inst.append(p)

        elif name != target_id:
            pass

    if N_spec>1:
        # After we're done with appending, we can figure out now which one is closes to the peak magnitude!
        sd = np.asanyarray(sd) # sd(spectra dates)
        dat_url = np.asanyarray(dat_url) # url to data
        dat_inst = np.asanyarray(dat_inst) # data instrument
        dat_inst = dat_inst.astype(str)

        # we should alter filter the instrument
        optical_inst = np.asanyarray(["P60", "P200", "Keck1", "NOT", "APO", "DCT", "LT", "Lick 3-m", "Lick", "GS", "FTN"]) # good instruments to run SNID on

        optical_check  = [set([dat_inst[i]]).issubset(optical_inst) for i in range(0,len(dat_inst))] # bool of instruments that belond to

        instruments_good = dat_inst[optical_check] # instruments that belong in the optical list
        url_good = dat_url[optical_check] # url's that belong in the optical list
        spectra_good_dates = sd[optical_check] # url to spectra that belong in the optical list

        near_peak = np.where((spectra_good_dates >= photometric_max_date) & (spectra_good_dates <= photometric_max_date + timedelta(days=14))) # +14 days of peak
        before_peak = np.where(spectra_good_dates >= photometric_max_date - timedelta(days=7)) # -7 days before peak
        before_peak_ext = np.where(spectra_good_dates >= photometric_max_date - timedelta(days=100)) # -50 days before peak
        near_peak_ext = np.where(spectra_good_dates <= photometric_max_date + timedelta(days=100))

        N_peak = len(near_peak) # number of sources near peak
        N_bpeak = len(before_peak) # number of sources before peak
        N_bpeak_ext = len(before_peak_ext)

        if np.shape(near_peak)[1] == 1: # found one spectrum +14 of peak magnitude
            print ("Found at least one spectrum within +14 days of peak magnitude. Now downloading...")
            dap = url_good[near_peak][0]
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(near_peak)[1] > 1:
            # If multiple spectra exist, take the closest one!
            mult_spec_14 = spectra_good_dates[near_peak] # dates to multiple spectra at +14 days after peak (N'xM')
            url_good_near_peak = url_good[near_peak] # url to the near peak (N'xM')
            sort_dates = np.argsort(mult_spec_14) # indx of those dates sorted!
            spec_near_peak = url_good_near_peak[sort_dates] # select those dates

            print ("Multiple spectra found within +14 days of peak magnitude. Now downloading...")
            dap = spec_near_peak[-1]
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(before_peak)[1] == 1: # look at -7 days before peak
            print ("Found at least one spectrum within -7 days before peak magnitude. Now downloading...")
            dap = url_good[before_peak][0] # datapath
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(before_peak)[1] > 1: # look at -7 days before peak for multiple ones
            print ("Multiple spectra found within -7 days before peak magnitude. Now downloading...")
            dap = url_good[before_peak][-1] # datapath choose most recent closest to the peak
            print ("DAP:%s"%dap)
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(near_peak_ext)[1] == 1: # look into the nebular phase
            print ("Found at least one spectrum within +100 days after peak magnitude. Now downloading...")
            dap = url_good[near_peak_ext][0] # datapath
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(near_peak_ext)[1] > 1:
            print ("Found multiple spectra within +100 days after peak magnitude. Now downloading...")
            dap = url_good[near_peak_ext][0] # datapath
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(before_peak_ext)[1] == 1:
            print ("Found spectra at least -100 days before peak magnitude. Now downloading...")
            dap = url_good[before_peak_ext][0] # datapath
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

        elif np.shape(before_peak_ext)[1] > 1:
            print ("Found multiple spectra at least -100 days before peak magnitude. Now downloading...")
            dap = url_good[before_peak_ext][-1] # datapath
            dap_ext = dap.split("/")[2].split(".")[1] # what is the file extension you're fetching?
            if dap_ext=="fits": # if it's fits make it into ascii
                print ("It was .fits and now making it into .ascii")
                dap = dap.replace(".fits", ".ascii")
                return (dap)
            elif dap_ext=="ascii":
                return (dap)

def fetch_single_ZTF_spectrum(target_id, user, passw, specfilter=True):
    """Returns single available spectra(data url paths) for a given ZTF target. Assumes all files are .ascii

    Input
    ------
    target_id: ZTF name of source (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)
    specfilter: True; Choose a single spectrum based of of the peak magnitude (PM) (+7 days after PM or -4 before PM) (bool)
    specfilter: False; Choose multiple spectra uploaded in the Marshall (bool)

    Output
    ------
    list_spec: List of download paths (from marshall) in format: spectra/data/ZTFid_date_inst_vn.ascii (list)
    """

    programidx = 0 # CLU
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', auth=(user, passw), data = {'programidx': str(programidx)})
    sources_clu = json.loads(r.text) # preliminary information on CLU objects
    list_spec, sd, dat_url, dat_inst = [], [], [], []

    for i in enumerate(sources_clu):
        name = i[1]['name'] # name generated from marshall

        if name == target_id:
            marsh_info = i[1] # fetch marshall information!
            s = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(user, passw), data = {'sourceid': str(marsh_info['id'])})
            sourceDict = json.loads(s.text)

            # Spectral and Photometric information --
            spec = sourceDict['uploaded_spectra'] # information on spectroscopy
            photometry_info = sourceDict['uploaded_photometry'] # information on photometry
            annotation = sourceDict['autoannotations']
            for ann in annotation:
                if ann['username']=="AndyTzanidakis":
                    print ("Found Autoannotation from @AndyTzanidakis... No need to coment on this source again")
                    return (None)

            completed = open("completed_list_07262019a.txt", 'r')
            names_ZTF = completed.read().split(" ") # names of the ZTF_ids that have been completed
            names_ZTF = np.asanyarray(names_ZTF)
            in_names = np.isin(str(target_id), names_ZTF)

            if in_names==True:
                print ("Source has already been commented on the Marshal.. Moving to next source...")
                return (None)
            elif in_names==False:
                pass

            N_spec = len(spec) # number of available spectra! -- check for uniqueness? [NOTE!]

            if specfilter==True:
                if N_spec==0:
                    # empty spectrum!
                    print ("Found empty spectrum!")
                    return (None)
            # Fetch the spectrum that will be the easisest to classify by +/- 7 days of max light
                if N_spec==1: # only one spectrum available
                    for spectra_path in spec:
                        print ("Only one spectrum found...")
                        datapath = str(spectra_path['datapath'])
                        ext = datapath.split("/")[2].split(".")[1]
                        if ext=="fits":
                            print ("It was .fits and now making it into .ascii")
                            datapath = datapath.replace(".fits", ".ascii")
                            return (datapath)
                        elif ext=="ascii":
                            return (datapath)
                else:
                    print ("Found Multiple spectra!")
                    return (None)

def download_spectra(data_path, target_id, date_dir_name, dtype='no_header'):
    """Downloads .ascii files from GROWTH Marshall for given ZTF object.

    Input
    -----
    data_path: Download URL from GROWTH Marshall -- usually called from Fetch_ZTF_Spectrum (str)
    target_id: ZTF name of source (str)
    date_dir_name: daving date directory you will be storing the downloaded files (str)
    dtype: data type you want for the ascii files. Usually no_header will not read the first row of the ascii file.

    Output
    ------
    Creates /spectra and /summary directories and downloads ascii file data path given.

    """
    # GROWTH Marshall url where data is available
    download_url = "http://skipper.caltech.edu:8080/growth-data/"

    # Complete path to url downloading data
    data_path_mod = data_path.split("/")[2].split('.')[0] #ZTFid_date_inst_vn

    try:
        # Stich final directories to the directory containing the data
        final_path_to_data = download_url + "spectra/data/" + data_path_mod + ".ascii"
        print ("this is the final path --> ")
        print (final_path_to_data)
    except:
        print ("This file you have requested is not available!")
        return (None)

    # Make a new directory with the ZTF name
    bash1 = subprocess.run("mkdir data/%s/%s"%(date_dir_name, target_id), shell=True)
    bash2 = subprocess.run("mkdir data/%s/%s/spectra"%(date_dir_name, target_id), shell=True) # create the spectra dir in the ZTF target dir
    bash3 = subprocess.run("mkdir data/%s/%s/summary"%(date_dir_name, target_id), shell=True) # create the summary dir where we save summary SNID plots

    # Fetch data from custom url
    data = ascii.read("%s"%final_path_to_data, format=dtype)
    # clean from nan values
    rmv_nan = np.where(np.isnan([data['col1'], data['col2']]))
    data.remove_rows(rmv_nan[1]) # will remove NaN identified values!
    min_lambda = min(data['col1']) # minimum wavelength
    max_lambda = max(data['col1'])# maximum wavelength

    if min_lambda>2500 and max_lambda<12000: # check that the downloaded ascii is within the optical range!
        print ("Looking in wavelength ranges: %s - %s"%(min_lambda, max_lambda))
        # Save to data path: ../n1/spectra
        # Generally in this directory we will have: ".asii", ".output"
        download_spectrum = ascii.write(data, "data/%s/%s/spectra/%s.ascii"%(date_dir_name, target_id, data_path_mod), format=dtype)
        stored_spectrum = "data/%s/%s/spectra/%s.ascii"%(date_dir_name, target_id, data_path_mod)
        return (stored_spectrum)
    else:
        return (False)

def fetch_ztf_z(target_id, user, passw):
    """ Fetch the associated ZTF redshift(z) value assigned from the GROWTH Marshall for a given ZTF target.

    Input
    -----
    target_id: ZTF name of source (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)

    Output
    ------
    z: Redhsift value assocaited with ZTF target (float)
    id_source: source id associated with ZTF target (float)
    """
    programidx = 0 # CLU
    r = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/list_program_sources.cgi', auth=(user, passw), data = {'programidx': str(programidx)})
    sources_clu = json.loads(r.text) # preliminary information on CLU objects

    for i in enumerate(sources_clu):
        name = i[1]['name'] # name generated from marshall
        if name == target_id:
            marsh_info = i[1] # fetch marshall information!
            id_source = marsh_info['id']
            z = marsh_info['redshift']

            if z==None or float(z)==None:
                s = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/source_summary.cgi', auth=(user, passw), data = {'sourceid': str(id_source)})
                sourceDict = json.loads(s.text)
                try:
                    cluz = sourceDict['CLU_z'] # CLU_z AUTO_COMMENT
                    cluz = cluz.split(" ") # split it and make it into a list
                except:
                    print ("assigned cluz to None")
                    cluz=None

                if cluz==None: # is there ever a case where there's no CLU_Z?
                    print ("Couldn't find reshift associated with this source on Marshal. Now setting z=0.04")
                    return (0.04, id_source)
                else:
                    print ("Found redshift CLU_z on the Marshal at z=%s"%cluz[0])
                    return (float(cluz[0]), id_source) # take the first redshift, since that should be the main SNe Galaxy from which we take the z value
            else:
                print ("Found reshift associated with it!")
                return (float(z), id_source)
        elif name != target_id:
            pass

def SNID_spectra(path_to_data, target_id, date_dir_name, user, passw, snid_args="rlapmin=5 fluxout=5 inter=0 plot=0"):
    """Will run SNID on .ascii files saved at ../date_dir_name/ZTF_id/spectra

    Input
    -----
    spec_list: list containing the download paths of the ZTF spectra (used for generating extentions) i.e /downloads/ZTFid_date_inst_v.ascii
    target_id: ZTF name of source (str)
    date_dir_name: daving date directory you will be storing the downloaded files (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)
    snid_args: Add custom arguments for SNID (i.e see snid-help). Default will save the TOP 5 templates from SNID.

    Output
    ------
    Will run SNID on .ascii files and extract those files (i.e .dat files) to the saving directory
    """
    SNID_prog = True # Snid progress, if true next steps will continue

    # File extension name (i.e ZTFid_date_inst_v.ascii)
    file_ext = path_to_data.split("/")[4]
    print ("file ext: %s"%file_ext)
    instrument = file_ext.split("_")[2] # instrument name

    print ("This is your instrument name: %s"%instrument)
    if str(instrument)=="Keck2":
        # We cannot run snid on IR spectra
        SNID_prog = False
        print ("Sorry I cannot process Keck2 spectra!")
        return (SNID_prog)

    # Final path to directory (i.e ../date_dir_name/ZTF_id/spectra/)
    path_to_save_dir = "data/%s/%s/spectra/"%(date_dir_name, target_id)
    print ("save path: %s"%path_to_save_dir)

    z = fetch_ztf_z(target_id, user, passw) # fetch Marshall redshift(z) and source-id

    print ('Now running SNID:  %s'%file_ext)
    Bac1 = subprocess.check_output("snid %s forcez=%s %s"%(snid_args, z[0], path_to_data), shell=True)
    Bac1_output = str(Bac1) # convert terminal output to string
    Bac1_output = np.asarray(Bac1_output.split('\\n')) # split

    snid_error_1 = np.where(Bac1_output==" Interactive mode off. No output written.")
    snid_error_2 = np.where(Bac1_output==" No peaks are good, setting z = 0.000.")

    if len(snid_error_1[0])==1 or len(snid_error_2[0])==1:
        SNID_prog = False
        print ("Problematic spectra identified, SNID had trouble running. Now adding to flagged spectra...")
        flag = open('failed_SNID_ia.txt', 'a')
        flag.write("%s,"%target_id)
        flag.close()
        return (SNID_prog)
    else:
        SNID_prog = True
        print ("SNID ran with no complications... now appending .output files to directories.")
        bac2 = subprocess.run("cp *.dat %s"%path_to_save_dir, shell=True) # copy-paste all .dat files move it into final directory
        bac3 = subprocess.run("cp *.output %s"%path_to_save_dir, shell=True) # copy-paste all .output files move into final directory
        bac4 = subprocess.run("rm *.dat", shell=True) # remove them to allow for the next iteration
        bac5 = subprocess.run("rm *.output", shell=True) # remove them locally to allow for the next itereation
        return (SNID_prog)

def load_snid_output(filename):
    """ Load a snid.output file and returns a pandas DataFrame of the best matches """
    f = open(filename).read().split("### rlap-ordered template listings ###")[-1].splitlines()
    return pandas.DataFrame([l.split() for l in f[2:]], columns=f[1][1:].split())

def Annotate_ZTF(source_id, comment_name, comment, user, passw, type="autoannotation"):
    """ Test...."""
    if type=="autoannotation":
        # set-up dictionary where we will post
        payload1 = {'action':'commit','id':-1,'sourceid':str(source_id),'datatype':'STRING','type':comment_name,'comment': "%s"%comment}
        request1 = requests.post('http://skipper.caltech.edu:8080/cgi-bin/growth/add_autoannotation.cgi', auth=(user, passw), data=payload1) # excecute Autocomment to Marshall

def upload_SNID_figure(spec_list, target_id, date_dir_name, user, passw, source_id, comment="SNID best fit"):
    """ Take the SNID fit, decide if it was a good fit, and upload to the marshall.

    Input
    ------
    target_id: ZTF name of source (str)
    date_dir_name: daving date directory you will be storing the downloaded files (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)
    source_id: Source id associated with target_id from the GROWTH Marshall (int)
    comment: Associated comment that will be alongside the SNID url (str)

    Output
    ------
    Attach SNID plots to the GROWTH Marshall.
    """

    data_extension = spec_list.split("/")[2] # select only the file extension
    data_extension = data_extension.split("_")
    extensions = data_extension[0:3] # Name, date, instrument
    v_s = data_extension[3].split(".")[0]
    extensions.append(v_s)

    path_to_output = "../%s/%s/summary/%s_%s_%s_%s_summary.jpg"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3])
    print (path_to_output) # note that the jpg file must be low DPI

    files = {'attachment' : open(path_to_output, 'rb')}

    payload3 = {'commit':'yes', "tablename":"sources", "tableid":str(source_id),'sourceid' : str(source_id),'comment' : "SNID best fit",
                'type':'comment', 'submit': 'Save Comment'}

    request1 = requests.post(url='http://skipper.caltech.edu:8080/cgi-bin/growth/edit_comment.cgi',
                             auth=(user, passw), files=files, data=payload3)
    print ("Done uploading")

def show_fits(spec_list, target_id, date_dir_name, user, passw, comp_n=6, plot=True):
    """ Show the TOP 3 SNID fits using Matplotlib.

    Input
    ------
    target_id: ZTF name of source (str)
    date_dir_name: daving date directory you will be storing the downloaded files (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)
    theta: SNID parameters (i.e sne_type, sne_name, z, z_err, age, rlap) (list)

    Output
    ------
    Will display the fits with Matplotlib
    """

    # Extract file extension
    data_ext = spec_list.split("/")[4] # ZTFid_date_inst_v...
    dat_rmv_acii = data_ext.split(".ascii")
    data_extension = data_ext.split("_")

    # Define theta:
    spec_output_data = load_snid_output("data/%s/%s/spectra/%s_snid.output"%(date_dir_name, target_id, dat_rmv_acii[0]))
    # First define local variables and store as theta!

    sne_name = spec_output_data['sn'] # load Sne template name SNID
    rlap = np.asarray(spec_output_data['rlap']) # load rlap score from SNID
    rlap = rlap.astype(float) # conver to float
    rlap = rlap[~np.isnan(rlap)]

    z = np.asarray(spec_output_data['z']) # load z score
    z = z.astype(float)
    z = z[~np.isnan(z)]

    z_err = spec_output_data['zerr'] # load z_err

    age = np.asarray(spec_output_data['age']) # load age
    age = age.astype(float)
    age = age[~np.isnan(age)]

    sne_type = spec_output_data['type'] # load sne_tpe
    sne_type = np.asanyarray(sne_type)

    theta = [sne_name, sne_type, z, z_err, rlap, age] # defined theta

    #print (theta[0], theta[1], theta[2], theta[3], theta[4])

    if data_extension[2]=="Gemini" or data_extension[2]=="Lick" :
        print ("Found Gemini!")
        extensions = data_extension[0:4] # Name, date, instrument
        v_s = data_extension[4].split(".")[0]
        extensions.append(v_s)
        comp = ["comp000%s"%i for i in range(1,comp_n)] # generate comp_list (for SNID templates)
        # Define path_to_output... where to find the output SNID file
        comp_path = ["data/%s/%s/spectra/%s_%s_%s_%s_%s_%s_snidflux.dat"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3],extensions[4], comp[i]) for i in range(0, len(comp))]
        data_path = "data/%s/%s/spectra/%s_%s_%s_%s_%s_snidflux.dat"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3], extensions[4])
    else:
        extensions = data_extension[0:3] # Name, date, instrument
        v_s = data_extension[3].split(".")[0]
        extensions.append(v_s)
        comp = ["comp000%s"%i for i in range(1,comp_n)] # generate comp_list (for SNID templates)
        print (comp)
        comp_path = ["data/%s/%s/spectra/%s_%s_%s_%s_%s_snidflux.dat"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3], comp[i]) for i in range(0, len(comp))]
        data_path = "data/%s/%s/spectra/%s_%s_%s_%s_snidflux.dat"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3])

    #print (data_path)
    data = ascii.read(data_path) # load SNID spectrum data
    data_lambda = data['wavelength[A]']
    data_flux = data['flux[arbitraty]']
    #plot = input("Display plot [y/n]: ")
    _, ax = plt.subplots(nrows=comp_n-3, ncols=1, figsize=(9,7))
    print ("Now going to display plot")

    for i in range (0, len(comp_path)-2):
        # Load tempalte data
        data_temp = ascii.read(comp_path[i])

        # Tweak headers
        if theta[0][i]=='rlap':
            theta[0][i]="None"

        if theta[1][i]=="cutoff" or theta[1][i]=="typename":
            theta[1][i]="None"

        ax[i].plot(data_lambda, data_flux, color='k', lw=2, label='data')
        ax[i].plot(data_temp['col1'], data_temp['col2'], color='red',
         label="SNe:%s (%s) z:%s+/-%s rlap:%s age:%s"%(theta[0][i], theta[1][i], theta[2][i], theta[3][i], theta[4][i], theta[5][i]))
        ax[i].set_xlabel(r'Wavelength [Å]', fontsize=7)
        ax[i].set_ylabel('Normalized Flux', fontsize=7)
        ax[i].legend(fontsize=8)
        ax[i].set_xlim(4000, 9000)
        ax[i].tick_params(axis = 'both', which = 'major', labelsize = 5, direction='in', length=5)
        plt.savefig('data/%s/%s/summary/%s_%s_%s_%s_summary.pdf'%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3]), bbox_inches='tight')
    #plt.show()
    #os.system("evince data/%s/%s/summary/%s_%s_%s_%s_summary.pdf"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3]))

def SNID_to_marshall(spec_list, target_id, date_dir_name, user, passw):
    """ Take the SNID fit, decide if it was a good fit, and upload to the marshall.

    Input
    ------
    target_id: ZTF name of source (str)
    date_dir_name: daving date directory you will be storing the downloaded files (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)

    Output
    ------
    Extrapolate data from SNID.output, view TOP3 SNID.output fits, view summary statistics, post to GROWTH Marshall
    """
    print ("Now reading SNID output!")
    # Extract file extension
    data_extension = spec_list.split("/")[2] # select only the file extension
    data_extension = data_extension.split("_")

    if data_extension[2]=="Gemini" or "Lick":
        print ("Found Gemini!")
        extensions = data_extension[0:4] # Name, date, instrument
        v_s = data_extension[4].split(".")[0]
        extensions.append(v_s)
        # Define path_to_output... where to find the output SNID file
        path_to_output = "../%s/%s/spectra/%s_%s_%s_%s_%s_snid.output"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3], extensions[4])
    else:
        extensions = data_extension[0:3] # Name, date, instrument
        v_s = data_extension[3].split(".")[0]
        extensions.append(v_s)
        path_to_output = "../%s/%s/spectra/%s_%s_%s_%s_snid.output"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3])

    # Load SNID output file
    snid_output = load_snid_output(path_to_output) # Pandas frame

    ########## Define Local Variables from the output file ##############
    sne_name = snid_output['sn'] # load Sne template name SNID
    rlap = np.asarray(snid_output['rlap']) # load rlap score from SNID
    rlap = rlap.astype(float) # conver to float
    rlap = rlap[~np.isnan(rlap)]

    z = np.asarray(snid_output['z']) # load z score
    z = z.astype(float)
    z = z[~np.isnan(z)]

    z_err = snid_output['zerr'] # load z_err

    age = np.asarray(snid_output['age']) # load age
    age = age.astype(float)
    age = age[~np.isnan(age)]

    sne_type = snid_output['type'] # load sne_tpe
    sne_type = np.asanyarray(sne_type)

    sne_type_trimmed = sne_type[0:3] # TOP 3 fits
    print (sne_type_trimmed)

    major_type = [sne_type_trimmed[i][0:2] for i in range(0,3)]
    major_type = np.asarray(major_type)
    major_type_uni = np.unique(major_type) # from the major types which are unique. i.e =1 (perfect match), =2 (good match), =3 no match!

    N_uni = np.shape(major_type_uni)[0] # how many unique pairs exist
    print ("Found N_unique pairs: %s"%N_uni)
    #############################################################################

    #if N_uni==3 or N_uni==None or N_uni==0:
        #print ("TOP 3 matches don't agree! Now added to flagged spectra...")
        #flag = open('Flagged_spectra_%s.txt'%date_dir_name, 'a')
        #flag.write("%s "%target_id)
        #flag.close()
        #return (None)

    #if major_type[0]=="AGN" or major_type[0]=="Gal" or major_type[0]=="LBV" or major_type[0]=="QSO" or major_type[0]=="typename" or major_type[0]=="M-star" or major_type[0]=="C-star" or major_type[0]=="cutoff":
        #print ("You cannot annotate source with Non-SNe type. Rejecting and adding to flagged spectra for now.")
        #flag = open('Flagged_spectra_%s.txt'%date_dir_name, 'a')
        #flag.write("%s "%target_id)
        #flag.close()
        #return (None)

    #ej_list = np.asarray(["AGN", "Gal", "LBV", "M-star", "QSO", "C-star", "typename", "cutoff"]) #what is the list of fits we want to avoid?
    #rej_list_data = set(sne_type_trimmed).issubset(rej_list)

    #if rej_list_data==True:
        #print ("Flagged as non-SNe spectrum...")
        #flag = open('Flagged_spectra_%s.txt'%date_dir_name, 'a')
        #flag.write("%s "%target_id)
        #flag.close()
        #return (None)

    med_10_z = np.median(z[0:11]) # look at the median of the top 10 fits

    theta = sne_type, sne_name, z, z_err, age, rlap # all paramters of interest

    source_id = fetch_ztf_z(str(extensions[0]), user, passw) # Fetch source id in case you want to use it for the Autoannotations!
    source_id = source_id[1] # these are the source id's!

    # Median rlap score correction... Generally if this score is above 10 the SNID fit will be good quality. Anything with rlap<10 is a bad fit
    #rlap_test = np.median(rlap[0:10])
    #if rlap_test< 10:
    #    print ("Median rlap score: %s is not good enough..."%rlap_test)
    #    flag = open('Flagged_spectra_%s.txt'%date_dir_name, 'a')
    #    flag.write("%s "%target_id)
    #    flag.close()
    #    return (None)
    #else:
        #print ("SNID median rlap score: %s"%rlap_test)

    # Here's how we're going to filter the fits
    #if med_10_z<=0.06 and N_uni==1 or N_uni==2: # if there's only 1 unique pair (Sne Classifications) & z is decent
        #print ("The SNID fit passed the filter (z_med_TOP10<0.06 and at least 2 identical pairs)")
        #d1 = input("Would you like to see the fits? [y/n]") # Show to user the fit with matplotlib
        #if d1=="y":
        #show_fits(spec_list, target_id, date_dir_name, user, passw, theta, plot=True)

        #d2 = input("Would you like to see Marshall Auto_Comments? [y/n]") # Display summaries
        #if d2=="y":
    snid_best_match = "ZTF_id:%s SNe:%s(%s), z:%s+/-%s age:%s, rlap:%s"%(extensions[0],sne_type[0], sne_name[0], z[0], z_err[0], age[0], rlap[0])
    show_fits(spec_list, target_id, date_dir_name, user, passw, theta, plot=False)
    snid_median_z = np.median(z[0:10]) # Show median of TOP 10
    snid_median_age = np.median(age[0:10]) # Show median of TOP 10
    snid_median_rlap = np.median(rlap[0:10]) # Show median of TOP 10
    snid_spectrum_date = extensions[1]

    print ("SNID_best_match: %s"%snid_best_match)
    print ("SNID_median_z: %s"%snid_median_z)
    print ("SNID_median_age: %s"%snid_median_age)
    print ("SNID_median_rlap: %s"%snid_median_rlap)
    print ("SNID_spectrum_date: %s"%snid_spectrum_date)
    #url_show = input("Show url[y/n]: ")
    #if url_show=="y":


    #post = input('Proceed to post [y/n]: ')
    #if post=="y":
        # source_id, comment_name, comment, user, passw, type="autoannotation"
    Annotate_ZTF(source_id, "SNID Best Match", snid_best_match, user, passw, type="autoannotation")# where are we going to find the auto_id?
    Annotate_ZTF(source_id, "SNID Median z", format(snid_median_z, '.4g'), user, passw, type="autoannotation")
    Annotate_ZTF(source_id, "SNID Median Age", format(snid_median_age, '.4g'), user, passw, type="autoannotation")
    Annotate_ZTF(source_id, "SNID Median rlap", snid_median_rlap, user, passw, type="autoannotation")
    Annotate_ZTF(source_id, "SNID Spectrum Date", snid_spectrum_date, user, passw, type="autoannotation")

        #post_img = input("Proceed to post fits [y/n]: ")
        #if post_img=="y":
    upload_SNID_figure(spec_list, target_id, date_dir_name, user, passw, source_id) # upload image to the marshall
    # Once compleated... make text file of all sources that have been added!
    comp_list = open('completed_list_07262019a.txt', 'a') # write to 07/26/2019
    comp_list.write("%s "%target_id)
    comp_list.close()

    #else:
        #print ("SNID had a bad fit... Now adding to flagged spectra...")
        #flag = open('Flagged_spectra_%s.txt'%date_dir_name, 'a')
        #flag.write("%s "%target_id)
        #flag.close()
        #return (None)

def SNID_fit_check(spec_list, target_id, date_dir_name, user, passw):
    """ Take the SNID fit, decide if it was a good fit, and upload to the marshall.

    Input
    ------
    target_id: ZTF name of source (str)
    date_dir_name: daving date directory you will be storing the downloaded files (str)
    user: Username for login cridentials for GROWTH Marshall (str)
    passw: Password for login cridentials for GROWTH Marshall (str)

    Output
    ------
    Extrapolate data from SNID.output, view TOP3 SNID.output fits, view summary statistics, post to GROWTH Marshall
    """
    print ("Now reading SNID output!")

    ########## Define Local Variables from the output file ##############
    data_ext = spec_list.split("/")[4] # data/date/ZTF/spectra/ZTFid_date_inst_v... : you are selecting just the name

    dat_rmv_acii = data_ext.split(".ascii") # this is just ZTF_date_instrument_v

    data_extension = data_ext.split("_") # list of all the extentions {ZTF, date, inst, v}
    print (data_extension)
    print ("This is the instrument I found: %s"%data_extension[2])

    if data_extension[2] == "Gemini" or data_extension[2] =="Lick":
        print ("Found Gemini or Lick!")
        extensions = data_extension[0:4] # Name, date, instrument
        v_s = data_extension[4].split(".")[0]
        extensions.append(v_s)
        # Define path_to_output... where to find the output SNID file
        path_to_output = "data/%s/%s/spectra/%s_%s_%s_%s_%s_snid.output"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3], extensions[4])
    else:
        print ("Normal case")
        extensions = data_extension[0:3] # Name, date, instrument
        v_s = data_extension[3].split(".")[0]
        extensions.append(v_s)
        path_to_output = "data/%s/%s/spectra/%s_%s_%s_%s_snid.output"%(date_dir_name, extensions[0], extensions[0], extensions[1], extensions[2], extensions[3])

    # Load SNID output file
    spec_output_data = load_snid_output(path_to_output) # Pandas frame

    sne_name = spec_output_data['sn'] # load Sne template name SNID
    rlap = np.asarray(spec_output_data['rlap']) # load rlap score from SNID
    rlap = rlap.astype(float) # conver to float
    rlap = rlap[~np.isnan(rlap)]

    z = np.asarray(spec_output_data['z']) # load z score
    z = z.astype(float)
    z = z[~np.isnan(z)]

    z_err = spec_output_data['zerr'] # load z_err

    age = np.asarray(spec_output_data['age']) # load age
    age = age.astype(float)
    age = age[~np.isnan(age)]

    sne_type = spec_output_data['type'] # load sne_tpe
    sne_type = np.asanyarray(sne_type)

    sne_type_trim = sne_type[0:3] # TOP 3 fits
    sne_name_trim = sne_name[0:3]
    age_trim = age[0:3]
    z_trim = z[0:3]
    rlap_trim = rlap[0:3]

    print ("#################################")
    print ("SN type summary: %s"%sne_type_trim)
    print ("SN rlap summary: %s"%rlap_trim)
    print ("SN z summary: %s"%z_trim)
    print ("SN age summary: %s"%age_trim)
    print ("#################################")

    return (sne_type_trim)
