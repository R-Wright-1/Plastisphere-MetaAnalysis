#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:58:42 2020

@author: robynwright
"""
from Bio import SeqIO
import cartopy.crs as ccrs
import csv
from itertools import chain
import math
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pandas as pd
import pickle
from scipy.cluster import hierarchy
from scipy.stats import pearsonr
import scipy.spatial.distance as ssd
import scipy.stats as stats
from skbio import DistanceMatrix
from skbio.stats.distance import anosim
from skbio.stats.distance import permanova
from skbio.stats.composition import ancom
from sklearn.cluster import AgglomerativeClustering
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

seed = 3 #random seed
fs_title, fs_main, fs_small = 14, 10, 8 #font sizes to use in figures
color_env = {'marine':'#03A498', 'freshwater':'#03B1FC', 'aquatic':'#9A75FC', 'terrestrial':'#FD9A64'} #colors for plotting each environment
label_env = ['Marine', 'Freshwater', 'Aquatic', 'Terrestrial'] #labels for each environment
ext, dpi = '.png', 600 #extension and dots per inch to save figures  with
color_source = {'aliphatic':'#5F9EA0', 'biofilm':'#FCA94A', 'other plastic':'#8B008B', 'unknown plastic':'#3593FC', 'planktonic':'#F9E79F', 'blank':'gray', 'early':'#FBC704', 'late':'#900C3F', 'collection':'gray'}
phylo_level_names = {0:'kingdoms', 1:'phyla', 2:'classes', 3:'orders', 4:'families', 5:'genera', 6:'species', 7:'ASVs'}
rename_plots = {'AmaralZettler':'Amaral Zettler $et$ $al$. 2015', 'AriasAndres':'Arias-Andres $et$ $al$. 2018', 'Canada':'Canada $et$ $al$. 2020', 'Curren':'Curren & Leong 2019', 'Delacuvellerie':'Delacuvellerie $et$ $al$. 2019', 'DeTender':'De Tender $et$ $al$. 2015', 'DeTenderB':'De Tender $et$ $al$. 2017', 'DussudHudec':'Dussud, Hudec $et$ $al$. 2018', 'DussudMeistertzheim':'Dussud, Meistertzheim $et$ $al$. 2018', 'ErniCassola':'Erni-Cassola $et$ $al$. 2019', 'Esan':'Esan $et$ $al$. 2019', 'Frere':'Frere $et$ $al$. 2018', 'Hoellein':'Hoellein $et$ $al$. 2014', 'HoelleinB':'Hoellein $et$ $al$. 2017', 'Jiang':'Jiang $et$ $al$. 2018', 'Kesy':'Kesy $et$ $al$. 2019', 'Kirstein':'Kirstein $et$ $al$. 2018', 'KirsteinB':'Kirstein $et$ $al$. 2019', 'McCormick':'McCormick $et$ $al$. 2014', 'McCormickB':'McCormick $et$ $al$. 2016', 'Oberbeckmann':'Oberbeckmann $et$ $al$. 2016', 'OberbeckmannB':'Oberbeckmann $et$ $al$. 2018', 'Ogonowski':'Ogonowski $et$ $al$. 2018', 'Parrish':'Parrish $et$ $al$. 2019', 'Pinto':'Pinto $et$ $al$. 2019', 'Pollet':'Pollet $et$ $al$. 2018', 'Rosato':'Rosato $et$ $al$. 2020', 'Syranidou':'Syranidou ', 'SyranidouPE':'Syranidou $et$ $al$. 2017a', 'SyranidouPS':'Syranidou $et$ $al$. 2017b', 'Tagg':'Tagg $et$ $al$. 2019', 'Woodall':'Woodall $et$ $al$. 2018', 'Wu':'Wu $et$ $al$. 2019', 'Xu':'Xu $et$ $al$. 2019', 'Zhang':'Zhang $et$ $al$. 2019', 'WaterOrSediment':'Water or Sediment', 'LabOrField':'Laboratory or Field', 'IncubationOrCollection':'Incubation or Collection', 'MaterialType':'Material type', 'PlasticTypeSpecific':'Plastic type (specific)', 'PlasticTypeGeneral':'Plastic type (general)', 'DEPTH':'Depth', 'IncubationTime':'Incubation time (specific)', 'IncubationGeneral':'Incubation time (general)', 'PrimerPair':'Primer pair', 'DNAExtraction':'DNA extraction method', 'lab':'Laboratory', 'not_plastic':'Not plastic', 'aged_oxope':'Aged Oxo-PE', 'freeliving':'Free living', 'particleassociated':'Particle associated', 'oxope':'Oxo-PE', 'rinse_pe':'PE rinse water', 'rinse_ps':'PS rinse water', 'rinse_wood':'Wood rinse water', 'bhet':'BHET', 'hdpe':'HDPE', 'ldpe':'LDPE', 'na':'NA', 'pa':'PA', 'pe':'PE', 'pes':'PES', 'pestur':'PESTUR', 'pet':'PET', 'phbv':'PHBV', 'pla':'PLA', 'pp':'PP', 'ps':'PS', 'pvc':'PVC', 'san':'SAN', 'w_pe':'Weathered PE', '10:14':'10:14 light:dark', '12:12':'12:12 light:dark', '16:08':'16:08 light:dark', '27F_519R':'27F-519R', '319F_806R':'319F-806R', '338F_806R':'338F-806R', '341F_785R':'341F-785R', '341F_802R':'341F-802R', '341F_806R':'341F-806R', '515F_806R':'515F-806R', '515FY_926R':'515FY-926R', '518F_926R':'518F-926R', '543F_783R':'543F-783R', '967F_1064R':'967F-1064R', 'B969F_BA1406R':'B969F-BA1406R', 'redextract_sigma':'REDExtract-$N$-AmpTM', 'gentra_puregene':'Gentra Puregene', 'purelink':'PureLink', 'powersoil':'PowerSoil', 'phenol_chloroform':'Phenol-Chloroform', 'powerbiofilm':'PowerBiofilm', 'ultraclean_soil':'UltraClean soil', 'fastdna_soil':'FastDNA soil', 'orders':'Order', 'classes':'Class', 'phyla':'Phylum', 'genera':'Genera', 'families':'Family', 'species':'Species', 'ASVs':'ASV', 'kingdoms':'Kingdom', 'PlasticOnly':'Plastic only', '534f_783r':'534F-783R', 'Phenol_chloroform':'Phenol-chloroform'}
meta_name_ind = {'Study':0, 'Latitude':1, 'Longitude':2, 'Environment':3, 'WaterOrSediment':4, 'LabOrField':5, 'IncubationOrCollection':6, 'Source':7, 'MaterialType':8, 'PlasticTypeSpecific':9, 'PlasticTypeGeneral':10, 'DEPTH':11, 'IncubationTime':12, 'IncubationGeneral':13, 'Temperature':14, 'Salinity':15, 'Light':16, 'Season':17, 'PrimerPair':18, 'DNAExtraction':19, 'PlasticOnly':20}
name_dict = {'La2020':'Latva\n$et$ $al$. 2020', 'AZ2015':'Amaral-Zettler\n$et$ $al$. 2015', 'AA2018':'Arias-Andres\n$et$ $al$. 2018', 'Ca2020':'Canada\n$et$ $al$. 2020', 'Cu2019':'Curren & Leong\n2019', 'De2019':'Delacuvellerie\n$et$ $al$. 2019', 'DT2015':'De Tender\n$et$ $al$. 2015', 'DT2017':'De Tender\n$et$ $al$. 2017', 'DH2018':'Dussud \n$et$ $al$. 2018a', 'DM2018':'Dussud\n$et$ $al$. 2018b', 'EC2019':'Erni-Cassola\n$et$ $al$. 2019', 'Es2019':'Esan\n$et$ $al$. 2019', 'Fr2018':'Frere\n$et$ $al$. 2018', 'Ho2014':'Hoellein\n$et$ $al$. 2014', 'Ho2017':'Hoellein\n$et$ $al$. 2017', 'Ji2018':'Jiang\n$et$ $al$. 2018', 'Ke2019':'Kesy\n$et$ $al$. 2019', 'Ki2018':'Kirstein\n$et$ $al$. 2018', 'Ki2019':'Kirstein\n$et$ $al$. 2019', 'MC2014':'McCormick\n$et$ $al$. 2014', 'MC2016':'McCormick\n$et$ $al$. 2016', 'Ob2016':'Oberbeckmann\n$et$ $al$. 2016', 'Ob2018':'Oberbeckmann\n$et$ $al$. 2018', 'Og2018':'Ogonowski\n$et$ $al$. 2018', 'Pi2019':'Pinto\n$et$ $al$. 2019', 'Po2018':'Pollet\n$et$ $al$. 2018', 'Ro2020':'Rosato\n$et$ $al$. 2020', 'Sy2019':'Syranidou\n$et$ $al$. 2019', 'SyPE20':'Syranidou\n$et$ $al$. 2017a', 'SyPS20':'Syranidou\n$et$ $al$. 2017b', 'Ta2019':'Tagg\n$et$ $al$. 2019', 'Wo2018':'Woodall\n$et$ $al$. 2018', 'Wr2019':'Wright\n$et$ $al$. 2020', 'Wu2019':'Wu\n$et$ $al$. 2019', 'Xu2019':'Xu\n$et$ $al$. 2019', 'Zh2019':'Zhang\n$et$ $al$. 2019', 'Br2016':'Bryant\n$et$ $al$. 2016', 'Pin201':'Pinnell\n$et$ $al$. 2019', 'Pa2019':'Parrish\n$et$ $al$. 2019'}
name_dict_2 = {'La2020':'Latva $et$ $al$. 2020', 'AZ2015':'Amaral-Zettler $et$ $al$. 2015', 'AA2018':'Arias-Andres $et$ $al$. 2018', 'Ca2020':'Canada $et$ $al$. 2020', 'Cu2019':'Curren & Leong 2019', 'De2019':'Delacuvellerie $et$ $al$. 2019', 'DT2015':'De Tender $et$ $al$. 2015', 'DT2017':'De Tender $et$ $al$. 2017', 'DH2018':'Dussud $et$ $al$.  2018a', 'DM2018':'Dussud $et$ $al$. 2018b', 'EC2019':'Erni-Cassola $et$ $al$. 2019', 'Es2019':'Esan $et$ $al$. 2019', 'Fr2018':'Frere $et$ $al$. 2018', 'Ho2014':'Hoellein $et$ $al$. 2014', 'Ho2017':'Hoellein $et$ $al$. 2017', 'Ji2018':'Jiang $et$ $al$. 2018', 'Ke2019':'Kesy $et$ $al$. 2019', 'Ki2018':'Kirstein $et$ $al$. 2018', 'Ki2019':'Kirstein $et$ $al$. 2019', 'MC2014':'McCormick $et$ $al$. 2014', 'MC2016':'McCormick $et$ $al$. 2016', 'Ob2016':'Oberbeckmann $et$ $al$. 2016', 'Ob2018':'Oberbeckmann $et$ $al$. 2018', 'Og2018':'Ogonowski $et$ $al$. 2018', 'Pi2019':'Pinto $et$ $al$. 2019', 'Po2018':'Pollet $et$ $al$. 2018', 'Ro2020':'Rosato $et$ $al$. 2020', 'Sy2019':'Syranidou $et$ $al$. 2019', 'SyPE20':'Syranidou $et$ $al$. 2017a', 'SyPS20':'Syranidou $et$ $al$. 2017b', 'Ta2019':'Tagg $et$ $al$. 2019', 'Wo2018':'Woodall $et$ $al$. 2018', 'Wr2019':'Wright $et$ $al$. 2020', 'Wu2019':'Wu $et$ $al$. 2019', 'Xu2019':'Xu $et$ $al$. 2019', 'Zh2019':'Zhang $et$ $al$. 2019', 'Br2016':'Bryant $et$ $al$. 2016', 'Pin201':'Pinnell $et$ $al$. 2019', 'Pa2019':'Parrish $et$ $al$. 2019'}
name_env = {'AZ2015':'marine', 'AA2018':'freshwater', 'Ca2020':'aquatic', 'Cu2019':'marine', 'De2019':'marine', 'DT2015':'marine', 'DT2017':'marine', 'DH2018':'marine', 'DM2018':'marine', 'EC2019':'marine', 'Es2019':'terrestrial', 'Fr2018':'marine', 'Ho2014':'freshwater', 'Ho2017':'freshwater', 'Ke2019':'aquatic', 'Ki2018':'marine', 'Ki2019':'marine', 'MC2014':'freshwater', 'MC2016':'freshwater', 'Ob2016':'marine', 'Ob2018':'aquatic', 'Og2018':'marine', 'Pi2019':'marine', 'Po2018':'marine', 'Ro2020':'marine', 'Sy2019':'marine', 'SyPE20':'marine', 'SyPS20':'marine', 'Ta2019':'aquatic', 'Wo2018':'marine', 'Wr2019':'marine', 'Wu2019':'freshwater', 'Xu2019':'marine', 'Zh2019':'terrestrial', 'Pa2019':'aquatic'}

def open_csv(fn): #open csv with file name fn, returning the data as rows
    '''
    Opens csv file and returns the data as rows
    Takes as input:
        - fn (.csv file name)
    Returns:
        - rows of the file as a list
    '''
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    return rows

def open_txt(fn):
    '''
    Opens text file and returns the data as rows
    Takes as input:
        - fn (.txt file name - assumes that data is tab delimited and lines terminated with \n)
    Returns:
        - rows of the file as a list
    '''
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f, delimiter='\t', lineterminator='\n'):
            rows.append(row)
    return rows

def open_pickle(fn):
    '''
    Opens python data object file and returns the data as a python object
    Takes as input:
        - fn (file name for python object)
    Returns:
        - python object
    '''
    with open(fn, 'rb') as pick: 
        pick = pickle.load(pick)
    return pick

def write_csv(fn, data):
    '''
    Save a .csv file of the data given
    Takes as input:
       - fn (the name with which to save the csv file)
       - data (a list of lists - the data to go in the csv file)
    '''
    with open(fn, 'w') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)
    return

def write_txt(fn, data):
    '''
    Save a .txt file (tab delimited with \n to break lines) of the data given
    Takes as input:
       - fn (the name with which to save the txt file)
       - data (a list of lists - the data to go in the txt file)
    '''
    with open(fn, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for row in data:
            writer.writerow(row)
    return

def write_pickle(fn, data): #write pickle python object with name fn and the data in data
    '''
    Write a pickle python object with the data given
    Takes as input:
        - fn (the name with which to save the data)
        - data (the python object to save)
    '''
    with open(fn, 'wb') as f:
        pickle.dump(data, f)
    return

def format_R(ft, basedir): 
    '''
    Function that transforms the feature table output by QIIME2 to a version that can be used by the R agglomerate script (phyloseq and unifrac)
    Takes as input:
        - feature table.txt with samples as columns and ASVs as rows
    Gives as output:
        - ft as a pandas dataframe with sorted sample names
        - tax_dict - a dictionary containing the full taxonomy for each ASV
        - saves various files: feature_table.csv, taxonomy.csv, taxonomy_name_only.csv, random_forest/taxonomy_name_only.csv, tax_dict.dictionary
    '''
    ft = open_txt(ft) #open the text file
    del ft[0] #delete the first line of the file that tells you it was created from a .biom file
    tax_dict = {} #create a dictionary that will have all of the taxon information
    tax_file = [['OTUID', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Species name']] #create a list that will have all of the taxon information
    tax_name_only = [['OTUID', 'Species name']]
    for a in range(len(ft)): #for each row of the feature table
        if a > 0: #if it's not the header row
            string = ft[a][-1].split(';') #split this taxon information based on semi-colons (the default format from QIIME2)
            new_string = []
            for b in range(len(string)):
                if b < 5 and 'Alteromonas' in string[b]:
                    break
                if 'uncultured' in string[b]:
                    continue
                elif 'Ambiguous' in string[b]:
                    continue
                elif 'Uncultured' in string[b]:
                    continue
                elif 'ambiguous' in string[b]:
                    continue
                else:
                    if string[b][:3] == ' D_':
                        string[b] = string[b][6:]
                    elif string[b][:2] == 'D_':
                        string[b] = string[b][5:]
                    new_string.append(string[b])
            sp_name = new_string[-1]
            while len(new_string) < 7:
                new_string.append('')
            new_string.append(sp_name)
            tax_name_only.append([ft[a][0], sp_name])
            tax_dict[ft[a][0]] = new_string #creates the taxon in the taxa dictionary
            tax_file.append([ft[a][0]]+new_string) #add this asv name and taxon information
        ft[a] = ft[a][:-1] #remove the taxon information from the feature table
    ft[0][0] = 'OTUID' #Rename the first cell to match the expected input in R
    write_csv(basedir+'feature_table.csv', ft) #write the new feature table to .csv
    ft = pd.read_csv(basedir+'feature_table.csv', header=0, index_col=0)
    ft.sort_index(inplace=True)
    ft.sort_index(axis=1, inplace=True)
    write_csv(basedir+'taxonomy.csv', tax_file) #write the taxonomy to .csv
    write_csv(basedir+'taxonomy_name_only.csv', tax_name_only) #write the taxonomy to .csv
    write_pickle(basedir+'tax_dict.dictionary', tax_dict) #write the taxon information to a python object
    return ft, tax_dict

def open_and_sort(fn):
    df = pd.read_csv(fn, header=0, index_col=0)
    df.sort_index(inplace=True)
    if 'unifrac' in fn:
        df.sort_index(axis=1, inplace=True)
    return df

def get_meta(meta): #get the information contained in the meta file as a dictionary
    '''
    Function to get the information contained in the metadata file as a dataframe and a dictionary
    Takes as input:
        - meta (.txt tab delimited file with sample names as rows and metadata categories as columns)
    Returns:
        - meta (the metadata file as a list)
        - meta_names (the names of the metadata categories)
        - meta_dict (a dictionary containing sample names as keys and all metadata stored)
    '''
    meta = open_txt(meta) #open the file
    meta_names = meta[0] #save the column names
    del meta[0] #delete the column names
    meta_dict = {} #create a dictionary of all information contained, with sample names as dictionary keys
    for a in range(len(meta)):
        meta_dict[meta[a][0]] = meta[a][1:]
        for b in range(len(meta[a])):
            if meta[a][b] not in rename_plots:
                try:
                    float(meta[a][b])
                except:
                    rename_plots[meta[a][b]] = meta[a][b].capitalize()
    return meta, meta_names, meta_dict

def get_meta_df(meta, meta_names, ft_samples):
    '''
    Function to check that only the samples that were included in the QIIME analysis (after filtering, rarefaction, etc) remain in the metadata dataframe
    Takes as input:
        - meta (the metadata file as a list (one of the outputs of the get_meta function)
        - meta_names (list of metadata categories (also output from the get_meta function))
        - ft_samples (the sample names from the feature table)
    Returns:
        - meta_df (a dataframe containing only the samples that are in the feature table)
    '''
    meta_df = pd.DataFrame(meta, columns=meta_names)
    meta_df.index = meta_df['#SampleID']
    meta_df.drop('#SampleID', axis=1, inplace=True)
    sn_meta = list(meta_df.index.values)
    for a in sn_meta:
        if a not in ft_samples:
            meta_df.drop(a, axis=0, inplace=True)
    return meta_df

def get_ASV_dict(ft, seqs):
    '''
    Function to rename the ASVs with a number (rather than the random sequence of letters/numbers that is standard for QIIME2 output)
    Takes as input:
        - ft (feature table pandas dataframe with sample names as columns and ASVs as rows)
        - seqs (name of the sequences fasta file that contains sequences for all of the ASVs)
    Returns:
        - ASV_dict (a dictionary containing the previous ASV names as keys and the new ASV names as data)
    '''
    asvs = list(ft.index.values)
    ASV_dict = {}
    for a in range(len(asvs)):
        ASV_dict[asvs[a]] = 'ASV'+str(a).zfill(4)
    seqs_rename = [] #create a list to add the sequence records to
    for record in SeqIO.parse(seqs, "fasta"): #open the sequences file and loop through it
        if record.id in asvs: #if the record id matches one found in the asv list
            record.id = ASV_dict[record.id]
            seqs_rename.append(record) #add the record to the list
    SeqIO.write(seqs_rename, "sequences_agglom_renamed.fasta", "fasta") #save the new list of sequence records
    return ASV_dict

def filter_seqs(ft, sf): 
    '''
    Function to filter the sequences file to contain only the sequences in the given feature table (used to reduce the PICRUSt analyses to only the important sequences)
    Takes as input:
        - ft (feature table pandas dataframe with sample anmes as columns and ASVs as rows)
        - sf (name of the sequences fasta file that contains sequences for all of the ASVs)
    Returns: 
        - the names of the new files saved - one of the new sequences fasta file and the other of the feature table saved in the correct format for running PICRUSt2
    '''
    asvs = list(ft.index.values) #get a list of the asv names
    seqs_agglom = [] #create a list to add the sequence records to
    for record in SeqIO.parse(sf, "fasta"): #open the sequences file and loop through it
        if record.id in asvs: #if the record id matches one found in the asv list
            seqs_agglom.append(record) #add the record to the list
    SeqIO.write(seqs_agglom, "picrust/sequences_agglom.fasta", "fasta") #save the new list of sequence records
    ft.to_csv('picrust/feature_table_agglom.txt', sep='\t')
    return 'picrust/sequences_agglom.fasta', 'picrust/feature_table_agglom.txt'

def make_KO_dict(fn):
    '''
    Function to make a dictionary with the kegg ortholog information in the file 
    Takes as input:
        - fn (the name of the file that contains the kegg ortholog information - if downloaded with the other data then this should be called kegg_list.csv)
    Returns:
        - two dictionaries, the first with only the kegg orthologs relating to xenobiotic degradation, pathogenesis or antimicrobial resistance and the second with all kegg orthologs
    '''
    KO_dict, KO_dict_full = {}, {} #set up the two dictionaries
    with open(fn, 'rU') as f: #open the .csv file
        for row in csv.reader(f): #for each row
            if row[0] == 'D': #check we are looking at kegg orthologs and not pathway names, which start with a C
                if len(row) > 3: #if the row has had something added to the end
                    if row[3] == 'Xenobiotics' or row[3] == 'Pathogen' or row[3] == 'Antimicrobial resistance': #if the KO is one of the ones we are interested in
                        KO_dict[row[1]] = row[2:] #add it to our reduced dictionary
            KO_dict_full[row[1]] = row[2:] #add it to the full dictionary
    write_pickle('KO_dict.dictionary', KO_dict) #write our shortened dictionary to file
    return KO_dict, KO_dict_full

def filter_picrust(picrust, KO_dict, KO_dict_full):
    '''
    Function to filter the PICRUSt2 output for only the kegg orthologs that we are interested in (i.e. those relating to xenobiotics biodegradation, pathogenesis and antimicrobial resistance)
    Takes as input:
        - picrust (file name of the unstratified picrust output)
        - KO_dict (dictionary with kegg orthologs as keys, containing only the orthologs of interest to us)
        - KO_dict_full (dictionary with kegg orthologs as keys, containing all orthologs)
    Returns:
        - picrust (dataframe containing only those kegg orthologs that we are interested in and that are present in at least one sample)
        - KO_dict (updated kegg ortholog dictionary, now containing any custom genes that we manually added to PICRUSt)
    '''
    picrust = pd.read_csv(picrust, sep='\t', header=0, index_col=0) #read in the picrust file
    keeping = [] #set up the list of kegg orthologs to keep
    ko = list(picrust.index.values) #get a list of all current kegg orthologs
    for k in range(len(ko)): #go through and check whether they are in the dictionary of kegg orthologs that we want to keep
        if ko[k] in KO_dict: #if they are
            keeping.append(True) #add 'True' to the list
        elif ko[k] not in KO_dict_full: #otherwise
            if ko[k][0] != 'K': #check whether we are looking at one of the genes that we added in with the HMM
                keeping.append(True) #keep it if we are
                KO_dict[ko[k]] = [ko[k], 'Xenobiotics'] #and update the dictionary
            else: #otherwise
                keeping.append(False) #add 'False' to the list
        else: #otherwise
            keeping.append(False) #add 'False' to the list
    picrust = picrust.loc[keeping, :] #only keep the rows with kegg orthologs relating to xenobiotics degradation, pathogenesis or antimicrobial resistance
    write_pickle('KO_dict.dictionary', KO_dict) #write the new dictionary to file
    picrust = picrust[picrust.max(axis=1) > 0] #only keep those lines where PICRUSt2 predicts that this gene was present
    write_pickle('picrust.dataframe', picrust) #write the new picrust dataframe to file
    return picrust, KO_dict

def study_map(dates, locs, basedir):
    '''
    Function to plot and save Figure 1A and 1B (cumulative number of studies at the beginning of 2020 and the world map showing all of these studies)
    It uses files that are not the metadata file in order to also include those studies with data that wasn't accessible
    Takes as input:
        - dates (file containing details of the numbers of studies published in each environment in each year, and whether the data for these studies are accessible, named Study_dates.csv in the data to be downloaded for this study)
        - locs (file containing details of study locations for all studies published up to 2020, including whether the data was accessible for each, named Study_location.csv in the data to be downloaded for this study)
        - basedir (location of the folder to save the figure to)
    Returns:
        - nothing, but saves figure 'study_map' to the figures folder
    '''
    dates = open_csv(dates) #open csv file with study dates
    locs = open_csv(locs) #open csv file with study locations
    plt.figure(figsize=(30,4)) #set up figure
    ax1 = plt.subplot2grid((40,9), (4,0), rowspan=32) #add axis for number of publications
    ax2 = plt.subplot2grid((1,9), (0,1), colspan=2, projection=ccrs.PlateCarree()) #add axis for map
    ax1.set_title('A', loc='left', fontsize=fs_main) #add axis label
    ax2.set_title('B', loc='left', fontsize=fs_main) #add axis label
    axins1 = plt.subplot2grid((2,9), (0,3), projection=ccrs.PlateCarree()) #add first inset map axis
    axins2 = plt.subplot2grid((3,13), (0,9), projection=ccrs.PlateCarree()) #add second inset map axis
    years, marine, freshwater, terrestrial, aquatic, marine_unav, freshwater_unav, terrestrial_unav, aquatic_unav = [], [], [], [], [], [], [], [], [] #set up lists to add data to
    mar, fw, terr, aqu, mar_unav, fw_unav, terr_unav, aqu_unav = 0, 0, 0, 0, 0, 0, 0, 0 #set all numbers to zero to start
    for a in range(1, len(dates)): #go through and see if there is data for each cell
        for b in range(len(dates[a])):
            if dates[a][b] == '': dates[a][b] = 0 #if there's nothing in the cell, make this equal to 0
            else: dates[a][b] = int(dates[a][b]) #if there is a number, make this into an integer
        years.append(dates[a][0]) #add year to list
        #add all cumulative numbers on
        mar += dates[a][1]
        mar_unav += dates[a][2]
        fw += dates[a][3]
        fw_unav += dates[a][4]
        aqu += dates[a][5]
        aqu_unav += dates[a][6]
        terr += dates[a][7]
        terr_unav += dates[a][8]
        #append the new cumulative numbers to the respective lists for each environment
        marine.append(mar), marine_unav.append(mar_unav), freshwater.append(fw), freshwater_unav.append(fw_unav), aquatic.append(aqu), aquatic_unav.append(aqu_unav), terrestrial.append(terr), terrestrial_unav.append(terr_unav)
    adding = [terrestrial, terrestrial_unav, aquatic, aquatic_unav, freshwater, freshwater_unav, marine, marine_unav]
    #now add all of the environments together
    for a in range(1, len(adding)):
        for b in range(len(adding[a])):
            adding[a][b] += adding[a-1][b]
    #set up the colors that we will use and the transparency for the studies with no data
    colors_2 = [color_env['terrestrial'], color_env['terrestrial'], color_env['aquatic'], color_env['aquatic'], color_env['freshwater'], color_env['freshwater'], color_env['marine'], color_env['marine']]
    a1, a2 = 1, 0.3
    alphas = [a1, a2, a1, a2, a1, a2, a1, a2]
    for a in range(len(adding)): #for each line that we have set up
        ax1.plot(years, adding[a], 'grey', lw=1) #plot the line in grey
        if a > 0: #and fill between that line and either the previous line or zero
            ax1.fill_between(years, adding[a-1], adding[a], facecolor=colors_2[a], alpha=alphas[a])
        else:
            ax1.fill_between(years, 0, adding[a], facecolor=colors_2[a], alpha=alphas[a])
    #show only the parts of the axis that we are interested in
    ax1.set_xlim([2012, 2020])
    ax1.set_ylim([0, 55])
    #add axis labels and change the font sizes to be consistent throughout
    ax1.set_ylabel('Cumulative number of\nPlastisphere publications', fontsize=fs_small)
    ax1.set_xlabel('Year', fontsize=fs_main)
    ax1.tick_params(axis='both', which='major', labelsize=fs_small)
    ax1.tick_params(axis='both', which='minor', labelsize=fs_small)
    #make a custom legend (this is so we have boxes for colors rather than lines)
    marine = mlines.Line2D([], [], color=color_env['marine'], marker='s', markersize=8, markeredgecolor='k', label=label_env[0], linestyle=' ')
    freshwater = mlines.Line2D([], [], color=color_env['freshwater'], marker='s', markersize=8, markeredgecolor='k', label=label_env[1], linestyle=' ')
    aquatic = mlines.Line2D([], [], color=color_env['aquatic'], marker='s', markersize=8, markeredgecolor='k', label=label_env[2], linestyle=' ')
    terrestrial = mlines.Line2D([], [], color=color_env['terrestrial'], marker='s', markersize=8, markeredgecolor='k', label=label_env[3], linestyle=' ')
    plt.sca(ax1)
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial], loc='upper left', fontsize=fs_small)  
    
    plt.sca(ax2)
    #make all locations into floats (rather than the default strings)
    for a in range(1, len(locs)):
        for b in range(len(locs[a])):
            if b == 1 or b == 2 or b == 3:
                locs[a][b] = float(locs[a][b])
    img = plt.imread("world_map.jpg") #get the background map picture and set the main and inset axis to show only the regions of interest
    ax2.imshow(img, extent=[-180, 180, 90, -90], alpha=0.6)
    axins1.imshow(img, extent=[-180, 180, 90, -90], alpha=0.6)
    axins2.imshow(img, extent=[-180, 180, 90, -90], alpha=0.6)
    extent0 = [-180, 180, -90, 90]
    ax2.set_extent(extent0, crs=ccrs.PlateCarree())
    extent1 = [-22, 30, 25, 65]
    lonmin1, lonmax1, latmin1, latmax1 = extent1
    axins1.set_extent(extent1, crs=ccrs.PlateCarree())
    extent2 = [105, 130, 15, 45]
    lonmin2, lonmax2, latmin2, latmax2 = extent2
    axins2.set_extent(extent2, crs=ccrs.PlateCarree())
    
    #add the inset axis to the main axis
    plt.draw()
    p1 = ax2.get_position()
    p2 = axins1.get_position()
    axins1.set_position([p1.x1-0.12, p1.y0-0.05, p2.width, p2.height])
    p2 = axins2.get_position()
    axins2.set_position([p1.x1-(p2.width)*1.1, p1.y0+(p2.height*0.1), p2.width, p2.height])
    
    #plot the lines around the inset areas and going to the box edges
    ax2.plot([lonmin1, lonmax1], [latmin1, latmin1], 'k', lw=0.5)
    ax2.plot([lonmin1, lonmax1], [latmax1, latmax1], 'k', lw=0.5)
    ax2.plot([lonmin1, lonmin1], [latmin1, latmax1], 'k', lw=0.5)
    ax2.plot([lonmax1,lonmax1], [latmin1, latmax1], 'k', lw=0.5)
    ax2.plot([lonmin1, -87.8], [latmin1, -2], 'k', lw=0.5)
    ax2.plot([lonmax1, 44], [latmin1, -2], 'k', lw=0.5)
    ax2.plot([lonmin2, lonmax2], [latmin2, latmin2], 'k', lw=0.5)
    ax2.plot([lonmin2, lonmax2], [latmax2, latmax2], 'k', lw=0.5)
    ax2.plot([lonmin2, lonmin2], [latmin2, latmax2], 'k', lw=0.5)
    ax2.plot([lonmax2,lonmax2], [latmin2, latmax2], 'k', lw=0.5)
    ax2.plot([lonmin2, 119], [latmin2, -17], 'k', lw=0.5)
    ax2.plot([lonmax2, 175], [latmin2, -17], 'k', lw=0.5)
    
    #plot each location, checking which color and marker these need based on the values in the other columns of the file
    for i in range(1, len(locs)):
        if locs[i][5] == 'Yes': alpha, edge = 1, 'k'
        else: alpha, edge = 1, 'w'
        if locs[i][6] == 'Lab': marker='^'
        elif locs[i][6] == 'Field': marker = 'o'
        elif locs[i][6] == 'Both': marker = '*'
        ax2.scatter(locs[i][2], locs[i][1], color=color_env[locs[i][4]], edgecolor=edge, marker=marker, s=15, linewidth=0.4, alpha=alpha)
        axins1.scatter(locs[i][2], locs[i][1], color=color_env[locs[i][4]], edgecolor=edge, marker=marker, s=25, linewidth=0.4, alpha=alpha)
        axins2.scatter(locs[i][2], locs[i][1], color=color_env[locs[i][4]], edgecolor=edge, marker=marker, s=20, linewidth=0.4, alpha=alpha)
    
    #again add a custom legend
    plt.sca(ax2)
    lab = mlines.Line2D([], [], color='w', marker='^', markersize=8, markeredgecolor='k', label='Laboratory', linestyle=' ')
    field = mlines.Line2D([], [], color='w', marker='o', markersize=8, markeredgecolor='k', label='Field', linestyle=' ')
    both = mlines.Line2D([], [], color='w', marker='*', markersize=8, markeredgecolor='k', label='Both', linestyle=' ')
    gap = mlines.Line2D([], [], color='w', marker='o', markersize=2, markeredgecolor='w', alpha=0, label='', linestyle=' ')
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial, gap, lab, field, both], loc='lower left', fontsize=fs_small)  
    
    #now save the figure (using the file extension and dpi specified at the top of the file) and close it
    plt.savefig(basedir+'/figures/study_map'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def study_metrics(meta_df, basedir):
    '''
    Function to plot basic metrics on the number of samples included for each environment, plastic type and whether these were from incubations or collections, i.e. Figure 1C, 1D, 1E
    Takes as input:
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - basedir (the directory to save the figure to)
    Returns:
        - nothing, but saves figure 'sample_metrics' to the figures folder
    '''
    #set up the figure and plots
    plt.figure(figsize=(10,3))
    ax1 = plt.subplot2grid((1,3), (0,0))
    ax2 = plt.subplot2grid((1,3), (0,1))
    ax3 = plt.subplot2grid((1,3), (0,2))
    samples = list(meta_df.index.values) #get all sample names
    envs, env_caps = ['terrestrial', 'aquatic', 'freshwater', 'marine'], ['Terrestrial', 'Aquatic', 'Freshwater', 'Marine'] #set up lists of the environment names
    for a in range(len(envs)): #for each of the environments
        bottom = 0 #reset the bottom
        for b in ['lab', 'field']: #for each of lab and field
            count = 0 #reset the count
            for c in range(len(samples)): #for each of the sample names
                if meta_df.loc[samples[c], 'Environment'] == envs[a] and meta_df.loc[samples[c], 'LabOrField'] == b: #if the metadata shows the sample is in the environment and in the lab/field category we are looking at
                    count += 1 #add it to the count
            if b == 'lab': #if we are looking at lab samples
                ax1.barh(a, count, left=bottom, color=color_env[envs[a]], hatch='//', edgecolor='k') #plot a bar with the correct colors and hatching
            else: #otherwise
                ax1.barh(a, count, left=bottom, color=color_env[envs[a]], edgecolor='k') #use the right color, but not hatched
            if count > 100: #if the count is above 100 (i.e. the bar is large enough to show text within)
                ax1.text((count/2)+bottom, a, str(count), color='w', ha='center', va='center', fontsize=fs_small) #add text saying the number of samples in this category
            bottom = count+bottom #now add this count to the bottom
            if b == 'field': #if we are looking at field samples (i.e. we have now plotted both bars for this environmet)
                ax1.text(bottom+20, a, str(bottom), color='k', ha='left', va='center', fontsize=fs_small) #print the total number of samples in this environment next to the bar
        pie, inc = [], [] #set up new lists
        for c in range(len(samples)): #for each of the samples
            if meta_df.loc[samples[c], 'Environment'] == envs[a]: #if this is the environment we are looking at
                pie.append(meta_df.loc[samples[c], 'PlasticTypeGeneral']) #add the plastic type to the list
                inc.append(meta_df.loc[samples[c], 'IncubationGeneral']) #add the incubation time to the other list
        unique = list(set(pie)) #get all of the unique values for plastic type in this environment
        pie_color = [color_source[x] for x in unique] #and get the colors that match each plastic type from the dictionary at the top
        pie_count = [pie.count(x) for x in unique] #and count how many samples of each plastic type
        pie_count = [(x/sum(pie_count))*100 for x in pie_count] #and get these counts as relative abundance
        left=0 #reset the left count
        for d in range(len(pie_count)): #for each sample type
            ax2.barh(a, pie_count[d], left=left, color=pie_color[d], label=unique[d], edgecolor='k') #plot the stacked bar corresponding to sample type
            fc = 'w' #set the fontcolor as white
            if unique[d] == 'planktonic' or unique[d] == 'biofilm': #unless it's planktonic or biofilm
                fc = 'k' #in which case, the fontcolor should be black
            if pie_count[d] > 5: #if this category has more than 5% of samples (i.e. is big enought to add text)
                ax2.text((pie_count[d]/2)+left, a, str(int(pie_count[d])), color=fc, ha='center', va='center', fontsize=fs_small) #add text showing the percentage of samples with this plastic type
            left += pie_count[d] #now add this number to the left count (to make the bar stacked)
        unique = list(set(inc)) #now we are doing the same for incubation time, so get the unique values for incubation time in this environment
        if len(unique) == 3: #if this list is equal to 3
            unique = ['early', 'late', 'collection'] #change the order to be this
        pie_color = [color_source[x] for x in unique] #and get the colors that match these incubation times
        pie_count = [inc.count(x) for x in unique] #count how many samples for each incubation time
        pie_count = [(x/sum(pie_count))*100 for x in pie_count] #and convert these counts to relative abundance
        left = 0 #reset the left count
        for d in range(len(pie_count)): #for each sample type
            ax3.barh(a, pie_count[d], left=left, color=pie_color[d], label=unique[d], edgecolor='k') #plot the stacked bar corresponding to sample type
            fc = 'w' #set the fontcolor as white
            if unique[d] == 'early': #unless it's an early incubation time
                fc = 'k' #in which case the fontcolor should be black
            if pie_count[d] > 5: #if this accounts for more than 5% of relative abundance (i.e. is big enough to add text)
                ax3.text((pie_count[d]/2)+left, a, str(int(pie_count[d])), color=fc, ha='center', va='center', fontsize=fs_small) #add text showing the percentage of samples with this incubation time
            left += pie_count[d] #now add the number to the left count 
    #set x and y labels, ticks, limits, etc.  and make patches and legends for each plot
    ax1.set_xlabel('Number of samples', fontsize=fs_small)
    plt.sca(ax1)
    plt.yticks([0, 1, 2, 3], env_caps, fontsize=fs_small)
    plt.xticks(fontsize=fs_small)
    plt.xlim([0, 1400])
    lab = patches.Patch(facecolor='w', edgecolor='k', hatch='//', label='Laboratory')
    field = patches.Patch(facecolor='w', edgecolor='k', label='Field')
    ax1.legend(handles=[lab, field], fontsize=fs_main, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=2)
    sources = ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic', 'blank']
    handles = [patches.Patch(facecolor=color_source[x], edgecolor='k', label=x.capitalize()) for x in sources]
    ax2.legend(handles=handles, fontsize=fs_small, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=2)
    plt.sca(ax2)
    plt.yticks([])
    plt.xticks(fontsize=fs_small)
    plt.xlabel('Relative abundance (%)', fontsize=fs_small)
    plt.xlim([0, 100])
    times = ['early', 'late', 'collection']
    handles = [patches.Patch(facecolor=color_source[x], edgecolor='k', label=x.capitalize()) for x in times]
    ax3.legend(handles=handles, fontsize=fs_small, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=3)
    plt.sca(ax3)
    plt.yticks([])
    plt.xticks(fontsize=fs_small)
    plt.xlabel('Relative abundance (%)', fontsize=fs_small)
    plt.xlim([0, 100])
    plt.subplots_adjust(wspace=0.1)
    ax1.set_title('C', loc='left', fontsize=fs_main)
    ax2.set_title('D', loc='left', fontsize=fs_main)
    ax3.set_title('E', loc='left', fontsize=fs_main)
    plt.savefig(basedir+'/figures/'+'sample_metrics'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def plot_box(ax, l, r, b, t, line_col):
    '''
    Function to plot a box on a given axes
    Takes as input:
        - ax (the axis to plot on)
        - l, r, b, t (left, right, bottom and top coordinates)
        - line_col (the color to use to plot the line)
    '''
    plt.sca(ax)
    plt.plot([l, r], [b, b], color=line_col, lw=2)
    plt.plot([l, r], [t, t], color=line_col, lw=2)
    plt.plot([l, l], [t, b], color=line_col, lw=2)
    plt.plot([r, r], [t, b], color=line_col, lw=2)
    return

def similarity_heatmap_combined(dist_matr_fn_w, dist_matr_fn_uw, basedir):
    '''
    Function to plot the heatmap showing both weighted and unweighted unifrac distances for all studies
    Takes as input:
        - dist_matr_fn_w (file name of the weighted unifrac - can be .csv or .txt - assumes this is a file with sample names as rows and columns with values showing their similarity in the matrix)
        - dist_matr_fn_uw (file name of the unweighted unifrac - as above)
        - basedir (name of the directory to save the figures to)
    Returns:
        - nothing, but saves figure 'unifrac_heatmap_combined' to the figures folder
    '''
    dist_matr_files = [dist_matr_fn_w, dist_matr_fn_uw]
    fig = plt.figure(figsize=(12,12))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    ax = [ax1, ax2]
    for z in range(len(dist_matr_files)):
        plt.sca(ax[z])
        dist_matr = dist_matr_files[z]
        dist_matr = dist_matr.reset_index().drop('index', axis=1)
        dist_matr = [list(dist_matr.columns)]+dist_matr.values.tolist()
        studies, first = [], [] #set up lists to store the names of all studies as well as the index for when these start within the distance matrices
        for a in range(len(dist_matr[0])): #loop through the study names
            if dist_matr[0][a] == '': #if the cell is empty, continue
                continue
            if dist_matr[0][a][:6] not in studies: #if the study part of the name (excluding sample number) isn't in the list, add it
                studies.append(dist_matr[0][a][:6])
                first.append(a) #and also add the index for where this study starts
        marine, freshwater, aquatic, terrestrial = [], [], [], [] #set up lists to store which samples are from which environment
        for a in range(len(studies)): #loop through all study names
            #for each, get the environment from the dictionary at the top
            #append the study name to the correct environment if it matches one of the environment types
            if name_env[studies[a]] == 'marine':
                marine.append(a)
            elif name_env[studies[a]] == 'freshwater':
                freshwater.append(a)
            elif name_env[studies[a]] == 'aquatic':
                aquatic.append(a)
            elif name_env[studies[a]] == 'terrestrial':
                terrestrial.append(a)
        order = marine+aquatic+freshwater+terrestrial #the order for the studies to be plotted (so that environments are grouped together)
        lens = [len(marine), len(aquatic), len(freshwater), len(terrestrial)] #get the length (/number of studies) of each of the environments
        text_cols = [] #set up a list for the text colors to be added to so these match the environment type
        for a in range(len(marine)):
            text_cols.append(color_env['marine'])
        for a in range(len(aquatic)):
            text_cols.append(color_env['aquatic'])
        for a in range(len(freshwater)):
            text_cols.append(color_env['freshwater'])
        for a in range(len(terrestrial)):
            text_cols.append(color_env['terrestrial'])
        new_text_cols = [] #now rearrange them so they are in the order they will be plotted
        for a in range(len(order)):
            for b in range(len(order)):
                if a == order[b]:
                    new_text_cols.append(text_cols[b])
        each_study, each_row = [], []
        #go  through the distance matrix and calculate means for the distance between samples within and between different studies
        for a in range(1, len(dist_matr)):
            this_num, this_row = [], []
            for b in range(1, len(dist_matr[a])):
                if b in first or b == len(dist_matr[a])-1:
                    if b == len(dist_matr[a])-1:
                        this_num.append(float(dist_matr[a][b]))
                    if this_num == []:
                        continue
                    this_row.append(np.mean(this_num))
                    this_num = []
                else:
                    if float(dist_matr[a][b]) != 0:
                        this_num.append(float(dist_matr[a][b]))
            if a in first or a == len(dist_matr)-1:
                if a == len(dist_matr)-1:
                    each_row.append(this_row)
                if each_row == []:
                    continue
                this_study = []
                for c in range(len(each_row[0])):
                    this_num = []
                    for d in range(len(each_row)):
                        this_num.append(each_row[d][c])
                    this_study.append(np.mean(this_num))
                each_study.append(this_study)
                each_row = []
            actual_row = []
            for d in range(len(order)):
                actual_row.append(this_row[order[d]])
            each_row.append(actual_row)
        #now change the order of the studies so that they are grouped by environment
        each_study_new, studies_new = [], []
        for a in range(len(order)):
            each_study_new.append(each_study[order[a]])
            studies_new.append(studies[order[a]])
        #reverse these as they will be plotted from the bottom not the top
        each_study_new.reverse()
        studies_new.reverse()
        each_study, studies = each_study_new, studies_new
        
        #get the minimum and maximum distances between studies so that colors can be normalised to this range
        colormap = mpl.cm.get_cmap('inferno_r', 256)                
        ma, mi = [], []
        for a in range(len(each_study)):
            ma.append(max(each_study[a]))
            mi.append(min(each_study[a]))
        norm = mpl.colors.Normalize(vmin=min(mi), vmax=max(ma))
        m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
        #set up figure and axes
        ally = []
        #now go through and get the equivalent color for each mean distance and plot these as a bar plot
        means = []
        for a in range(len(each_study)):
            colors = []
            x, y, bottom = [], [], []
            mean = []
            for b in range(len(each_study[a])):
                if b != a:
                    mean.append(each_study[a][b])
                color = m.to_rgba(each_study[a][b])
                colors.append(color)
                x.append(b+1)
                y.append(1)
                bottom.append(a)
                plt.bar(b+1, 1, bottom=a, color=color, edgecolor='k', width=1)
            means.append(np.mean(mean))
            ally.append(a+0.5)
        cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colormap),ax=ax[z], shrink=0.4, pad=0.02, orientation='vertical')
        cb.ax.tick_params(labelsize=fs_small)
        if z == 0:
            cb.set_label('Weighted\nunifrac distance', fontsize=fs_main)
            il = 'Weighted unifrac distance'
        else:
            cb.set_label('Unweighted\nunifrac distance', fontsize=fs_main)
            il = 'Unweighted unifrac distance'
        #change x and y limits
        plt.xlim([0.5,ally[-1]+1])
        plt.ylim([0,ally[-1]+0.5])
        colors = ['#03A498', '#03B1FC', '#9A75FC', '#FD9A64']
        line_col = 'w'
    
        #plot white boxes around samples of the same environments
        l, r, b, t = 0.54, lens[0]+0.43, ally[-1]-lens[0]+0.57, ally[-1]+0.48#left, right, bottom, top
        plot_box(ax[z], l, r, b, t, line_col)
        l, r, b, t = l+lens[0], r+lens[1], b-lens[1], t-lens[0]
        plot_box(ax[z], l, r, b, t, line_col)
        l, r, b, t = l+lens[1], r+lens[2], b-lens[2], t-lens[1]
        plot_box(ax[z], l, r, b, t, line_col)
        l, r, b, t = l+lens[2], r+lens[3], b-lens[3], t-lens[2]
        plot_box(ax[z], l, r, b, t, line_col)
        for a in range(len(studies)): #rename the studies to a full study name rather than the short code they were given
            studies[a] = name_dict_2[studies[a]]
        new_means = pd.DataFrame([means], columns=studies, index=[il])
        if z == 0:
            means_df = new_means
        else:
            means_df = pd.concat([means_df, new_means])
            write_pickle(basedir+'/mean_unifrac.dataframe', means_df)
        plt.sca(ax[z])
        empty = []
        text_cols.reverse()
        for a in range(len(ally)): #plot the study names with the correct colors  on the y axis
            plt.text(0, ally[a], studies[a], color=text_cols[a], va='center', ha='right', fontsize=fs_small)
            empty.append('')
        plt.yticks(ally, empty, fontsize=fs_small) #remove the y ticks
        if z == 0:
            studies.reverse()
            for a in range(len(x)): #not plot the study names on the x axis (the list is reversed as we plotted from the bottom up)
                plt.text(x[a], ally[-1]+1.5, studies[a], rotation=90, color=text_cols[-(a+1)], va='bottom', ha='center', fontsize=fs_small)
            plt.xticks(x, empty, fontsize=fs_small) #remove the x ticks
            ax[z].xaxis.tick_top()
        else:
            plt.xticks([])
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(basedir+'/figures/unifrac_heatmap_combined'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def transform_for_NMDS(similarities, n_jobs=1): #transform the similarity matrix to 2D space (n_components=2)
    '''
    For a similarity matrix, calculate the best fit of points in 2D 
    Takes as input:
        - similarities (a distance matrix with no row or sample names)
        - n_jobs (number of processors to use for constraining - doesn't seem to work so we leave this as the default 1)
    '''
    X_true = similarities.iloc[0:].values
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity="precomputed", n_jobs=n_jobs)
    pos = mds.fit(similarities).embedding_
    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12, dissimilarity="precomputed", random_state=seed, n_jobs=n_jobs, n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)
    npos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((npos ** 2).sum())
    clf = PCA()
    npos = clf.fit_transform(npos)
    return npos, nmds.stress_

def get_single_nmds(rows, filter_on, filt_ind, color_on, color_ind, ax, leg, colors, names, meta_dict, second_filter='', second_filter_ind='', npos='', n_jobs=1, get_stats=False):
    '''
    Function to get a single nmds plot on an axis
    Takes as input:
        - rows (a distance matrix with no row or column names, i.e. containing only numbers)
        - filter_on (if we want to filter out samples from the matrix, then this contains what we are filtering on)
        - filter_ind (and the index in the metadata that we are filtering on)
        - color_on (the category that we are coloring samples on)
        - color_ind (the index for this category in the metadata)
        - ax (the axis we are plotting on)
        - leg (the position of the index we are plotting)
        - colors (the colors to use for plotting categories)
        - names (the sample names, in the same order as the samples occur in rows, the distance matrix)
        - meta_dict (a dictionary containing all metadata for each sample, with sample names as keys)
        - second_filter (if we want to filter samples based on a second category, then this contains what we are filtering on. Default is '')
        - second_filter_ind (and the index in the metadata that we are filtering on. Default is '')
        - npos (the nmds positions - this saves calculating the nmds points for a second time if we are plotting the same data again, e.g. just with different colored points)
        - n_jobs (the number of processors to use for calculating the nmds plot - this doesn't seem to work if we put more than 1, so isn't actually used)
        - get_stats (this tells us whether to get stats based on our categories - e.g. permanova and anosim - default is False)
    '''
    plt.sca(ax)
    color = []
    #If we are filtering, go through and save only the samples that are included into new lists
    if filter_on != 'none':
        new_names, new_ind = [], []
        for a in range(len(names)):
            if meta_dict[names[a]][filt_ind] == filter_on:
                if second_filter != '':
                    if meta_dict[names[a]][second_filter_ind] != second_filter:
                        continue
                new_names.append(names[a])
                new_ind.append(a)
        new_rows = []
        if new_names == []:
            return
        for b in range(len(rows)):
            if b in new_ind:
                this_row = []
                for c in range(len(rows[b])):
                    if c in new_ind:
                        this_row.append(rows[b][c])
                new_rows.append(this_row)
        names = new_names
        rows = new_rows
    #get the appropriate colors for all included samples
    groups = []
    for a in range(len(names)):
        added = False
        groups.append(meta_dict[names[a]][color_ind])
        for b in range(len(color_on)):
            if meta_dict[names[a]][color_ind] == color_on[b]:
                color.append(colors[b])
                added = True
        if not added:
            print('Didnt find ', color_on[b], 'in the metadata file so no color for this or in the list of colors')
    #if we don't already have the values for npos, then transform the similarity matrix for this
    if npos == '':
        s = pd.DataFrame(rows)
        npos, stress = transform_for_NMDS(s, n_jobs=n_jobs)
    else:
        s = pd.DataFrame(rows)
    if get_stats:
        dm = DistanceMatrix(s)
        ans = anosim(dm, np.array(groups))
        perm = permanova(dm, np.array(groups))
        string = ans[0]+': '+ans[1]+'='+str(round(ans[4], 3))+r', $p$='+str(round(ans[5], 3))
        string += '\n'+perm[0]+': '+perm[1]+'='+str(round(perm[4], 3))+r', $p$='+str(round(perm[5], 3))
        t = plt.text(0.02, 0.98, string, ha='left', va='top', transform = ax.transAxes, fontsize=fs_main, fontweight='bold')
        t.set_bbox(dict(facecolor='white', alpha=0.8))
        
    size = 20
    #plot the values for nmds 1 (npos[a,0]) and nmds 2 (npos[a,1]) for all samples, giving the color as determined by sample type
    for a in range(len(rows[0])):
        plt.scatter(npos[a,0], npos[a,1], color=color[a], marker='o', s=size, edgecolor='k')
    #get the legend handles incase these are being plotted
    
    handles = []
    for a in range(len(color_on)):
        if color_on[a] in rename_plots:
            this_color = rename_plots[color_on[a]]
        else:
            this_color = color_on[a]
        handles.append(mlines.Line2D([], [], color=colors[a], marker='o', markersize=fs_main, markeredgecolor='k', label=this_color, linestyle=' '))    
    if filter_on != 'none':
        return handles
    else:
        return npos, handles

def nmds_plot_study_env(dist_matr_fn_w, dist_matr_fn_uw, meta_dict, basedir, n_jobs=1):
    '''
    Function to make nmds plots for weighted and unweighted unifrac distances for all samples, either colored by environment or by study
    Takes as input:
        - dist_matr_fn_w (file name of the weighted unifrac - can be .csv or .txt - assumes this is a file with sample names as rows and columns with values showing their similarity in the matrix)
        - dist_matr_fn_uw (file name of the unweighted unifrac - as above)
        - meta_dict (dictionary containing metadata for all samples, with sample names as keys)
        - basedir (name of the directory to save the figures to)
        - n_jobs (number of processors to use for calculating the best ordination - doesn't seem to work so we leave with the default of 1)
    Returns:
        - nothing, but saves figure 'nmds_overall' to the figures folder
    '''
    plt.figure(figsize=(15,15))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (0,1))
    ax3 = plt.subplot2grid((2,2), (1,0))
    ax4 = plt.subplot2grid((2,2), (1,1))
    
    dist = [dist_matr_fn_w, dist_matr_fn_uw]
    ax = [[ax1, ax2], [ax3, ax4]]
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    color_env = ['#03A498', '#03B1FC', '#9A75FC', '#FD9A64']
    study, study_index = ['AmaralZettler','AriasAndres','Canada','Curren','Delacuvellerie','DeTender','DeTenderB','DussudHudec','DussudMeistertzheim','ErniCassola','Esan','Frere','Hoellein','HoelleinB','Jiang','Kesy','Kirstein','KirsteinB','McCormick','McCormickB','Oberbeckmann','OberbeckmannB','Ogonowski','Parrish', 'Pinto','Pollet','Rosato','Syranidou','SyranidouPE','SyranidouPS','Tagg','Woodall','Wu','Xu','Zhang'], 0
    colormap_40b, colormap_40c, colormap_40a = mpl.cm.get_cmap('tab20b', 256), mpl.cm.get_cmap('tab20c', 256), mpl.cm.get_cmap('tab20', 256)
    norm, norm2, norm3 = mpl.colors.Normalize(vmin=0, vmax=19), mpl.colors.Normalize(vmin=20, vmax=39), mpl.colors.Normalize(vmin=40, vmax=59)
    m, m2, m3 = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40b), mpl.cm.ScalarMappable(norm=norm2, cmap=colormap_40c), mpl.cm.ScalarMappable(norm=norm3, cmap=colormap_40a)
               
    for z in range(len(dist)):
        dist_matr = pd.read_csv(basedir+dist[z], header=0, index_col=0) #read in the distance matrix
        dist_matr = dist_matr.astype('float') #turn all values into floats
        names = list(dist_matr.columns) #get a list of column names
        dist_matr = dist_matr.rename_axis('ID').values #and turn it into a matrix with no row or column names
        color_study = []
        for a in range(len(study)):
            if a < 20:
                color_study.append(m.to_rgba(a))
            elif a < 40:
                color_study.append(m2.to_rgba(a))
            else:
                color_study.append(m3.to_rgba(a))
        filter_on, filter_index = 'none', 'none'
        second_filter, second_filter_ind = '', ''
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, envs, env_index, ax[z][0], 'upper left', color_env, names, meta_dict, second_filter, second_filter_ind, '', n_jobs=n_jobs, get_stats=True)
        if z == 0:
            plt.sca(ax[z][0])
            plt.legend(handles=handles, loc='upper right', fontsize=fs_main, edgecolor='k')
        else:
            plt.sca(ax[z][0])
            plt.legend(handles=handles, loc='lower left', fontsize=fs_main, edgecolor='k')
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, study, study_index, ax[z][1], 'upper right', color_study, names, meta_dict, second_filter, second_filter_ind, npos, n_jobs=n_jobs, get_stats=True)
        if z == 0:
            plt.sca(ax[z][1])
            plt.legend(handles=handles, bbox_to_anchor=(1.05,1.025), fontsize=fs_main, edgecolor='k')
    titles = [r'$\bf{Environment}$', r'$\bf{Study}$', '', '']
    title_letter = ['A', 'B', 'C', 'D']
    axes = [ax1, ax2, ax3, ax4]
    ylabs, xlabs = [r'$\bf{Weighted}$ $\bf{unifrac}$'+'\nnMDS2', '', r'$\bf{Unweighted}$ $\bf{unifrac}$'+'\nnMDS2', ''], ['', '', 'nMDS1', 'nMDS1']
    for a in range(len(titles)):
        axes[a].set_title(titles[a], fontsize=fs_title)
        axes[a].set_title(title_letter[a], fontsize=fs_title, loc='left')
        axes[a].set_xlabel(xlabs[a], fontsize=fs_title)
        axes[a].set_ylabel(ylabs[a], fontsize=fs_title)
    plt.savefig(basedir+'/figures/nmds_overall'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def get_venn_labels(data): 
    '''
    Calculate the number of overlapping ASVs between different samples
    It takes as input:
        - data (which should be a list of lists, containing the taxa present in different samples)
    '''
    N = len(data)
    sets_data = [set(data[i]) for i in range(N)]
    s_all = set(chain(*data))
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if  key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if  key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value
    labels = []
    for k in set_collections:
        labels.append(str(len(set_collections[k])))
    return labels

def get_venn(labels, ax, colors, names=['A', 'B', 'C', 'D', 'E'], plotting=[True, True, True, True, True]):
    '''
    This is a function to get a venn diagram with up to 5 ellipses (only plotted if that sample type is present, but the layout remains the same)
    It takes as input:
        - labels (a list of labels as output by the get_venn_labels function)
        - ax (axis to plot the venn diagram on)
        - colors (list of colors to use for each ellipse, expected that this will have length 5)
        - names (list of names to use for each ellipse, expected that this will have length 5. Default=['A', 'B', 'C', 'D', 'E'])
        - plotting (list of True/False for whether to plot each ellipse, expected that this will have length 5. Default is plotting=[True, True, True, True, True])
    '''
    #list of locations for each point of the ellipses
    el1, el2, el3, el4, el5 = [0.428, 0.469, 0.558, 0.578, 0.489], [0.449, 0.543, 0.523, 0.432, 0.383], [0.87, 0.87, 0.87, 0.87, 0.87], [0.50, 0.50, 0.50, 0.50, 0.50], [155.0, 82.0, 10.0, 118.0, 46.0]
    for a in range(len(el1)): #plot all ellipses 
        if plotting[a]: #only if this value is true for each
            e = patches.Ellipse(xy=(el1[a], el2[a]), width=el3[a], height=el4[a], angle=el5[a], color=colors[a], alpha=0.4)
            ax.add_patch(e)
    #list of text locations
    x = [0.27, 0.72, 0.55, 0.91, 0.78, 0.84, 0.76, 0.51, 0.39, 0.42, 0.50, 0.67, 0.70, 0.51, 0.64, 0.10, 0.20, 0.76, 0.65, 0.18, 0.21, 0.81, 0.74, 0.27, 0.34, 0.33, 0.51, 0.25, 0.28, 0.36, 0.51]
    y = [0.11, 0.11, 0.13, 0.58, 0.64, 0.41, 0.55, 0.90, 0.15, 0.78, 0.15, 0.76, 0.71, 0.74, 0.67, 0.61, 0.31, 0.25, 0.23, 0.50, 0.37, 0.37, 0.40, 0.70, 0.25, 0.72, 0.22, 0.58, 0.39, 0.66, 0.47]
    for a in range(len(x)):
        if labels[a] != '0': #only plot the overlap number if it is bigger than 0 (it will probably be 0 only if that sample type wasn't present and we aren't plotting the ellipse anyway)
            ax.text(x[a], y[a], labels[a], horizontalalignment='center',verticalalignment='center',fontsize=fs_small,color="black")
    #list of text locations
    x, y = [0.02, 0.72, 0.97, 0.88, 0.12], [0.72, 0.94, 0.74, 0.05, 0.05]
    for a in range(len(x)):
        if plotting[a]: #add the sample type name if the value for plotting is true
            ax.text(x[a], y[a], names[a], horizontalalignment='center',verticalalignment='center',fontsize=fs_main,color="black")
    return


def get_diversity(diversity, sample):
    '''
    function to calculate a range of different diversity metrics
    It takes as input:
        - diversity (the name of the diversity metric we want, can be 'Simpsons', 'Shannon', 'Richness', 'Evenness', 'Maximum' (Maximum is not a diversity metric, the function will just return the maximum abundance value given in sample)
        - sample (a list of abundance values that should correspond to one sample)
    Returns:
        - The diversity index for the individual sample
    '''
    for a in range(len(sample)):
        sample[a] = float(sample[a])
    total = sum(sample)
    if diversity == 'Simpsons':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)**2
        simpsons = 1-(sum(sample))
        return simpsons
    elif diversity == 'Shannon':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)
            if sample[b] != 0:
                sample[b] = -(sample[b] * (np.log(sample[b])))
        shannon = sum(sample)
        return shannon
    elif diversity == 'Richness':
        rich = 0
        for b in range(len(sample)):
            if sample[b] != 0:
                rich += 1
        return rich
    elif diversity == 'Evenness':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)
            if sample[b] != 0:
                sample[b] = -(sample[b] * (np.log(sample[b])))
        shannon = sum(sample)
        rich = 0
        for b in range(len(sample)):
            if sample[b] != 0:
                rich += 1
        even = shannon/(np.log(rich))
        return even
    elif diversity == 'Maximum':
        ma = (max(sample)/total)*100
        return ma
    return

def get_cols(num):
    '''
    This is a function to get a list of either 20 or 40 colors depending on the input 'num' (assumed that this is a number between 0 and 40)
    '''
    colormap_20 = mpl.cm.get_cmap('tab20', 256)
    norm = mpl.colors.Normalize(vmin=0, vmax=19)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_20)
    colors_20 = []
    for a in range(20):
        colors_20.append(m.to_rgba(a))
    colormap_40b = mpl.cm.get_cmap('tab20b', 256)
    norm = mpl.colors.Normalize(vmin=0, vmax=19)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40b)
    colors_40 = []
    for a in range(20):
        colors_40.append(m.to_rgba(a))
    colormap_40c = mpl.cm.get_cmap('tab20c', 256)
    norm = mpl.colors.Normalize(vmin=20, vmax=39)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40c)
    for a in range(20):
        a += 20
        colors_40.append(m.to_rgba(a))
    if num < 21:
        return colors_20
    else:
        return colors_40

def bar_dendro_venn(ft, ft_full, meta_dict, basedir, tax_dict):
    '''
    Function to make the dendrogram, stacked bar, heatmap, simpsons diversity and venn diagrams of shared ASVs
    Takes as input:
        - ft (dataframe with samples as columns and ASVs as rows)
        - meta_dict (dictionary containing metadata for all samples, with sample names as keys)
        - basedir (name of the directory to save the figures to)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
    Returns:
        - nothing, but saves figure 'dendro_venn' to the figures folder
    '''
    #set up the figure and all axes
    fig = plt.figure(figsize=(12,19))
    ax1 = plt.subplot2grid((2,7), (0,0), frameon=False)
    ax2 = plt.subplot2grid((2,7), (0,2), colspan=3)
    ax3 = plt.subplot2grid((2,7), (0,5))
    ax4 = plt.subplot2grid((2,7), (0,6))
    ax3_colbar = plt.subplot2grid((50, 7), (24,6))
    ax5 = plt.subplot2grid((27,2), (16,0), rowspan=4, frameon=False)
    ax6 = plt.subplot2grid((27,2), (16,1), rowspan=4, frameon=False)
    ax7 = plt.subplot2grid((27,2), (21,0), rowspan=4, frameon=False)
    ax8 = plt.subplot2grid((27,2), (21,1), rowspan=4, frameon=False)
    ft_full = pd.read_csv(ft_full, header=0, index_col=0)
    ft_full = ft_full/20
    ft = ft*100 #convert the feature table to % relative abundance
    colnames = list(ft.columns) #get the column names (sample names)
    env, source, env_source_dict, env_source = [], [], {}, []
    for a in range(len(colnames)): #use the sample names to get the metadata categories that these samples belong to (environment and source)
        env.append(meta_dict[colnames[a]][3])
        source.append(meta_dict[colnames[a]][10])
        if meta_dict[colnames[a]][5] == 'lab':
            meta_dict[colnames[a]][5] = 'laboratory'
        env_source_dict[colnames[a]] = meta_dict[colnames[a]][3]+' '+'\n'+meta_dict[colnames[a]][10]
        env_source.append(meta_dict[colnames[a]][3]+' '+'\n'+meta_dict[colnames[a]][10])
    envs, sources = ['marine', 'freshwater', 'aquatic', 'terrestrial'], ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic']
    colors = [color_source[sources[0]], color_source[sources[1]], color_source[sources[2]], color_source[sources[3]], color_source[sources[4]]]
    venn_ax = [ax5, ax6, ax7, ax8]
    #sort the data and make venn diagrams with ASVs that overlap between treatments
    for a in range(len(envs)):
        asv_present, all_unique_asv, plotting = [], [], []
        names = []
        for b in range(len(sources)):
            keeping = []
            for c in range(len(env)):
                if env[c] == envs[a] and source[c] == sources[b]:
                    keeping.append(True)
                else:
                    keeping.append(False)
            source_ft = ft.loc[:, keeping] #this should now only keep the samples that belong to envs[a] and sources[b]
            ma = list(source_ft.max(axis=1))
            all_asv = list(source_ft.index.values) #get all ASVs
            asvs_here = []
            #only keep the ASVs if they are present (i.e. have abundance above 0)
            for d in range(len(all_asv)):
                if ma[d] == 'nan':
                    continue
                if ma[d] > 0:
                    asvs_here.append(all_asv[d])
            asv_present.append(asvs_here)
            all_unique_asv += asvs_here
            #tell it to only plot if this treatment exists
            if len(asvs_here) == 0:
                names.append(sources[b].capitalize())
                plotting.append(False)
            else:
                names.append(sources[b].capitalize())
                plotting.append(True)
        all_unique_asv = set(all_unique_asv) #get all ASVs that are present in this environment
        labels = get_venn_labels(asv_present) #get the labels based on which ASVs are present in each treatment (source)
        get_venn(labels, venn_ax[a], colors, names=names, plotting=plotting) #get the venn diagram for this environment
        venn_ax[a].set_title(envs[a].capitalize()+'\nTotal ASVs = '+str(len(all_unique_asv)), fontweight="bold", fontsize=fs_main) #change the title for this axes
        plt.sca(venn_ax[a])
        plt.xticks([]), plt.yticks([]) #remove the x and y ticks
    #get simpsons diversity for each sample
    indiv_div = ['Simpsons']
    for cn in colnames:
        indiv_div.append(get_diversity('Simpsons', list(ft_full.loc[:, cn])))
    #rename all sample names by a combination of the environment and source
    ft_full.rename(columns=env_source_dict, inplace=True)
    env_source = list(set(env_source)) #now get a list of all unique environments and sources
    ft_full = ft_full.reset_index()
    ft_full.rename(columns={'OTUID':'index'}, inplace=True)
    ft_full = ft_full.append(pd.Series(indiv_div, index=ft_full.columns), ignore_index=True) #add the simpsons diversity to the overall dataframe
    ft_full.index = ft_full['index']
    ft_full.drop('index', axis=1, inplace=True)
    ft_dendro = ft_full.drop('Simpsons', axis=0) #make a new dataframe without the simpsons index of diversity
    ft_dendro = ft_dendro.groupby(by=ft_dendro.columns, axis=1).mean() #group this by the environment/source name, taking a mean rather than a sum for each ASV
    ft_dendro.to_csv('grouped_samples_for_unifrac.csv')
    os.rename(basedir+'qiime_output/tree.nwk', 'tree.nwk')
    os.system("/usr/local/bin/Rscript unifrac_grouped_samples.R")
    os.rename('tree.nwk', basedir+'qiime_output/tree.nwk')
    os.rename('grouped_samples_for_unifrac.csv', basedir+'grouped_samples_for_unifrac.csv')
    os.rename('weighted_unifrac.csv', basedir+'weighted_unifrac_grouped_samples.csv')
    w_uf = pd.read_csv(basedir+'weighted_unifrac_grouped_samples.csv', header=0, index_col=0)
    snames = list(w_uf.columns)
    r_name_dict = {}
    for s in snames:
        r_name_dict[s] = s.replace('..', ' \n').replace('.', ' ')
    w_uf.rename(columns=r_name_dict, index=r_name_dict, inplace=True)
    plt.sca(ax1)
    #ft_dendro_T = ft_dendro.transpose() #transpose the dataframe
    Z = ssd.squareform(w_uf)
    Z = hierarchy.linkage(Z, "ward")
    mpl.rcParams['lines.linewidth'] = 2
    hierarchy.set_link_color_palette(['k'])
    dn = hierarchy.dendrogram(Z, above_threshold_color='k', orientation='left') #plot this dendrogram of sample groupings
    y_labels, locs, xlocs, labels = list(ax1.get_yticklabels()), list(ax1.get_yticks()), list(ax1.get_xticks()), [] #get all labels so that we can work out which order samples are plotted in (and apply this to the other plots)
    for y in y_labels:
        labels.append(y.get_text())
    plot_labels, pl, colors = [], [], []
    grouping = list(w_uf.columns)
    sources.append('blank')
    #get the appropriate color and label for the order of samples
    for a in range(len(labels)):
        for b in range(len(envs)):
            if envs[b] in grouping[int(labels[a])-1]:
                colors.append(color_env[envs[b]])
        plot_labels.append(grouping[int(labels[a])-1].capitalize())
        pl.append(grouping[int(labels[a])-1])
    #add the yticks (with color corresponding to the environment that they came from)
    plt.yticks([])
    for a in range(len(locs)):
        ax1.text(xlocs[0]-0.15, locs[a], plot_labels[a], fontsize=fs_main, bbox=dict(facecolor=colors[a], alpha=0.4, pad=0.5, edgecolor='w'), ha='center', va='center')
    plt.sca(ax1)
    plt.xticks([])
    asv = list(ft_dendro.index.values)
    phyla, clss, spec = {}, {}, {}
    #for all ASVs, get the phylum, class and species that they belong to (although we don't actually use the species information)
    for a in range(len(asv)):
        tax = tax_dict[asv[a]]
        phyla[asv[a]] = tax[1]
        clss[asv[a]] = tax[3]
        spec[asv[a]] = tax[-1]
    #make new feature tables where ASVs are named by phyla/class and then group all of these phylums/classes (summing rather than taking means now)
    ft_T = ft_dendro.transpose()
    ft_phyla = ft_T.rename(columns=phyla)
    ft_phyla = ft_phyla.groupby(by=ft_phyla.columns, axis=1).sum()
    ft_class = ft_T.rename(columns=clss)
    ft_class = ft_class.groupby(by=ft_class.columns, axis=1).sum()
    ft_phyla = ft_phyla[ft_phyla.columns[ft_phyla.max() > 1]] #only keep those phyla that have a maximum abundance above 1%
    yt = []
    plot_names = [['Alteromonadales', 'Alteromonas macleodii'], 'Bacillales', 'Flavobacteriales', 'Oceanospirillales', 'Rhodobacterales', 'Rhodospirillales', 'Sphingomonadales', 'Vibrionales']
    colormap = mpl.cm.get_cmap('hot', 256)
    norm = mpl.colors.Normalize(vmin=0, vmax=15)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
    for a in range(len(pl)): #for every environment/source combination
        lst = list(ft_phyla.loc[pl[a], :]) #get a list of the relative abundance values for each phylum
        left = 0
        colors = get_cols(len(lst)) #get the colors for plotting (based on number of phyla above 1%)
        if len(lst) > 40: #if we have more than 40 phyla, print a warning
            print('Warning: there are not enough distinct colors for the stacked bar of the dendro and venn plot')
        for b in range(len(lst)): #plot the stacked bar chart for each phyla
            ax2.barh(a, lst[b], left=left, height=0.8, edgecolor='k', color=colors[b])
            left += lst[b]
        yt.append(a)
        ax2.barh(a, 100-left, left=left, height=0.8, edgecolor='k', color='k') #make the sample up to 100% (accounting for those phyla removed due to being <1% either above here or in the R script that werent agglomerated due to low abundance)
        left = 0
        x3 = []
        for c in range(len(plot_names)): #for each of the families that we are specifically interested in (due to them having been reported in a large number of Plastisphere studies)
            try:
                num = list(ft_class.loc[pl[a], plot_names[c]]) #get the abundance in the sample for this family only
                if len(num) > 1: #this accounts for if it is Alteromonadales (where for some reason QIIME has Alteromonas macleodii separately)
                    num = sum(num) #if this was the case, just add Alteromonas macleodii to the other Alteromonadales
                ax3.barh(a, 1, left=left, height=1, edgecolor='k', color=m.to_rgba(num)) #now make the heatmap square, and scaled to this abundance
            except:
                ax3.barh(a, 1, left=left, height=1, edgecolor='k', color=m.to_rgba(ft_class.loc[pl[a], plot_names[c]])) #if it wasn't more than one value, then just make the heatmap square and scaled to this abundance
            x3.append(left+0.5) #add the y value
            left += 1
    #make the legend for the stacked bar of phyla, with the correct colors for each phylum
    handles = []
    phy = list(ft_phyla.columns)
    for a in range(len(phy)):
        handles.append(mlines.Line2D([], [], color=colors[a], marker='o', markersize=8, markeredgecolor='k', label=phy[a], linestyle=' '))
    handles.append(mlines.Line2D([], [], color='k', marker='o', markersize=8, markeredgecolor='k', label='Other', linestyle=' '))
    ax2.legend(handles=handles, bbox_to_anchor=(1.03, -0.06), ncol=round(len(lst)/4), fontsize=fs_small)
    #change the x and y limits and ticks for these plots
    ax2.set_ylim([-0.5, len(pl)-0.5])
    ax3.set_ylim([-0.5, len(pl)-0.5])
    ax2.set_xlim([0, 100])
    ax3.set_xlim([0,8])
    plt.sca(ax2)
    plt.yticks(yt, [])
    plt.sca(ax3)
    plot_names[0] = plot_names[0][0]
    #add the ticks for each of the families of interest
    for a in range(len(plot_names)):
        plot_names[a] = r'$'+plot_names[a]+'$'
    plt.xticks(x3, plot_names, rotation=90, fontsize=fs_main)
    plt.yticks(yt, [])
    #now get only the simpsons diversity index for each (kept separately because we don't want to group it by sample as we want to be able to plot the outliers etc in the boxplot)
    ft_simps = ft_full.loc[['Simpsons'], :]
    for a in range(len(pl)): #for each sample type, get all samples and then make a box plot
        this_simp = ft_simps.loc[:, [pl[a]]]
        ax4.boxplot(this_simp, positions=[a], vert=False, showfliers=False)
    #add x ticks etc and some titles to our plots
    plt.sca(ax4)
    plt.yticks(yt, [])
    ax2.set_xlabel('Relative abundance (%)', fontsize=fs_main)
    ax2.set_title('Mean relative abundance of phyla', fontweight='bold', fontsize=fs_title-2)
    ax3.set_title('Mean relative\nabundance of\nkey orders', fontweight='bold', fontsize=fs_title-2)
    ax4.set_title('Simpsons\nindex\nof diversity', fontweight='bold', fontsize=fs_title-2)
    cb1 = mpl.colorbar.ColorbarBase(ax3_colbar, cmap=colormap, norm=norm, orientation='horizontal') #get the colorbar for the heatmap
    cb1.set_label('Relative abundance (%)', fontsize=fs_main)
    cb1.set_ticks([0, 5, 10, 15])
    cb1.ax.tick_params(labelsize=fs_small)
    plt.savefig(basedir+'/figures/dendro_venn'+ext, dpi=dpi, bbox_inches='tight') #save the figure
    plt.close()
    return

def make_colorbar_fig(basedir):
    '''
    Function to plot a colorbar with the colors used for the early vs late ANCOM and metacoder plots
    Takes as input:
        - basedir (name of the directory to save the figures to)
    '''
    #Set up list of names to plot and the figure
    cols = ['early', 'late', 'aliphatic', 'other plastic', 'unknown plastic', 'biofilm']
    fig = plt.figure(figsize=(4, 2))
    ax1 = plt.subplot(311)
    for a in range(len(cols)): #for each of the treatments that we want to plot
        ax1.bar(a, 1, color=color_source[cols[a]], edgecolor='k', width=1) #make a bar using the color in the dictionary
        if cols[a] in rename_plots: #now rename the treatments in our list based on the rename_plots dictionary
            cols[a] = rename_plots[cols[a]]
        else:
            cols[a] = cols[a].capitalize() #If it wasn't there, then just capitalize the name
    #set x and y limits, ticks and save figure
    plt.xlim([-0.5, 5.5])
    plt.ylim([0,1])
    plt.xticks([0, 1, 2, 3, 4, 5], cols, fontsize=fs_main, rotation=90)
    plt.yticks([])
    plt.savefig(basedir+'/figures/colorbars'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def group_ft_level(ft, level, tax_dict, basedir, rename=False, saving=False): #This script was written specifically for grouping before ancom, but the rename option is there to make it compatible with other functions that might want to regroup but not rename to an ASV (this renaming is only to have a representative sequence for each group)
    '''
    Function to group a feature table to a particular phylogenetic level
    Takes as input:
        - ft (dataframe with sample names as columns and ASVs as rows)
        - level (integer indicating the phylogenetic level to group to, where 1 is phylum and 7 is ASV)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - basedir (name of the directory to save the figures to)
        - rename (boolean indicating whether to rename the grouped ASVs or not, default is False)
        - saving (boolean indicating whether to save the resulting taxonomy file (this is used to rename to the level that we have used for plots made in R), default is False)
    '''
    asv = list(ft.index.values) #get a list of ASVs
    asv_tax_dict, asv_tax = {}, [] #set up dictionary and list
    for a in range(len(asv)): #for each ASV
        asv_tax_dict[asv[a]] = tax_dict[asv[a]][level] #get the taxonomy at the level that we want, and add this to a dictionary
        asv_tax.append(tax_dict[asv[a]][level]) #also add it to a list
    unique_tax = list(set(asv_tax)) #get a list of the unique taxa at this level
    representative_asv, rep_asv_dict = [], {} #set up a dictionary and list
    for a in range(len(unique_tax)): #for each of the unique taxa
        for b in range(len(asv_tax)): #for each ASV
            if unique_tax[a] == asv_tax[b]: #if tthe ASV is the unique taxa that we are looking at
                representative_asv.append([asv[b], unique_tax[a]]) #add this ASV (i.e. the first) as the representative ASV for this taxa
                rep_asv_dict[unique_tax[a]] = asv[b] #and add it to the dictionary
                break #now break this loop, as we only need to do this once
    ft = ft.rename(index=asv_tax_dict) #now rename the feature table at the level that we want
    ft = ft.groupby(by=ft.index, axis=0).sum() #and group all ASVs that now have the same name together (taking a sum of the abundance for each)
    if rename: #if we're renaming the taxa
        ft = ft.rename(index=rep_asv_dict) #rename the feature table again, with the representative ASV for each taxa
        if saving: #if we are saving this file
            with open('taxonomy_name_only.csv', 'w') as f: #open the csv file and add the rows to it
                writer = csv.writer(f)
                writer.writerow(['OTUID', 'Species name'])
                for asv in representative_asv:
                    writer.writerow(asv)
    return ft

def tree_heatmap(ft, meta_dict, basedir, tax_dict, level=7):
    '''
    Function to perform ANCOM tests for statistical significance at a given phylogenetic level between treatments at different time points, and plot the results of this as a tree and heatmap using an R script
    Takes as input:
        - ft (dataframe containing sample names as columns and ASVs as rows)
        - meta_dict (dictionary containing metadata for all samples, with sample names as keys)
        - basedir (name of the directory to save the figures to)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - level (phylogenetic level at which to perform the comparison, where 1 is phylum and 7 is ASV)
    Returns:
        - plots for each environment at the given phylogenetic level. Figure names will take the form environment_ancom_significance_level.pdf
    '''
    ft = ft*100 #multiply the relative abundance for each value in the feature table by 100 
    if level < 7: #if this isn't at the ASV level, then group by taxa
        ft = group_ft_level(ft, level, tax_dict, basedir+'/ancom', rename=True, saving=True)
    samples = list(ft.columns) #get a list of sample names
    env, source, inc_time, source_inc, source_inc_dict = [], [], [], [], {} #set up lists for environment, incubation time, source and a combination of these
    envs = ['marine', 'freshwater', 'aquatic', 'terrestrial']
    os.rename(basedir+'agglom/reduced_tree.tree', 'reduced_tree.tree') #rename the tree to save it to the directory that contains the script (this means that we don't need to put anything additional into the R script, rather than needing to change the file names etc.)
    for a in range(len(samples)): #for each sample
        env.append(meta_dict[samples[a]][3]) #get the environment that it came from, using the metadata dictionary
        source.append(meta_dict[samples[a]][10]) #get the sample type/source for this sample
        inc_time.append(meta_dict[samples[a]][13]) #get the incubation time
        source_inc.append(meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13]) #add a combination of the source and incubation time to another list
        source_inc_dict[samples[a]] = meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13] #and also to a dictionary
    #set up lists of the comparisons to be carried out, as well as the names to use
    comparisons = [['aliphatic early', 'aliphatic late'], ['other plastic early', 'other plastic late'], ['biofilm early', 'biofilm late'],
                   ['aliphatic early', 'other plastic early'], ['aliphatic early', 'unknown plastic early'], ['aliphatic early', 'biofilm early'],
                   ['aliphatic late', 'other plastic late'], ['aliphatic late', 'unknown plastic late'], ['aliphatic late', 'biofilm late'],
                   ['other plastic early', 'unknown plastic early'], ['other plastic early', 'biofilm early'],
                   ['other plastic late', 'unknown plastic late'], ['other plastic late', 'biofilm late'],
                   ['unknown plastic early', 'biofilm early'], ['unknown plastic late', 'biofilm late']]
    comp_names = ['Aliphatic: early vs late', 'Other plastic: early vs late', 'Biofilm: early vs late',
                  'Early: aliphatic vs other plastic', 'Early: aliphatic vs unknown plastic', 'Early: aliphatic vs biofilm', 
                  'Late: aliphatic vs other plastic', 'Late: aliphatic vs unknown plastic', 'Late: aliphatic vs biofilm', 
                  'Early: other plastic vs unknown plastic', 'Early: other plastic vs biofilm', 
                  'Late: other plastic vs unknown plastic', 'Late: other plastic vs biofilm',
                  'Early: unknown plastic vs biofilm', 'Late: unknown plastic vs biofilm']
    for a in range(len(envs)): #for each environment
        this_env, keeping = [], []
        for b in range(len(env)): #for each environment in the list of environments that the samples belong to
            if env[b] == envs[a]: #if it is the same as the environment that we are looking at
                keeping.append(True) #add True to the keeping list
            else:
                keeping.append(False) #otherwise, add False to the keeping list
        env_ft = ft.loc[:, keeping] #get a new dataframe with only those samples that are in envs[a]
        env_ft.rename(columns=source_inc_dict, inplace=True) #rename the remaining samples using the dictionary of sample type and incubation time
        all_comps, colors_list = [], [[], []] #set up lists
        for c in range(len(comparisons)): #for each comparison in the list of comparisons
            try: #try means that the script won't break if this doesn't work - this means that if we don't have the right samples for the comparison in this environment, the script will just continue to the next comparison
                comp_env_ft = env_ft.loc[:, comparisons[c]] #only keep the samples with names matching this comparison
                if len(comp_env_ft.loc[:, comparisons[c][0]].columns) < 2 or len(comp_env_ft.loc[:, comparisons[c][1]].columns) < 2: #if we don't have enough samples remaining
                    all_comps.append(False) #add False to all_comps
                    continue #and continue to the next comparison
                else:
                    all_comps.append(True) #otherwise, add True to all_comps
            except: #if this broke the script
                all_comps.append(False) #add False to all_comps
                continue #and continue to the next comparison
            comp_env_ft = comp_env_ft.transpose() #transpose the feature table
            comp_env_ft = comp_env_ft[comp_env_ft.columns[comp_env_ft.max() > 0]] #only keep taxa if they are present in at least one sample
            comp_env_ft = comp_env_ft.replace(to_replace=0,value = 0.0001) #if we have any zeroes, replace these with 
            comp_env_ft = comp_env_ft.fillna(value=0.0001) #if we have any na values, fill these with 0.0001
            #perform the ANCOM comparison with this transposed feature table, a list of the sample names and holm-bonferroni false discovery rate correction:
            ancom_df, percentile_df = ancom(comp_env_ft, pd.Series(list(comp_env_ft.index.values), index=list(comp_env_ft.index.values)), multiple_comparisons_correction='holm-bonferroni')
            ancom_results, percentile_results, tax_names = ancom_df.values.tolist(), percentile_df.values.tolist(), ancom_df.index.values #get lists of the parts of the ancom data that we are interested in
            significant, medians, not_sig, not_sig_medians = [], [], [], [] #set up some lists
            for d in range(len(ancom_results)): #for each of the ancom results (i.e. each taxa)
                if ancom_results[d][1] == True: #if the ancom result is significant
                    if percentile_results[d][2] == 0.0001 and percentile_results[d][7] == 0.0001: #if the percentile results aren't both 0.0001 (i.e. added in artificially)
                        not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]]) #add the medians to the list of non-significant results
                        not_sig.append(tax_names[d]) #add the name of this taxa to the list of non-significant results
                    else:
                        significant.append(tax_names[d]) #otherwise, add the name of this taxa to the list of significant results
                        medians.append([percentile_results[d][2], percentile_results[d][7]]) #add the medians to the list of significant results
                else: #if the ancom result is not significant
                    not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]]) #add the medians to the list of non-significant results
                    not_sig.append(tax_names[d]) #add the name of this taxa to the list of non-significant results
            differences = []
            for e in range(len(medians)): #for each of the significant sets of medians
                if 0 in medians[e] or 1 in medians[e]: #if either value is 0 or 1
                    medians[e][0] += 0.0001 #add 0.0001 to both numbers
                    medians[e][1] += 0.0001
                if medians[e][0] == medians[e][1]: #if the medians are equal to one another
                    diff = 0 #make the difference be 0
                else: #otherwise
                    diff = math.log2(medians[e][1])/math.log2(medians[e][0]) #calculate the log2 fold change
                    diff = math.pow(2, diff) #and then convert the number back to not log2
                    if diff < 1 and diff > 0: #if the difference is between 0 and 1, then make this a negative number and divide by 1
                        diff = -(1/diff)
                differences.append(diff) #add the difference to the list
            if len(significant) < 1: #if we only have one significant result
                all_comps[c] = False #don't bother plotting this competition
                continue #and continue
            if isinstance(this_env, list): #if this_env is a list (i.e. we haven't made  a dataframe with the results yet)
                this_env = pd.DataFrame(data=differences, index=significant, columns=[comp_names[c]]) #get a dataframe with the results of this competition
            else: #otherwise
                new_df = pd.DataFrame(data=differences, index=significant, columns=[comp_names[c]]) #make a new dataframe with the results of this competition
                this_env = pd.concat([this_env, new_df], axis=1, sort=False) #and add it to the dataframe with the results of the other competitions
            str1, str2 = comparisons[c][0].split(" "), comparisons[c][1].split(" ") #split the comparisons based on spaces, so we can look at the first and second parts of each
            if str1[-1] == 'early' and str2[-1] == 'late': #if the comparison was with the same sample type, between early and late samples
                col1, col2 = color_source['early'], color_source['late'] #set the color to be that of early/late in the dictionary
            else: #otherwise
                if len(str1) > 2: #if the length of the comparison is more than two (the sample name was longer/too long after splitting based on spaces)
                    str1[0] = str1[0]+' '+str1[1] #then add the names back together
                    str1[1] = str1[2]
                    str1 = str1[:-1]
                if len(str2) > 2: #otherwise if this is true for the second part, do as above
                    str2[0] = str2[0]+' '+str2[1]
                    str2[1] = str2[2]
                    str2 = str2[:-1]
                col1 = color_source[str1[0]] #get the color for the first sample type
                col2 = color_source[str2[0]] #and the second sample type 
            colors_list[0].append(str(col1)+", "+str(col2)), colors_list[1].append(comp_names[c]) #add a string to the color list
        if sum(all_comps) > 0: #if we had results for some competitions
            cols_df = pd.DataFrame(data=[colors_list[0]], index=['Colors'], columns=colors_list[1]) #get the colors as a dataframe
            this_env = this_env.append(cols_df) #add this dataframe to the one of all comparisons for this environment
            this_env = this_env.fillna(value=0) #and fill in any na values
            lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'specific_species']
            this_env.to_csv(basedir+'/ancom/'+envs[a]+'_'+lvls[level]+'_ancom_significant.csv') #save this dataframe with the differences for each taxa, for each comparison
            this_env.to_csv('ancom_significant.csv') #and save it as a generic file, that the R script can get to
            os.system("/usr/local/bin/Rscript plot_ancom_tree_heatmap.R") #now call the R script
            os.rename('tree_and_heatmap.pdf', basedir+'/figures/ancom/'+envs[a]+'_ancom_significance_'+lvls[level]+'.pdf') #rename the output of the R script
    os.remove('ancom_significant.csv') #remove this from the folder with the scripts (each environment will be saved in the ancom folder anyway)
    os.remove('taxonomy_name_only.csv') #remove this, too
    make_colorbar_fig(basedir) #make the colorbar
    os.rename('reduced_tree.tree', basedir+'agglom/reduced_tree.tree') #and put the reduced tree back where it came from
    return

def metacoder(ft, tax_dict, meta_dict, basedir):
    '''
    This will separate the feature table to environments and perform pairwise comparisons between the same plastic types at early and late incubation times
    as well as different samples at the same incubation times.
    It takes as input:
        - ft (a feature table containing the abundance of ASVs as rows and samples as columns)
        - tax_dict (a dictionary containing ASV names as keys and taxonomy as values)
        - meta_dict (a dictionary containing all metadata, with sample names as keys)
        - basedir (the base directory to use for saving figures and files to)
    '''
    ft = ft*100 #get the abundances in the feature table as relative abundance
    env, source, inc_time, source_inc, source_inc_dict, source_dict = [], [], [], [], {}, {} #set up some lists
    envs = ['marine', 'freshwater', 'aquatic', 'terrestrial'] #add list of the environments that we want to plot
    samples = list(ft.columns) #get a list of sample names
    for a in range(len(samples)): #for each of these sample names
        env.append(meta_dict[samples[a]][3]) #get the environment this sample was in
        source.append(meta_dict[samples[a]][10]) #get the plastic type/control type of the sample
        source_dict[samples[a]] = meta_dict[samples[a]][10] #add this to a new dictionary containing only this information
        inc_time.append(meta_dict[samples[a]][13]) #get the incubation time (early/late/collection)
        source_inc.append(meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13]) #get the name of the sample type and incubation time together
        source_inc_dict[samples[a]] = meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13] #and also add this to a new dictionary
    comparisons1 = [['aliphatic early', 'aliphatic late'], ['other plastic early', 'other plastic late'], ['biofilm early', 'biofilm late'],
                   ['aliphatic early', 'other plastic early'], ['aliphatic early', 'unknown plastic early'], ['aliphatic early', 'biofilm early'],
                   ['aliphatic late', 'other plastic late'], ['aliphatic late', 'unknown plastic late'], ['aliphatic late', 'biofilm late'],
                   ['other plastic early', 'unknown plastic early'], ['other plastic early', 'biofilm early'],
                   ['other plastic late', 'unknown plastic late'], ['other plastic late', 'biofilm late'],
                   ['unknown plastic early', 'biofilm early'], ['unknown plastic late', 'biofilm late']] #set up a list of all comparisons we will make
    comparisons2 = [['aliphatic', 'other plastic'], ['aliphatic', 'unknown plastic'], ['aliphatic', 'biofilm'], ['aliphatic', 'planktonic'], 
                    ['other plastic', 'unknown plastic'], ['other plastic', 'biofilm'], ['other plastic', 'planktonic'], 
                    ['unknown plastic', 'biofilm'], ['unknown plastic', 'planktonic'], 
                    ['biofilm', 'planktonic']] #and also a second lot of comparisons
    comps = [comparisons1, comparisons2] #add these two sets of comparisons to a list
    asvs = list(ft.index.values) #get a list of all ASVs
    lineage = [] #and an empty list - this will contain the full lineage of each taxon, as required by the R metacoder package for these plots
    for a in range(len(asvs)): #for each asv
        lineage_single = tax_dict[asvs[a]] #get all taxonomic levels for this ASV
        if lineage_single[0] == 'Archaea': #if it is archaea
            ft.drop(asvs[a], axis=0, inplace=True) #remove it from the feature table (there were too few archaea for useful comparisons to be made)
            continue #and continue to the next ASV
        string = '' #add an empty string
        for b in range(len(lineage_single)): #for each level in the taxonomy
            if b > 0 and lineage_single[b] == lineage_single[b-1]: #if this level is the same as the previous level (i.e. it is unclassified at lower taxonomic levels)
                length = b #then this is as many levels as we will want to add
        for b in range(length): #for each of these levels
            if lineage_single[b] != '': #if they're not empty strings
                if string != '': #and we have already add a level to our overall string
                    string += ';' #add a semi-colon to separate levels
                string += lineage_single[b] #and add this level to the string
        lineage.append(string) #and then add it to our overall list
    ft.insert(0, "lineage", lineage, True) #insert this list at the beginning of our feature table
    for a in range(len(envs)): #for each environment
        for z in range(len(comps)): #and each set of comparisons
            comparisons = comps[z] #set the competition we are using
            if comps[z] == comparisons1: #if its the first set of competitions
                treat = source_inc #we're renaming all of our samples using the source and incubation time
                treat_dict = source_inc_dict #and using the dictionary with source and incubation time
            else: #otherwise
                treat = source #we're only renaming them based on the source
                treat_dict = source_dict #and using the source dictionary
            for b in range(len(comparisons)): #for each of our comparisons
                if a < 1 or b > 3:
                    continue
                if len(comparisons[b]) < 3: #if we haven't already
                    comparisons[b].append('lineage') #add the 'lineage' to our list
                keeping = [True] #and tell it that we're keeping the first column (with out lineage information)
                for c in range(len(env)): #for each name in our list of environments
                    if envs[a] == env[c] and treat[c] in comparisons[b]: #if it is the same as the environment that we're currently looking at, and the name of the sample matches the current comparison
                        keeping.append(True) #then we'll be keeping this column
                    else: #otherwise
                        keeping.append(False) #we won't be keeping this column
                env_comp_ft = ft.loc[:, keeping] #now make a new feature table with only those columns that correspond to samples in this comparison
                env_comp_ft.rename(columns=treat_dict, inplace=True) #rename all of the samples to just the general sample type
                names = list(set(env_comp_ft.columns)) #now get a list of those names
                if names == ['lineage'] or len(names) == 2: #if we didn't keep any samples, or we only have one sample and the lineage column
                    continue #then carry on to the next comparison as we don't have enough samples for this
                names_dict = {} #make a new name dictionary
                for d in range(len(names)): #for each name in our list
                    if names[d] == comparisons[b][0]: #if it's equal to the first name in our comparison list
                        names_dict[names[d]] = 'Treat1' #then rename it Treat1
                        names[d] = names[d].replace(' early', '') #and in the list (not in the column) remove the 'early' so we can get the right color for this sample type
                        names[d] = names[d].replace(' late', '') #same for late
                        col1 = color_source[names[d]] #get the color
                    elif names[d] == comparisons[b][1]: #otherwise, if it's equal to the second name in our comparison list
                        names_dict[names[d]] = 'Treat2' #rename it Treat2
                        names[d] = names[d].replace(' early', '') #and again remove the 'early'
                        names[d] = names[d].replace(' late', '') #or late
                        col2 = color_source[names[d]] #and get the color
                if comparisons == comparisons1 and b < 3: #if its one of the comparisons between early and late samples in the same sample type
                    col1, col2 = color_source['early'], color_source['late'] #then get the early and late colors
                env_comp_ft.rename(columns=names_dict, inplace=True) #rename the columns with the dictionary we have created
                command = "/usr/local/bin/Rscript metacoder.R '"+col1+"' '"+col2+"'" #make a string of the command that will be used to run metacoder in R, with the color names
                env_comp_ft.to_csv(basedir+'/metacoder/metacoder_'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.csv') #save the .csv file with the name of this comparison, incase we need to look in the future
                env_comp_ft.to_csv('metacoder.csv') #save the .csv file with the information for running this comparison in metacoder
                os.system(command) #run metacoder in R (through the terminal)
                os.rename('metacoder.pdf', basedir+'/figures/metacoder/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.pdf') #rename the first file, that doesn't have labels on the metacoder plot
                os.rename('metacoder_labels.pdf', basedir+'/figures/metacoder/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'_labels.pdf') #and then rename the second one, that does have labels on the metacoders plot
    return

def generate_rf(X_train, y_train, X_test, y_test, rc='cls', est=10000, n_jobs=None): #generate a random forest based on the data being split to train and test
    '''
    Function to generate a random forest model for one dataset, depending on whether we want to use the classification or regression model
    It takes as input:
        - X_train, y_train, X_test, y_test (this should be the output from the sklearn train_test_split function)
        - rc (a string of either 'cls' or 'reg', depending on whether we want to perform classification or regression, respectively. Default='cls' as this will work on numeric or non-numeric data, where regression requires numeric)
        - est (the number of estimators to use. In general, the higher this number the better/more accurate the model can fit. Default=10000)
        - n_jobs (the number of processors to use. Default=None, meaning 1 is used)
    Returns:
        - rf (the random forest model)
        - rf.score (the classification accuracy of this model)
        - rf.feature_importances_ (the importance values for each included feature/taxa)
        - rf.oob_score_ (the out of bag score - we don't really use this as we have the training dataset, but this may be useful if we didn't have enough sample to have a test/train dataset, and we return it anyway for comparison sake)
    '''
    if rc == 'cls': #if we are using a classification (i.e. discrete categories)
        rf = RandomForestClassifier(n_estimators=est, min_samples_leaf=3, n_jobs=n_jobs, random_state=seed, oob_score=True)
    else: #if we are using a regression (i.e. continuous categories)
        rf = RandomForestRegressor(n_estimators=est, min_samples_leaf=3, n_jobs=n_jobs, random_state=seed, oob_score=True)
    rf.fit(X_train, y_train) #fit out data to either the regressor or classifier
    return rf, rf.score(X_test, y_test), rf.feature_importances_, rf.oob_score_ #return the forest (i.e. features), score (how well it classifies the test data) and feature importances (ASV importances)

def get_single_forest(ft, meta_df, sn, level, est=10000, n_jobs=None):
    '''
    Function to calculate random forest models for each metadata category for a feature table at one phylogenetic level
    Takes as input:
        - ft (feature table/dataframe with sample names as columns and ASVs/taxa names as rows)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - sn (name/string to use to save the results of these random forests)
        - level (phylogenetic level that these comparisons are being performed at, where 1 is phylum and 7 is ASV)
        - est (integer, number of estimators to use to calculate the random forests. Default=10000)
        - n_jobs (integer, number of processors to use for calculations. Default is None, meaning that 1 will be used)
    Returns:
        - scores (scores for how well the random forest models performed for each metadata category)
        - dfs (dataframe containing the importance values for each taxon across all metadata categories)
    '''
    cls_reg = {'Study':'cls', 'Latitude':'reg', 'Longitude':'reg', 'Environment':'cls', 'WaterOrSediment':'cls',
       'LabOrField':'cls', 'IncubationOrCollection':'cls', 'Source':'cls', 'MaterialType':'cls',
       'PlasticTypeSpecific':'cls', 'PlasticTypeGeneral':'cls', 'DEPTH':'reg', 'IncubationTime':'reg',
       'IncubationGeneral':'cls', 'Temperature':'reg', 'Salinity':'reg', 'Light':'cls', 'Season':'cls',
       'PrimerPair':'cls', 'DNAExtraction':'cls', 'CollectionDate':'cls', 'PlasticOnly':'cls'} #dictionary telling it whether each category should use a classification (cls) or regression (reg) algorithm for the random forests
    ft_in = ft #save the original feature table
    ft = ft.rename_axis('ID').values #rename the feature table axis (i.e. remove the indice and column names)
    max_abs_scaler = preprocessing.MaxAbsScaler() #scale all ASVs to the maximum absolute value
    ft_scale = max_abs_scaler.fit_transform(ft) #apply this to our feature table
    meta_split, reg_cls, col_names = [], [], []
    for col in meta_df.columns: #make lists of all of the sample names based on the sample type within each meta category
        meta_split.append(meta_df[col])
        reg_cls.append(cls_reg[col])
        col_names.append(col)
    scores, importances, oob_scores = [], [], []
    count = 0
    for a in range(len(meta_split)): #for each meta-category
        meta = meta_split[a].rename_axis('ID').values #split the metadata dataframe to be only the category we are currently looking at
        if len(set(meta)) == 2 and 'unknown' in meta:
            continue
        keeping = []
        for z in range(len(meta)):
            if meta[z] == 'unknown':
                keeping.append(False)
            else:
                keeping.append(True)
        meta_reduced = meta[keeping]
        ft_scale_reduced = ft_scale[keeping]
        X_train, X_test, y_train, y_test = train_test_split(ft_scale_reduced, meta_reduced, test_size=0.2) #split the data to training and test data (with 80% of samples being used for training and 20% for testing)
        RF, RF_score, RF_importances, RF_oob = generate_rf(X_train, y_train, X_test, y_test, rc=reg_cls[a], est=est, n_jobs=None) #generate the random forest for this meta category
        #print(col_names[a], RF_score) #print the category and random forest score/classification accuracy
        scores.append(RF_score), importances.append(RF_importances), oob_scores.append(RF_oob) #add the scores and importances for this category to the overall lists
        this_rf = {'ASV':ft_in.columns, 'Importance':RF_importances} #make a dataframe with the information on how important each ASV/taxa is
        df_rf = pd.DataFrame(this_rf, columns = ['ASV', 'Importance'])
        df_rf.sort_values(by=['Importance'], ascending=False, inplace=True) #sort the values by importance
        df_rf_score = df_rf.set_index('ASV')
        df_rf_score = df_rf_score
        df_rf = df_rf[df_rf.max(axis=1) > 0.00]
        df_rf = df_rf[:50] #now only keep the top 50 features
        df_rf = df_rf.set_index('ASV') #if we didn't use get_heatmap_random_forest, then we need to uncomment this line
        df_rf = df_rf_score.rename(columns={'Importance':col_names[a]}) #rename the column to be this meta category
        if count == 0:
            dfs = df_rf
            count += 1
        else:
            dfs = dfs.merge(right=df_rf, on='ASV')
    if scores == []:
        return [], []
    scores = pd.DataFrame([scores, oob_scores], index=['Score', 'OOB_score'], columns=list(dfs.columns))
    dfs = pd.concat([dfs, scores], sort=True)
    dfs.to_csv(sn+'.csv')
    return scores, dfs

def get_random_forests(ft, tax_dict, meta_df, basedir, est=10000, n_jobs=None, sn=None):
    '''
    Function to call a function to get random forests for each phylogenetic level, across all metadata categories
    Takes as input:
        - ft (feature table/dataframe with sample names as columns and ASVs/taxa names as rows)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - basedir (name of the directory to save the files to)
        - est (integer, number of estimators to use to calculate the random forests. Default=10000)
        - n_jobs (integer, number of processors to use for calculations. Default is None, meaning that 1 will be used)
    Returns:
        - nothing, but saves a .csv file for each phylogenetic level, containing information on taxon importance at that level
    '''
    this_sn = str(sn)
    asv = list(ft.index.values)
    for a in [1, 2, 3, 4, 5, 6, 7]:
        this_asv = list(asv)
        new_names, other_levels = {}, {}
        for b in range(len(this_asv)):
            this_tax = tax_dict[this_asv[b]]
            prev = ''
            for c in range(len(this_tax)):
                if this_tax[c] != '':
                    prev = this_tax[c]
                else:
                    this_tax[c] = prev
            if a < 7:
                new_names[this_asv[b]] = this_tax[a]
                this_asv[b] = this_tax[a]
            else:
                new_names[this_asv[b]] = this_asv[b]
            other_levels[this_asv[b]] = this_tax[:a]
        this_ft = ft.rename(index=new_names) #now rename the feature table using the new names that we have just defined
        this_ft = this_ft.groupby(by=this_ft.index, axis=0).sum() #group all of the names that are together, taking a sum for each phyla
        this_ft = this_ft.transpose()
        if sn == None:
            sn = basedir+'/random_forest/'+phylo_level_names[a]+'_overall'
        else:
            sn = this_sn+phylo_level_names[a]
        get_single_forest(this_ft, meta_df, sn, a, est=est, n_jobs=None) #get all of the individual forests for this phylogenetic level
    return

def get_single_forest_leave_one_dataset_out(ft, meta_df, sn, level, study, meta_dict, est=10000, n_jobs=None):
    '''
    Function to calculate random forest models for each metadata category for a feature table at one phylogenetic level
    Takes as input:
        - ft (feature table/dataframe with sample names as columns and ASVs/taxa names as rows)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - sn (name/string to use to save the results of these random forests)
        - level (phylogenetic level that these comparisons are being performed at, where 1 is phylum and 7 is ASV)
        - est (integer, number of estimators to use to calculate the random forests. Default=10000)
        - n_jobs (integer, number of processors to use for calculations. Default is None, meaning that 1 will be used)
    Returns:
        - scores (scores for how well the random forest models performed for each metadata category)
        - dfs (dataframe containing the importance values for each taxon across all metadata categories)
    '''
    cls_reg = {'Study':'cls', 'Latitude':'reg', 'Longitude':'reg', 'Environment':'cls', 'WaterOrSediment':'cls',
       'LabOrField':'cls', 'IncubationOrCollection':'cls', 'Source':'cls', 'MaterialType':'cls',
       'PlasticTypeSpecific':'cls', 'PlasticTypeGeneral':'cls', 'DEPTH':'reg', 'IncubationTime':'reg',
       'IncubationGeneral':'cls', 'Temperature':'reg', 'Salinity':'reg', 'Light':'cls', 'Season':'cls',
       'PrimerPair':'cls', 'DNAExtraction':'cls', 'CollectionDate':'cls', 'PlasticOnly':'cls'} #dictionary telling it whether each category should use a classification (cls) or regression (reg) algorithm for the random forests
    ft_in = ft #save the original feature table
    ft = ft.rename_axis('ID').values #rename the feature table axis (i.e. remove the indice and column names)
    max_abs_scaler = preprocessing.MaxAbsScaler() #scale all ASVs to the maximum absolute value
    ft_scale = max_abs_scaler.fit_transform(ft) #apply this to our feature table
    meta_split, reg_cls, col_names = [], [], []
    for col in meta_df.columns: #make lists of all of the sample names based on the sample type within each meta category
        meta_split.append(meta_df[col])
        reg_cls.append(cls_reg[col])
        col_names.append(col)
    scores, importances, oob_scores = [], [], []
    count = 0
    samples = meta_df.index.values
    for a in range(len(meta_split)): #for each meta-category
        meta = meta_split[a].rename_axis('ID').values #split the metadata dataframe to be only the category we are currently looking at
        if len(set(meta)) == 2 and 'unknown' in meta:
            continue
        keeping, keep_train, keep_test = [], [], []
        for z in range(len(meta)):
            if meta[z] == 'unknown':
                keeping.append(False)
            else:
                keeping.append(True)
        meta_reduced = meta[keeping]
        ft_scale_reduced = ft_scale[keeping]
        samples_reduced = samples[keeping]
        for samn in samples_reduced:
            if meta_dict[samn][0] == study:
                keep_test.append(True)
                keep_train.append(False)
            else:
                keep_train.append(True)
                keep_test.append(False)
        X_train = ft_scale_reduced[keep_train]
        y_train = meta_reduced[keep_train]
        X_test = ft_scale_reduced[keep_test]
        y_test = meta_reduced[keep_test]
        if len(y_test) < 2:
            continue
        #X_train, X_test, y_train, y_test = train_test_split(ft_scale_reduced, meta_reduced, test_size=0.2) #split the data to training and test data (with 80% of samples being used for training and 20% for testing)
        RF, RF_score, RF_importances, RF_oob = generate_rf(X_train, y_train, X_test, y_test, rc=reg_cls[a], est=est, n_jobs=None) #generate the random forest for this meta category
        #print(col_names[a], RF_score) #print the category and random forest score/classification accuracy
        scores.append(RF_score), importances.append(RF_importances), oob_scores.append(RF_oob) #add the scores and importances for this category to the overall lists
        this_rf = {'ASV':ft_in.columns, 'Importance':RF_importances} #make a dataframe with the information on how important each ASV/taxa is
        df_rf = pd.DataFrame(this_rf, columns = ['ASV', 'Importance'])
        df_rf.sort_values(by=['Importance'], ascending=False, inplace=True) #sort the values by importance
        df_rf_score = df_rf.set_index('ASV')
        df_rf_score = df_rf_score
        df_rf = df_rf[df_rf.max(axis=1) > 0.00]
        df_rf = df_rf[:50] #now only keep the top 50 features
        df_rf = df_rf.set_index('ASV') #if we didn't use get_heatmap_random_forest, then we need to uncomment this line
        df_rf = df_rf_score.rename(columns={'Importance':col_names[a]}) #rename the column to be this meta category
        if count == 0:
            dfs = df_rf
            count += 1
        else:
            dfs = dfs.merge(right=df_rf, on='ASV')
    if scores == []:
        return [], []
    scores = pd.DataFrame([scores], index=['Score'], columns=list(dfs.columns))
    #dfs.to_csv(sn+'.csv')
    return scores


def get_random_forests_leave_one_dataset_out(ft, tax_dict, meta_df, basedir, meta_dict, est=10000, n_jobs=None, sn=None):
    '''
    Function to call a function to get random forests for each phylogenetic level, across all metadata categories
    Takes as input:
        - ft (feature table/dataframe with sample names as columns and ASVs/taxa names as rows)
        - tax_dict (dictionary containing taxonomy information with ASV names as keys)
        - meta_df (dataframe containing samples as rows and metadata categories as columns)
        - basedir (name of the directory to save the files to)
        - est (integer, number of estimators to use to calculate the random forests. Default=10000)
        - n_jobs (integer, number of processors to use for calculations. Default is None, meaning that 1 will be used)
    Returns:
        - nothing, but saves a .csv file for each phylogenetic level, containing information on taxon importance at that level
    '''
    this_sn = str(sn)
    asv = list(ft.index.values)
    for a in [4, 5, 6, 7]:
        this_asv = list(asv)
        new_names, other_levels = {}, {}
        for b in range(len(this_asv)):
            this_tax = tax_dict[this_asv[b]]
            prev = ''
            for c in range(len(this_tax)):
                if this_tax[c] != '':
                    prev = this_tax[c]
                else:
                    this_tax[c] = prev
            if a < 7:
                new_names[this_asv[b]] = this_tax[a]
                this_asv[b] = this_tax[a]
            else:
                new_names[this_asv[b]] = this_asv[b]
            other_levels[this_asv[b]] = this_tax[:a]
        this_ft = ft.rename(index=new_names) #now rename the feature table using the new names that we have just defined
        this_ft = this_ft.groupby(by=this_ft.index, axis=0).sum() #group all of the names that are together, taking a sum for each phyla
        this_ft = this_ft.transpose()
        studies = list(set(meta_df.loc[:, 'Study']))
        this_df = meta_df.drop('Study', axis=1)
        overall_score = pd.read_csv(basedir+'/random_forest/'+phylo_level_names[a]+'_overall.csv', header=0, index_col=0)
        dfs = overall_score[-2:-1]
        dfs.drop('Study', axis=1, inplace=True)
        for study in studies:
            if sn == None:
                sn = basedir+'/random_forest/leave_one_dataset_out/'+phylo_level_names[a]+'_'+study
            else:
                sn = this_sn+phylo_level_names[a]
            scores = get_single_forest_leave_one_dataset_out(this_ft, this_df, sn, a, study, meta_dict, est=est, n_jobs=None) #get all of the individual forests for this phylogenetic level
            name = rename_plots[study]
            name = name.replace('\n', ' ')
            scores = scores.rename(index={'Score':name})
            sn = None
            dfs = pd.concat([dfs, scores], sort=True)
        dfs["Mean"] = dfs.mean(axis=1)
        dfs.sort_values(by=['Mean'], ascending=False, inplace=True)
        dfs.to_csv(basedir+'/random_forest/leave_one_dataset_out/'+phylo_level_names[a]+'.csv')
        dfs = dfs*100
        dfs.drop('Mean', axis=1, inplace=True)
        dfs = dfs.rename(columns=rename_plots).fillna(value=0).astype('int32').transpose()
        dfs["Mean"] = dfs.mean(axis=1)
        dfs.sort_values(by=['Mean'], ascending=False, inplace=True)
        dfs.drop('Mean', axis=1, inplace=True)
        fig = plt.figure(figsize=(15, 10))
        ax1 = plt.subplot(111)
        annotate_heatmap(ax1, dfs, cmap='inferno', yticks=True, xticks=True, rnd=0, annotate=True)
        plt.savefig(basedir+'/figures/random_forest/leave_one_dataset_out/'+phylo_level_names[a]+'_'+str(est)+ext, dpi=dpi, bbox_inches='tight')
        plt.close()
    return

def sort_by_tax(df, level, other_levels, tax_dict=False, names=False):
    '''
    Function to sort a dataframe based on the higher-level taxonomy (i.e. if the dataframe contains Oceanospirillales, Alteromonadales, Vibrionales, Rhodobacterales and Rhizobiales, then these will be grouped alphabetically to Alpha (Rhizobiales & Rhodobacterales) and Gamma (Alteromonadales, Oceanospirillales & Vibrionales) proteobacteria)
    Takes as input:
        - df (dataframe containing taxa as rows and sample names as columns)
        - level (numeric value of the current level, where 1 is phylum and 7 is ASV)
        - other_levels (dictionary containing taxa at this level as the key and the other higher levels for each as a list)
        - tax_dict (the dictionary containing all taxa - this is False as default and only useful is the level is 7)
        - names (this is redundant but left in so as not to mess up other functions. It is False as default, but is over-written within the function)
    Returns:
        - df_tax (an ordered dataframe containing taxa as rows and sample names as columns)
        - names (a dictionary containing the current taxa as keys and the lowest classification also if the level we are working at is ASV)
    '''
    #function to sort the dataframe df (with index values being the taxon/ASV names) by the higher phylogeny
    current_tax = df.index.values #get the current values
    df_tax = df
    sort_levels, names = [], {}
    for a in range(level): #for each phylogenetic level thaqt is higher
        new_col = []
        for b in range(len(current_tax)):
            if level < 6 and current_tax[b] == 'Alteromonas macleodii':
                current_tax[b] = 'Marinimicrobia (SAR406 clade)'
            if current_tax[b][0] == '_':
                current_tax[b] = current_tax[b].replace('_', '')
            new_col.append([other_levels[current_tax[b]][a], current_tax[b]]) #append the appropriate value for the current higher level as well as the current taxon name
            if tax_dict != False: #if tax_dict isn't False (i.e. if we are looking at the level of ASV - we use tax_dict as here we already worked out the lowest classification for each ASV)
                names[current_tax[b]] = tax_dict[current_tax[b]][-1] #get the lowest classification for each ASV
        new_col = pd.DataFrame(new_col, columns=[str(a), 'ASV']) #turn the classifications at the current level into a dataframe
        new_col = new_col.set_index(['ASV']) #set the ASV (note: this does not necessarily mean this is an ASV, rather it is the taxon names that we are using for this level) column as the index column
        df_tax = new_col.merge(right=df_tax, on='ASV') #and add it to the precious dataframe
        sort_levels.append(str(a)) #add the name of the column that we made with the classifications
    df_tax = df_tax.sort_values(sort_levels, ascending=False) #and now sort the dataframe based on these classifications (this will go in order, so it will keep all proteobacteria together, with all e.g. alpha-, gamma-proteobacteria)
    df_tax.drop(sort_levels, axis=1, inplace=True) #now remove these columns from the sorted dataframe
    return df_tax, names #return them along with the dataframe

def plot_colorbar(ax, cmap, df, name, fs=fs_main, orientation='horizontal'):
    '''
    Function to add a colorbar to an axis
    Takes as input:
        - ax (the axis to plot the colorbar on - note that this will take up the whole axis, rather than just adding a padded version on the side)
        - cmap (the colormap to plot on the bar, as a string not a matplotlib colormap object)
        - df (the data that the colorbar corresponds to - this can be a list of values or a dataframe)
        - name (a string of the name that should me added as a label)
        - fs (the fontsize to use for the axis tick labels - default is the main size set at the top of the script)
        - orientation (string of whether to make this a horizontal or vertical colorbar, default is 'horizontal')
    '''
    #function to plot a colorbar on the axis ax using the colormap cmap
    #it also takes the dataframe of values that it should be based on, a name and the fontsize (using fs_main as default)
    mins, maxs = [], []
    if isinstance(df, list):
        mins.append(df[0])
        maxs.append(df[1])
    else:
        for a in df.columns: #get the minimum and maximum values for each column
            mins.append(df[a].min())
            maxs.append(df[a].max())
    mi = min(mins) #now get the minimum and maximum values overall
    ma = max(maxs)
    if ma < 0:
        mi = 0
        ma = 0.1
    cmap = mpl.cm.get_cmap(cmap, 256) #tell it which colormap to use
    norm = mpl.colors.Normalize(vmin=0, vmax=ma) #the minimum and maximum values to normalize to
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation=orientation) #get the colorbar for the heatmap
    plot_names = []
    for a in range(len(name)):
        if name[a] in rename_plots:
            plot_names.append(rename_plots[name[a]])
        else:
            plot_names.append(name[a])
    cb1.set_label(name, fontsize=fs) #set the label names
    cb1.ax.tick_params(labelsize=fs_small) #and the size of the ticks
    return

def annotate_heatmap(ax, df, cmap='inferno', yticks=True, xticks=True, rnd=1, annotate=True, italics=False, vmax=False, annotate_only=False):
    '''
    Function to plot a heatmap and to optionally annotate each cell of the heatmap
    Takes as input:
        - ax (the axis to plot the heatmap on)
        - df (a dataframe containing the data to plot as floats/integers in the cells. Y labels will be the index values and X labels the column headers)
        - cmap (string of which colormap to use for plotting, default is 'inferno')
        - yticks (optional boolean or string for whether to have ytick labels. Can be True (default), False or 'right' if you want the yticks to show on the right of the plot)
        - xticks (optional boolean or string for whether to have xtick labels. Can be True (default), False or 'top' if you want the xticks to show on the top of the plot)
        - rnd (number of decimal places to round the numbers to when the heatmap is annotated, default is 1)
        - annotate (boolean that decides whether to annotate the heatmap or not. Default is True)
        - italics (boolean of whether to make y values italic or not - i.e. if they are species names)
        - vmax (maximum value to scale color values to in the heatmap, default is False but this should be a number to work)
        - annotate_only (boolean as to whether to only annotate an existing heatmap)
    '''
    #get a heatmap plot using the dataframe df on the axis ax, with the colormap cmap, 
    #with additional options for whether we are to annotate the boxes of the heatmap, 
    #add x or y tick labels and how many decimal places to round each text label to
    plt.sca(ax) #set the axis 
    maxs = [] #get the maximum value of the dataframe
    for a in  df.columns:
        maxs.append(df[a].max())
    ma = max(maxs)
    if annotate: #if we are adding text labels
        for a in range(len(list(df.columns))): #go through each column
            for b in range(len(list(df.index.values))): #and now each row of that column
                num = round(df.iloc[b, a], rnd) #get the value corresponding to the row and columns
                if ma > num:
                    if num < (0.6)*ma: col = 'w' #if the value is below 60% of the maximum value, make the text color white
                    else: col = 'k' #otherwise, make the text color black
                    if cmap == 'PuBu': #if we are using the colormap 'PuBu', do the opposite in terms of text color
                        if num < (0.6)*ma: col = 'k'
                        else: col = 'w'
                else:
                    col = 'k'
                    if num == ma and cmap == 'PuBu': col = 'w'
                plt.text(a+0.5, b+0.5, str(num), ha='center', va='center', color=col, fontsize=fs_small) #and plot this text label
    if annotate_only:
        return
    if vmax != False:
        plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0, vmax=vmax) #make the heatmap
    else:
        plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0) #make the heatmap
    rownames = list(df.index.values)
    for a in range(len(rownames)):
        if len(rownames[a]) > 35:
            row = rownames[a].split(' ')
            if len(row) == 1:
                row = rownames[a].split('-')
            string = ''
            for b in range(len(row)):
                if b > 0 and b < 2:
                    string += '\n'+row[b]+' '
                else:
                    string += row[b]+' '
            rownames[a] = string
    if yticks == 'right':
        plt.yticks(np.arange(0.5, len(df.index), 1), rownames, fontsize=fs_small, linespacing=0.8) #add them based in the row values of the dataframe
        ax.yaxis.tick_right()
    elif yticks: #if we are adding y tick labels
        plt.yticks(np.arange(0.5, len(df.index), 1), rownames, fontsize=fs_main, linespacing=0.8) #add them based in the row values of the dataframe
    else:
        plt.yticks([]) #otherwise, add none
    if xticks == 'top':
        colnames = list(df.columns)
        for col in range(len(colnames)):
            ind = col
            col = colnames[col].split(' ')
            string = ''
            for a in range(len(col)):
                if a > 0:
                    string += '\n'+col[a]
                else:
                    string += col[a]+' '
            colnames[ind] = string
        plt.xticks(np.arange(0.5, len(df.columns), 1), colnames, fontsize=fs_main)
        ax.xaxis.tick_top()
    elif xticks: #if we are adding x tick labels (do the same as for the y tick labels)
        plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90, fontsize=fs_main)
    else:
        plt.xticks([])
    return

def get_heatmap_random_forest(df_rf, ft_in, colnames, basedir, level, other_levels, meta_name, predict, tax_dict, ASV_dict, other_folder=False):
    '''
    This function will get an individual heatmap plot for one category of one random forest, by sorting out all of the data and ordering taxa etc
    It takes as input:
        - df_rf (dataframe containing list of ASVs in one column with a list of their importances in a second column)
        - ft_in (the original feature table containing a list of ASVs as rows and samples as columns)
        - colnames (the groupings within the metadata category)
        - basedir (the base directory to save everything to)
        - level (taxonomic level that we are plotting)
        - meta_name (the metadata category that we are plotting)
        - predict (the accuracy with which this category performed)  
    '''
    #ft = ft_in.transpose() #transpose the feature table
    #df_rf = df_rf.set_index(['ASV']) #in our list of
    ft = ft_in.reset_index() #reset the indices in our feature table
    ft = ft.rename(columns={'index':'ASV'}) #rename the index column
    ft = ft.set_index(['ASV']) #set the index column as ASV again
    numeric, non_numeric = [], []
    for a in colnames: #for each grouping name
        try: #if it is numeric
            colnames[a] = float(colnames[a]) #get it as a float
            colnames[a] = round(colnames[a]) #and then round it to the nearest whole number (the random forest regression was still based on the full values, but we do this to simplify plotting)
            numeric.append(colnames[a])
        except:
            colnames[a] = colnames[a] #otherwise, just leave it as it is (but the try: function doesn' work without an except: option)
            non_numeric.append(colnames[a])
    new_col_order = sorted(list(set(numeric))) + sorted(list(set(non_numeric))) #add together the numeric and non-numeric lists of column names, removing duplicates and sorting them in ascending order
    ft = ft.rename(columns=colnames) #rename the columns
    ft = ft.groupby(by=ft.columns, axis=1).mean() #if we have the same names now, group them together (taking a mean for each abundance value)
    ft = ft.reindex(new_col_order, axis=1) #re-index based on the order of our unique sorted columns
    df_rf.sort_values(by=['Importance'], ascending=False, inplace=True)
    df_rf = df_rf[df_rf.max(axis=1) > 0.00]
    df_rf = df_rf[:50]
    df_rf_tax, n = sort_by_tax(df_rf, level, other_levels) #sort the taxa/ASVs in the random forest data frame by their higher phylogenetic grouping
    merge = df_rf_tax.merge(right=ft, on='ASV') #merge the importance dataframe with the feature table (this means we only keep those ASVs that are important)
    merge.drop('Importance', axis=1, inplace=True) #now remove the importance column again
    merge = merge.div(merge.max(axis=1), axis=0) #normalise within each ASV (to make the colors visible for all even if there are large differences in abundance)
    #set up the figure size based on the number of metadata groupings we have, but checking that it won't be too small if we only have a couple of categories
    width = len(merge.columns)+3
    if width < 15:
        width = 15
    start_col2 = int(width-(width/5))+2
    width_col = int(width/5)
    plt.figure(figsize=(max([width, 25])/3, 20))
    try:
        ax1 = plt.subplot2grid((merge.shape[0], width+3), (3,0), colspan=2, rowspan=merge.shape[0]-2)
    except:
        print(merge, meta_name)
        return
    ax2 = plt.subplot2grid((merge.shape[0], width+3), (3,2), colspan=width, rowspan=merge.shape[0]-2)
    ax3 = plt.subplot2grid((merge.shape[0], width+3), (0,0), colspan=width_col)
    ax4 = plt.subplot2grid((merge.shape[0], width+3), (0,start_col2), colspan=width_col)
    #get the heatmaps for the taxon importances (df_rf_tax) and taxon abundances (merge) separately and with different color schemes
    if level == 7:
        asv_plot = {}
        prev_asv = list(df_rf_tax.index.values)
        for asv in prev_asv:
            if asv in ASV_dict:
                asv_plot[asv] = ASV_dict[asv]+': '+tax_dict[asv][-1]
            else:
                asv_plot[asv] = tax_dict[asv][-1]
        df_rf_tax_plot = df_rf_tax.rename(columns=rename_plots, index=asv_plot)
    else:
        df_rf_tax_plot = df_rf_tax.rename(columns=rename_plots, index=rename_plots)
    merge_plot = merge.rename(columns=rename_plots, index=rename_plots)
    annotate_heatmap(ax1, df_rf_tax_plot, cmap='inferno', rnd=3)
    annotate_heatmap(ax2, merge_plot, cmap='viridis', yticks=False, rnd=1, annotate=False)
    #now get the appropriate color bars for each dataframe
    plot_colorbar(ax3, 'inferno', df_rf_tax_plot, name='Taxon importance')
    plot_colorbar(ax4, 'viridis', merge_plot, name='Normalized abundance')
    if meta_name in rename_plots:
        meta_name_plot = rename_plots[meta_name]
    else:
        meta_name_plot = meta_name
    plt.subplots_adjust(hspace=0.5)
    if other_folder == False:
        ax2.set_title(meta_name_plot+'\nInformative '+phylo_level_names[level]+'\nClassification accuracy = '+str(predict)+'%')
        plt.savefig(basedir+'/figures/random_forest/'+phylo_level_names[level]+'_'+meta_name+'_single_forest'+ext, dpi=dpi, bbox_inches='tight')
    else: 
        of = other_folder.split('_')
        if of[0] in rename_plots:
            of[0] = rename_plots[of[0]]
        if len(of) > 1:
            if of[1] in rename_plots:
                of[1] = rename_plots[of[1]]
            of = of[0]+' '+of[1]
        else:
            of = of[0]
        ax2.set_title(of+'\n'+meta_name_plot+'\nInformative '+phylo_level_names[level]+'\nClassification accuracy = '+str(predict)+'%')
        plt.savefig(basedir+'/figures/random_forest/single_environment/'+other_folder+'_'+phylo_level_names[level]+'_'+meta_name+'_single_forest'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return df_rf_tax

def get_random_forest_plots(ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=False, skip_individual=False):
    '''
    Function to get all random forest plots, i.e. at all Phylogenetic levels and including the overall heatmap showing how accurate all random forests were for all metadata categories at all phylogenetic levels
    Takes as input:
        - ft (feature table as a dataframe with samples as columns and taxa as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to use to open random forest information as well as to save figures to)
        - other_folder (optional input string (file path) to tell the function to save to a folder that is not the default directory for random forests, default is False)
        - skip_individual (optional boolean, default is False. If True no individual heatmaps (i.e. for each metadata category and each phylogenetic level) will be made)
    '''
    asv = list(ft.index.values)
    count = 0
    for a in [1, 2, 3, 4, 5, 6, 7]:
        this_asv = list(asv)
        new_names, other_levels = {}, {}
        for b in range(len(this_asv)):
            this_tax = tax_dict[this_asv[b]]
            prev = ''
            for c in range(len(this_tax)):
                if this_tax[c] != '':
                    prev = this_tax[c]
                else:
                    this_tax[c] = prev
            if a < 7:
                new_names[this_asv[b]] = this_tax[a]
                this_asv[b] = this_tax[a]
            else:
                new_names[this_asv[b]] = this_asv[b]
            other_levels[this_asv[b]] = this_tax[:a]
        this_ft = ft.rename(index=new_names) #now rename the feature table using the new names that we have just defined
        this_ft = this_ft.groupby(by=this_ft.index, axis=0).sum() #group all of the names that are together, taking a sum for each phyla
        if other_folder == False:
            df_rf = pd.read_csv(basedir+'/random_forest/'+phylo_level_names[a]+'_overall.csv', header=0, index_col=0)
        else:
            df_rf = pd.read_csv(basedir+'/random_forest/single_environment/'+other_folder+'_'+phylo_level_names[a]+'.csv', header=0, index_col=0)
        scores, oob_scores = df_rf.loc[['Score'], :], df_rf.loc[['OOB_score'], :]
        meta_names = list(df_rf.columns)
        df_rf.drop(['Score', 'OOB_score'], axis=0, inplace=True)
        df_rf = df_rf.reset_index()
        df_rf = df_rf.rename(columns={'index':'ASV'})
        df_rf.index = df_rf['ASV']
        df_rf.drop(['ASV'], axis=1, inplace=True)
        df_rf_full = df_rf.copy()
        if a < 7:
            df_rf.rename(index={'Alteromonas macleodii':'Marinimicrobia (SAR406 clade)'}, inplace=True)
            df_rf_sorted, names = sort_by_tax(df_rf, a, other_levels, tax_dict=False, names=False)
        else:
            df_rf_sorted, names = sort_by_tax(df_rf, a, other_levels, tax_dict=tax_dict, names=True)
        df_rf_sorted["Mean"] = df_rf_sorted.mean(axis=1)
        df_rf_sorted.sort_values(by=['Mean'], ascending=False, inplace=True)
        df_rf_sorted = df_rf_sorted[df_rf_sorted.max(axis=1) > 0.00]
        df_rf_sorted = df_rf_sorted[:50]
        df_mean = df_rf_sorted.drop(meta_names, axis=1, inplace=False)
        df_rf_sorted.drop(['Mean'], axis=1, inplace=True)
        if other_folder == False:
            fig = plt.figure(figsize=(len(df_rf_sorted.columns)/3, 0.22*df_rf_sorted.shape[0]))
            ax1 = plt.subplot2grid((df_rf_sorted.shape[0], len(df_rf_sorted.columns)+3), (3,2), colspan=len(df_rf_sorted.columns), rowspan=2)
            ax2 = plt.subplot2grid((df_rf_sorted.shape[0], len(df_rf_sorted.columns)+3), (5,0), colspan=2, rowspan=df_rf_sorted.shape[0]-2)
            ax3 = plt.subplot2grid((df_rf_sorted.shape[0], len(df_rf_sorted.columns)+3), (5,2), colspan=len(df_rf_sorted.columns), rowspan=df_rf_sorted.shape[0]-2)
            ax4 = plt.subplot2grid((df_rf_sorted.shape[0], len(df_rf_sorted.columns)+3), (0,0), colspan=int(len(df_rf_sorted.columns)/5))
            n1 = int(len(df_rf_sorted.columns)-len(df_rf_sorted.columns)/5)+3
            ax5 = plt.subplot2grid((df_rf_sorted.shape[0], len(df_rf_sorted.columns)+3), (0,n1), colspan=int(len(df_rf_sorted.columns)/5))
            n2 = int(n1-len(df_rf_sorted.columns)/5)-1
            ax6 = plt.subplot2grid((df_rf_sorted.shape[0], len(df_rf_sorted.columns)+3), (0,n2), colspan=int(len(df_rf_sorted.columns)/5))
            
        else:
            fig = plt.figure(figsize=(5, 0.22*df_rf_sorted.shape[0]))
            ax1 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (7,1), colspan=5, rowspan=1)
            ax2 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (9,0), colspan=1, rowspan=df_rf_sorted.shape[0]-3)
            ax3 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (9,1), colspan=5, rowspan=df_rf_sorted.shape[0]-3)
            ax4 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (0,6), colspan=3)
            ax5 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (3,6), colspan=3)
            ax6 = plt.subplot2grid((df_rf_sorted.shape[0], 9), (6,6), colspan=3)
            
        #plot colorbars for each of the three sets of values (the sums, the importances within categories and the scores, respectively)
        scores = scores*100
        scores = scores.astype('int32')
        plot_colorbar(ax4, 'summer', df_mean, name='Mean taxon importance', fs=fs_small)
        plot_colorbar(ax5, 'inferno', df_rf_sorted, name='Taxon importance', fs=fs_small)
        plot_colorbar(ax6, 'PuBu', scores, name='Classification accuracy (%)', fs=fs_small)
        if other_folder == False:
            ax1.set_title('Top informative                              \n'+phylo_level_names[a]+'                              ', fontsize=fs_main)
        else:
            of = other_folder.split('_')
            if of[0] in rename_plots:
                of[0] = rename_plots[of[0]]
            if len(of) > 1:
                if of[1] in rename_plots:
                    of[1] = rename_plots[of[1]]
                of = of[0]+' '+of[1]
            else:
                of = of[0]
            ax1.set_title('Top informative '+phylo_level_names[a]+'\n'+of, fontsize=fs_main) 
        #now get the heatmaps for each dataframe (same as for the colorbars)
        df_rf_sorted_names = df_rf_sorted.rename(columns=rename_plots, index=rename_plots)
        df_rf_sorted_names = df_rf_sorted_names.div(df_rf_sorted_names.max(axis=1), axis=0)
        
        if a == 7:
            mean_names = list(df_rf_sorted_names.index.values)
            new_names = {}
            for b in range(len(mean_names)):
                if mean_names[b] in tax_dict and mean_names[b] in ASV_dict:
                    new_names[mean_names[b]] = ASV_dict[mean_names[b]]+': '+tax_dict[mean_names[b]][-1]
            df_rf_sorted_names.rename(index=new_names, inplace=True)   
        annotate_heatmap(ax1, scores, cmap='PuBu', yticks=False, xticks=False, rnd=0, annotate=True)
        annotate_heatmap(ax2, df_mean, cmap='summer', yticks=False, xticks=True, rnd=3, annotate=True)
        annotate_heatmap(ax3, df_rf_sorted_names, cmap='inferno', yticks='right', xticks=True, rnd=3, annotate=False, italics=True)
        if other_folder == False:
            plt.savefig(basedir+'/figures/random_forest/'+phylo_level_names[a]+'_overall_forest'+ext, dpi=dpi, bbox_inches='tight') #save the figure
        #else:
        #    plt.savefig(basedir+'/figures/random_forest/single_environment/'+other_folder+'_'+phylo_level_names[a]+ext, dpi=dpi, bbox_inches='tight')
        plt.close()
        for n in range(len(meta_names)):
            if skip_individual:
                continue
            mn = list(meta_names)
            mn.remove(meta_names[n])
            new_names = {}
            colnames = list(this_ft.columns)
            for b in range(len(colnames)):
                if other_folder == False:
                    new_names[colnames[b]] = meta_dict[colnames[b]][n]
                else:
                    if meta_names[n] in ['Source', 'MaterialType', 'PlasticTypeSpecific', 'PlasticTypeGeneral']:
                        new_names[colnames[b]] = meta_dict[colnames[b]][n+7]
                    elif meta_names[n] in ['PlasticOnly']:
                        new_names[colnames[b]] = meta_dict[colnames[b]][n+16]
            score = scores.loc['Score', meta_names[n]]
            this_df = df_rf_full.drop(mn, axis=1, inplace=False)
            try:
                this_df.drop(['Score', 'OOB_score'], axis=0, inplace=True)
            except:
                do_nothing = True
            this_df = this_df.rename(columns={meta_names[n]:'Importance'})
            get_heatmap_random_forest(this_df, this_ft, new_names, basedir, a, other_levels, meta_names[n], score, tax_dict, ASV_dict, other_folder)
        scores = scores.rename(index={'Score':phylo_level_names[a]})
        if count == 0:
            all_scores = scores
            count += 1
        else:
            all_scores = pd.concat([scores, all_scores], sort=True)
    fig = plt.figure(figsize=(5, len(list(all_scores.columns))/5))
    ax1 = plt.subplot(111)
    all_scores = all_scores.transpose()
    all_scores = all_scores.rename(columns=rename_plots, index=rename_plots)
    all_scores = all_scores[all_scores.columns[::-1]]
    all_scores['Mean'] = all_scores.mean(axis=1)
    all_scores = all_scores.sort_values(by=['Mean'], ascending=True)
    all_scores.drop('Mean', axis=1, inplace=True)
    annotate_heatmap(ax1, all_scores, cmap='inferno', yticks=True, xticks=True, rnd=0, annotate=True)
    if other_folder == False:
        plt.savefig(basedir+'/figures/overall_forest_scores'+ext, dpi=dpi, bbox_inches='tight')
    else:
        of = other_folder.split('_')
        if of[0] in rename_plots:
            of[0] = rename_plots[of[0]]
        if len(of) > 1:
            if of[1] in rename_plots:
                of[1] = rename_plots[of[1]]
            if of[1] != 'Laboratory' and of[0] != 'Freshwater':
                of = of[0]+' '+of[1]
            else:
                of = of[0]+'\n'+of[1]
        else:
            of = of[0]
        plt.ylabel(of, fontsize=fs_main, fontweight='bold')
        plt.savefig(basedir+'/figures/random_forest/single_environment/'+other_folder+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return all_scores

def leave_one_out_plots(basedir, meta_df):
    mean_uf = open_pickle(basedir+'/mean_unifrac.dataframe')
    levels = []
    for a in [1, 2, 3]:
        level = pd.read_csv(basedir+'/random_forest/leave_one_dataset_out/'+phylo_level_names[a]+'.csv', header=0, index_col=0)
        level = level*100
        level.drop(['Mean'], axis=1, inplace=True)
        score = pd.DataFrame(level.loc['Score', :])
        level = level - level.loc['Score'].values.squeeze()
        level.drop(['Score'], axis=0, inplace=True)
        levels.append(level)
        score.rename(columns={'Score':phylo_level_names[a]}, inplace=True)
        if a == 1:
            scores = score
        else:
            scores = pd.concat([scores, score], axis=1)
    studies = list(sorted(levels[0].index.values))
    studies_uf = list(sorted(mean_uf.columns))
    scores['Mean'] = scores.mean(axis=1)
    scores = scores.sort_values('Mean', ascending=False)
    meta = list(scores.index.values)
    fig, [[ax1, ax2, ax3, ax4], [ax5, ax6, ax7, ax8], [ax9, ax10, ax11, ax12], [ax13, ax14, ax15, ax16], [ax17, ax18, ax19, ax20]] = plt.subplots(5,4, figsize=(20,20))
    axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18, ax19, ax20]
    studies_meta = list(sorted(set(meta_df['Study'].values)))
    list_studies = list(meta_df.Study)
    
    colormap_40b, colormap_40c, colormap_40a = mpl.cm.get_cmap('tab20b', 256), mpl.cm.get_cmap('tab20c', 256), mpl.cm.get_cmap('tab20', 256)
    norm, norm2, norm3 = mpl.colors.Normalize(vmin=0, vmax=19), mpl.colors.Normalize(vmin=20, vmax=39), mpl.colors.Normalize(vmin=40, vmax=59)
    m, m2, m3 = mpl.cm.ScalarMappable(norm=norm, cmap=colormap_40b), mpl.cm.ScalarMappable(norm=norm2, cmap=colormap_40c), mpl.cm.ScalarMappable(norm=norm3, cmap=colormap_40a)
    color_study = []
    for a in range(len(studies_meta)):
        if a < 20:
            color_study.append(m.to_rgba(a))
        elif a < 40:
            color_study.append(m2.to_rgba(a))
        else:
            color_study.append(m3.to_rgba(a))
    for a in range(len(axs)):
        if a == 19:
            axs[a].axis('off')
            continue
        scr = scores.loc[meta[a], 'Mean']
        meta_name = meta[a]
        if meta[a] in rename_plots: meta[a] = rename_plots[meta[a]]
        meta[a] += ': '+str(int(scr))+'%'
        axs[a].set_title(meta[a], fontsize=fs_title)
        x, y = [], []
        
        for b in range(len(studies)):
            this_study = []
            for c in range(len(levels)):
                num = levels[c].loc[studies[b], meta_name]
                if num < -100:
                    num = -100
                this_study.append(num)
            this_study = np.mean(this_study)
            uf = mean_uf.loc['Weighted unifrac distance', studies_uf[b]]
            ns = list_studies.count(studies_meta[b])
            axs[a].scatter(ns, this_study, color=color_study[b], label=studies[b])
            x.append(ns)
            y.append(this_study)
        x = [0 if math.isnan(n) else n for n in x]
        x = [0 if math.isinf(n) else n for n in x]
        y = [0 if math.isnan(n) else n for n in y]
        y = [0 if math.isinf(n) else n for n in y]
        
        corr, _ = pearsonr(x, y)
        anchored_text = AnchoredText(r'r$^{2}$='+str(round(corr,3)), loc='upper right')
        axs[a].add_artist(anchored_text)
        axs[a].set_ylim([-110, 50])
        if a == 3:
            axs[a].legend(bbox_to_anchor=(1.05, 1), fontsize=fs_main)
        if a > 14:
            axs[a].set_xlabel('Sample number')
        if a % 4 == 0:
            axs[a].set_ylabel('Change in classification accuracy (%)')
    plt.savefig(basedir+'/figures/leave_one_dataset_out'+ext, dpi=dpi, bbox_inches='tight')
    
    return

def get_environment_random_forest(ft, tax_dict, meta_df, meta_dict, basedir, est=10000, n_jobs=1):
    '''
    Function to calculate random forest models for specific metadata categories in different environments.
    It takes as input:
        - ft (feature table data frame containing taxa as rows and sample names as columns)
        - tax_dict (dictionary containing the full taxonomy for each ASV)
        - meta_df (dataframe containing all metadata information, with samples as rows and metadata categories as columns)
        - meta_dict (dictionary containing all metadata with sample names as keys)
        - basedir (directory to save .csv output files to)
        - est (number of estimators to use for calculating random forests, default is 10000. Note that lower numbers will lead to lower accuracy)
        - n_jobs (number of processors to use for calculating random forests, default is 1)
    It uses the other general random forest function but with filtered versions of the metadata and feature tables, and could easily be edited to do this separately for lab/field studies and to do more than just 'PlasticTypeGeneral'
    '''
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    labfield, labindex = [['lab', 'field']], 5
    
    all_samples = list(ft.columns)
    all_samples_meta = list(meta_df.index.values)
    for a in range(len(envs)):
        for b in range(len(labfield)):
            meta_names = list(meta_df.columns)
            keeping = []
            for c in range(len(all_samples)):
                if meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] == labfield[b]:
                    keeping.append(True)
                elif meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] in labfield[b]:
                    keeping.append(True)
                else:
                    keeping.append(False)
            single_ft = ft.loc[:, keeping]
            if len(single_ft.index.values) == 0:
                continue
            single_ft = single_ft[single_ft.max(axis=1) > 0]
            keeping_meta = []
            for d in range(len(all_samples_meta)):
                if meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] == labfield[b]:
                    keeping_meta.append(True)
                elif meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] in labfield[b]:
                    keeping_meta.append(True)
                else:
                    keeping_meta.append(False)
            single_meta_df = meta_df.loc[keeping_meta]
            meta_names.remove('PlasticTypeGeneral')
            single_meta_df.drop(meta_names, axis=1, inplace=True)
            if len(single_meta_df.index.values) == 0:
                continue
            sn = basedir+'/random_forest/single_environment/'+envs[a]+'_'
            get_random_forests(single_ft, tax_dict, single_meta_df, basedir, est=est, n_jobs=n_jobs, sn=sn) 
    return

def get_environment_random_forest_plots(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir):
    '''
    Function to make figures for all random forests calculated on PlasticTypeGeneral in different environments
    Takes as input:
        - ft (feature table as a dataframe with samples as columns and taxa as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to use to open random forest information as well as to save figures to)
    Uses the get_random_forest_plots function from above but using reduced feature and metadata tables 
    This could easily be edited to take lab and field separately as well as together
    '''
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    labfield, labindex = [['lab', 'field']], 5
    all_samples = list(ft.columns)
    all_samples_meta = list(meta_df.index.values)
    for a in range(len(envs)):
        for b in range(len(labfield)):
            if envs[a] == 'freshwater' and labfield[b] == 'lab':
                continue
            meta_names = list(meta_df.columns)
            keeping = []
            for c in range(len(all_samples)):
                if meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] == labfield[b]:
                    keeping.append(True)
                elif meta_dict[all_samples[c]][env_index] == envs[a] and meta_dict[all_samples[c]][labindex] in labfield[b]:
                    keeping.append(True)
                else:
                    keeping.append(False)
            single_ft = ft.loc[:, keeping]
            if len(single_ft.index.values) == 0:
                continue
            single_ft = single_ft[single_ft.max(axis=1) > 0]
            keeping_meta = []
            for d in range(len(all_samples_meta)):
                if meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] == labfield[b]:
                    keeping_meta.append(True)
                elif meta_dict[all_samples_meta[d]][env_index] == envs[a] and meta_dict[all_samples_meta[d]][labindex] in labfield[b]:
                    keeping_meta.append(True)
                else:
                    keeping_meta.append(False)
            single_meta_df = meta_df.loc[keeping_meta]
            single_meta_df.drop(meta_names, axis=1, inplace=True)
            if len(single_meta_df.index.values) == 0:
                continue
            name = envs[a]
            scores = get_random_forest_plots(single_ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=name, skip_individual=True)
            scores.rename(index={'Plastic type (general)':envs[a].capitalize()}, inplace=True)
            if a == 0:
                all_scores = scores
            else:
                all_scores = pd.concat([all_scores, scores], sort=False)
    return all_scores

def get_overall_random_forest_plot(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir):
    '''
    Function to make figures for all random forests calculated on PlasticTypeGeneral in different environments and different metadata categories in all environments
    Takes as input:
        - ft (feature table as a dataframe with samples as columns and taxa as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to use to open random forest information as well as to save figures to)
    '''
    metadata_scores = get_random_forest_plots(ft, tax_dict, ASV_dict, meta_dict, basedir, skip_individual=True)
    env_scores = get_environment_random_forest_plots(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir)
    env_scores['Mean'] = env_scores.mean(axis=1)
    env_scores = env_scores.sort_values(by=['Mean'], ascending=True)
    env_scores.drop('Mean', axis=1, inplace=True)
    
    fig = plt.figure(figsize=(5, 8))
    ax1 = plt.subplot2grid((26, 1), (0, 0), rowspan=20)
    ax2 = plt.subplot2grid((26, 1), (22, 0), rowspan=4)
    
    annotate_heatmap(ax1, metadata_scores, cmap='inferno', yticks=True, xticks=False, rnd=0, annotate=True)
    annotate_heatmap(ax2, env_scores, cmap='inferno', yticks=True, xticks=True, rnd=0, annotate=True)
    ax1.set_title('A', loc='left', fontsize=fs_title, fontweight='bold')
    ax2.set_title('B', loc='left', fontsize=fs_title, fontweight='bold')
    ax1.set_title('Metadata categories', fontsize=fs_title, fontweight='bold')
    ax2.set_title('Plastic type (general)', fontsize=fs_title, fontweight='bold')
    plt.savefig(basedir+'/figures/random_forest_metadata_env'+ext, dpi=dpi, bbox_inches='tight')
    return

def get_ancom(ft_abun):
    '''
    Function to perform ANCOM tests for statistical significance on a given feature table
    Takes as input:
        - ft_abun (a feature table dataframe with sample names as rows and taxon names as columns)
    Returns:
        - significant (a list of taxon names corresponding to all taxa that were determined as being statistically significantly differentially abundant between treatments)
    '''
    other_ft = ft_abun.copy()
    other_ft = other_ft.replace(to_replace=0,value = 0.0001)
    other_ft = other_ft.fillna(value=0.0001)
    other_ft = other_ft.transpose()
    ancom_df, percentile_df = ancom(other_ft, pd.Series(list(other_ft.index.values), index=list(other_ft.index.values)), multiple_comparisons_correction='holm-bonferroni')
    ancom_results, percentile_results, tax_names = ancom_df.values.tolist(), percentile_df.values.tolist(), ancom_df.index.values#, list(percentile_df.columns.values)
    significant, medians, not_sig, not_sig_medians = [], [], [], []
    for d in range(len(ancom_results)):
        if ancom_results[d][1] == True:
            if percentile_results[d][2] == 0.0001 and percentile_results[d][7] == 0.0001:
                not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]])
                not_sig.append(tax_names[d])
                significant.append(tax_names[d])
                medians.append([percentile_results[d][2], percentile_results[d][7]])
            else:
                not_sig_medians.append([percentile_results[d][2], percentile_results[d][7]])
                not_sig.append(tax_names[d])
    return significant

def make_env_rf_plot(ft, tax_dict, basedir, ASV_dict, meta_dict):
    '''
    Function to make heatmaps (of abundance and including whether they were significantly differentially abundant according to ANCOM) of all taxa at all phylogenetic levels with a feature importance value of above 0.01 in the random forest models constructed for separate environments
    Takes as input:
        - ft (feature table dataframe containing sample names as columns and ASV names as rows)
        - tax_dict (dictionary with taxon names as keys and full taxonomic information to be returned as a list)
        - basedir (directory to use to open random forest information as well as to save the resulting figure to)
        - ASV_dict (dictionary containing all ASV assigned readable name with QIIME2 ASV name as keys)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
    Returns:
        - nothing, but saves the figure environment_random_forest to the figures directory
    '''
    envs, ei = ['marine', 'aquatic', 'freshwater', 'terrestrial'], 3
    ft = ft*100
    envs_rf_imp = []
    for a in range(len(envs)):
        if os.path.exists(basedir+'/random_forest/single_environment/'+envs[a]+'.csv'):
            envs_rf_imp.append(pd.read_csv(basedir+'/random_forest/single_environment/'+envs[a]+'.csv', header=0, index_col=0))
            continue
        samples = list(ft.columns)
        keeping, count, rename = [], 0, {}
        for b in range(len(samples)):
            if meta_dict[samples[b]][ei] == envs[a]:
                keeping.append(True)
                rename[samples[b]] = meta_dict[samples[b]][10]
                count += 1 
            else:
                keeping.append(False)
        env_ft = ft.loc[:, keeping]
        env_ft = env_ft[env_ft.max(axis=1) > 0]
        env_ft.rename(columns=rename, inplace=True)
        scores = []
        asv = list(env_ft.index.values)
        for b in [7, 6, 5, 4, 3, 2, 1]:
            lvl = phylo_level_names[b]
            rf = pd.read_csv(basedir+'random_forest/single_environment/'+envs[a]+'_'+lvl+'.csv', header=0, index_col=0)
            score = rf.loc[['Score']].values
            scores.append(score[0][0])
            rf.drop(['Score', 'OOB_score'], axis=0, inplace=True)
            rf = rf[rf.max(axis=1) > 0.005]
            lvl_rename, other_levels = {}, {}
            for c in range(len(asv)):
                tax = tax_dict[asv[c]]
                for d in range(len(tax)):
                    if tax[d] == 'Alteromonas macleodii' and d < 6:
                        for e in range(d, len(tax)):
                            tax[e] = 'Marinimicrobia (SAR406 clade)'
                            break
                    if tax[d] == '':
                        if d < 7:
                            if tax[-1] == tax[5]:
                                tax[d] = tax[-1]
                            else:
                                tax[d] = tax[-1]
                        else:
                            tax[d] = tax[-1]
                lvl_rename[asv[c]] = tax[b]
                adding = tax[:b]
                while len(adding) < 7:
                    adding.append(adding[-1])
                if b == 7:
                    other_levels[asv[c]] = adding
                else:
                    other_levels[tax[b]] = adding
            if b < 7:
                env_ft_level = env_ft.rename(index=lvl_rename)
                env_ft_level = env_ft_level.groupby(by=env_ft_level.index, axis=0).sum()
            else:
                env_ft_level = env_ft
            if 'blank' in env_ft_level.columns:
                env_ft_level.drop('blank', axis=1, inplace=True)
            significant = get_ancom(env_ft_level)
            env_ft_level = env_ft_level.groupby(by=env_ft_level.columns, axis=1).mean()
            rows = [['Taxa', 'kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera', 'species']+list(env_ft_level.columns)+['Importance', 'ANCOM']]
            rf_id = list(rf.index.values)
            for e in range(len(rf_id)):
                if b == 7:
                    this_name = ASV_dict[rf_id[e]]+': '+tax_dict[rf_id[e]][-1]
                else:
                    this_name = rename_plots[lvl]+': '+rf_id[e]
                if b < 6 and rf_id[e] == 'Alteromonas macleodii':
                    print(rf_id[e], this_name)
                    rf_id[e] = 'Marinimicrobia (SAR406 clade)'
                    this_name = 'Marinimicrobia (SAR406 clade)'
                try:
                    this_row = [this_name]+other_levels[rf_id[e]]+list(env_ft_level.loc[rf_id[e], :].values)+list(rf.loc[rf_id[e], :].values)
                except:
                    print(rf_id[e], this_name, envs[a], lvl)
                if rf_id[e] in significant:
                    this_row.append('S')
                else:
                    this_row.append('NS')
                if b == 7:
                    rows.append(this_row)
                else:
                    rows.append(this_row)
            if len(this_row) < 2:
                continue
            new_df = pd.DataFrame(rows[1:], columns=rows[0])
            if b == 7:
                all_rows = new_df
            else:
                all_rows = all_rows.append(new_df, ignore_index=True)
        all_rows = all_rows.sort_values(['kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera', 'species'], ascending=False)
        all_rows.to_csv('Test_rf_plots.csv')
        all_rows = all_rows.set_index('Taxa')
        envs_rf_imp.append(all_rows)
        all_rows.to_csv(basedir+'/random_forest/single_environment/'+envs[a]+'.csv')
    length, width = [], []
    [(length.append(i.shape[0]), width.append(i.shape[1]-9)) for i in envs_rf_imp]
    ma = max(length)
    fig = plt.figure(figsize=(25, int(ma/5)))
    start = 0
    for a in range(len(envs_rf_imp)):
        env_rf = envs_rf_imp[a]
        env_rf = env_rf.drop(['kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera', 'species'], axis=1)
        treats1, treats2 = list(env_rf.columns)[:-2]+['ANCOM'], list(env_rf.columns)[:-2]+['Importance']
        rf = env_rf.drop(treats1, axis=1)
        sig = env_rf.drop(treats2, axis=1)
        env_rf.drop(['Importance', 'ANCOM'], axis=1, inplace=True)
        ax = plt.subplot2grid((ma, 60), (0, start), colspan=env_rf.shape[1], rowspan=env_rf.shape[0])
        start += env_rf.shape[1]
        ax_imp = plt.subplot2grid((ma, 60), (0, start), rowspan=env_rf.shape[0])
        start += 10
        env_rf.rename(columns=rename_plots, inplace=True)
        if 'Other plastic' in env_rf.columns:
            if 'Unknown plastic' in env_rf.columns:
                env_rf = env_rf.reindex(['Aliphatic', 'Other plastic', 'Unknown plastic', 'Biofilm', 'Planktonic'], axis=1)
            else:
                env_rf = env_rf.reindex(['Aliphatic', 'Other plastic', 'Biofilm', 'Planktonic'], axis=1)
        annotate_heatmap(ax_imp, rf, cmap='inferno', yticks=False, xticks=True, rnd=2, annotate=True)
        env_rf = env_rf.div(env_rf.max(axis=1), axis=0)
        annotate_heatmap(ax, env_rf, cmap='viridis', yticks=False, xticks=True, annotate=False)
        nums = np.arange(0.5, env_rf.shape[0]+0.5, 1)
        tax = list(env_rf.index.values)
        empty = []
        for b in range(len(nums)):
            color = 'k'
            if list(sig.loc[tax[b], :])[0][0] == 'S':
                color = '#900C3F'
            ax_imp.text(1.5, nums[b], tax[b], fontsize=fs_main, color=color, va='center')
            empty.append('')
        plt.sca(ax_imp)
        plt.yticks(nums, empty)
        ax_imp.yaxis.tick_right()
        ax.set_title(envs[a].capitalize(), fontsize=fs_title, fontweight='bold')
    plt.savefig(basedir+'/figures/environment_random_forest'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def get_picrust_individual_plot(this_env, axes):
    '''
    Function to perform a T-test and plot a volcano plot for PICRUSt metagenome predictions in one particular environment between plastic and biofilm samples
    Takes as input:
        - this_env (dataframe containing sample names - either 'plastic' or 'biofilm' - as columns and kegg orthologs as rows)
        - axes (the axis object to plot the volcano plot)
    Returns:
        - this_env (the original dataframe with the significance level of the T-test for that kegg ortholog)
        - sig (a list of all kegg orthologs with p values below 0.05 and fold change values either above 2 or below -2)
    '''
    ko = list(this_env.index.values)
    sig = []
    all_p = []
    plt.sca(axes[0])
    for a in range(len(ko)):
        group1 = this_env.loc[ko[a], 'plastic'].values
        group2 = this_env.loc[ko[a], 'biofilm'].values
        t, p = stats.ttest_ind(group1, group2)
        all_p.append(p)
        if p < 0.05:
            if this_env.loc[ko[a], 'FC'] > 2:
                color = '#B90404'
                sig.append(ko[a])
            elif this_env.loc[ko[a], 'FC'] < -2:
                color = '#0498B9'
                sig.append(ko[a])
            else:
                color = '#B1B1B1'
        else:
            color = '#B1B1B1'
        axes[0].scatter(this_env.loc[ko[a], 'FC'], -math.log10(p), marker='o', color=color, s=20)
    this_env['p'] = all_p
    axes[0].set_xlabel(r'log$_{2}$ mean fold change', fontsize=fs_main)
    plt.xscale("symlog")
    plt.yticks(fontsize=fs_main)
    plt.xticks(fontsize=fs_main)
    plt.xlim([-10, 10])
    return this_env, sig

def get_picrust_individual_bar(this_env, KO_dict, axes):
    '''
    This is a function to get a barplot of the top 20 up- and down-regulated kegg orthologs from the picrust2 predicted metagenome in one environment
    Takes as input:
        - this_env (dataframe containing kegg orthologs as rows and 'FC' (fold change) and 'p' (p value from T-test) for the current comparison)
        - KO_dict (a dictionary containing kegg ortholog numbers as keys and returning the description as well as whether the ortholog is classified as xenobiotics degradation, pathogenesis or antimicrobial resistance)
        - axes (the axis object to plot the barplot on)
    '''
    top = this_env.sort_values(['FC'], ascending=False)
    bottom = this_env.sort_values(['FC'], ascending=True)
    top = top[:20]
    bottom = bottom[:20]
    plotting = pd.concat([top, bottom], sort=True)
    plotting = plotting.sort_values(['FC'], ascending=False)
    plotting = plotting[::-1]
    ko = list(plotting.index.values)
    names = []
    y = []
    for a in range(40):
        alpha = 1
        if plotting['FC'].iloc[a] > 2:
            color = '#B90404'
            if plotting['p'].iloc[a] > 0.05:
                alpha = 0.5
        elif plotting['FC'].iloc[a] < -2 and plotting['p'].iloc[a] < 0.05:
            color = '#0498B9'
            if plotting['p'].iloc[a] > 0.05:
                alpha = 0.5
        else:
            color = '#B1B1B1'
            alpha = 0.5
        colormap = mpl.cm.get_cmap('viridis', 256)                
        norm = mpl.colors.Normalize(vmin=-10, vmax=3)
        m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
        col1, col2 = m.to_rgba(plotting['plastic mean'].iloc[a]), m.to_rgba(plotting['biofilm mean'].iloc[a])
        if ko[a] in KO_dict:
            name = ko[a]+': '+KO_dict[ko[a]][0]
            xenpath = KO_dict[ko[a]][1]
            if xenpath == 'Xenobiotics':
                color = '#02A3B6'
            elif xenpath == 'Pathogen':
                color = '#DD6801'
            elif xenpath == 'Antimicrobial resistance':
                color = '#F9E79F'
            adding, spaces = True, 0
            new_name = ''
            for b in range(len(name)):
                if name[b] == '[':
                    adding = False
                    continue
                elif name[b] == ' ':
                    spaces += 1
                else:
                    spaces = 0
                if spaces > 2:
                    adding = False
                if adding:
                    new_name += name[b]
            name = new_name
        else:
            name = ko[a]
        names.append(name)
        y.append(a+1)
        axes[1].barh(a+1, plotting['FC'].iloc[a], color=color, edgecolor='k', alpha=alpha)
    plt.sca(axes[1])
    plt.yticks(y, names, fontsize=fs_small, linespacing=0.9)
    plt.ylim([0.5, 40.5])
    plt.xscale("symlog")
    axes[1].yaxis.tick_right()
    plt.xlim([-10, 10])
    axes[1].set_xlabel(r'log$_{2}$ mean fold change'+'\n'+r'Control biofilms $vs$ Plastic samples', fontsize=fs_main)
    plt.xticks(fontsize=fs_main)
    return names

def picrust_plots(picrust, KO_dict, meta_dict, basedir):
    '''
    Function to make a volcano and barplot of the kegg orthologs from the picrust2 predicted metagenome relating to xenobiotics degradation, pathogenesis or antimicrobial resistance in each environment
    Takes as input:
        - picrust (a dataframe containing sample names as columns and kegg orthologs as rows)
        - KO_dict (a dictionary containing kegg ortholog numbers as keys and returning the description as well as whether the ortholog is classified as xenobiotics degradation, pathogenesis or antimicrobial resistance)
        - meta_dict (dictionary with all metadata  for each sample, with sample names as keys)
        - basedir (directory to save the resulting figure to)
    Returns:
        - nothing, but saves a figure picrust_plot to the base directory
    '''
    envs, env_index = ['marine', 'aquatic', 'freshwater', 'terrestrial'], 3
    source, source_index = [['aliphatic', 'other plastic', 'unknown plastic'], ['biofilm']], 10
    samples = list(picrust.columns)
    fig = plt.figure(figsize=(30,30))
    ax1_volc = plt.subplot2grid((4, 3), (0,0))
    ax2_volc = plt.subplot2grid((4, 3), (1,0))
    ax3_volc = plt.subplot2grid((4, 3), (2,0))
    ax4_volc = plt.subplot2grid((4, 3), (3,0))
    ax1_bar = plt.subplot2grid((4, 6), (0,2))
    ax2_bar = plt.subplot2grid((4, 6), (1,2))
    ax3_bar = plt.subplot2grid((4, 6), (2,2))
    ax4_bar = plt.subplot2grid((4, 6), (3,2))
    all_ax = [ax1_volc, ax1_bar, ax2_volc, ax2_bar, ax3_volc, ax3_bar, ax4_volc, ax4_bar]
    axes = [[ax1_volc, ax1_bar], [ax2_volc, ax2_bar], [ax3_volc, ax3_bar], [ax4_volc, ax4_bar]]
    titles = ['Marine', 'Marine', 'Aquatic', 'Aquatic', 'Freshwater', 'Freshwater', 'Terrestrial', 'Terrestrial']
    ax1_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    ax2_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    ax3_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    ax4_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    for a in range(len(all_ax)):
        all_ax[a].set_title(titles[a], fontsize=fs_title, fontweight='bold')
    for a in range(len(envs)):
        keeping = []
        rename = {}
        for b in range(len(samples)):
            meta = meta_dict[samples[b]]
            if meta[env_index] == envs[a]:
                if meta[source_index] in source[0]:
                    keeping.append(True)
                    rename[samples[b]] = 'plastic'
                elif meta[source_index] in source[1]:
                    keeping.append(True)
                    rename[samples[b]] = 'biofilm'
                else:
                    keeping.append(False)
            else:
                keeping.append(False)
        this_env = picrust.loc[:, keeping]
        this_env.rename(columns=rename, inplace=True)
        this_env = this_env[this_env.max(axis=1) > 0]
        this_env = this_env.replace(to_replace=0,value = 0.0001)
        this_env = this_env.applymap(math.log2)
        this_env['plastic mean'] = this_env[['plastic']].mean(axis=1)
        this_env['biofilm mean'] = this_env[['biofilm']].mean(axis=1)
        this_env['FC'] = this_env['plastic mean']-this_env['biofilm mean']
        this_env, sig = get_picrust_individual_plot(this_env, axes[a])
        if len(this_env[this_env['p'] < 0.05].index.values) > 40:
            this_env = this_env[this_env['p'] < 0.05]
        this_env['plastic std'] = this_env[['plastic']].std(axis=1)
        this_env['biofilm std'] = this_env[['biofilm']].std(axis=1)
        this_env.drop(['plastic', 'biofilm'], axis=1, inplace=True)
        this_env = this_env.sort_values(['plastic mean'], ascending=False)
        plotted_ko = get_picrust_individual_bar(this_env, KO_dict, axes[a])
        for c in range(len(plotted_ko)):
            plotted_ko[c] = [plotted_ko[c]]
        write_csv(basedir+'/picrust/'+envs[a]+'_plotted.csv', plotted_ko)
    plastic = mlines.Line2D([], [], color='#B90404', marker='o', markersize=8, markeredgecolor='k', label='More abundant in\nplastic samples', linestyle=' ')
    biofilm = mlines.Line2D([], [], color='#0498B9', marker='o', markersize=8, markeredgecolor='k', label='More abundant in\ncontrol biofilms', linestyle=' ')
    xenobiotic = mlines.Line2D([], [], color='#02A3B6', marker='s', markersize=8, markeredgecolor='k', label='Xenobiotic\ndegradation', linestyle=' ')
    pathogen = mlines.Line2D([], [], color='#DD6801', marker='s', markersize=8, markeredgecolor='k', label='Pathogenesis', linestyle=' ')
    amr = mlines.Line2D([], [], color='#F9E79F', marker='s', markersize=8, markeredgecolor='k', label='Antimicrobial\nresistance', linestyle=' ')
    ax1_volc.legend(handles=[plastic, biofilm], loc='upper left', fontsize=fs_main)
    ax1_bar.legend(handles=[xenobiotic, pathogen, amr], loc='lower right', fontsize=fs_main)
    #plt.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.subplots_adjust(wspace=0.1, hspace=0.3)
    plt.savefig(basedir+'/figures/picrust_plot'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return
