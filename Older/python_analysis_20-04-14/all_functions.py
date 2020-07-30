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
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
import pandas as pd
import pickle
from scipy.cluster import hierarchy
import scipy.stats as stats
from skbio import DistanceMatrix
from skbio.stats.distance import anosim
from skbio.stats.distance import permanova
from skbio.stats.composition import ancom
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

seed = 3 #random seed
fs_title, fs_main, fs_small = 10, 8, 6 #font sizes to use in figures
color_env = {'marine':'#03A498', 'freshwater':'#03B1FC', 'aquatic':'#9A75FC', 'terrestrial':'#FD9A64'} #colors for plotting each environment
label_env = ['Marine', 'Freshwater', 'Aquatic', 'Terrestrial'] #labels for each environment
ext, dpi = '.png', 600 #extension and dots per inch to save figures  with
color_source = {'aliphatic':'#5F9EA0', 'biofilm':'#FCA94A', 'other plastic':'#8B008B', 'unknown plastic':'#3593FC', 'planktonic':'#F9E79F', 'blank':'gray', 'early':'#FBC704', 'late':'#900C3F', 'collection':'gray'}
phylo_level_names = {0:'kingdoms', 1:'phyla', 2:'classes', 3:'orders', 4:'families', 5:'genera', 6:'species', 7:'ASVs'}
rename_plots = {'AmaralZettler':'Amaral Zettler $et$ $al$. 2015', 'AriasAndres':'Arias-Andres $et$ $al$. 2018', 'Canada':'Canada $et$ $al$. 2020', 'Curren':'Curren & Leong 2019', 'Delacuvellerie':'Delacuvellerie $et$ $al$. 2019', 'DeTender':'De Tender $et$ $al$. 2015', 'DeTenderB':'De Tender $et$ $al$. 2017', 'DussudHudec':'Dussud, Hudec $et$ $al$. 2018', 'DussudMeistertzheim':'Dussud, Meistertzheim $et$ $al$. 2018', 'ErniCassola':'Erni-Cassola $et$ $al$. 2019', 'Esan':'Esan $et$ $al$. 2019', 'Frere':'Frere $et$ $al$. 2018', 'Hoellein':'Hoellein $et$ $al$. 2014', 'HoelleinB':'Hoellein $et$ $al$. 2017', 'Jiang':'Jiang $et$ $al$. 2018', 'Kesy':'Kesy $et$ $al$. 2019', 'Kirstein':'Kirstein $et$ $al$. 2018', 'KirsteinB':'Kirstein $et$ $al$. 2019', 'McCormick':'McCormick $et$ $al$. 2014', 'McCormickB':'McCormick $et$ $al$. 2016', 'Oberbeckmann':'Oberbeckmann $et$ $al$. 2016', 'OberbeckmannB':'Oberbeckmann $et$ $al$. 2018', 'Ogonowski':'Ogonowski $et$ $al$. 2018', 'Parrish':'Parrish $et$ $al$. 2019', 'Pinto':'Pinto $et$ $al$. 2019', 'Pollet':'Pollet $et$ $al$. 2018', 'Rosato':'Rosato $et$ $al$. 2020', 'Syranidou':'Syranidou ', 'SyranidouPE':'Syranidou $et$ $al$. 2017a', 'SyranidouPS':'Syranidou $et$ $al$. 2017b', 'Tagg':'Tagg $et$ $al$. 2019', 'Woodall':'Woodall $et$ $al$. 2018', 'Wu':'Wu $et$ $al$. 2019', 'Xu':'Xu $et$ $al$. 2019', 'Zhang':'Zhang $et$ $al$. 2019', 'WaterOrSediment':'Water or Sediment', 'LabOrField':'Laboratory or Field', 'IncubationOrCollection':'Incubation or Collection', 'MaterialType':'Material type', 'PlasticTypeSpecific':'Plastic type (specific)', 'PlasticTypeGeneral':'Plastic type (general)', 'DEPTH':'Depth', 'IncubationTime':'Incubation time (specific)', 'IncubationGeneral':'Incubation time (general)', 'PrimerPair':'Primer pair', 'DNAExtraction':'DNA extraction method', 'lab':'Laboratory', 'not_plastic':'Not plastic', 'aged_oxope':'Aged Oxo-PE', 'freeliving':'Free living', 'particleassociated':'Particle associated', 'oxope':'Oxo-PE', 'rinse_pe':'PE rinse water', 'rinse_ps':'PS rinse water', 'rinse_wood':'Wood rinse water', 'bhet':'BHET', 'hdpe':'HDPE', 'ldpe':'LDPE', 'na':'NA', 'pa':'PA', 'pe':'PE', 'pes':'PES', 'pestur':'PESTUR', 'pet':'PET', 'phbv':'PHBV', 'pla':'PLA', 'pp':'PP', 'ps':'PS', 'pvc':'PVC', 'san':'SAN', 'w_pe':'Weathered PE', '10:14':'10:14 light:dark', '12:12':'12:12 light:dark', '16:08':'16:08 light:dark', '27F_519R':'27F-519R', '319F_806R':'319F-806R', '338F_806R':'338F-806R', '341F_785R':'341F-785R', '341F_802R':'341F-802R', '341F_806R':'341F-806R', '515F_806R':'515F-806R', '515FY_926R':'515FY-926R', '518F_926R':'518F-926R', '543F_783R':'543F-783R', '967F_1064R':'967F-1064R', 'B969F_BA1406R':'B969F-BA1406R', 'redextract_sigma':'REDExtract-$N$-AmpTM', 'gentra_puregene':'Gentra Puregene', 'purelink':'PureLink', 'powersoil':'PowerSoil', 'phenol_chloroform':'Phenol-Chloroform', 'powerbiofilm':'PowerBiofilm', 'ultraclean_soil':'UltraClean soil', 'fastdna_soil':'FastDNA soil', 'orders':'Order', 'classes':'Class', 'phyla':'Phylum', 'genera':'Genera', 'families':'Family', 'species':'Species', 'ASVs':'ASV', 'kingdoms':'Kingdom', 'PlasticOnly':'Plastic only', '534f_783r':'534F-783R', 'Phenol_chloroform':'Phenol-chloroform'}
meta_name_ind = {'Study':0, 'Latitude':1, 'Longitude':2, 'Environment':3, 'WaterOrSediment':4, 'LabOrField':5, 'IncubationOrCollection':6, 'Source':7, 'MaterialType':8, 'PlasticTypeSpecific':9, 'PlasticTypeGeneral':10, 'DEPTH':11, 'IncubationTime':12, 'IncubationGeneral':13, 'Temperature':14, 'Salinity':15, 'Light':16, 'Season':17, 'PrimerPair':18, 'DNAExtraction':19, 'PlasticOnly':20}
"""
Automate the creation of these dictionaries for each study
"""
name_dict = {'La2020':'Latva\n$et$ $al$. 2020', 'AZ2015':'Amaral-Zettler\n$et$ $al$. 2015', 'AA2018':'Arias-Andres\n$et$ $al$. 2018', 'Ca2020':'Canada\n$et$ $al$. 2020', 'Cu2019':'Curren & Leong\n2019', 'De2019':'Delacuvellerie\n$et$ $al$. 2019', 'DT2015':'De Tender\n$et$ $al$. 2015', 'DT2017':'De Tender\n$et$ $al$. 2017', 'DH2018':'Dussud \n$et$ $al$. 2018a', 'DM2018':'Dussud\n$et$ $al$. 2018b', 'EC2019':'Erni-Cassola\n$et$ $al$. 2019', 'Es2019':'Esan\n$et$ $al$. 2019', 'Fr2018':'Frere\n$et$ $al$. 2018', 'Ho2014':'Hoellein\n$et$ $al$. 2014', 'Ho2017':'Hoellein\n$et$ $al$. 2017', 'Ji2018':'Jiang\n$et$ $al$. 2018', 'Ke2019':'Kesy\n$et$ $al$. 2019', 'Ki2018':'Kirstein\n$et$ $al$. 2018', 'Ki2019':'Kirstein\n$et$ $al$. 2019', 'MC2014':'McCormick\n$et$ $al$. 2014', 'MC2016':'McCormick\n$et$ $al$. 2016', 'Ob2016':'Oberbeckmann\n$et$ $al$. 2016', 'Ob2018':'Oberbeckmann\n$et$ $al$. 2018', 'Og2018':'Ogonowski\n$et$ $al$. 2018', 'Pi2019':'Pinto\n$et$ $al$. 2019', 'Po2018':'Pollet\n$et$ $al$. 2018', 'Ro2020':'Rosato\n$et$ $al$. 2020', 'Sy2019':'Syranidou\n$et$ $al$. 2019', 'SyPE20':'Syranidou\n$et$ $al$. 2017a', 'SyPS20':'Syranidou\n$et$ $al$. 2017b', 'Ta2019':'Tagg\n$et$ $al$. 2019', 'Wo2018':'Woodall\n$et$ $al$. 2018', 'Wr2019':'Wright\n$et$ $al$. 2020', 'Wu2019':'Wu\n$et$ $al$. 2019', 'Xu2019':'Xu\n$et$ $al$. 2019', 'Zh2019':'Zhang\n$et$ $al$. 2019', 'Br2016':'Bryant\n$et$ $al$. 2016', 'Pin201':'Pinnell\n$et$ $al$. 2019', 'Pa2019':'Parrish\n$et$ $al$. 2019'}
name_dict_2 = {'La2020':'Latva $et$ $al$. 2020', 'AZ2015':'Amaral-Zettler $et$ $al$. 2015', 'AA2018':'Arias-Andres $et$ $al$. 2018', 'Ca2020':'Canada $et$ $al$. 2020', 'Cu2019':'Curren & Leong 2019', 'De2019':'Delacuvellerie $et$ $al$. 2019', 'DT2015':'De Tender $et$ $al$. 2015', 'DT2017':'De Tender $et$ $al$. 2017', 'DH2018':'Dussud $et$ $al$.  2018a', 'DM2018':'Dussud $et$ $al$. 2018b', 'EC2019':'Erni-Cassola $et$ $al$. 2019', 'Es2019':'Esan $et$ $al$. 2019', 'Fr2018':'Frere $et$ $al$. 2018', 'Ho2014':'Hoellein $et$ $al$. 2014', 'Ho2017':'Hoellein $et$ $al$. 2017', 'Ji2018':'Jiang $et$ $al$. 2018', 'Ke2019':'Kesy $et$ $al$. 2019', 'Ki2018':'Kirstein $et$ $al$. 2018', 'Ki2019':'Kirstein $et$ $al$. 2019', 'MC2014':'McCormick $et$ $al$. 2014', 'MC2016':'McCormick $et$ $al$. 2016', 'Ob2016':'Oberbeckmann $et$ $al$. 2016', 'Ob2018':'Oberbeckmann $et$ $al$. 2018', 'Og2018':'Ogonowski $et$ $al$. 2018', 'Pi2019':'Pinto $et$ $al$. 2019', 'Po2018':'Pollet $et$ $al$. 2018', 'Ro2020':'Rosato $et$ $al$. 2020', 'Sy2019':'Syranidou $et$ $al$. 2019', 'SyPE20':'Syranidou $et$ $al$. 2017a', 'SyPS20':'Syranidou $et$ $al$. 2017b', 'Ta2019':'Tagg $et$ $al$. 2019', 'Wo2018':'Woodall $et$ $al$. 2018', 'Wr2019':'Wright $et$ $al$. 2020', 'Wu2019':'Wu $et$ $al$. 2019', 'Xu2019':'Xu $et$ $al$. 2019', 'Zh2019':'Zhang $et$ $al$. 2019', 'Br2016':'Bryant $et$ $al$. 2016', 'Pin201':'Pinnell $et$ $al$. 2019', 'Pa2019':'Parrish $et$ $al$. 2019'}
name_env = {'AZ2015':'marine', 'AA2018':'freshwater', 'Ca2020':'aquatic', 'Cu2019':'marine', 'De2019':'marine', 'DT2015':'marine', 'DT2017':'marine', 'DH2018':'marine', 'DM2018':'marine', 'EC2019':'marine', 'Es2019':'terrestrial', 'Fr2018':'marine', 'Ho2014':'freshwater', 'Ho2017':'freshwater', 'Ke2019':'aquatic', 'Ki2018':'marine', 'Ki2019':'marine', 'MC2014':'freshwater', 'MC2016':'freshwater', 'Ob2016':'marine', 'Ob2018':'aquatic', 'Og2018':'marine', 'Pi2019':'marine', 'Po2018':'marine', 'Ro2020':'marine', 'Sy2019':'marine', 'SyPE20':'marine', 'SyPS20':'marine', 'Ta2019':'aquatic', 'Wo2018':'marine', 'Wr2019':'marine', 'Wu2019':'freshwater', 'Xu2019':'marine', 'Zh2019':'terrestrial', 'Pa2019':'aquatic'}

def write_csv(fn, data): #Save a csv file with file name fn and data in a list
    '''
    add doc string here
    '''
    with open(fn, 'w') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)
    return

def write_txt(fn, data): #Save a csv file with file name fn and data in a list
    '''
    add doc string here
    '''
    with open(fn, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for row in data:
            writer.writerow(row)
    return

def open_csv(fn): #open csv with file name fn, returning the data as rows
    '''
    add doc string here
    '''
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    return rows

def open_txt(fn): #open text file with file name fn, returning the data as rows
    '''
    add doc string here
    '''
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f, delimiter='\t', lineterminator='\n'):
            rows.append(row)
    return rows

def write_pickle(fn, data): #write pickle python object with name fn and the data in data
    '''
    add doc string here
    '''
    with open(fn, 'wb') as f:
        pickle.dump(data, f)
    return

def make_KO_dict(fn):
    KO_dict = {}
    KO_dict_full ={}
    with open(fn, 'rU') as f:
        for row in csv.reader(f):
            if row[0] == 'D':
                if len(row) > 3:
                    if row[3] == 'Xenobiotics' or row[3] == 'Pathogen' or row[3] == 'Antimicrobial resistance':
                        if row[3] == 'Xenobiotics':
                            KO_dict[row[1]] = row[2:]
                        else:
                            KO_dict[row[1]] = row[2:]
            KO_dict_full[row[1]] = row[2:]
    write_pickle('KO_dict.dictionary', KO_dict)
    return KO_dict, KO_dict_full

def filter_picrust(picrust, KO_dict, KO_dict_full):
    keeping = []
    ko = list(picrust.index.values)
    for k in range(len(ko)):
        if ko[k] in KO_dict:
            keeping.append(True)
        elif ko[k] not in KO_dict_full:
            if ko[k][0] != 'K':
                keeping.append(True)
                print(ko[k])
                KO_dict[ko[k]] = [ko[k], 'Xenobiotics']
            else:
                keeping.append(False)
        else:
            keeping.append(False)
    picrust = picrust.loc[keeping, :]
    write_pickle('KO_dict.dictionary', KO_dict)
    picrust = picrust[picrust.max(axis=1) > 0]
    write_pickle('picrust.dataframe', picrust)
    return picrust, KO_dict

def format_R(ft): #format the raw data files for phyloseq and unifrac analyses in R
    '''
    Function that transforms the feature table output by QIIME2 to a version that can be used by the R agglomerate script
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
    write_csv('feature_table.csv', ft) #write the new feature table to .csv
    ft = pd.read_csv('feature_table.csv', header=0, index_col=0)
    ft.sort_index(inplace=True)
    ft.sort_index(axis=1, inplace=True)
    write_csv('taxonomy.csv', tax_file) #write the taxonomy to .csv
    write_csv('taxonomy_name_only.csv', tax_name_only) #write the taxonomy to .csv
    write_csv('random_forest/taxonomy_name_only.csv', tax_name_only) #write the taxonomy to .csv
    write_pickle('tax_dict.dictionary', tax_dict) #write the taxon information to a python object
    return ft, tax_dict

def get_meta(meta): #get the information contained in the meta file as a dictionary
    '''
    add doc string here
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
    add doc string here
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

def filter_seqs(ft, sf): #filter the sequences file to contain only the sequences in the given feature table
    '''
    add doc string here
    '''
    ft = pd.read_csv(ft, header=0, index_col=0) #open the feature table using pandas
    asvs = list(ft.index.values) #get a list of the asv names
    seqs_agglom = [] #create a list to add the sequence records to
    for record in SeqIO.parse(sf, "fasta"): #open the sequences file and loop through it
        if record.id in asvs: #if the record id matches one found in the asv list
            seqs_agglom.append(record) #add the record to the list
    SeqIO.write(seqs_agglom, "picrust/sequences_agglom.fasta", "fasta") #save the new list of sequence records
    ft.to_csv('picrust/feature_table_agglom.txt', sep='\t')
    return 'picrust/sequences_agglom.fasta', 'picrust/feature_table_agglom.txt', ft

def study_map(dates, locs, basedir):
    '''
    add doc string here
    '''
    dates = open_csv(dates) #open csv file with study dates
    locs = open_csv(locs) #open csv file with study locations
    plt.figure(figsize=(30,4)) #set up figure
    ax1 = plt.subplot2grid((40,9), (4,0), rowspan=32) #add axis for number of publications
    ax2 = plt.subplot2grid((1,9), (0,1), colspan=2, projection=ccrs.PlateCarree()) #add axis for map
    ax1.set_title('A', loc='left', fontsize=fs_title) #add axis label
    ax2.set_title('B', loc='left', fontsize=fs_title) #add axis label
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
    ax1.set_ylabel('Cumulative number of\nPlastisphere publications', fontsize=fs_main)
    ax1.set_xlabel('Year', fontsize=fs_main)
    ax1.tick_params(axis='both', which='major', labelsize=fs_main)
    ax1.tick_params(axis='both', which='minor', labelsize=fs_main)
    #make a custom legend (this is so we have boxes for colors rather than lines)
    marine = mlines.Line2D([], [], color=color_env['marine'], marker='s', markersize=8, markeredgecolor='k', label=label_env[0], linestyle=' ')
    freshwater = mlines.Line2D([], [], color=color_env['freshwater'], marker='s', markersize=8, markeredgecolor='k', label=label_env[1], linestyle=' ')
    aquatic = mlines.Line2D([], [], color=color_env['aquatic'], marker='s', markersize=8, markeredgecolor='k', label=label_env[2], linestyle=' ')
    terrestrial = mlines.Line2D([], [], color=color_env['terrestrial'], marker='s', markersize=8, markeredgecolor='k', label=label_env[3], linestyle=' ')
    plt.sca(ax1)
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial], loc='upper left', fontsize=fs_main)  
    
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
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial, gap, lab, field, both], loc='lower left', fontsize=fs_main)  
    
    #now save the figure (using the file extension and dpi specified at the top of the file) and close it
    plt.savefig(basedir+'/figures/study_map'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def study_metrics(meta_df, basedir):
    plt.figure(figsize=(10,3))
    ax1 = plt.subplot2grid((1,3), (0,0))
    samples = list(meta_df.index.values)
    envs = ['terrestrial', 'aquatic', 'freshwater', 'marine']
    env_caps = ['Terrestrial', 'Aquatic', 'Freshwater', 'Marine']
    ax2 = plt.subplot2grid((1,3), (0,1))
    ax3 = plt.subplot2grid((1,3), (0,2))
    for a in range(len(envs)):
        bottom = 0
        for b in ['lab', 'field']:
            count = 0
            for c in range(len(samples)):
                if meta_df.loc[samples[c], 'Environment'] == envs[a] and meta_df.loc[samples[c], 'LabOrField'] == b:
                    count += 1
            if b == 'lab':
                ax1.barh(a, count, left=bottom, color=color_env[envs[a]], hatch='//', edgecolor='k')
            else:
                ax1.barh(a, count, left=bottom, color=color_env[envs[a]], edgecolor='k')
            if count > 100:
                ax1.text((count/2)+bottom, a, str(count), color='w', ha='center', va='center', fontsize=fs_main)
            bottom = count+bottom
            if b == 'field':
                ax1.text(bottom+20, a, str(bottom), color='k', ha='left', va='center', fontsize=fs_main)
        pie, inc = [], []
        for c in range(len(samples)):
            if meta_df.loc[samples[c], 'Environment'] == envs[a]:
                pie.append(meta_df.loc[samples[c], 'PlasticTypeGeneral'])
            if meta_df.loc[samples[c], 'Environment'] == envs[a]:
                inc.append(meta_df.loc[samples[c], 'IncubationGeneral'])
        unique = list(set(pie))
        pie_color = [color_source[x] for x in unique]
        pie_count = [pie.count(x) for x in unique]
        pie_count = [(x/sum(pie_count))*100 for x in pie_count]
        left=0
        for d in range(len(pie_count)):
            ax2.barh(a, pie_count[d], left=left, color=pie_color[d], label=unique[d], edgecolor='k')
            fc = 'w'
            if unique[d] == 'planktonic' or unique[d] == 'biofilm':
                fc = 'k'
            if pie_count[d] > 5:
                ax2.text((pie_count[d]/2)+left, a, str(int(pie_count[d])), color=fc, ha='center', va='center', fontsize=fs_main)
            left += pie_count[d]
        unique = list(set(inc))
        if len(unique) == 3:
            unique = ['early', 'late', 'collection']
        pie_color = [color_source[x] for x in unique]
        pie_count = [inc.count(x) for x in unique]
        pie_count = [(x/sum(pie_count))*100 for x in pie_count]
        left = 0
        for d in range(len(pie_count)):
            ax3.barh(a, pie_count[d], left=left, color=pie_color[d], label=unique[d], edgecolor='k')
            fc = 'w'
            if unique[d] == 'early':
                fc = 'k'
            if pie_count[d] > 5:
                ax3.text((pie_count[d]/2)+left, a, str(int(pie_count[d])), color=fc, ha='center', va='center', fontsize=fs_main)
            left += pie_count[d]
    ax1.set_xlabel('Number of samples', fontsize=fs_main)
    plt.sca(ax1)
    plt.yticks([0, 1, 2, 3], env_caps, fontsize=fs_main)
    plt.xticks(fontsize=fs_main)
    plt.xlim([0, 1400])
    lab = patches.Patch(facecolor='w', edgecolor='k', hatch='//', label='Laboratory')
    field = patches.Patch(facecolor='w', edgecolor='k', label='Field')
    ax1.legend(handles=[lab, field], fontsize=fs_main, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=2)
    sources = ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic', 'blank']
    handles = [patches.Patch(facecolor=color_source[x], edgecolor='k', label=x.capitalize()) for x in sources]
    ax2.legend(handles=handles, fontsize=fs_main, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=2)
    plt.sca(ax2)
    plt.yticks([])
    plt.xticks(fontsize=fs_main)
    plt.xlabel('Relative abundance (%)', fontsize=fs_main)
    plt.xlim([0, 100])
    times = ['early', 'late', 'collection']
    handles = [patches.Patch(facecolor=color_source[x], edgecolor='k', label=x.capitalize()) for x in times]
    ax3.legend(handles=handles, fontsize=fs_main, bbox_to_anchor=(0., 0, 1., -.25), loc='upper left', borderaxespad=0., mode='expand', ncol=3)
    plt.sca(ax3)
    plt.yticks([])
    plt.xticks(fontsize=fs_main)
    plt.xlabel('Relative abundance (%)', fontsize=fs_main)
    plt.xlim([0, 100])
    plt.subplots_adjust(wspace=0.1)
    ax1.set_title('C', loc='left', fontsize=fs_title)
    ax2.set_title('D', loc='left', fontsize=fs_title)
    ax3.set_title('E', loc='left', fontsize=fs_title)
    plt.savefig(basedir+'/figures/'+'sample_metrics'+ext, dpi=dpi, bbox_inches='tight')
    #plt.close()
    return

def transform_for_NMDS(similarities, n_jobs=1): #transform the similarity matrix to 2D space (n_components=2)
    '''
    add doc string here
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
    add doc string here
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
        string = ans[0]+'\n'+ans[1]+'='+str(round(ans[4], 3))+r', $p$='+str(round(ans[5], 3))
        string += '\n'+perm[0]+'\n'+perm[1]+'='+str(round(perm[4], 3))+r', $p$='+str(round(perm[5], 3))
        plt.text(0.05, 0.98, string, ha='left', va='top', transform = ax.transAxes, fontsize=fs_main, fontweight='bold')
        
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

def nmds_plot_study_env(dist_matr_fn_w, dist_matr_fn_uw, meta_dict, basedir, n_jobs):
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
        dist_matr = pd.read_csv(dist[z], header=0, index_col=0) #read in the distance matrix
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
            plt.legend(handles=handles, loc='upper right', fontsize=fs_main)
        else:
            plt.sca(ax[z][0])
            plt.legend(handles=handles, loc='lower left', fontsize=fs_main)
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, study, study_index, ax[z][1], 'upper right', color_study, names, meta_dict, second_filter, second_filter_ind, npos, n_jobs=n_jobs, get_stats=True)
        if z == 0:
            plt.sca(ax[z][1])
            plt.legend(handles=handles, bbox_to_anchor=(1.05,1.025), fontsize=fs_main)
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
    return

def plot_box(ax, l, r, b, t, line_col):
    '''
    add doc string here
    '''
    plt.sca(ax)
    plt.plot([l, r], [b, b], color=line_col, lw=2)
    plt.plot([l, r], [t, t], color=line_col, lw=2)
    plt.plot([l, l], [t, b], color=line_col, lw=2)
    plt.plot([r, r], [t, b], color=line_col, lw=2)
    return

def similarity_heatmap_combined(dist_matr_fn_w, dist_matr_fn_uw, basedir):
    '''
    add doc string here
    '''
    dist_matr_files = [dist_matr_fn_w, dist_matr_fn_uw]
    fig = plt.figure(figsize=(12,12))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
    ax = [ax1, ax2]
    for z in range(len(dist_matr_files)):
        plt.sca(ax[z])
        dist_matr_fn = dist_matr_files[z]
        if dist_matr_fn[-4:] == '.csv':
            dist_matr = open_csv(dist_matr_fn) #get distance matrix of csv file
        else:
            dist_matr = open_txt(dist_matr_fn) #get distance matrix of txt file
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
        for a in range(len(each_study)):
            colors = []
            x, y, bottom = [], [], []
            for b in range(len(each_study[a])):
                color = m.to_rgba(each_study[a][b])
                colors.append(color)
                x.append(b+1)
                y.append(1)
                bottom.append(a)
                plt.bar(b+1, 1, bottom=a, color=color, edgecolor='k', width=1)
            ally.append(a+0.5)
        cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colormap),ax=ax[z], shrink=0.4, pad=0.02, orientation='vertical')
        cb.ax.tick_params(labelsize=fs_small)
        if z == 0:
            cb.set_label('Weighted\nunifrac distance', fontsize=fs_main)
        else:
            cb.set_label('Unweighted\nunifrac distance', fontsize=fs_main)
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
    #plt.close()
    return

def generate_rf(X_train, y_train, X_test, y_test, rc='cls', est=10000, n_jobs=None): #generate a random forest based on the data being split to train and test
    '''
    add doc string here
    '''
    if rc == 'cls': #if we are using a classification (i.e. discrete categories)
        rf = RandomForestClassifier(n_estimators=est, min_samples_leaf=3, n_jobs=n_jobs, random_state=seed, oob_score=True)
    else: #if we are using a regression (i.e. continuous categories)
        rf = RandomForestRegressor(n_estimators=est, min_samples_leaf=3, n_jobs=n_jobs, random_state=seed, oob_score=True)
    rf.fit(X_train, y_train) #fit out data to either the regressor or classifier
    return rf, rf.score(X_test, y_test), rf.feature_importances_, rf.oob_score_ #return the forest (i.e. features), score (how well it classifies the test data) and feature importances (ASV importances)

def annotate_heatmap(ax, df, cmap='inferno', yticks=True, xticks=True, rnd=1, annotate=True, italics=False, vmax=False):
    '''
    add doc string here
    '''
    #get a heatmap plot using the dataframe df on the axis ax, with the colormap cmap, 
    #with additional options for whether we are to annotate the boxes of the heatmap, 
    #add x or y tick labels and how many decimal places to round each text label to
    plt.sca(ax) #set the axis 
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
    return

def sort_by_tax(df, level, other_levels, tax_dict=False, names=False):
    '''
    add doc string here
    '''
    #function to sort the dataframe df (with index values being the taxon/ASV names) by the higher phylogeny
    current_tax = df.index.values #get the current values
    df_tax = df
    sort_levels, names = [], {}
    for a in range(level): #for each phylogenetic level thaqt is higher
        new_col = []
        for b in range(len(current_tax)):
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
    add doc string here
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

def get_heatmap_random_forest(df_rf, ft_in, colnames, basedir, level, other_levels, meta_name, predict, tax_dict, ASV_dict, other_folder=False):
    '''
    This function will get an individual heatmap plot for one category of one random forest
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
            colnames[a] = floa(colnames[a]) #get it as a float
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
    """
    if level == 7: #if at ASV level
        #get the phylogenetic tree and heatmap
        ft_R = merge.reset_index()
        ft_R.to_csv(basedir+'/random_forest_R/random_forest.csv', index=False)
        os.system("/usr/local/bin/Rscript plot_RF_tree_heatmap.R") #run the R script that plots the phylogenetic tree (of only the important ASVs) and heatmap for this meta category
        os.rename(basedir+'/random_forest_R/tree_and_heatmap.pdf', basedir+'/figures/random_forest/'+phylo_level_names[level]+'_'+meta_name+'_tree_and_heatmap.pdf') #rename the resulting figure and move it to the figures directory
    """
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

def get_single_forest(ft, meta_df, sn, level, est=10000, n_jobs=None):
    '''
    add doc string here
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
    dfs = pd.concat([dfs, scores])
    dfs.to_csv(sn+'.csv')
    return scores, dfs

def get_random_forests(ft, tax_dict, meta_df, basedir, est=10000, n_jobs=None, sn=None):
    '''
    add doc string here
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

def get_environment_random_forest(ft, tax_dict, meta_df, meta_dict, basedir, est=10000, n_jobs=1):
    '''
    Add doc string here
    '''
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    labfield, labindex = ['lab', 'field', ['lab', 'field']], 5
    
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
            meta_names.remove('Source')
            meta_names.remove('PlasticTypeSpecific')
            meta_names.remove('MaterialType')
            meta_names.remove('PlasticTypeGeneral')
            meta_names.remove('PlasticOnly')
            single_meta_df.drop(meta_names, axis=1, inplace=True)
            if len(single_meta_df.index.values) == 0:
                continue
            if b < 2:
                sn = basedir+'/random_forest/single_environment/'+envs[a]+'_'+labfield[b]+'_'
            else:
                sn = basedir+'/random_forest/single_environment/'+envs[a]+'_'
            get_random_forests(single_ft, tax_dict, single_meta_df, basedir, est=est, n_jobs=n_jobs, sn=sn)
    return

def get_random_forest_plots(ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=False, skip_individual=False):
    '''
    add doc string here
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
        else:
            plt.savefig(basedir+'/figures/random_forest/single_environment/'+other_folder+'_'+phylo_level_names[a]+ext, dpi=dpi, bbox_inches='tight')
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
                print('Scores already dropped for '+phylo_level_names[a]+' '+other_folder)
            this_df = this_df.rename(columns={meta_names[n]:'Importance'})
            get_heatmap_random_forest(this_df, this_ft, new_names, basedir, a, other_levels, meta_names[n], score, tax_dict, ASV_dict, other_folder)
        scores = scores.rename(index={'Score':phylo_level_names[a]})
        if count == 0:
            all_scores = scores
            count += 1
        else:
            all_scores = pd.concat([scores, all_scores])
    fig = plt.figure(figsize=(5, len(list(all_scores.columns))/5))
    ax1 = plt.subplot(111)
    all_scores = all_scores.transpose()
    all_scores = all_scores.rename(columns=rename_plots, index=rename_plots)
    all_scores = all_scores[all_scores.columns[::-1]]
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
    return

def get_environment_random_forest_plots(ft, meta_df, tax_dict, ASV_dict, meta_dict, basedir):
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    labfield, labindex = ['lab', 'field', ['lab', 'field']], 5
    all_samples = list(ft.columns)
    all_samples_meta = list(meta_df.index.values)
    for a in range(len(envs)):
        if envs[a] == 'marine':
            continue
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
            #meta_names.remove('Source')
            #meta_names.remove('MaterialType')
            #meta_names.remove('PlasticTypeSpecific')
            meta_names.remove('PlasticOnly')
            single_meta_df.drop(meta_names, axis=1, inplace=True)
            if len(single_meta_df.index.values) == 0:
                continue
            if b == 2:
                name = envs[a]
            else:
                name = envs[a]+'_'+labfield[b]
            get_random_forest_plots(single_ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=name)
            """
            try:
                if name != 'marine_field':
                    get_random_forest_plots(single_ft, tax_dict, ASV_dict, meta_dict, basedir, other_folder=name)
            except:
                print('Couldnt get '+name+' plots')
            """
    return

def group_ft_level(ft, level, tax_dict, basedir, rename=False, saving=False): #This script was written specifically for grouping before ancom, but the rename option is there to make it compatible with other functions that might want to regroup but not rename to an ASV (this renaming is only to have a representative sequence for each group)
    '''
    add doc string here
    '''
    asv = list(ft.index.values)
    asv_tax_dict, asv_tax = {}, []
    for a in range(len(asv)):
        asv_tax_dict[asv[a]] = tax_dict[asv[a]][level]
        asv_tax.append(tax_dict[asv[a]][level])
    unique_tax = list(set(asv_tax))
    representative_asv, rep_asv_dict = [], {}
    for a in range(len(unique_tax)):
        for b in range(len(asv_tax)):
            if unique_tax[a] == asv_tax[b]:
                representative_asv.append([asv[b], unique_tax[a]])
                rep_asv_dict[unique_tax[a]] = asv[b]
                break
    ft = ft.rename(index=asv_tax_dict)
    ft = ft.groupby(by=ft.index, axis=0).sum()
    if rename:
        ft = ft.rename(index=rep_asv_dict)
        if saving:
            with open(basedir+'/taxonomy_name_only.csv', 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['OTUID', 'Species name'])
                for asv in representative_asv:
                    writer.writerow(asv)
    return ft

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
    add doc string here
    #get a venn diagram with specified colors and names, with up to 5 ellipses (these are not plotted if that sample type is not present)
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
    This is a function to get a list of either 20 or 40 colors depending on the input 'num'
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

def bar_dendro_venn(ft, meta_dict, basedir, tax_dict):
    '''
    Add doc string here
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
    ft = pd.read_csv(ft, header=0, index_col=0) #get the feature table
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
        indiv_div.append(get_diversity('Simpsons', list(ft.loc[:, cn])))
    #rename all sample names by a combination of the environment and source
    ft.rename(columns=env_source_dict, inplace=True)
    env_source = list(set(env_source)) #now get a list of all unique environments and sources
    ft = ft.reset_index()
    ft = ft.append(pd.Series(indiv_div, index=ft.columns), ignore_index=True) #add the simpsons diversity to the overall dataframe
    ft.index = ft['index']
    ft.drop('index', axis=1, inplace=True)
    ft_dendro = ft.drop('Simpsons', axis=0) #make a new dataframe without the simpsons index of diversity
    ft_dendro = ft_dendro.groupby(by=ft_dendro.columns, axis=1).mean() #group this by the environment/source name, taking a mean rather than a sum for each ASV
    plt.sca(ax1)
    ft_dendro_T = ft_dendro.transpose() #transpose the dataframe
    Z = hierarchy.linkage(ft_dendro_T, 'average', metric='braycurtis') #calculate the bray-curtis distance between the means of sample groupings
    mpl.rcParams['lines.linewidth'] = 2
    hierarchy.set_link_color_palette(['k'])
    dn = hierarchy.dendrogram(Z, above_threshold_color='k', orientation='left') #plot this dendrogram of sample groupings
    y_labels, locs, xlocs, labels = list(ax1.get_yticklabels()), list(ax1.get_yticks()), list(ax1.get_xticks()), [] #get all labels so that we can work out which order samples are plotted in (and apply this to the other plots)
    for y in y_labels:
        labels.append(y.get_text())
    plot_labels, pl, colors = [], [], []
    grouping = list(ft_dendro.columns)
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
        ax1.text(xlocs[0]-0.75, locs[a], plot_labels[a], fontsize=fs_main, bbox=dict(facecolor=colors[a], alpha=0.4, pad=0.5, edgecolor='w'), ha='center', va='center')
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
    ft_phyla = ft_dendro_T.rename(columns=phyla)
    ft_phyla = ft_phyla.groupby(by=ft_phyla.columns, axis=1).sum()
    ft_class = ft_dendro_T.rename(columns=clss)
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
    ax2.legend(handles=handles, bbox_to_anchor=(1.03, -0.06), ncol=round(len(lst)/4), fontsize=fs_main)
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
    ft_simps = ft.loc[['Simpsons'], :]
    for a in range(len(pl)): #for each sample type, get all samples and then make a box plot
        this_simp = ft_simps.loc[:, [pl[a]]]
        ax4.boxplot(this_simp, positions=[a], vert=False, showfliers=False)
    #add x ticks etc and some titles to our plots
    plt.sca(ax4)
    plt.yticks(yt, [])
    ax2.set_xlabel('Relative abundance (%)', fontsize=fs_main)
    ax2.set_title('Mean relative abundance of phyla', fontweight='bold', fontsize=fs_title)
    ax3.set_title('Mean relative\nabundance of\nkey families', fontweight='bold', fontsize=fs_title)
    ax4.set_title('Simpsons index\nof diversity', fontweight='bold', fontsize=fs_title)
    cb1 = mpl.colorbar.ColorbarBase(ax3_colbar, cmap=colormap, norm=norm, orientation='horizontal') #get the colorbar for the heatmap
    cb1.set_label('Relative abundance (%)', fontsize=fs_main)
    cb1.set_ticks([0, 5, 10, 15])
    cb1.ax.tick_params(labelsize=fs_small)
    plt.savefig(basedir+'/figures/dendro_venn'+ext, dpi=dpi, bbox_inches='tight') #save the figure
    plt.close()
    return

def make_colorbar_fig(basedir):
    cols = ['early', 'late', 'aliphatic', 'other plastic', 'unknown plastic', 'biofilm']
    fig = plt.figure(figsize=(4, 2))
    ax1 = plt.subplot(311)
    for a in range(len(cols)):
        ax1.bar(a, 1, color=color_source[cols[a]], edgecolor='k', width=1)
        if cols[a] in rename_plots:
            cols[a] = rename_plots[cols[a]]
        else:
            cols[a] = cols[a].capitalize()
    plt.xlim([-0.5, 5.5])
    plt.ylim([0,1])
    plt.xticks([0, 1, 2, 3, 4, 5], cols, fontsize=fs_main, rotation=90)
    plt.yticks([])
    plt.savefig(basedir+'/figures/colorbars'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def tree_heatmap(ft, meta_dict, basedir, tax_dict, level=7):
    '''
    add doc string here
    '''
    ft = pd.read_csv(ft, header=0, index_col=0)
    ft = ft*100
    if level < 7:
        ft = group_ft_level(ft, level, tax_dict, basedir+'/ancom', rename=True, saving=True)
    samples = list(ft.columns)
    env, source, inc_time, source_inc, source_inc_dict = [], [], [], [], {}
    envs = ['marine', 'freshwater', 'aquatic', 'terrestrial']
    for a in range(len(samples)):
        env.append(meta_dict[samples[a]][3])
        source.append(meta_dict[samples[a]][10])
        inc_time.append(meta_dict[samples[a]][13])
        source_inc.append(meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13])
        source_inc_dict[samples[a]] = meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13]
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
    for a in range(len(envs)):
        this_env = []
        keeping = []
        for b in range(len(env)):
            if env[b] == envs[a]:
                keeping.append(True)
            else:
                keeping.append(False)
        env_ft = ft.loc[:, keeping]
        env_ft.rename(columns=source_inc_dict, inplace=True)
        all_comps, colors_list = [], [[], []]
        for c in range(len(comparisons)):
            try:
                comp_env_ft = env_ft.loc[:, comparisons[c]]
                if len(comp_env_ft.loc[:, comparisons[c][0]].columns) < 2 or len(comp_env_ft.loc[:, comparisons[c][1]].columns) < 2:
                    all_comps.append(False)
                    continue
                else:
                    all_comps.append(True)
            except:
                all_comps.append(False)
                continue
            comp_env_ft = comp_env_ft.transpose()
            comp_env_ft = comp_env_ft[comp_env_ft.columns[comp_env_ft.max() > 0]]
            comp_env_ft = comp_env_ft.replace(to_replace=0,value = 0.0001)
            comp_env_ft = comp_env_ft.fillna(value=0.0001)
            ancom_df, percentile_df = ancom(comp_env_ft, pd.Series(list(comp_env_ft.index.values), index=list(comp_env_ft.index.values)), multiple_comparisons_correction='holm-bonferroni')
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
            differences = []
            for e in range(len(medians)):
                if medians[e][0] == 0:
                    medians[e][0] = 0.0001
                elif medians[e][1] == 0:
                    medians[e][1] = 0.0001
                if medians[e][0] == 1:
                    medians[e][0] = 1.0001
                elif medians[e][1] == 1:
                    medians[e][1] = 1.0001
                if medians[e][0] == medians[e][1]:
                    diff = 0
                else:
                    diff = math.log2(medians[e][1])/math.log2(medians[e][0])
                    diff = math.pow(2, diff)
                    if diff < 1 and diff > 0:
                        diff = -(1/diff)
                differences.append(diff)
            if len(significant) < 1:
                all_comps[c] = False
                continue
            if isinstance(this_env, list):
                this_env = pd.DataFrame(data=differences, index=significant, columns=[comp_names[c]])
            else:
                new_df = pd.DataFrame(data=differences, index=significant, columns=[comp_names[c]])
                this_env = pd.concat([this_env, new_df], axis=1, sort=False)
            str1, str2 = comparisons[c][0].split(" "), comparisons[c][1].split(" ")
            if str1[-1] == 'early' and str2[-1] == 'late':
                col1, col2 = '#ff0000', '#25df50'
            else:
                if len(str1) > 2:
                    str1[0] = str1[0]+' '+str1[1]
                    str1[1] = str1[2]
                    str1 = str1[:-1]
                if len(str2) > 2:
                    str2[0] = str2[0]+' '+str2[1]
                    str2[1] = str2[2]
                    str2 = str2[:-1]
                col1 = color_source[str1[0]]
                col2 = color_source[str2[0]]
            colors_list[0].append(str(col1)+", "+str(col2)), colors_list[1].append(comp_names[c])
        if sum(all_comps) > 0:
            cols_df = pd.DataFrame(data=[colors_list[0]], index=['Colors'], columns=colors_list[1])
            this_env = this_env.append(cols_df)
            this_env = this_env.fillna(value=0)
            lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'specific_species']
            this_env.to_csv(basedir+'/ancom/'+envs[a]+'_'+lvls[level]+'_ancom_significant.csv')
            this_env.to_csv(basedir+'/ancom/ancom_significant.csv')
            os.system("/usr/local/bin/Rscript plot_ancom_tree_heatmap.R")
            os.rename(basedir+'/ancom/tree_and_heatmap.pdf', basedir+'/figures/ancom/'+envs[a]+'_ancom_significance_'+lvls[level]+'.pdf')
    make_colorbar_fig(basedir)
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
                env_comp_ft.to_csv(basedir+'/metacoder/metacoder.csv') #save the .csv file with the information for running this comparison in metacoder
                os.system(command) #run metacoder in R (through the terminal)
                os.rename(basedir+'/metacoder/metacoder.pdf', basedir+'/figures/metacoder/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.pdf') #rename the first file, that doesn't have labels on the metacoder plot
                os.rename(basedir+'/metacoder/metacoder_labels.pdf', basedir+'/figures/metacoder/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'_labels.pdf') #and then rename the second one, that does have labels on the metacoders plot
    return

def get_ancom_phylo(ft_abun):
    '''
    Get doc string
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

def group_phylo_plot(ft, tax_dict, plot_on, sn, title, scores, rfs, ASV_dict, abun_fts):
    '''
    This will make the phylogenetic plots with all taxa grouped into the phylogenetic level above (a bit like a stacked bar chart)
    It takes as input:
        - ft (a feature table - this may contain abundance, but ASV names are all that is necessary)
        - tax_dict (a dictionary containing ASV names as keys and taxonomy as values)
        - plot_on (a list of lists, each containing the important features at one phylogenetic level)
        - sn (name to save the resulting figure under)
        - title (plot title)
        - scores (the classification accuracy of each phylogenetic level)
        - rfs (a list of the results of each random forest containing one column with taxa and one with feature importance values - this will be used to colour the cells)
        - ASV_dict (a dictionary containing ASV names as keys and ASV numbers are returned as values (rather than a random collection of numbers and letters))
    Note that it is assumed that the plot will start at the kingdom level, but plot_on, scores and rfs will start at the phylum level
    '''
    asvs = list(ft.index.values) #get a list of the ASVs in the feature table
    asv_tax, names = [], []
    rng = len(plot_on) #we will only go to the level that the data in plot_on goes to
    for a in range(len(asvs)): #for each of these ASVs
        this_tax = tax_dict[asvs[a]] #get the full taxonomy
        last = ''
        new_tax = []
        for b in range(len(this_tax)-1): #for each of the taxonomic levels, excluding the last value that is just the lowest level available
            if this_tax[b] != '': #if there is a value
                new_tax.append(this_tax[b]) #append it to the new list
                last = this_tax[b] #and change the last value to be it
            else:
                new_tax.append(last) #otherwise, take the last value that existed (i.e. the lowest classification level for this ASV)
        new_tax.append(asvs[a]) #and add the ASV name
        asv_tax.append(new_tax[:rng+1]) #append down to the level that plot_on goes to
    for a in range(rng+1): #for each of the levels that plot_on goes to
        names.append(phylo_level_names[a]) #get the name of this taxonomic level
    tax_df = pd.DataFrame(asv_tax, columns=names) #turn the taxa only into a dataframe
    tax_df = tax_df.sort_values(names, ascending=True) #sort them by name - this automatically does kingdom first, followed by all of the taxonomic levels present in plot_on
    new_cols, col_names = [], [] #set up some lists
    already_added = [[], [], [], [], [], [], [], []]
    for a in range(rng): #for all of the levels that plot_on goes to
        a = rng-a-1 #reverse the direction, so we add the lowest taxonomic level first
        this_col = []
        col_names.append(names[a+1]+'_adding') #and make a new name, to differentiate this from the name of the taxonomic level
        for b in range(len(asvs)): #for each of the ASVs present
            if tax_df[names[a+1]].iloc[b] in plot_on[a]: #see if the name of this level is in the list of important features
                if tax_df[names[a+1]].iloc[b] not in already_added[a+1]: #if it is, check if it has already been added or not (this ensure that we do not have multiple values for e.g. proteobacteria that are not important features at any lower classification level)
                    this_col.append(1) #if it's not already been added, add a 1 to the respective new column
                    this_row = tax_df.iloc[[b]].values[0] #get the whole row
                    for c in range(len(this_row)): #and add the names from this row to the already added lists
                        already_added[c].append(this_row[c])
                else: #if it wasn't an important feature or was already added, add a 0 to this column
                    this_col.append(0)
            else:
                this_col.append(0)
        new_cols.append(this_col) #add the column to the list of new columns
    for a in range(len(new_cols)): #now add the new columns to the dataframe
        tax_df[col_names[a]] = new_cols[a]
    tax_df['Max'] = tax_df[col_names].max(axis=1) #add a column with the maximum values
    tax_df = tax_df[tax_df['Max'] > 0] #keep only the rows where we want to add one of the taxonomic levels (i.e. get rid of all that aren't important and those already added)
    tax_df.drop('Max', axis=1, inplace=True) #and get rid of the maximum column again
    #set up the figure, scaling the size to the number of taxa present
    fig = plt.figure(figsize=(14, int(len(tax_df.index.values)/5)))
    ttl = fig.suptitle(title, fontsize=fs_title, fontweight='bold')
    ttl.set_position([.5, 0.9])
    ax1 = plt.subplot2grid((int(len(tax_df.index.values)/1.5), 27), (0, 6), rowspan=int(len(tax_df.index.values)/1.5), colspan=10) #the main axes 
    ax2 = plt.subplot2grid((int(len(tax_df.index.values)/1.5), 27), (0, 0), colspan=3) #an axes for the colorbar
    ma = 0
    cmap = 'viridis'
    colormap = mpl.cm.get_cmap(cmap, 256)                
    norm = mpl.colors.Normalize(vmin=0, vmax=0.05)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
    count_all = 1
    names_legend = []
    if os.path.exists(os.getcwd()+'/random_forest/env_phylo/'+title+'_plot_order.csv'):
        order = []
        with open(os.getcwd()+'/random_forest/env_phylo/'+title+'_plot_order.csv', 'rU') as f:
            for row in csv.reader(f):
                if row[0] in tax_df.ASVs.values:
                    order.append(row[0])
        order.reverse()
        cols = list(tax_df.columns)
        tax_df = tax_df.set_index('ASVs')
        tax_df = tax_df.reindex(order)
        tax_df = tax_df.reset_index()
        tax_df = tax_df.reindex(cols, axis=1)
        print(tax_df)
    tax_df.to_csv(os.getcwd()+'/random_forest/env_phylo/'+title+'_tax_df_asv_test.csv')
    
    for a in range(rng+1): #for all of the levels in plot_on, and also for one extra (the kingdom level)
        if a < 4:
            names_legend.append(names[a].capitalize())
        all_names = list(tax_df.loc[:, names[a]]) #get a list of all taxa names
        unique = []
        for b in range(len(all_names)): #and filter it for only the unique names
            if all_names[b] not in unique:
                unique.append(all_names[b])
        unique.reverse() #reverse this list, as we are plotting from the bottom
        bottom = 0
        y = [] 
        for b in range(len(unique)): #now for each of the unique names
            color, edgecolor = 'w', 'k'
            count = tax_df.groupby(names[a]).size()[unique[b]] #count how many times this value occurs in the list
            text_col = 'k'
            if a > 0: #if we are below the kingdom level, check whether this was an important feature
                lis = tax_df.index[tax_df[names[a]] == unique[b]].tolist() #get a list of index values for where this value occurs in the dataframe
                for c in range(len(lis)): #for each of these values
                    if tax_df[names[a]+'_adding'].loc[lis[c]] > 0: #if it is important (i.e. if we added a 1 in the loop above to the corresponding column)
                        color = 'r' #give this cell a color
                        num = float(rfs[a-1].loc[[tax_df[names[a]].loc[lis[c]]]].iloc[0]) #and find the feature importance value from the other dataframe
                        tax_df[names[a]+'_adding'].loc[lis[c]] = float(rfs[a-1].loc[[tax_df[names[a]].loc[lis[c]]]].iloc[0])
                        color = m.to_rgba(num) #get the corresponding value from the above colormap
                        tax_df[names[a]+'_adding'].loc[lis[c]] = mpl.colors.to_hex(color)
                        if num > ma: #if the number if above the maximum already seen (defined as 0 initially)
                            ma = num #then replace the maximum already seen with this number
                        if num < 0.025:
                            text_col = 'w'
                    else: #otherwise (it wasn't an important feature)
                        if names[a-1] in rename_plots: #if we need to rename the previous name
                            this_name = rename_plots[names[a-1]] #rename from the dictionary
                        else:
                            this_name = names[a-1] #otherwise just get the previous name
                        n = this_name.capitalize()+': ' #capitalize the first letter
                        if ':' in tax_df[names[a-1]].loc[lis[c]]: #if this name was already from the previous taxonomic level
                            tax_df[names[a]].loc[lis[c]] = tax_df[names[a-1]].loc[lis[c]] #then keep it as it is
                        else: 
                            tax_df[names[a]].loc[lis[c]] = n+tax_df[names[a-1]].loc[lis[c]] #otherwise, change the name in this cell to indicate the name of the last level that was an important feature (this is for plotting the y ticks and giving names of all important features)
            if a < 4:
                if count >= 10:
                    ax1.text(a, bottom+(count/2), str(count_all)+': '+unique[b], rotation=90, ha='center', va='center', color=text_col, fontsize=fs_main)
                    names_legend.append(str(count_all)+': '+unique[b])
                else:
                    names_legend.append(str(count_all)+': '+unique[b])
                    ax1.text(a, bottom+(count/2), str(count_all), ha='center', va='center', color=text_col, fontsize=fs_main)
                count_all += 1
            ax1.bar(a, count, bottom=bottom, color=color, edgecolor=edgecolor, width=1) #make the plot for this taxon
            bottom += count #change the bottom value, so we will stack on top of it next time
            y.append(b+0.5) #and add the position for the y value
    y.reverse() #reverse the y ticks
    plot_colorbar(ax2, cmap, [0, 0.05], 'Taxon importance', fs=fs_main) #get the colorbar plot - if you want to use it, replace 0.05 with ma
    plt.sca(ax1) #set plt. to be ax1
    plt.xlim([-0.5, rng+0.5]) #change the x and y limits so we have no additional space around the plots
    plt.ylim([0, len(tax_df.index.values)])
    nums = [0, 1, 2, 3, 4, 5, 6, 7]
    yplot = tax_df[names[rng]].values
    for a in range(len(yplot)):
        if yplot[a] in tax_dict and yplot[a] in ASV_dict:
            yplot[a] = ASV_dict[yplot[a]]+': '+tax_dict[yplot[a]][-1]
    plotted = list(tax_df[names[rng]].values)
    for a in range(len(names)): #for each of the names
        try:
            this_score = int(scores[names[a]]*100) #try to get the classification accuracy score (and multiply by 100)
        except:
            this_score = -200 #otherwise, set the score as -200 (this should only be for the kingdom level, where there is no score)
        if names[a] in rename_plots: #if we need to, rename the name
            names[a] = rename_plots[names[a]]
        if this_score > -200: #if we had a score
            names[a] = names[a] +'\n'+str(this_score)+'%' #add the score to the name so it will be plotted as an x tick label
    plt.xticks(nums[:rng+1], names, rotation=90) #plot the x ticks
    plotted.reverse()
    list_plotting = []
    count = 0
    ancom_sigs = []
    for a in range(len(abun_fts)):
        drop = []
        these_taxa = list(abun_fts[a].index.values)
        for b in range(len(these_taxa)):
            if these_taxa[b] not in tax_df.iloc[:, a+2].values:
                drop.append(these_taxa[b])
        ancom_df = abun_fts[a].drop(drop, axis=0)
        if title == 'Incubation time (specific)':
            ancom_df.drop('unknown', axis=1, inplace=True)
        if len(ancom_df.index.values) > 0:
            sigs = get_ancom_phylo(ancom_df)
            ancom_sigs.append(sigs)
        else:
            ancom_sigs.append([])
        abun_fts[a] = abun_fts[a].groupby(by=abun_fts[a].columns, axis=1).mean()
    
    for a in range(len(plotted)):
        for b in range(len(abun_fts)):
            if plotted[a] in abun_fts[b].index.values:
                list_plotting.append(abun_fts[b].loc[plotted[a], :])
                if count == 0:
                    df_plotting = abun_fts[b].loc[plotted[a], :].to_frame().transpose()
                    count += 1
                else:
                    new_row = abun_fts[b].loc[plotted[a], :].to_frame().transpose()
                    df_plotting = pd.concat([df_plotting, new_row])
    df_plotting = df_plotting.rename(columns=rename_plots)
    if 'Blank' in df_plotting.columns:
        df_plotting.drop('Blank', axis=1, inplace=True)
    if 'Other plastic' in df_plotting.columns:
        if 'Unknown plastic' in df_plotting.columns:
            df_plotting = df_plotting.reindex(['Aliphatic', 'Other plastic', 'Unknown plastic', 'Biofilm', 'Planktonic'], axis=1)
        else:
            df_plotting = df_plotting.reindex(['Aliphatic', 'Other plastic', 'Biofilm', 'Planktonic'], axis=1)
    elif 'Unknown plastic' in df_plotting.columns:
        df_plotting = df_plotting.reindex(['Unknown plastic', 'Biofilm'], axis=1)
    mc = len(df_plotting.columns)
    mc = mc*2
    if mc > 5:
        mc = 11
    elif mc < 4:
        mc = 4
    ax3 = plt.subplot2grid((int(len(tax_df.index.values)/1.5), 27), (0, 16), rowspan=int(len(tax_df.index.values)/1.5), colspan=mc)
    ax4 = plt.subplot2grid((int(len(tax_df.index.values)/1.5), 27), (0, 3), colspan=3) #an axes for the colorbar
    df_plotting = df_plotting.div(df_plotting.max(axis=1), axis=0)
    annotate_heatmap(ax3, df_plotting, cmap='inferno', yticks=False, xticks=True, annotate=False)
    plot_colorbar(ax4, 'inferno', df_plotting, 'Normalized\ntaxon abundance', fs=fs_main)
    plt.sca(ax3)
    plotted.reverse()
    #plt.yticks(y, plotted) #add the y ticks (using the important features for which we changed names above)
    ax5 = plt.subplot2grid((int(len(tax_df.index.values)/1.5), 27), (3, 0), colspan=6, rowspan=int(len(names_legend))-15, frameon=False)
    names_legend.reverse()
    for a in range(len(names_legend)):
        if names_legend[a][0] in ['K', 'P', 'C', 'O']:
            ax5.text(1, a+1, names_legend[a], fontsize=fs_main, fontweight='bold')
        else:
            ax5.text(1, a+1, names_legend[a], fontsize=fs_main)
    ax5.set_xlim([0.95, 3])
    ax5.set_ylim(0, len(names_legend)+1)
    plt.sca(ax5)
    plt.xticks([]), plt.yticks([])
    x = len(df_plotting.columns)+0.15
    empty = []
    for a in range(len(plotted)):
        color = 'k'
        for b in range(len(ancom_sigs)):
            if plotted[a] in ancom_sigs[b]:
                color = '#A70341'
        ax3.text(x, y[a], plotted[a], color=color, fontsize=fs_main, ha='left', va='center')
        empty.append(' ')
    #tax_df.to_csv(os.getcwd()+'/random_forest/env_phylo/'+title+'_tax_df.csv')
    #df_plotting.to_csv(os.getcwd()+'/random_forest/env_phylo/'+title+'_df_plotting.csv')
    plt.sca(ax3)
    plt.yticks(y, empty)
    ax3.yaxis.tick_right() #change the y ticks to be on the right hand side (i.e. with the lowest taxonomic labels)
    plt.sca(ax1)
    plt.yticks([])
    #plt.yticks(y, plotted) #add the y ticks (using the important features for which we changed names above)
    

    #ax1.set_title(title) #set the title
    plt.subplots_adjust(wspace=1)
    plt.savefig(sn+'_0.005'+ext, dpi=dpi, bbox_inches='tight') #save the figure using the savename, file extension and dpi options already defined
    plt.close()
    return

def group_phylo_plot_level(this_ft, level, tax_dict, ASV_dict):
    asv = list(this_ft.index.values)
    rename = {}
    for a in range(len(asv)):
        tax = tax_dict[asv[a]]
        if tax[level] == '':
            using = tax[-1]
        else:
            using = tax[level]
        if level == 7:
            rename[asv[a]] = ASV_dict[asv[a]]+': '+using
        else:
            rename[asv[a]] = rename_plots[phylo_level_names[level]]+': '+using
    new_ft = this_ft.rename(index=rename)
    new_ft = new_ft.groupby(by=new_ft.index, axis=0).sum()
    #sig_ancom = get_ancom_phylo(new_ft)
    return new_ft#, sig_ancom

def myround(x, base=5):
    return base * round(x/base)

def make_phylo_plot(ft, tax_dict, basedir, ASV_dict, meta_dict, env=False, def_sn=False, rf_fol = False):
    ''' 
    This function will make a phylogenetic plot for each metadata category for which you carried out random forest models.
    The plots are like stacked bar plots, but show all phylogenetic levels, with each lower level at the same height as the lower levels.
    Taxa will be coloured only if the feature was above 0.0001 (i.e. the proportion by which the accuracy of that model decreased without the feature)
    It takes as input:
        - ft (a feature table - this may contain abundance, but ASV names are all that is necessary)
        - tax_dict (a dictionary containing ASV names as keys and taxonomy as values)
        - basedir (file path to the base directory that you are working from)
        - ASV_dict (a dictionary containing ASV names as keys and ASV numbers are returned as values (rather than a random collection of numbers and letters))
    '''
    levels = []
    for a in [1, 2, 3, 4, 5, 6, 7]: #for each of the phylogenetic levels we used for random forests
        #get the scores and feature importances for this level
        if rf_fol == False:
            rf_level = pd.read_csv(basedir+'/random_forest/'+phylo_level_names[a]+'_overall.csv', header=0, index_col=0)
        else:
            rf_level = pd.read_csv(rf_fol+phylo_level_names[a]+'.csv', header=0, index_col=0)
        scores = rf_level.loc[['Score']]
        scores = scores.rename(index={'Score':phylo_level_names[a]})
        if a == 1: #if we don't have a dataframe of scores yet
            all_scores = scores #make one with these scores
        else:
            all_scores = pd.concat([all_scores, scores]) #otherwise, add them to the others
        rf_level.drop(['Score', 'OOB_score'], axis=0, inplace=True) #now get rid of the scores from the main dataframe
        levels.append(rf_level) #and add this level to the list of dataframes
    names = list(levels[0].columns) #get the metadata category names
    ft = ft*100
    for a in range(len(names)): #for each of the metadata categories
        if rf_fol != False:
            if names[a] != 'PlasticTypeGeneral':
                continue
        ft_abun = ft.copy()
        samplenames = list(ft_abun.columns)
        rename = {}
        for b in range(len(samplenames)):
            rename[samplenames[b]] = meta_dict[samplenames[b]][meta_name_ind[names[a]]]
        ft_abun.rename(columns=rename, inplace=True)
        colnames = list(ft_abun.columns)
        numeric = {}
        for b in range(len(colnames)): #for each grouping name
            try: #if it is numeric
                num = float(colnames[b]) #get it as a float
                num = myround(num) #and then round it to the nearest whole number (the random forest regression was still based on the full values, but we do this to simplify plotting)
                numeric[colnames[b]] = num
            except:
                do_nothing=True
        ft_abun.rename(columns=numeric, inplace=True)
        #abun_fts, ancom_sigs = [], []
        abun_fts = []
        for c in [1, 2, 3, 4, 5, 6, 7]:
            level_ft = group_phylo_plot_level(ft_abun, c, tax_dict, ASV_dict)
            #level_ft, sig_ancom = group_phylo_plot_level(ft_abun, c, tax_dict, ASV_dict)
            #level_ft = level_ft.groupby(by=level_ft.columns, axis=1).mean()
            abun_fts.append(level_ft)
            #ancom_sigs.append(sig_ancom)
        on = list(names) #get a list of the names
        on.remove(names[a])  #without the current metadata category names
        names_to_add = [] #set up a list for the important features at each phylogenetic level
        rfs = [] #and for the individual random forests for each level
        for b in range(len(levels)): #for all levels
            this_set = levels[b].copy() #copy the dataframe for this level
            this_set = this_set.drop(on, axis=1, inplace=False) #drop the other metadata categories
            this_set = this_set[this_set.max(axis=1) > 0.005] #keep only the features that are above 0.0001 importance
            names_to_add.append(list(this_set.index.values)) #add all of the remaining feature names to the list
            rfs.append(this_set) #and add the dataframe to the list
        if def_sn == False:
            sn = basedir+'/figures/random_forest/phylo_plot_'+names[a] #set up the save names
        else:
            sn = def_sn+names[a]
        
        if names[a] in rename_plots: #change the name of the metadata category if necessary (this just removes spaces and capitalizes, etc)
            title = rename_plots[names[a]]
        else:
            title = names[a]
        if env != False:
            title = rename_plots[env]
        group_phylo_plot(ft, tax_dict, names_to_add, sn, title, all_scores[names[a]], rfs, ASV_dict, abun_fts) #get the plot
        #group_phylo_plot(ft, tax_dict, names_to_add, sn, title, all_scores[names[a]], rfs, ASV_dict, abun_fts, ancom_sigs) #get the plot
    return

def make_env_phylo_plot(ft, tax_dict, basedir, ASV_dict, meta_dict):
    envs, ei = ['marine'], 3#, 'aquatic', 'freshwater', 'terrestrial'], 3
    for a in range(len(envs)):
        samples = list(ft.columns)
        keeping = []
        count = 0
        for b in range(len(samples)):
            if meta_dict[samples[b]][ei] == envs[a]:
                keeping.append(True)
                count += 1 
            else:
                keeping.append(False)
        env_ft = ft.loc[:, keeping]
        env_ft = env_ft[env_ft.max(axis=1) > 0]
        rf_fol = basedir+'/random_forest/single_environment/'+envs[a]+'_'
        sn = basedir+'/figures/random_forest/'+envs[a]+'_phylo_plot_'
        make_phylo_plot(env_ft, tax_dict, basedir, ASV_dict, meta_dict, env=envs[a], def_sn=sn, rf_fol=rf_fol)
    return

def plot_individual_study(study_name, study_df, study_meta_df, study_w_uf, study_uw_uf, meta_dict, tax_dict, ASV_dict, basedir, est=10000, n_jobs=1):
    '''
    Get a figure for an individual study, including weighted and unweighted unifrac nMDS plots 
    and heatmaps to show the classification accuracy of random forest models as well as the top ASVs
    Input:
        - study_name (the name of the study that we are plotting, a string)
        - study_df (a feature table dataframe containing only samples for this study, with samples as columns and ASVs as rows)
        - study_meta_df (a dataframe containing all metadata for the samples in this study, with samples as rows and categories as columns. It should only contain the categories (columns) that you want to carry out random forest models for)
        - study_w_uf (dataframe containing the weighted unifrac distance matrix with only samples for this study)
        - study_uw_uf (dataframe containing the unweighted unifrac distance matrix with only samples for this study)
        - meta_dict (a dictionary with sample names as keys and all metadata information returned)
        - tax_dict (a dictionary containing ASVs as keys that returns all other taxonomic levels)
        - basedir (path to the base directory for analyses as a string)
        - est (number of estimators to use for random forests; default is 10,000)
        - n_jobs (number of threads to use for calculating random forests; default is 1)
    '''
    this_fol = basedir+'/separate_studies/'+study_name+'/' #define the folder name for this study
    if os.path.exists(this_fol): #if the folder already exists, do nothing
        do_nothing = True
    else: #otherwise, make the folder
        os.system("mkdir "+this_fol)
    #set up the figure and subplots, using 
    fig = plt.figure(figsize=(20,10))
    fig.suptitle(rename_plots[study_name], fontsize=fs_title+4, fontweight='bold')
    ax_nmds_w = plt.subplot2grid((2,3), (0,0))
    ax_nmds_uw = plt.subplot2grid((2,3), (1,0))
    ax_nmds_w.set_title('Weighted Unifrac', fontsize=fs_title, fontweight='bold')
    ax_nmds_uw.set_title('Unweighted Unifrac', fontsize=fs_title, fontweight='bold')
    uf, ax_uf = [study_w_uf, study_uw_uf], [ax_nmds_w, ax_nmds_uw]
    filter_on, filter_index, second_filter, second_filter_ind = 'none', 'none', '', ''
    source, source_index, color_src = ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic', 'blank'], 10, ['#5F9EA0', '#8B008B', '#3593FC', '#FCA94A', 'yellow', '#C3C3C3']
    scores, dfs = get_single_forest(study_df.transpose(), study_meta_df, this_fol+'/random_forest', 7, est=est, n_jobs=n_jobs)
    for a in range(len(uf)):
        names = list(uf[a].columns) #get a list of column names
        dist_matr = uf[a].rename_axis('ID').values
        npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, source, source_index, ax_uf[a], 'best', color_src, names, meta_dict, second_filter, second_filter_ind, '', n_jobs=1)
        plt.sca(ax_uf[a])
        plt.legend(handles=handles, loc='best', fontsize=fs_main)
        plt.xlabel('nMDS1', fontsize=fs_small), plt.ylabel('nMDS2', fontsize=fs_small)
        plt.yticks(fontsize=fs_small), plt.xticks(fontsize=fs_small)
    if len(list(study_meta_df.columns)) == 0 or (isinstance(scores, list) and scores == []) or np.isnan(scores.loc['Score', :].max()) or scores.loc['Score', :].max() <= 0:
        cs = len(study_df.columns)
        if cs <= 8:
            cs = 16
        elif cs > 30:
            cs = 34
        ax1 = plt.subplot2grid((20,66), (0,30), rowspan=20, colspan=cs)
        ax_colbar1 = plt.subplot2grid((3,66), (0,30+cs+1))
        asv_tax = {}
        asv = list(study_df.index.values)
        for a in range(len(asv)):
            asv_tax[asv[a]] = ASV_dict[asv[a]]+': '+r'$'+tax_dict[asv[a]][-1]+'$'
        study_df.rename(index=asv_tax, inplace=True)
        study_df = study_df*100
        study_df["Mean"] = study_df.mean(axis=1)
        study_df.sort_values(by=['Mean'], ascending=False, inplace=True)
        study_df = study_df[:30]
        study_df = study_df.drop(['Mean'], axis=1, inplace=False)
        study_df = study_df.reindex(sorted(study_df.columns), axis=1)
        if cs > 20:
            rename = {}
            nums = []
            for a in range(len(study_df.columns)):
                rename[study_df.columns[a]] = a+1
                if a+1 % 10 == 0:
                    nums.append(a+1)
            study_df.rename(columns=rename, inplace=True)
            annotate_heatmap(ax1, study_df, cmap='viridis', yticks=True, xticks=False, rnd=3, annotate=False, italics=True)
            plt.sca(ax1)
            plt.xticks(nums)
        else:
            annotate_heatmap(ax1, study_df, cmap='viridis', yticks=True, xticks=True, rnd=3, annotate=True, italics=True)
        plot_colorbar(ax_colbar1, 'viridis', [0, max(study_df.max())], 'ASV relative abundance (%)', fs=fs_main, orientation='vertical')
        string = r'$\bf{Relative}$ $\bf{abundance}$ $\bf{of}$ $\bf{top}$ $\bf{30}$ $\bf{ASVs}$'
        ax1.set_title(string, fontsize=fs_title)
        plt.savefig(this_fol+'overall_plot'+ext, dpi=dpi, bbox_inches='tight')
        plt.subplots_adjust(hspace=0.2, wspace=0.2)
        plt.close()
        return
    ax_rf_scores = plt.subplot2grid((20,33), (0,13), colspan=19)
    ax_mean_importance = plt.subplot2grid((20,33), (2,13), rowspan=18)
    ax_rf_1_importance = plt.subplot2grid((20,33), (2,14), rowspan=18)
    ax_rf_1 = plt.subplot2grid((20,33), (2,15), rowspan=18, colspan=5)
    ax_rf_2_importance = plt.subplot2grid((20,33), (2,20), rowspan=18)
    ax_rf_2 = plt.subplot2grid((20,33), (2,21), rowspan=18, colspan=5)
    ax_rf_3_importance = plt.subplot2grid((20,33), (2,26), rowspan=18)
    ax_rf_3 = plt.subplot2grid((20,33), (2,27), rowspan=18, colspan=5)
    ax_colbar1 = plt.subplot2grid((3,66), (0,65))
    ax_colbar2 = plt.subplot2grid((3,66), (1,65))
    ax_colbar3 = plt.subplot2grid((3,66), (2,65))
    ax1, ax2, ax3 = [ax_rf_1_importance, ax_rf_1], [ax_rf_2_importance, ax_rf_2], [ax_rf_3_importance, ax_rf_3]
    axes = [ax1, ax2, ax3]
    scores.drop('OOB_score', axis=0, inplace=True)
    dfs.drop(['Score', 'OOB_score'], axis=0, inplace=True)
    scores = scores*100
    scores_rename = scores.rename(columns=rename_plots)
    ax_rf_scores.set_title('Random forest classification accuracy (%)', fontsize=fs_title, fontweight='bold')
    annotate_heatmap(ax_rf_scores, scores_rename, cmap='PuBu', yticks=False, xticks='top', rnd=1, annotate=True, italics=False)
    meta_cats = list(scores.columns)
    scores_list = []
    for cat in meta_cats:
        scores_list.append(scores.loc['Score', cat])
    plot_colorbar(ax_colbar1, 'PuBu', [0, max(scores_list)], 'Classification accuracy (%)', fs=fs_main, orientation='vertical')
    scores_list = [x for _,x in sorted(zip(scores_list,meta_cats))]
    scores_list.reverse()
    scores_list, drop = scores_list[:3], scores_list[3:]
    count, count1 = 0, 0
    for a in range(len(scores_list)):
        if 'incubation' in scores_list[a] or 'Incubation' in scores_list[a]:
            if count > 0:
                rem = str(scores_list[a])
                scores_list[a] = drop[0]
                drop.remove(drop[0])
                drop.append(rem)
            else:
                count += 1
        if 'type' in scores_list[a] or 'Type' in scores_list[a]:
            if count1 > 0:
                rem = str(scores_list[a])
                scores_list[a] = drop[0]
                drop.remove(drop[0])
                drop.append(rem)
            else:
                count1 += 1
    dfs.drop(drop, axis=1, inplace=True)
    dfs["Mean"] = dfs.mean(axis=1)
    dfs.sort_values(by=['Mean'], ascending=False, inplace=True)
    dfs = dfs[:30]
    df_mean = dfs.drop(scores_list, axis=1, inplace=False)
    dfs.drop(['Mean'], axis=1, inplace=True)
    asv_tax = {}
    asv = list(dfs.index.values)
    for a in range(len(asv)):
        asv_tax[asv[a]] = ASV_dict[asv[a]]+': '+r'$'+tax_dict[asv[a]][-1]+'$'
    df_mean = df_mean.reset_index()
    df_mean.rename(columns={'index':'ASV'}, inplace=True)
    df_mean = df_mean.set_index('ASV')
    study_df = study_df.merge(right=df_mean, on='ASV')
    study_df.sort_values(by=['Mean'], ascending=False, inplace=True)
    study_df.drop('Mean', axis=1, inplace=True)
    df_mean.rename(columns=rename_plots, index=asv_tax, inplace=True)
    annotate_heatmap(ax_mean_importance, df_mean, cmap='inferno', yticks=True, xticks=True, rnd=3, annotate=True, italics=True)
    plot_colorbar(ax_colbar2, 'inferno', [0, df_mean.iloc[0,0]], 'ASV importance', fs=fs_main, orientation='vertical')
    mas = []
    this_dfs = []
    for a in range(len(scores_list)):
        ax = axes[a]
        this_df = dfs.copy()
        mn = list(this_df.columns)
        mn.remove(scores_list[a])
        this_df.drop(mn, axis=1, inplace=True)
        this_df.rename(columns={scores_list[a]:'Importance'}, inplace=True)
        annotate_heatmap(ax[0], this_df, cmap='inferno', yticks=False, xticks=True, rnd=3, annotate=True)
        rename_cols = {}
        samples = list(study_df.columns)
        for s in range(len(samples)):
            group = study_meta_df.loc[samples[s], scores_list[a]]
            rename_cols[samples[s]] = group
        this_study_df = study_df.rename(columns=rename_cols, inplace=False)
        this_study_df.rename(columns=rename_plots, inplace=True)
        this_study_df = this_study_df.groupby(by=this_study_df.columns, axis=1).mean()
        this_study_df = this_study_df*100
        mas.append(max(this_study_df.max()))
        if scores_list[a] in rename_plots:
            scores_list[a] = rename_plots[scores_list[a]]
        if a == 1:
            string = r'$\bf{Top}$ $\bf{3}$ $\bf{random}$ $\bf{forest}$ $\bf{models}$ $\bf{and}$ $\bf{important}$ $\bf{ASVs}$'
            ax[1].set_title(string+'\n'+scores_list[a], fontsize=fs_title)
        else:
            ax[1].set_title(scores_list[a], fontsize=fs_title)
        this_dfs.append(this_study_df)
    annotate_heatmap(ax_rf_1, this_dfs[0], cmap='viridis', yticks=False, xticks=True, rnd=1, annotate=True, vmax=max(mas))
    if len(this_dfs) > 1:
        annotate_heatmap(ax_rf_2, this_dfs[1], cmap='viridis', yticks=False, xticks=True, rnd=1, annotate=True, vmax=max(mas))
        if len(this_dfs) > 2:
            annotate_heatmap(ax_rf_3, this_dfs[2], cmap='viridis', yticks=False, xticks=True, rnd=1, annotate=True, vmax=max(mas))
        else:
            ax_rf_3_importance = plt.subplot2grid((20,33), (2,26), rowspan=18, frameon=False)
            ax_rf_3 = plt.subplot2grid((20,33), (2,27), rowspan=18, colspan=5, frameon=False)
            plt.sca(ax_rf_3_importance)
            plt.yticks([])
            plt.xticks([])
            plt.sca(ax_rf_3)
            plt.yticks([])
            plt.xticks([])
    else:
        ax_rf_2_importance = plt.subplot2grid((20,33), (2,20), rowspan=18, frameon=False)
        ax_rf_2 = plt.subplot2grid((20,33), (2,21), rowspan=18, colspan=5, frameon=False)
        plt.sca(ax_rf_2_importance)
        plt.yticks([])
        plt.xticks([])
        plt.sca(ax_rf_2)
        plt.yticks([])
        plt.xticks([])
        ax_rf_3_importance = plt.subplot2grid((20,33), (2,26), rowspan=18, frameon=False)
        ax_rf_3 = plt.subplot2grid((20,33), (2,27), rowspan=18, colspan=5, frameon=False)
        plt.sca(ax_rf_3_importance)
        plt.yticks([])
        plt.xticks([])
        plt.sca(ax_rf_3)
        plt.yticks([])
        plt.xticks([])
    plot_colorbar(ax_colbar3, 'viridis', [0, max(mas)], 'ASV relative abundance (%)', fs=fs_main, orientation='vertical')
    plt.subplots_adjust(hspace=0.4, wspace=0.5)
    plt.savefig(this_fol+'overall_plot'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def plots_per_study(ft_df, meta_df, meta_dict, w_uf, uw_uf, tax_dict, ASV_dict, basedir, est=10000, n_jobs=1):
    '''
    Get figures for each study, showing weighted and unweighted unifrac distances between samples, random forest models based on ASVs and
    plots of the top 3 random forest models and 30 most informative ASVs
    Input:
        - ft_df (a dataframe of a feature table with all samples across all studies, with ASVs as rows and samples as columns)
        - meta_df (a dataframe with samples as rows and metadata categories as columns)
        - meta_dict (a dictionary with sample names as keys and returning all sample metadata)
        - w_uf (the file name of a weighted unifrac distance matrix for all samples)
        - uw_uf (the file name of an unweighted unifrac distance matrix for all samples)
        - tax_dict (a dictionary containing ASVs as keys that returns all other taxonomic levels)
        - basedir (path to the base directory for analyses as a string)
        - est (number of estimators to use for random forests; default is 10,000)
        - n_jobs (number of threads to use for calculating random forests; default is 1)
    '''
    studies = sorted(list(set(meta_df['Study']))) #get the unique study names from the metadata file
    samples = list(ft_df.columns) #get a list of all sample names
    samples_meta = list(meta_df.index.values) #get a list of all samples names in the metadata file
    w_uf = pd.read_csv(w_uf, header=0, index_col=0) #read in the weighted unifrac distance matrix
    uw_uf = pd.read_csv(uw_uf, header=0, index_col=0) #read in the unweighted unifrac distance matrix
    cols_uf = list(w_uf.columns) #get a list of sample names from the distance matrix
    ft_df = ft_df.reset_index() #rest the index column
    ft_df.rename(columns={'index':'ASV'}, inplace=True) #rename the previous index column
    ft_df = ft_df.set_index('ASV') #now reset the original index column, now named 'ASV'
    adding = False
    for a in range(len(studies)): #for each study
        if a > 23:
            continue
        print(studies[a])
        keeping, keeping_meta, keeping_uf = [], [], [] #set up lists
        #check each of the feature table, metadata dataframe and unifrac sample names to see whether they are in this study. 
        #if they are, add True to the appropriate list, if not, add False
        for b in range(len(samples)): 
            if meta_dict[samples[b]][0] == studies[a]:
                keeping.append(True)
            else:
                keeping.append(False)
        for c in range(len(samples_meta)):
            if meta_dict[samples_meta[c]][0] == studies[a]:
                keeping_meta.append(True)
            else:
                keeping_meta.append(False)
        for d in range(len(cols_uf)):
            if meta_dict[cols_uf[d]][0] == studies[a]:
                keeping_uf.append(True)
            else:
                keeping_uf.append(False)
        study_df = ft_df.loc[:, keeping] #now get only those that we are keeping as a new dataframe for the feature table
        study_df = study_df[study_df.max(axis=1) > 0] #and only keep ASVs that are present in this study
        study_w_uf = w_uf.loc[keeping_uf, keeping_uf] #now get only those we are keeping as a new dataframe for the weighted unifrac
        study_uw_uf = uw_uf.loc[keeping_uf, keeping_uf] #now get only those we are keeping as a new dataframe for the unweighted unifrac
        study_meta_df = meta_df.loc[keeping_meta, :] #now get only those we are keeping as a new dataframe for the metadata dataframe
        meta_cols = list(study_meta_df.columns) #get a list of the metadata categories
        not_adding = []
        for col in meta_cols: #for each of those columns
            vals = list(study_meta_df[col]) #get a list of unique values
            if len(list(set(vals))) < 2: #if there are fewer than two unique values in the category, then don't use this for random forests
                not_adding.append(col) #so add it to not adding
        study_meta_df.drop(not_adding, axis=1, inplace=True) #and drop those columns that we don't want from the dataframe
        plot_individual_study(studies[a], study_df, study_meta_df, study_w_uf, study_uw_uf, meta_dict, tax_dict, ASV_dict, basedir, est=est, n_jobs=n_jobs) #now get the plot for this study
    return

def get_picrust_individual_plot(this_env, axes):
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
    top = this_env.sort_values(['FC'], ascending=False)
    bottom = this_env.sort_values(['FC'], ascending=True)
    top = top[:20]
    bottom = bottom[:20]
    plotting = pd.concat([top, bottom])
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
        #axes[2].bar([0], [1], bottom=[a], color=col1, width=1, edgecolor='k')
        #axes[2].bar([1], [1], bottom=[a], color=col2, width=1, edgecolor='k')
        #axes[2].text(0, a+0.5, int(plotting['plastic mean'].iloc[a]), ha='center', va='center')
        #axes[2].text(1, a+0.5, int(plotting['biofilm mean'].iloc[a]), ha='center', va='center')
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
            """
            name, count = '', 0
            for c in range(len(new_name)):
                if count > 35:
                    if new_name[c] == ' ' or new_name[c] == ';' or new_name[c] == '-' or new_name[c] == '_':
                        name += '\n'
                        count = 0
                    else:
                        name += new_name[c]
                        count += 1
                        continue
                name += new_name[c]
                count += 1
            """
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
    #plt.sca(axes[2])
    #plt.yticks([])
    #plt.xticks([0, 1], ['Plastic', 'Biofilm'], rotation=90)
    #plt.ylim([0, 40])
    axes[1].set_xlabel(r'log$_{2}$ mean fold change'+'\n'+r'Control biofilms $vs$ Plastic samples', fontsize=fs_main)
    plt.xticks(fontsize=fs_main)
    #axes[2].set_title('Abundance')
    return names

def picrust_plots(picrust, KO_dict, meta_dict, basedir):
    envs, env_index = ['marine', 'aquatic', 'freshwater', 'terrestrial'], 3
    source, source_index = [['aliphatic', 'other plastic', 'unknown plastic'], ['biofilm']], 10
    samples = list(picrust.columns)
    fig = plt.figure(figsize=(20,20))
    ax1_volc = plt.subplot2grid((4, 3), (0,0))
    ax2_volc = plt.subplot2grid((4, 3), (1,0))
    ax3_volc = plt.subplot2grid((4, 3), (2,0))
    ax4_volc = plt.subplot2grid((4, 3), (3,0))
    ax1_bar = plt.subplot2grid((4, 6), (0,2))
    ax2_bar = plt.subplot2grid((4, 6), (1,2))
    ax3_bar = plt.subplot2grid((4, 6), (2,2))
    ax4_bar = plt.subplot2grid((4, 6), (3,2))
    """
    fig = plt.figure(figsize=(20, 15))
    ax1_volc = plt.subplot2grid((4,16), (0,0), colspan=2)
    ax2_volc = plt.subplot2grid((4,16), (0,4), colspan=2)
    ax3_volc = plt.subplot2grid((4,16), (0,8), colspan=2)
    ax4_volc = plt.subplot2grid((4,16), (0,12), colspan=2)
    ax1_bar = plt.subplot2grid((4,16), (1,0), rowspan=3, colspan=2)
    ax2_bar = plt.subplot2grid((4,16), (1,4), rowspan=3, colspan=2)
    ax3_bar = plt.subplot2grid((4,16), (1,8), rowspan=3, colspan=2)
    ax4_bar = plt.subplot2grid((4,16), (1,12), rowspan=3, colspan=2)
    #ax1_heat = plt.subplot2grid((4,16), (1,0), rowspan=3, colspan=1)
    #ax2_heat = plt.subplot2grid((4,16), (1,4), rowspan=3, colspan=1)
    #ax3_heat = plt.subplot2grid((4,16), (1,8), rowspan=3, colspan=1)
    #ax4_heat = plt.subplot2grid((4,16), (1,12), rowspan=3, colspan=1)
    """
    all_ax = [ax1_volc, ax1_bar, ax2_volc, ax2_bar, ax3_volc, ax3_bar, ax4_volc, ax4_bar]
    axes = [[ax1_volc, ax1_bar], [ax2_volc, ax2_bar], [ax3_volc, ax3_bar], [ax4_volc, ax4_bar]]
    titles = ['Marine', 'Marine', 'Aquatic', 'Aquatic', 'Freshwater', 'Freshwater', 'Terrestrial', 'Terrestrial']
    ax1_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    ax2_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    ax3_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    ax4_volc.set_ylabel(r'-log$_{10}$(P value)', fontsize=fs_main)
    for a in range(len(all_ax)):
        all_ax[a].set_title(titles[a])
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
    plastic = mlines.Line2D([], [], color='#B90404', marker='o', markersize=8, markeredgecolor='k', label='Upregulated in\nplastic samples', linestyle=' ')
    biofilm = mlines.Line2D([], [], color='#0498B9', marker='o', markersize=8, markeredgecolor='k', label='Upregulated in\ncontrol biofilms', linestyle=' ')
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





"""""""""
Not using
"""""""""

def add_map_features():
    plt.figure(figsize=(30,4)) #set up figure
    ax1 = plt.subplot2grid((40,9), (4,0), rowspan=32) #add axis for number of publications
    ax2 = plt.subplot2grid((1,9), (0,1), colspan=2, projection=ccrs.PlateCarree()) #add axis for map
    ax1.set_title('A', loc='left', fontsize=fs_title) #add axis label
    ax2.set_title('B', loc='left', fontsize=fs_title) #add axis label
    axins1 = plt.subplot2grid((2,9), (0,3), projection=ccrs.PlateCarree()) #add first inset map axis
    axins2 = plt.subplot2grid((3,13), (0,9), projection=ccrs.PlateCarree()) #add second inset map axis
    
    plt.sca(ax2)
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
    
    
    return

def get_important_forests_environment(ft, tax_dict, ASV_dict, meta_dict, basedir):
    envs, ei, si = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3, 10
    level = 5
    asv = list(ft.index.values)
    tax_rename = {}
    for a in range(len(asv)):
        tax = tax_dict[asv[a]]
        this_level = ''
        if tax[level] == '':
            this_level = tax[-1]
        else:
            this_level = tax[level]
        tax_rename[asv[a]] = this_level
    ft.rename(index=tax_rename, inplace=True)
    ft = ft.groupby(by=ft.index.values, axis=0).sum()
    ft = ft.reset_index() #reset the indices in our feature table
    ft = ft.rename(columns={'index':'ASV'}) #rename the index column
    ft = ft.set_index('ASV')
    ft = ft*100
    for a in range(len(envs)):
        keeping = []
        colnames = list(ft.columns)
        rename = {}
        for b in range(len(colnames)):
            if meta_dict[colnames[b]][ei] == envs[a]:
                if meta_dict[colnames[b]][si] != 'planktonic' and meta_dict[colnames[b]][si] != 'blank':
                    keeping.append(True)
                    rename[colnames[b]] = meta_dict[colnames[b]][si]
                else:
                    keeping.append(False)
            else:
                keeping.append(False)
        single_ft = ft.loc[:, keeping]
        single_ft = single_ft.rename(columns=rename)
        rf = pd.read_csv(basedir+'/random_forest/single_environment/'+envs[a]+'_'+phylo_level_names[level]+'.csv', header=0, index_col=0)
        cols = list(rf.columns)
        cols.remove('PlasticTypeGeneral')
        rf.drop(cols, axis=1, inplace=True), rf.drop(['Score', 'OOB_score'], axis=0, inplace=True)
        rf = rf.sort_values(['PlasticTypeGeneral'], ascending=False)
        """
        #This is to make box plots for the top 10 significant taxa (random forest) if they have a significant ANOVA test
        fig = plt.figure(figsize=(20,10))
        ax1, ax2, ax3, ax4, ax5 = plt.subplot2grid((2,5), (0,0)), plt.subplot2grid((2,5), (0,1)), plt.subplot2grid((2,5), (0,2)), plt.subplot2grid((2,5), (0,3)), plt.subplot2grid((2,5), (0,4))
        ax6, ax7, ax8, ax9, ax10 = plt.subplot2grid((2,5), (1,0)), plt.subplot2grid((2,5), (1,1)), plt.subplot2grid((2,5), (1,2)), plt.subplot2grid((2,5), (1,3)), plt.subplot2grid((2,5), (1,4))
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]
        rf = rf.reset_index() #reset the indices in our feature table
        rf = rf.rename(columns={'index':'ASV'}) #rename the index column
        rf = rf.set_index('ASV')
        single_ft = single_ft.merge(right=rf, on='ASV')
        single_ft = single_ft.sort_values(['PlasticTypeGeneral'], ascending=False)
        single_ft.drop(['PlasticTypeGeneral'], axis=1, inplace=True)
        taxa = list(single_ft.index.values)
        x = []
        count = 0
        for c in range(len(taxa)):
            if count > 9:
                continue
            colnames = list(set(single_ft.columns))
            vals = []
            maxs = []
            for d in range(len(colnames)):
                if c == 0:
                    x.append(d+1)
                new_vals = single_ft.loc[taxa[c], colnames[d]]
                vals.append(new_vals)
                maxs.append(max(new_vals))
            f, p = stats.f_oneway(*[list(vals[f]) for f in range(len(vals))])
            if p < 0.05:
                axes[count].boxplot(vals, positions=x, showfliers=False, showmeans=True)
                axes[count].set_title(r'$'+taxa[c]+'$')
                sigs = []
                adding = 0
                for e in range(len(colnames)-1):
                    this_sig = []
                    for f in range(e+1, len(colnames)):
                        ttest, ttest_p = stats.ttest_ind(vals[e], vals[f])
                        if ttest_p > 0.05:
                            continue
                        if adding == 0:
                            adding = max(maxs)+max(maxs)*0.1
                        else:
                            adding += max(maxs)*0.1
                        print(adding)
                        #axes
                        this_sig.append(ttest_p)
                        #axes[count].text(x[e], adding, '*')
                        
                    sigs.append(this_sig)
                plt.sca(axes[count])
                for f in range(len(colnames)):
                    if colnames[f] in rename_plots:
                        colnames[f] = rename_plots[colnames[f]]
                plt.xticks(x, colnames, rotation=90)
                count += 1
            else:
                do_nothing=True
            
        ax1.set_ylabel('Relative abundance (%)'), ax6.set_ylabel('Relative abundance (%)')
        removex = [ax1, ax2, ax3, ax4, ax5]
        for ax in removex:
            plt.sca(ax)
            plt.xticks([])
        plt.savefig(basedir+'/figures/random_forest_taxa/'+envs[a]+'_'+phylo_level_names[level]+ext, dpi=dpi, bbox_inches='tight')
        plt.close()
        """
    return
