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
from skbio.stats.composition import ancom
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

fs_title, fs_main, fs_small = 10, 8, 6
fig_h, fig_w = 11.69, 8.27
color_env = {'marine':'#03A498', 'freshwater':'#03B1FC', 'aquatic':'#9A75FC', 'terrestrial':'#FD9A64'} #colors for plotting each environment
label_env = ['Marine', 'Freshwater', 'Aquatic', 'Terrestrial'] #labels for each environment
ext, dpi = '.png', 600 #extension and dots per inch to save figures  with
color_source = {'aliphatic':'#5F9EA0', 'biofilm':'#FCA94A', 'other plastic':'#8B008B', 'unknown plastic':'#3593FC', 'planktonic':'#F9E79F', 'blank':'gray'}
"""
Automate the creation of these dictionaries for each study
"""
name_dict = {'La2020':'Latva\n$et$ $al$. 2020', 'AZ2015':'Amaral-Zettler\n$et$ $al$. 2015', 'AA2018':'Arias-Andres\n$et$ $al$. 2018', 'Ca2020':'Canada\n$et$ $al$. 2020', 'Cu2019':'Curren & Leong\n2019', 'De2019':'Delacuvellerie\n$et$ $al$. 2019', 'DT2015':'De Tender\n$et$ $al$. 2015', 'DT2017':'De Tender\n$et$ $al$. 2017', 'DH2018':'Dussud $et$ $al$. \n2018a', 'DM2018':'Dussud $et$ $al$. \n2018b', 'EC2019':'Erni-Cassola\n$et$ $al$. 2019', 'Es2019':'Esan\n$et$ $al$. 2019', 'Fr2018':'Frere\n$et$ $al$. 2018', 'Ho2014':'Hoellein\n$et$ $al$. 2014', 'Ho2017':'Hoellein\n$et$ $al$. 2017', 'Ji2018':'Jiang\n$et$ $al$. 2018', 'Ke2019':'Kesy\n$et$ $al$. 2019', 'Ki2018':'Kirstein\n$et$ $al$. 2018', 'Ki2019':'Kirstein\n$et$ $al$. 2019', 'MC2014':'McCormick\n$et$ $al$. 2014', 'MC2016':'McCormick\n$et$ $al$. 2016', 'Ob2016':'Oberbeckmann\n$et$ $al$. 2016', 'Ob2018':'Oberbeckmann\n$et$ $al$. 2018', 'Og2018':'Ogonowski\n$et$ $al$. 2018', 'Pi2019':'Pinto\n$et$ $al$. 2019', 'Po2018':'Pollet\n$et$ $al$. 2018', 'Ro2020':'Rosato\n$et$ $al$. 2020', 'Sy2019':'Syranidou\n$et$ $al$. 2019', 'SyPE20':'Syranidou\n$et$ $al$. 2017a', 'SyPS20':'Syranidou\n$et$ $al$. 2017b', 'Ta2019':'Tagg\n$et$ $al$. 2019', 'Wo2018':'Woodall\n$et$ $al$. 2018', 'Wr2019':'Wright\n$et$ $al$. 2020', 'Wu2019':'Wu\n$et$ $al$. 2019', 'Xu2019':'Xu\n$et$ $al$. 2019', 'Zh2019':'Zhang\n$et$ $al$. 2019', 'Br2016':'Bryant\n$et$ $al$. 2016', 'Pin201':'Pinnell\n$et$ $al$. 2019'}
name_env = {'AZ2015':'marine', 'AA2018':'freshwater', 'Ca2020':'aquatic', 'Cu2019':'marine', 'De2019':'marine', 'DT2015':'marine', 'DT2017':'marine', 'DH2018':'marine', 'DM2018':'marine', 'EC2019':'marine', 'Es2019':'terrestrial', 'Fr2018':'marine', 'Ho2014':'freshwater', 'Ho2017':'freshwater', 'Ke2019':'aquatic', 'Ki2018':'marine', 'Ki2019':'marine', 'MC2014':'freshwater', 'MC2016':'freshwater', 'Ob2016':'marine', 'Ob2018':'aquatic', 'Og2018':'marine', 'Pi2019':'marine', 'Po2018':'marine', 'Ro2020':'marine', 'Sy2019':'marine', 'SyPE20':'marine', 'SyPS20':'marine', 'Ta2019':'aquatic', 'Wo2018':'marine', 'Wr2019':'marine', 'Wu2019':'freshwater', 'Xu2019':'marine', 'Zh2019':'terrestrial'}

def write_csv(fn, data): #Save a csv file with file name fn and data in a list
    with open(fn, 'w') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)
    return

def write_txt(fn, data): #Save a csv file with file name fn and data in a list
    with open(fn, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        for row in data:
            writer.writerow(row)
    return

def open_csv(fn): #open csv with file name fn, returning the data as rows
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    return rows

def open_txt(fn): #open text file with file name fn, returning the data as rows
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f, delimiter='\t', lineterminator='\n'):
            rows.append(row)
    return rows

def write_pickle(fn, data): #write pickle python object with name fn and the data in data
    with open(fn, 'wb') as f:
        pickle.dump(data, f)
    return

def format_R(ft): #format the raw data files for phyloseq and unifrac analyses in R
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
    write_csv('taxonomy.csv', tax_file) #write the taxonomy to .csv
    write_csv('taxonomy_name_only.csv', tax_name_only) #write the taxonomy to .csv
    write_pickle('tax_dict.dictionary', tax_dict) #write the taxon information to a python object
    return ft, tax_dict

def get_meta(meta): #get the information contained in the meta file as a dictionary
    meta = open_txt(meta) #open the file
    meta_names = meta[0] #save the column names
    del meta[0] #delete the column names
    meta_dict = {} #create a dictionary of all information contained, with sample names as dictionary keys
    for a in range(len(meta)):
        meta_dict[meta[a][0]] = meta[a][1:]
    return meta, meta_names, meta_dict

def filter_seqs(ft, sf): #filter the sequences file to contain only the sequences in the given feature table
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
    colors_2 = [color_env['marine'], color_env['marine'], color_env['freshwater'], color_env['freshwater'], color_env['aquatic'], color_env['aquatic'], color_env['terrestrial'], color_env['terrestrial']]
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
    lab = mlines.Line2D([], [], color='w', marker='^', markersize=8, markeredgecolor='k', label='Lab', linestyle=' ')
    field = mlines.Line2D([], [], color='w', marker='o', markersize=8, markeredgecolor='k', label='Field', linestyle=' ')
    both = mlines.Line2D([], [], color='w', marker='*', markersize=8, markeredgecolor='k', label='Both', linestyle=' ')
    gap = mlines.Line2D([], [], color='w', marker='o', markersize=2, markeredgecolor='w', alpha=0, label='', linestyle=' ')
    plt.legend(handles=[marine,freshwater,aquatic, terrestrial, gap, lab, field, both], loc='lower left', fontsize=fs_main)  
    
    #now save the figure (using the file extension and dpi specified at the top of the file) and close it
    plt.savefig(basedir+'/figures/study_map'+ext, dpi=dpi, bbox_inches='tight')
    plt.close()
    return

def transform_for_NMDS(similarities): #transform the similarity matrix to 2D space (n_components=2)
    seed = np.random.RandomState(seed=3)
    X_true = similarities.iloc[0:].values
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed, dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_
    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12, dissimilarity="precomputed", random_state=seed, n_jobs=1, n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)
    npos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((npos ** 2).sum())
    clf = PCA()
    npos = clf.fit_transform(npos)
    return npos, nmds.stress_

def get_single_nmds(rows, filter_on, filt_ind, color_on, color_ind, ax, leg, colors, names, meta_dict, second_filter='', second_filter_ind='', npos=''):
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
    for a in range(len(names)):
        for b in range(len(color_on)):
            if meta_dict[names[a]][color_ind] == color_on[b]:
                color.append(colors[b])
    #if we don't already have the values for npos, then transform the similarity matrix for this
    if npos == '':
        s = pd.DataFrame(rows)
        npos, stress = transform_for_NMDS(s)
    size = 20
    #plot the values for nmds 1 (npos[a,0]) and nmds 2 (npos[a,1]) for all samples, giving the color as determined by sample type
    print(len(npos), len(color))
    for a in range(len(rows[0])):
        plt.scatter(npos[a,0], npos[a,1], color=color[a], marker='o', s=size, edgecolor='k')
    #get the legend handles incase these are being plotted
    handles = []
    for a in range(len(color_on)):
        handles.append(mlines.Line2D([], [], color=colors[a], marker='o', markersize=fs_main, markeredgecolor='k', label=color_on[a], linestyle=' '))
    return npos, handles

def nmds_plot(dist_matr_fn, meta_dict, basedir):
    #set up figure and axes
    plt.figure(figsize=(15,15))
    ax1 = plt.subplot2grid((2,2), (0,0))
    ax2 = plt.subplot2grid((2,2), (0,1))
    ax3 = plt.subplot2grid((4,4), (2,0))
    ax4 = plt.subplot2grid((4,4), (2,1))
    ax5 = plt.subplot2grid((4,4), (2,2))
    ax6 = plt.subplot2grid((4,4), (2,3))
    ax7 = plt.subplot2grid((4,4), (3,0))
    ax8 = plt.subplot2grid((4,4), (3,2))
    #these can be uncommented if there will now be lab-based freshwater (ax9) or terrestrial (ax10) studies
    #ax9 = plt.subplot2grid((4,4), (3,1))
    #ax10 = plt.subplot2grid((4,4), (3,3))
    
    if dist_matr_fn[-4:] == '.csv':
        dist_matr = open_csv(dist_matr_fn) #open the csv file of the distance matrix
    else:
        dist_matr = open_txt(dist_matr_fn) #open the txt file of the distance matrix
    for a in range(len(dist_matr)):
        dist_matr[a] = dist_matr[a][1:] #take all numbers but not the row name
    names = dist_matr[0] #take the column names (sample names) from the first row 
    del dist_matr[0] #delete this row from the distance matrix
    for a in range(len(dist_matr)): #turn all values in the matrix into floats
        for b in range(len(dist_matr)):
            dist_matr[a][b] = float(dist_matr[a][b])
    
    #set up all color options for the nmds plots - for coloring based on environment, source and study
    envs, env_index = ['marine', 'freshwater', 'aquatic', 'terrestrial'], 3
    color_env = ['#03A498', '#03B1FC', '#9A75FC', '#FD9A64']
    
    source, source_index = ['aliphatic', 'other plastic', 'unknown plastic', 'biofilm', 'planktonic', 'blank'], 10
    color_source = ['#5F9EA0', '#8B008B', '#3593FC', '#FCA94A', 'yellow', '#C3C3C3']
    
    study, study_index = ['AmaralZettler2015','AriasAndres2018','Canada2020','Curren2019','Delacuvellerie2019','DeTender2015','DeTender2017','DussudHudec2018','DussudMeistertzheim2018','ErniCassola2019','Esan2019','Frere2018','Hoellein2014','Hoellein2017','Jiang2018','Kesy2019','Kirstein2018','Kirstein2019','McCormick2014','McCormick2016','Oberbeckmann2016','Oberbeckmann2018','Ogonowski2018','Pinto2019','Pollet2018','Rosato2020','Syranidou2019','SyranidouPE2017','SyranidouPS2017','Tagg2019','Woodall2018','Wright2020','Wu2019','Xu2019','Zhang2019'], 0
    color_study = [(0.2150767543107187, 0.6494256446006359, 0.9751873123180738), (0.19871051089612168, 0.9949001161033785, 0.5171863529790246), (0.8005217032861965, 0.9808524817393911, 0.07919858947341929), (0.7800432594341148, 0.05614030366701428, 0.9610189983758901), (0.9730904559406481, 0.8504292507337813, 0.11446201949258006), (0.08991025975344336, 0.9047578028038412, 0.9811497599648167), (0.9813929893132025, 0.05434493047499922, 0.6900350279640525), (0.11741250827908623, 0.16664119434192998, 0.9789145143788549), (0.9775787324850631, 0.042833970592817794, 0.3633178889558728), (0.9932441025402244, 0.4181904516868864, 0.11816245993731866), (0.9527604880433751, 0.2050145933864841, 0.05030854621609282), (0.47226467042905784, 0.9792001373201127, 0.04537164567869634), (0.9631274626909789, 0.21165590086011288, 0.21165590086011288), (0.1493280952250614, 0.9739436651258512, 0.9032623305629265), (0.06488974861545382, 0.9575112138086248, 0.7279799799018097), (0.9717950954037526, 0.09835331685778903, 0.5475519458242837), (0.6648633784092389, 0.9590343952439288, 0.16703550376591825), (0.5484424509141032, 0.19589379540263863, 0.9670939793339678), (0.212011506314699, 0.9835147571548518, 0.1124626997546797), (0.24615647818088732, 0.1509195356286156, 0.9842427829609864), (0.9588337609318365, 0.07464252425362417, 0.22621816482703205), (0.1652313333209956, 0.34649769281916903, 0.9582716561255071), (0.9533751865974686, 0.675410329349355, 0.06894154989892509), (0.09233059229333751, 0.44545868121504567, 0.9751508145976095), (0.32613834523101703, 0.06268882221557415, 0.9847621527696233), (0.9655508415393406, 0.5497733629079831, 0.10953838553360429), (0.1903099435917901, 0.9721944317625262, 0.3690263980308156), (0.9261814242997817, 0.0492633867541552, 0.9519731312864173), (0.9855382669820902, 0.026275180087142136, 0.8485006831399537), (0.9614059316276589, 0.9868945709577983, 0.09479219440291908), (0.06578923981775453, 0.9968090810748765, 0.11899037360387621), (0.3685840351292384, 0.9593160822362674, 0.13229121628642715), (0.056776251268179534, 0.9622398986535152, 0.5741840497740859), (0.6403409835890199, 0.07858523328440081, 0.9722875633144761), (0.14424182754131876, 0.7505049729717633, 0.960365292543841)]
    #for the first two plots, we are only filtering/coloring based on environment/study (and not also separating to environment and lab/field as in subsequent plots)
    filter_on, filter_index = 'none', 'none'
    second_filter, second_filter_ind = '', ''
    npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, envs, env_index, ax1, 'upper left', color_env, names, meta_dict, second_filter, second_filter_ind, '')
    plt.sca(ax1)
    plt.legend(handles=handles, loc='upper left', fontsize=fs_main)
    npos, handles = get_single_nmds(dist_matr, filter_on, filter_index, study, study_index, ax2, 'upper right', color_study, names, meta_dict, second_filter, second_filter_ind, npos)
    plt.sca(ax2)
    plt.legend(handles=handles, bbox_to_anchor=(1.05,1.025), fontsize=fs_main)
    
    #for the second row of plots, we are only plotting those samples that are taken from the field, and are separating each plot to be only one environment
    second_filter, second_filter_ind = 'field', 5
    npos, handles = get_single_nmds(dist_matr, envs[0], env_index, source, source_index, ax3, '', color_source, names, meta_dict, second_filter, second_filter_ind)
    npos, handles = get_single_nmds(dist_matr, envs[1], env_index, source, source_index, ax5, '', color_source, names, meta_dict, second_filter, second_filter_ind)
    npos, handles = get_single_nmds(dist_matr, envs[2], env_index, source, source_index, ax4, '', color_source, names, meta_dict, second_filter, second_filter_ind)
    npos, handles = get_single_nmds(dist_matr, envs[3], env_index, source, source_index, ax6, 'upper right', color_source, names, meta_dict, second_filter, second_filter_ind)
    plt.sca(ax6)
    plt.legend(handles=handles, bbox_to_anchor=(1.101,1.025), fontsize=fs_main)
    
    #for the third row of plots, we are only plotting those samples that are taken from the lab, and are separating each plot to be only one environment
    second_filter, second_filter_ind = 'lab', 5
    npos, handles = get_single_nmds(dist_matr, envs[0], env_index, source, source_index, ax7, '', color_source, names, meta_dict, second_filter, second_filter_ind)
    npos, handles = get_single_nmds(dist_matr, envs[1], env_index, source, source_index, ax8, '', color_source, names, meta_dict, second_filter, second_filter_ind)
    #uncomment these if there will now be lab-based freshwater (ax9) or terrestrial (ax10) studies
    #npos, handles = get_single_nmds(rows, envs[2], env_index, source, source_index, ax9, '', color_source, names, second_filter, second_filter_ind)
    #npos, handles = get_single_nmds(rows, envs[3], env_index, source, source_index, ax10, '', color_source, names, second_filter, second_filter_ind)
    
    #now add titles and axis labels to all plots
    all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    titles = ['All environments', 'All environments', 'Marine', 'Aquatic', 'Freshwater', 'Terrestrial', '', '', '', '']
    yax = ['nMDS2', '', 'Field\nnMDS2', '', '', '', 'Lab\nnMDS2', '', '', 'nMDS2']
    xax = ['nMDS1', 'nMDS1', '', 'nMDS1', '', 'nMDS1', 'nMDS1', 'nMDS1']
    title_left = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    for a in range(len(all_ax)):
        plt.sca(all_ax[a])
        plt.title(titles[a], fontsize=fs_title, fontweight='bold'), plt.ylabel(yax[a], fontsize=fs_main, fontweight='bold'), plt.xlabel(xax[a], fontsize=fs_main, fontweight='bold')
        plt.title(title_left[a], loc='left', fontsize=fs_title, fontweight='bold')
    new_fn = dist_matr_fn.split('/')
    plt.savefig(basedir+'/figures/'+new_fn[1][:-4]+'_nmds'+ext, dpi=600, bbox_inches='tight')
    plt.close()
    return

def plot_box(ax, l, r, b, t, line_col):
    plt.sca(ax)
    plt.plot([l, r], [b, b], color=line_col, lw=2)
    plt.plot([l, r], [t, t], color=line_col, lw=2)
    plt.plot([l, l], [t, b], color=line_col, lw=2)
    plt.plot([r, r], [t, b], color=line_col, lw=2)
    return

def similarity_heatmap(dist_matr_fn, basedir):
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
    order = marine+freshwater+aquatic+terrestrial #the order for the studies to be plotted (so that environments are grouped together)
    lens = [len(marine), len(freshwater), len(aquatic), len(terrestrial)] #get the length (/number of studies) of each of the environments
    text_cols = [] #set up a list for the text colors to be added to so these match the environment type
    for a in range(len(marine)):
        text_cols.append(color_env['marine'])
    for a in range(len(freshwater)):
        text_cols.append(color_env['freshwater'])
    for a in range(len(aquatic)):
        text_cols.append(color_env['aquatic'])
    for a in range(len(terrestrial)):
        text_cols.append(color_env['terrestrial'])
    new_text_cols = [] #now rearrange them so they are in the order they will be plotted
    for a in range(len(order)):
        for b in range(len(order)):
            if a == order[b]:
                new_text_cols.append(text_cols[b])
    each_study, each_row = [], []
    #go  through the distance matrix and calculate means for the distance between samples within and between different studies
    """
    go through this and comment more thoroughly?
    """
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
    plt.figure(figsize=(10,12))
    ax1 = plt.subplot(111)
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
    cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=colormap),ax=ax1, shrink=0.4, pad=0.02, orientation='horizontal')
    cb.ax.tick_params(labelsize=fs_small)
    cb.set_label('Unifrac distance', fontsize=fs_main)
    #change x and y limits
    plt.xlim([0.5,ally[-1]+1])
    plt.ylim([0,ally[-1]+0.5])
    colors = ['#03A498', '#03B1FC', '#9A75FC', '#FD9A64']
    line_col = 'w'

    #plot white boxes around samples of the same environments
    l, r, b, t = 0.54, lens[0]+0.43, ally[-1]-lens[0]+0.57, ally[-1]+0.48#left, right, bottom, top
    plot_box(ax1, l, r, b, t, line_col)
    l, r, b, t = l+lens[0], r+lens[1], b-lens[1], t-lens[0]
    plot_box(ax1, l, r, b, t, line_col)
    l, r, b, t = l+lens[1], r+lens[2], b-lens[2], t-lens[1]
    plot_box(ax1, l, r, b, t, line_col)
    l, r, b, t = l+lens[2], r+lens[3], b-lens[3], t-lens[2]
    plot_box(ax1, l, r, b, t, line_col)
    for a in range(len(studies)): #rename the studies to a full study name rather than the short code they were given
        studies[a] = name_dict[studies[a]]
    empty = []
    text_cols.reverse()
    for a in range(len(ally)): #plot the study names with the correct colors  on the y axis
        plt.text(0.2, ally[a], studies[a], color=text_cols[a], va='center', ha='right', fontsize=fs_small)
        empty.append('')
    plt.yticks(ally, empty, fontsize=fs_small) #remove the y ticks
    studies.reverse()
    for a in range(len(x)): #not plot the study names on the x axis (the list is reversed as we plotted from the bottom up)
        plt.text(x[a], ally[-1]+0.8, studies[a], rotation=90, color=text_cols[-(a+1)], va='bottom', ha='center', fontsize=fs_small)
    plt.xticks(x, empty, fontsize=fs_small) #remove the x ticks
    ax1.xaxis.tick_top()
    new_fn = dist_matr_fn.split('/') #get only the second part of the file name (allows for the same function for agglomerated and not agglomerated data)
    plt.savefig(basedir+'/figures/'+new_fn[1][:-4]+'_heatmap'+ext, dpi=600, bbox_inches='tight') #save the figure
    plt.close()
    return

def generate_rf(X_train, y_train, X_test, y_test, rc='cls', est=10000): #generate a random forest based on the data being split to train and test
    if rc == 'cls': #if we are using a classification (i.e. discrete categories)
        rf = RandomForestClassifier(n_estimators=est, min_samples_leaf=3)
    else: #if we are using a regression (i.e. continuous categories)
        rf = RandomForestRegressor(n_estimators=est, min_samples_leaf=3)
    rf.fit(X_train, y_train) #fit out data to either the regressor or classifier
    return rf, rf.score(X_test, y_test), rf.feature_importances_ #return the forest (i.e. features), score (how well it classifies the test data) and feature importances (ASV importances)

def get_rf_tree(df_rf, fn, basedir, meta, ft): #for each meta category, get a phylogenetic tree and heatmap from R
    ft = pd.read_csv(ft, header=0, index_col=0) #read in the csv file
    ft.columns = meta #rename the columns
    ft = ft.reset_index()
    ft = ft.rename(columns={'index':'ASV'})
    ft = ft.set_index(['ASV'])
    df_rf = df_rf.set_index(['ASV'])
    cols = list(ft.columns)
    #if the columns are numeric, then turn them into a float rather than a string
    for a in range(len(cols)):
        try:
            cols[a] = float(cols[a])
            cols[a] = round(cols[a])
        except:
            cols[a] = cols[a]
    ft.columns=cols
    ft = ft.groupby(by=ft.columns, axis=1).mean() #group all of the columns that are the same
    ft = ft*100 #get relative abundance (%) - not really necessary 
    ft = ft.reindex(sorted(ft.columns), axis=1) #sort the columns either alphabetically or numerically
    merge = df_rf.merge(right=ft, on='ASV') #merge the importance dataframe with the feature table (this means we only keep those ASVs that are important)
    merge.drop('Importance', axis=1, inplace=True) #now remove the importance column again
    merge = merge.div(merge.max(axis=1), axis=0) #normalise within each ASV (to make the colors visible for all even if there are large differences in abundance)
    merge = merge.reset_index()
    merge.to_csv(basedir+'/random_forest_R/random_forest.csv', index=False) #save the csv to the random forest R directory (this means that we don't need to worry about telling the R script the file name as it is the same for all)
    os.system("/usr/local/bin/Rscript plot_RF_tree_heatmap.R") #run the R script that plots the phylogenetic tree (of only the important ASVs) and heatmap for this meta category
    os.rename(basedir+'/random_forest_R/tree_and_heatmap.pdf', basedir+'/figures/random_forest/'+fn[:-4]+'.pdf') #rename the resulting figure and move it to the figures directory
    return

def get_overall_rf_tree(df_rf, top_features_plot, basedir, names): #get the overall phylogenetic tree and heatmap for all meta categories
    ft = pd.DataFrame(top_features_plot, columns=['ASV']+names) #open the previous .csv file
    ft = ft.set_index(['ASV'])
    ft = ft.div(ft.max(axis=1), axis=0) #normalise within each ASV (to make it so we can see all plot colors)
    ft = ft.reset_index()
    ft.to_csv(basedir+'/random_forest_R/random_forest.csv', index=False) #save this csv to the directory that R uses
    os.system("/usr/local/bin/Rscript plot_RF_tree_heatmap.R") #run the R script
    os.rename(basedir+'/random_forest_R/tree_and_heatmap.pdf', basedir+'/figures/top_features_random_forest.pdf') #move the resulting plot to the figures directory
    return

def get_single_forest(asv, meta, asv_names, basedir, ft, est=10000): #get forests for each of the metadata categories
    #set up a dictionary saying which meta categories are continuous (reg) or discrete (cls) and will use regressors or classifiers, respectively
    cls_reg = {'Study':'cls', 'Latitude':'reg', 'Longitude':'reg', 'Environment':'cls', 'WaterOrSediment':'cls',
       'LabOrField':'cls', 'IncubationOrCollection':'cls', 'Source':'cls', 'MaterialType':'cls',
       'PlasticTypeSpecific':'cls', 'PlasticTypeGeneral':'cls', 'DEPTH':'reg', 'IncubationTime':'reg',
       'IncubationGeneral':'cls', 'Temperature':'reg', 'Salinity':'reg', 'Light':'cls', 'Season':'cls',
       'PrimerPair':'cls', 'DNAExtraction':'cls', 'CollectionDate':'cls'}
    #normalise all of the data
    asv = asv.rename_axis('ID').values
    max_abs_scaler = preprocessing.MaxAbsScaler()
    asv_scale = max_abs_scaler.fit_transform(asv)
    meta_split, reg_cls, col_names = [], [], []
    for col in meta.columns: #make lists of all of the sample names based on the sample type within each meta category
        meta_split.append(meta[col])
        reg_cls.append(cls_reg[col])
        col_names.append(col)
    scores, importances = [], []
    for a in range(len(meta_split)): #for each meta category
        meta = meta_split[a].rename_axis('ID').values #split the metadata based in this category
        X_train, X_test, y_train, y_test = train_test_split(asv_scale, meta, test_size=0.2) #split the data to training and test data (with 80% of samples being used for training and 20% for testing)
        RF, RF_score, RF_importances = generate_rf(X_train, y_train, X_test, y_test, rc=reg_cls[a], est=est) #generate the random forest for this meta category
        #add the scores to the overall list
        scores.append(RF_score)
        importances.append(RF_importances)
        print(col_names[a], RF_score)
        #make a dataframe with this forest
        this_rf = {'ASV':asv_names, 'Importance':RF_importances}
        df_rf = pd.DataFrame(this_rf, columns = ['ASV', 'Importance'])
        df_rf.sort_values(by=['Importance'], ascending=False, inplace=True) #sort the values by importance
        df_rf = df_rf[:100] #now only keep the top 100 features
        df_rf.to_csv(basedir+'/random_forest/'+col_names[a]+'_random_forest.csv', index=False) #save this to a .csv file
        get_rf_tree(df_rf, col_names[a]+'_random_forest.csv', basedir, meta, ft) #get the phylogenetic tree and heatmap plot
    return col_names, scores, importances

def get_meta_rf(fn, sample_names):
    with open('metadata.txt', 'rU') as f:
        meta = []
        for row in csv.reader(f, delimiter='\t', lineterminator='\n'):
            meta.append(row)
    new_meta, row_names = [meta[0][1:]], []
    for a in range(len(sample_names)):
        for b in range(len(meta)):
            if sample_names[a] == meta[b][0]:
                new_meta.append(meta[b][1:])
                row_names.append(meta[b][0])
    meta = pd.DataFrame(new_meta[1:], columns=new_meta[0], index=row_names)
    return meta

def random_forests(ft, tax, basedir, est=10000):
    #number of estimators to use for random forests (this should be 10,000 for robustness, but can be reduced to 1000 (or even less) to save time and check that everything is running OK before doing a more thorough analysis)
    asv = pd.read_csv(ft, header=0, index_col=0) #fead feature table
    sample_names = asv.columns #get sample names
    asv = asv.T #transpose
    asv_names = list(asv.columns.values) #get list of ASVs
    rn_meta = get_meta_rf('meta_dropped_samples.dictionary', sample_names) #get the meta table in the correct format
    meta = rn_meta
    col_names, scores, importances = get_single_forest(asv, meta, asv_names, basedir, ft, est) #calculate the scores and importances for all meta categories
    feature_importance = []
    for a in range(len(importances)): #calculate the importances for all ASVs, weighted by the importance within each meta category and the score for that category
        for b in range(len(importances[a])):
            importances[a][b] *= scores[a]
    trans_importances = []
    for a in range(len(importances[0])):
        this_importance = []
        for b in range(len(importances)):
            this_importance.append(importances[b][a])
        feature_importance.append(sum(this_importance))
        trans_importances.append(this_importance)
    #sort the ASVs (features) based on how important they are 
    feature_names, features = list(asv.columns.values), []
    for a in range(len(feature_importance)):
        names = [feature_names[a], tax[feature_names[a]]]
        features.append(names+trans_importances[a])
    features = sorted(zip(feature_importance, features))
    features.reverse()
    top_features, top_features_full, asv_scores = [], [], []
    top_features_plot = []
    #get full details of the taxonomy for only the top 100 important ASVs (features)
    for a in range(100):
        top_features.append(features[a][1][0])
        this_tax = features[a][1][1]
        while len(this_tax) < 7:
            this_tax.append('')
        asv_scores.append(features[a][0])
        top_features_full.append([features[a][0]]+[features[a][1][0]]+this_tax+features[a][1][2:])
        top_features_plot.append([features[a][1][0]]+features[a][1][2:])
    #get a dataframe with the importance of all ASVs
    new_df = {'ASV':top_features, 'Importance':asv_scores}
    df_rf = pd.DataFrame(new_df, columns=['ASV', 'Importance'])
    #save a csv file with full information on the top 100 ASVs as well as their contribution to each category and the importance of each meta category
    with open('random_forest/top_features_full.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['Score', 'ASV', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Species name']+col_names)
        writer.writerow(['Single scores', '', '', '', '', '', '', '', '', '']+scores)
        for feat in top_features_full:
            writer.writerow(feat)
    get_overall_rf_tree(df_rf, top_features_plot, basedir, col_names) #get the phylogenetic tree and heatmap for all features
    return

def get_venn_labels(data, fill=["number"]): #calculate the number of overlapping ASVs between different sample types
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
    #get a venn diagram with specified colors and names, with up to 5 ellipses (these are not plotted if that sample type is not present)
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
    #function to calculate a range of different diversity metrics depending on the diversity that we put in, as well as sample which should be a list of abundance values
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
    #get a list of either 20 or 40 colors depending on how many unique phyla we have
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
    #set up the figure and all axes
    fig = plt.figure(figsize=(12,20))
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
        env_source_dict[colnames[a]] = meta_dict[colnames[a]][3]+' '+meta_dict[colnames[a]][5]+'\n'+meta_dict[colnames[a]][10]
        env_source.append(meta_dict[colnames[a]][3]+' '+meta_dict[colnames[a]][5]+'\n'+meta_dict[colnames[a]][10])
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
    plt.savefig(basedir+'/figures/dendro_venn'+ext, dpi=dpi, bbox_inches='tight') #save the figure
    plt.close()
    return

def group_ft_level(ft, level, tax_dict, basedir, rename=False): #This script was written specifically for grouping before ancom, but the rename option is there to make it compatible with other functions that might want to regroup but not rename to an ASV (this renaming is only to have a representative sequence for each group)
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
        with open(basedir+'/ancom/taxonomy_name_only.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['OTUID', 'Species name'])
            for asv in representative_asv:
                writer.writerow(asv)
    return ft

def tree_heatmap(ft, meta_dict, basedir, tax_dict, level=7):
    ft = pd.read_csv(ft, header=0, index_col=0)
    ft = ft*100
    if level < 7:
        ft = group_ft_level(ft, level, tax_dict, basedir, rename=True)
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
    return

def metacoder(ft, tax_dict, meta_dict, basedir):
    ft = pd.read_csv(ft, header=0, index_col=0)
    ft = ft*100
    env, source, inc_time, source_inc, source_inc_dict, source_dict = [], [], [], [], {}, {}
    envs = ['marine', 'freshwater', 'aquatic', 'terrestrial']
    samples = list(ft.columns)
    for a in range(len(samples)):
        env.append(meta_dict[samples[a]][3])
        source.append(meta_dict[samples[a]][10])
        source_dict[samples[a]] = meta_dict[samples[a]][10]
        inc_time.append(meta_dict[samples[a]][13])
        source_inc.append(meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13])
        source_inc_dict[samples[a]] = meta_dict[samples[a]][10]+' '+meta_dict[samples[a]][13]
    comparisons1 = [['aliphatic early', 'aliphatic late'], ['other plastic early', 'other plastic late'], ['biofilm early', 'biofilm late'],
                   ['aliphatic early', 'other plastic early'], ['aliphatic early', 'unknown plastic early'], ['aliphatic early', 'biofilm early'],
                   ['aliphatic late', 'other plastic late'], ['aliphatic late', 'unknown plastic late'], ['aliphatic late', 'biofilm late'],
                   ['other plastic early', 'unknown plastic early'], ['other plastic early', 'biofilm early'],
                   ['other plastic late', 'unknown plastic late'], ['other plastic late', 'biofilm late'],
                   ['unknown plastic early', 'biofilm early'], ['unknown plastic late', 'biofilm late']]
    comparisons2 = [['aliphatic', 'other plastic'], ['aliphatic', 'unknown plastic'], ['aliphatic', 'biofilm'], ['aliphatic', 'planktonic'], 
                    ['other plastic', 'unknown plastic'], ['other plastic', 'biofilm'], ['other plastic', 'planktonic'], 
                    ['unknown plastic', 'biofilm'], ['unknown plastic', 'planktonic'], 
                    ['biofilm', 'planktonic']]
    comps = [comparisons1, comparisons2]
    asvs = list(ft.index.values)
    lineage = []
    for a in range(len(asvs)):
        lineage_single = tax_dict[asvs[a]]
        if lineage_single[0] == 'Archaea':
            ft.drop(asvs[a], axis=0, inplace=True)
            continue
        string = ''
        for b in range(len(lineage_single)):
            if b > 0 and lineage_single[b] == lineage_single[b-1]:
                length = b
        for b in range(length):
            if lineage_single[b] != '':
                if string != '':
                    string += ';'
                string += lineage_single[b]
        lineage.append(string)
    ft.insert(0, "lineage", lineage, True)
    for a in range(len(envs)):
        for z in range(len(comps)):
            comparisons = comps[z]
            if comps[z] == comparisons1:
                treat = source_inc
                treat_dict = source_inc_dict
            else:
                treat = source
                treat_dict = source_dict
            for b in range(len(comparisons)):
                if len(comparisons[b]) < 3:
                    comparisons[b].append('lineage')
                keeping = [True]
                for c in range(len(env)):
                    if envs[a] == env[c] and treat[c] in comparisons[b]:
                        keeping.append(True)
                    else:
                        keeping.append(False)
                env_comp_ft = ft.loc[:, keeping]
                env_comp_ft.rename(columns=treat_dict, inplace=True)
                names = list(set(env_comp_ft.columns))
                if names == ['lineage'] or len(names) == 2:
                    break
                names_dict = {}
                for d in range(len(names)):
                    if names[d] == comparisons[b][0]:
                        names_dict[names[d]] = 'Treat1'
                        names[d] = names[d].replace(' early', '')
                        names[d] = names[d].replace(' late', '')
                        col1 = color_source[names[d]]
                    elif names[d] == comparisons[b][1]:
                        names_dict[names[d]] = 'Treat2'
                        names[d] = names[d].replace(' early', '')
                        names[d] = names[d].replace(' late', '')
                        col2 = color_source[names[d]]
                    else:
                        names_dict['lineage'] = 'lineage'
                if comparisons == comparisons1 and b < 3:
                    col1, col2 = '#ff0000', '#25df50'
                env_comp_ft.rename(columns=names_dict, inplace=True)
                command = "/usr/local/bin/Rscript metacoder.R '"+col1+"' '"+col2+"'"
                env_comp_ft.to_csv(basedir+'/metacoder/metacoder.csv')
                os.system(command)
                os.rename(basedir+'/metacoder/metacoder.pdf', basedir+'/figures/metacoder/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'.pdf')
                os.rename(basedir+'/metacoder/metacoder_labels.pdf', basedir+'/figures/metacoder/'+envs[a]+'_'+comparisons[b][0]+'_'+comparisons[b][1]+'_labels.pdf')
                
    return

def map_abundance():
    
    return

def picrust_plots():
    
    return

def plots_per_study():
    
    return

