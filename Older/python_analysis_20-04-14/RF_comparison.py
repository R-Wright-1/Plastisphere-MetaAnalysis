#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 14:21:08 2020

@author: robynwright
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from lifelines.utils import concordance_index

fs_title, fs_main, fs_small = 14, 10, 8 #font sizes to use in figures
basedir = '/Users/robynwright/Documents/OneDrive/Papers_writing/Plastisphere_Meta-analysis/test_recreate_analyses/paper_data/'
rename_plots = {'AmaralZettler':'Amaral Zettler $et$ $al$. 2015', 'AriasAndres':'Arias-Andres $et$ $al$. 2018', 'Canada':'Canada $et$ $al$. 2020', 'Curren':'Curren & Leong 2019', 'Delacuvellerie':'Delacuvellerie $et$ $al$. 2019', 'DeTender':'De Tender $et$ $al$. 2015', 'DeTenderB':'De Tender $et$ $al$. 2017', 'DussudHudec':'Dussud, Hudec $et$ $al$. 2018', 'DussudMeistertzheim':'Dussud, Meistertzheim $et$ $al$. 2018', 'ErniCassola':'Erni-Cassola $et$ $al$. 2019', 'Esan':'Esan $et$ $al$. 2019', 'Frere':'Frere $et$ $al$. 2018', 'Hoellein':'Hoellein $et$ $al$. 2014', 'HoelleinB':'Hoellein $et$ $al$. 2017', 'Jiang':'Jiang $et$ $al$. 2018', 'Kesy':'Kesy $et$ $al$. 2019', 'Kirstein':'Kirstein $et$ $al$. 2018', 'KirsteinB':'Kirstein $et$ $al$. 2019', 'McCormick':'McCormick $et$ $al$. 2014', 'McCormickB':'McCormick $et$ $al$. 2016', 'Oberbeckmann':'Oberbeckmann $et$ $al$. 2016', 'OberbeckmannB':'Oberbeckmann $et$ $al$. 2018', 'Ogonowski':'Ogonowski $et$ $al$. 2018', 'Parrish':'Parrish $et$ $al$. 2019', 'Pinto':'Pinto $et$ $al$. 2019', 'Pollet':'Pollet $et$ $al$. 2018', 'Rosato':'Rosato $et$ $al$. 2020', 'Syranidou':'Syranidou ', 'SyranidouPE':'Syranidou $et$ $al$. 2017a', 'SyranidouPS':'Syranidou $et$ $al$. 2017b', 'Tagg':'Tagg $et$ $al$. 2019', 'Woodall':'Woodall $et$ $al$. 2018', 'Wu':'Wu $et$ $al$. 2019', 'Xu':'Xu $et$ $al$. 2019', 'Zhang':'Zhang $et$ $al$. 2019', 'WaterOrSediment':'Water or Sediment', 'LabOrField':'Laboratory or Field', 'IncubationOrCollection':'Incubation or Collection', 'MaterialType':'Material type', 'PlasticTypeSpecific':'Plastic type (specific)', 'PlasticTypeGeneral':'Plastic type (general)', 'DEPTH':'Depth', 'IncubationTime':'Incubation time (specific)', 'IncubationGeneral':'Incubation time (general)', 'PrimerPair':'Primer pair', 'DNAExtraction':'DNA extraction method', 'lab':'Laboratory', 'not_plastic':'Not plastic', 'aged_oxope':'Aged Oxo-PE', 'freeliving':'Free living', 'particleassociated':'Particle associated', 'oxope':'Oxo-PE', 'rinse_pe':'PE rinse water', 'rinse_ps':'PS rinse water', 'rinse_wood':'Wood rinse water', 'bhet':'BHET', 'hdpe':'HDPE', 'ldpe':'LDPE', 'na':'NA', 'pa':'PA', 'pe':'PE', 'pes':'PES', 'pestur':'PESTUR', 'pet':'PET', 'phbv':'PHBV', 'pla':'PLA', 'pp':'PP', 'ps':'PS', 'pvc':'PVC', 'san':'SAN', 'w_pe':'Weathered PE', '10:14':'10:14 light:dark', '12:12':'12:12 light:dark', '16:08':'16:08 light:dark', '27F_519R':'27F-519R', '319F_806R':'319F-806R', '338F_806R':'338F-806R', '341F_785R':'341F-785R', '341F_802R':'341F-802R', '341F_806R':'341F-806R', '515F_806R':'515F-806R', '515FY_926R':'515FY-926R', '518F_926R':'518F-926R', '543F_783R':'543F-783R', '967F_1064R':'967F-1064R', 'B969F_BA1406R':'B969F-BA1406R', 'redextract_sigma':'REDExtract-$N$-AmpTM', 'gentra_puregene':'Gentra Puregene', 'purelink':'PureLink', 'powersoil':'PowerSoil', 'phenol_chloroform':'Phenol-Chloroform', 'powerbiofilm':'PowerBiofilm', 'ultraclean_soil':'UltraClean soil', 'fastdna_soil':'FastDNA soil', 'orders':'Order', 'classes':'Class', 'phyla':'Phylum', 'genera':'Genera', 'families':'Family', 'species':'Species', 'ASVs':'ASV', 'kingdoms':'Kingdom', 'PlasticOnly':'Plastic only', '534f_783r':'534F-783R', 'Phenol_chloroform':'Phenol-chloroform'}
ext = '.png'
norm_names = {'rare':'Rarefied', 'rel_abun':'Relative\nabundance', 'log':'Log', 'clr':'CLR'}

def annotate_heatmap(ax, df, cmap='inferno', yticks=True, xticks=True, rnd=1, annotate=True, italics=False, vmax=False, vmin=False, annotate_only=False):
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
                if cmap == 'Purples' or cmap == 'Blues':
                    if num > 0.9: col = 'w'
                    else: col = 'k'
                plt.text(a+0.5, b+0.5, str(num), ha='center', va='center', color=col, fontsize=fs_small) #and plot this text label
    if annotate_only:
        return
    if vmax != False and vmin != False:
        plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=vmin, vmax=vmax) #make the heatmap
    elif vmax != False:
        plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=0, vmax=vmax) #make the heatmap
    elif vmin != False:
        plt.pcolor(df, edgecolor='k', cmap=cmap, vmin=vmin) #make the heatmap
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

def get_rf_comparison(basedir):
    norm = ['rare', 'rel_abun', 'log', 'clr']
    levels = ['phyla', 'classes', 'orders', 'families', 'genera', 'species']  
    order_cols, order_rows = [], []
    fig = plt.figure(figsize=(13, 10)) 
    for l in range(len(levels)):
        level = levels[l]
        dfs, dfs_tax = [], []
        for n in norm:
            df = pd.read_csv(basedir+'random_forest/'+n+'/'+level+'_overall.csv', header=0, index_col=0)
            try:
                df_tax = df.drop(['Score', 'OOB_score'], axis=0)
            except:
                df_tax = df.drop(['Score'], axis=0)
            df_tax['Mean'] = df_tax.mean(axis=1)
            df_tax = pd.DataFrame(df_tax.loc[:, 'Mean'])
            df = pd.DataFrame(df.loc['Score', :])
            df = df.transpose()
            df.rename(index={'Score':n}, inplace=True)
            dfs.append(df)
            df_tax = df_tax.transpose()
            df_tax.rename(index={'Mean':n}, inplace=True)
            dfs_tax.append(df_tax)
        dfs = pd.concat(dfs)*100
        dfs_tax = pd.concat(dfs_tax).fillna(value=0)
        for a in range(2):
            dfs['Mean'] = dfs.mean(axis=1)
            dfs = dfs.sort_values(by='Mean', axis=0, ascending=True)
            if a == 1: 
                dfs = dfs.sort_values(by='Mean', axis=0, ascending=False)
            else:
                df_mean = pd.DataFrame(dfs.loc[:, 'Mean'])
            dfs.drop('Mean', axis=1, inplace=True)
            dfs = dfs.transpose()
        dfs.rename(index=norm_names, columns=rename_plots, inplace=True)
        df_mean.rename(index=norm_names, inplace=True)
        order_rows = ['Rarefied', 'Relative\nabundance', 'Log', 'CLR']
        if order_cols == []:
            order_cols = list(dfs.columns)
        else:
            dfs = dfs.loc[:, order_cols]
        dfs = dfs.loc[order_rows, :]
        df_mean = df_mean.loc[order_rows, :]
        ax = plt.subplot2grid((6,120), (l, 5), colspan=76)
        ax_mean = plt.subplot2grid((6, 120), (l, 0), colspan=4)
        ax_con = plt.subplot2grid((6,6), (l, 4))
        ax_con_tax = plt.subplot2grid((6,6), (l, 5))
        xtcks = False
        if l == 5: xtcks = True
        annotate_heatmap(ax, dfs, xticks=xtcks, yticks=False)
        
        annotate_heatmap(ax_mean, df_mean, xticks=xtcks)
        ax_mean.set_ylabel(level.capitalize(), fontsize=fs_main, fontweight='bold')
        concs = []
        for a in range(4):
            conc = []
            for b in range(4):
                l1 = dfs.iloc[a, :].values
                l2 = dfs.iloc[b, :].values
                conc.append(concordance_index(l1, l2))
            concs.append(conc)
        concs = pd.DataFrame(concs, index=dfs.index.values, columns=dfs.index.values)
        annotate_heatmap(ax_con, concs, cmap='Blues', rnd=2, yticks=False, xticks=xtcks, vmin=0.75)
        concs_tax = []
        for a in range(4):
            conc_tax = []
            for b in range(4):
                l1 = dfs_tax.iloc[a, :].values
                l2 = dfs_tax.iloc[b, :].values
                conc_tax.append(concordance_index(l1, l2))
            concs_tax.append(conc_tax)
        concs_tax = pd.DataFrame(concs_tax, index=dfs_tax.index.values, columns=dfs_tax.index.values)
        concs_tax.rename(columns=norm_names, inplace=True)
        annotate_heatmap(ax_con_tax, concs_tax, cmap='Purples', rnd=2, yticks=False, xticks=xtcks, vmin=0.75)
        
        if l == 0:
            ax.set_title('Classification accuracy (%)', fontsize=fs_title, fontweight='bold')
            ax_con.set_title('Concordance in\nclassification accuracy', fontsize=fs_main, fontweight='bold')
            ax_con_tax.set_title('Concordance in\nfeature importance', fontsize=fs_main, fontweight='bold')
    
    plt.savefig(basedir+'/figures/RF_compare'+ext, dpi=600, bbox_inches='tight')
    plt.close()
    return 

get_rf_comparison(basedir)
            