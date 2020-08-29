#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 11:39:06 2020

@author: robynwright
"""

import pandas as pd
import os 

dir1 = '/Users/robynwright/Documents/OneDrive/Papers_writing/Plastisphere_Meta-analysis/test_recreate_analyses/recreate_analyses/random_forest/'
dir2 = '/Users/robynwright/Documents/OneDrive/Papers_writing/Plastisphere_Meta-analysis/test_recreate_analyses/recreate_analyses/random_forest/rel_abun/'

levels = ['phyla', 'classes', 'orders', 'families', 'genera', 'species']
cut = 50

for level in levels:
    fi1 = pd.read_csv(dir1+level+'_overall.csv', header=0, index_col=0)
    fi2 = pd.read_csv(dir2+level+'_overall.csv', header=0, index_col=0)
    fi1.drop(['Score'], axis=0, inplace=True)
    fi2.drop(['Score'], axis=0, inplace=True)
    """
    fi1['Mean'] = fi1.mean(axis=1)
    fi2['Mean'] = fi2.mean(axis=1)
    fi1 = pd.DataFrame(fi1.loc[:, 'Mean'])
    fi2 = pd.DataFrame(fi2.loc[:, 'Mean'])
    fi1 = fi1.sort_values(by='Mean', axis=0, ascending=False)[:cut]
    fi2 = fi2.sort_values(by='Mean', axis=0, ascending=False)[:cut]
    """
    
    fi1 = fi1[fi1.max(axis=1) > 0.1]
    fi2 = fi1[fi1.max(axis=1) > 0.1]
    
    vals1 = list(fi1.index.values)
    vals2 = list(fi2.index.values)
    
    vals1 = [x for x in vals1 if x != 'Score']
    vals1 = [x for x in vals1 if x != 'OOB_score']
    
    vals2 = [x for x in vals2 if x != 'Score']
    vals2 = [x for x in vals2 if x != 'OOB_score']
    
    overlap = [x for x in vals1 if x in vals2]
    total = list(set(vals1+vals2))
    
    if len(overlap) > 0:
        print(len(overlap)/len(total))
    else:
        print('None here')
