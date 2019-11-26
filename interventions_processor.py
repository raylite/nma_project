#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:27:18 2019

@author: ja18581
"""
import pandas as pd
import os
from pathlib import Path

def load_interventions():
    interventions = pd.read_csv('intervention.csv', index_col=False).sort_values(by='name')
    interventions = interventions.drop('Unnamed: 0', axis = 1)
    interventions = interventions.applymap(lambda x: str(x).replace('"', '').lower())
    return interventions

def break_by_type(interventions):
    drugs = pd.DataFrame(interventions[interventions['type'] == 'drug'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['drugs'])
    behaviours = pd.DataFrame(interventions[interventions['type'] == 'behavioral'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['behaviours'])
    procedures = pd.DataFrame(interventions[interventions['type'] == 'procedure'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['procedures'])
    devices = pd.DataFrame(interventions[interventions['type'] == 'device'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['devices'])
    dietary_supps = pd.DataFrame(interventions[interventions['type'] == 'dietary supplement'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['dietary supplements'])
    tests = pd.DataFrame(interventions[interventions['type'] == 'diagnostic test'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['diagnostic tests'])
    radiations = pd.DataFrame(interventions[interventions['type'] == 'radiation'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['radiation'])
    biologicals = pd.DataFrame(interventions[interventions['type'] == 'biological'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['biological'])
    genetics = pd.DataFrame(interventions[interventions['type'] == 'genetic'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['genetic'])
    comb_prods = pd.DataFrame(interventions[interventions['type'] == 'combination product'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['combination product'])
    others = pd.DataFrame(interventions[interventions['type'] == 'other'].sort_values(by= 'name').reset_index(drop = True)['name'].unique().tolist(), columns=['others'])
    
    return [drugs, behaviours, procedures, devices, dietary_supps, tests,radiations, biologicals, genetics, comb_prods, others]

def write_to_disk(content, name):
    dir_path = 'ctgov_sep_interventions'
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    content.to_csv(Path('.', dir_path, f'{name}.csv'), index = False)
    
    
if __name__=='__main__':
    ctgov_intervention = load_interventions()
    interventions_df_list = break_by_type(ctgov_intervention)
    names_list = ['drugs', 'behaviours', 'procedures', 'devices', 'dietary_supps', 'diagnostic_tests','radiations', 
                  'biologicals', 'genetics', 'combined_prods', 'others']
    for intervention_type, name in zip(interventions_df_list, names_list):
        write_to_disk(intervention_type, name)
    
    