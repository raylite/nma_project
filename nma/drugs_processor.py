#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:07:36 2019

@author: ja18581
"""
import pandas as pd
import re
from pathlib import Path

class trial_reg_intervention_processor:
    def __init__(self, path=None, inter_type=None):
        self.path = path
        self.type = inter_type
        
    def load_drugs_list(self, folder, inter_type):
        self.path = folder
        self.type = inter_type
        drugs_list = pd.read_csv(Path('.', self.path, f'{self.type}.csv'), )
    
        return drugs_list

    def clean_drug_details(self,drugs_list):
        """
        extract drugs names, dosage and concentration
        
        :input: name of drug from trial register
        
        :output: dataframe of drug name, dosage and concentration
        """
        dosage_pattern = re.compile(r'[0-9]*(\,|\.)?[0-9]+\s*(mg|ml|grams|mcg)|[0-9]*(\,|\.)?[0-9]+\s*(mg|ml|grams|mcg)(/|per)*\s*[a-z0-9]*')
        conc_pattern = re.compile(r'^%*\s*[0-9]*\.*[0-9]+%*')
        name_pattern = re.compile('|'.join([r'(mg|ml|grams|mcg|l|\(*\+*)(/)\s*([a-z0-9]|\-*\)*)*', r'^[a-z]+\s*[a-z0-9]\s*:|^[a-z]+\s*[0-9]']))
        
        drugs_split_info = []
        
            
        for drug in drugs_list['drugs']:
            try:
                dosage = dosage_pattern.search(drug).group(0).strip()
            except:
                dosage = None
                
            name = name_pattern.sub('', drug).strip()
            
            try:
                concentration = conc_pattern.search(drug).group(0).strip()
            except:
                concentration = None
            
            drugs_split_info.append({'dose':dosage,
                                      'name': name,
                                      'concentration': concentration
                                      })
        return pd.DataFrame(drugs_split_info)

#if __name__=='__main__':
    #dir_path ='ctgov_sep_interventions'
    #intervention_type = 'drugs'
    
    
    #drug_names_df = load_drugs_list(dir_path, intervention_type)
    
    #drug_names = clean_drug_details(drug_names_df)
    