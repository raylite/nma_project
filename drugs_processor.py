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
        self.common_terms = ['day', 'days','drug', 'drugs', 'drug:', 'complete','dose', 'doses','dosage','capsule','tablet'
                             'syrup','solution','spray', 'nasa', 'iv', 'injection', 'year', 'day', 'days', 'years','capsules'
                             'month', 'months', 'subject', 'sujects', 'patient', 'patients', 'parent', 'parents', 'male',
                             'female', 'man', 'men', 'woman', 'women', 'child', 'children', 'adult', 'old', 'age', 'once',
                             'twice', 'thrice', 'one', 'two', 'three', 'time', 'times', 'morning', 'afternoon', 'night',
                             'administartion', 'administered','placebo', 'week', 'weeks', 'cycles', 'of', 'at', 'with',
                             'in', 'treatment', 'treatments', 'per', 'hour', 'hours', 'the', 'gel', '(', ')','phase',
                             'stage', 'line','art','period', 'skin', 'nasal', 'over', 'minutes', 'weekly', 'week', 
                             'follow', 'followed', 'arm', 'cream', 'inhibitor', 'intravenous', 'daily', 'body', '[', ']',
                             'group', 'groups', 'grouped', 'w/v', 'prescribed', 'post-operatively', 'prescribed', 'blocking',
                             'blinded', 'blind', 'blinding', 'to', 'type', 'single', 'high', 'single', 'topical', 'drop',
                             'take', 'taken','in', 'soak', 'cycle', 'cycles', 'supplement', 'extract', 'other', 'than',
                             'agent', 'primary', 'agent', 'agents', 'factor', 'factors', 'continued', 'continued','therapy',
                             'therapeutic', 'fix', 'fixed', 'combination', 'will', 'be', 'assign', 'assigned', 'on','level',
                             'levels', 'use', 'using', 'used', 'base', 'based']
    
        return drugs_list

    def clean_drug_details(self,drugs_list):
        """
        extract drugs names, dosage and concentration
        
        :input: name of drug from trial register
        
        :output: dataframe of drug name, dosage and concentration
        """
        dosage_pattern = re.compile(r'(((\d*\.*\d\s*-)|\.)*\d+((\,|\.)\d+)?\s*\-*(l|min|mg|ml|grams|mcg|ug|ng|gm|cc|g|iu|kg|mcgr|units|µg|ppm|milliliters))\s*(/|per)*\s*\d*(day|hour|m2|kg|mg|h|ml|g|qd|μl)*')
        conc_pattern = re.compile(r'(\d*((\,|\.)?\d+)\s*%)|%\s*(\d*((\,|\.)?\d+))')
        misc = re.compile(r'\(((\+|-|/)*|(\d*[a-z]{1}-*))\)-*')
        drug_prefix = re.compile(r'[a-z]+\s*[a-z0-9]\s*:')
        all_patterns = [dosage_pattern,conc_pattern, misc,drug_prefix]
        name_pattern = re.compile('|'.join(p.pattern for p in all_patterns))
        separator = re.compile(r'\b\s*(&|plus|and|or|;|/|\+|\,|vs)[^/(day|hour|m2|kg|mg|h|ml|g|0-9|/|\-|m)]')
        
        drugs_split_info = []
        separators_list = ['&', 'plus', 'and', 'or', ';', '/', '+', ',', 'vs']
        
        for drugs in drugs_list['drugs']:
            #print(f'DRUGS: {drugs}')
            split_list = separator.split(drugs)
            #put if statemetn to ignore the separators from list
            for drug in split_list:
                if drug.strip() in separators_list:
                    continue
                #print(f'DRUG: {drug}')
                try:
                    dosage = dosage_pattern.search(drug).group(0).strip()
                except:
                    dosage = None                
                try:
                    concentration = conc_pattern.search(drug).group(0).strip()
                except:
                    concentration = None
                drugl = drug.split()
                drugl = ([term for term in drugl if term.lower() not in self.common_terms])
                drugl = ' '.join(drugl)
                name = name_pattern.sub('', drugl).strip()
                
                #print(f'NAME: {name}; CONC: {concentration}; DOSE: {dosage}')
                
                if name and len(name)>1:
                    drugs_split_info.append({'dose':dosage,
                                          'cleaned_name': name,
                                          'source_name': drugs,
                                          'split_name': drug,
                                          'concentration': concentration
                                          })
        return pd.DataFrame(drugs_split_info)

# =============================================================================
# if __name__=='__main__':
#     dir_path ='ctgov_sep_interventions'
#     intervention_type = 'drugs'
#     
#     drug = trial_reg_intervention_processor()
#     
#     drug_names_df = drug.load_drugs_list(dir_path, intervention_type)
#     
#     drug_names = drug.clean_drug_details(drug_names_df)
# =============================================================================
    