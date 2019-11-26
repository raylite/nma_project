#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:59:57 2019

@author: ja18581
"""

import pandas as pd
from pathlib import Path
import re
from datetime import datetime
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from nltk.text import TextCollection
import string


#local import
from pubmed_term_search import pubmed_term_search

class searchIntervention(pubmed_term_search):
    def __init__(self, main_intervention=None, drugs_filename=None):
        self.__drugs_filename = drugs_filename
        #self.__drugnames = None
        self.__pubmed_search = pubmed_term_search()
        
    def load_drugs_list(self, data_dir, filename):
        self.__drugs_filename = filename
        drugnames = pd.read_csv(Path('.', data_dir, self.__drugs_filename))
        return drugnames
        
    def search_pubmed(self, drug_term, start_year= None, end_year=None, yrange=True):
        cur_year = datetime.today().year
        if not start_year and not end_year:
            search_string = f'{drug_term} AND randomized controlled trial [pt] OR controlled clinical trial [pt] OR randomized [tiab] \
            OR placebo [tiab] OR drug therapy [sh] OR randomly [tiab] OR trial [tiab] OR groups [tiab]'
        elif start_year and end_year:
            search_string = f'{drug_term} AND randomized controlled trial [pt] OR controlled clinical trial [pt] OR randomized [tiab] \
            OR placebo [tiab] OR drug therapy [sh] OR randomly [tiab] OR trial [tiab] OR groups [tiab] AND {start_year}[pdat]: {end_year} [pdat]'
        elif start_year and not end_year and yrange:
            search_string = f'{drug_term} AND randomized controlled trial [pt] OR controlled clinical trial [pt] OR randomized [tiab] \
            OR placebo [tiab] OR drug therapy [sh] OR randomly [tiab] OR trial [tiab] OR groups [tiab] AND {start_year}[pdat]:{cur_year}'
        elif start_year and not end_year and not yrange:
            search_string = f'{drug_term} AND randomized controlled trial [pt] OR controlled clinical trial [pt] OR randomized [tiab] \
            OR placebo [tiab] OR drug therapy [sh] OR randomly [tiab] OR trial [tiab] OR groups [tiab] AND {start_year} [pdat]' 
        self.__pubmed_search.terms_search(search_string)
        search_results = self.__pubmed_search.fetch_details()
        
        return search_results
    
    def check_variant(self, main_interv_term, drugs_list, publications):
        #compared = {}
        drugs_list.fillna('', inplace=True)
        publications.fillna('', inplace=True)
                            
        for index, study in publications.iterrows():
            if main_interv_term.lower() in study['Abstract'].lower():
                publications.at[index, 'main'] = main_interv_term
            else:
                publications.at[index, 'main'] = ''
        return publications
        
    
    def search_codrug(self, main_interv_term, drugs_list, publications, terms_driven = True):
        drugs_list.fillna('', inplace=True)
        #publications.fillna('', inplace=True)
        
        
        if terms_driven:
            new_list = drugs_list[drugs_list['Resolved'].str.lower() != main_interv_term.lower()]
        else:
            new_list = drugs_list
        #make unique
        new_list = new_list.groupby('Resolved', as_index = False).agg(lambda x: ','.join(x.unique()))#may have to group  y rxcui in future
        for index, study in publications.iterrows():
            d_list = []
            for drug in new_list['Resolved']:
                #print(f'Current drug: {drug}')
                try:
                    if study['Objective'] and drug != '':
                        if re.search(r'\b{}\b'.format(drug.lower()), study['Objective'].lower()):
                            print (f'Drug: {drug} found in Objective')
                            d_list.append(drug)
                        elif re.search(r'\b{}\b'.format(drug.lower()), ', '.join([study['Objective'].lower(), study['Background'].lower()])):
                             print (f'Drug: {drug} found in Objective + Background')
                             d_list.append(drug)
                    elif drug !='' and study['Method']:
                        if re.search(r'\b{}\b'.format(drug.lower()), study['Method'].lower()):
                            print (f'Drug: {drug} found in Method')
                            d_list.append(drug)     
                    elif re.search(r'\b{}\b'.format(drug.lower()), re.split(r'RESULTS', study['Abstract'])[0].lower()) and drug != '':
                        print (f'Drug: {drug} found in Abstract0')
                        d_list.append(drug)
                    elif re.search(r'\b{}\b'.format(drug.lower()), study['Abstract'].lower()) and drug != '' and 'RESULTS' not in study['Abstract']:
                        print (f'Drug: {drug} found in Abstract1')
                        d_list.append(drug)
                except:
                    continue
            
            if terms_driven:
                publications.at[index, 'drug'] = main_interv_term
                publications.at[index, 'co_drug'] = ','.join(d_list)
            else:
                publications.at[index, 'co_drug'] = ','.join(d_list)
            
        return publications
    
    def __tokenize(self, doc):
        common_terms = ['randomized', 'randomize', 'randomised', 'randomise', 'random', 'clinical',
                      'research', 'trial', 'trials', 'affect', 'also', 'control', 'controlled', 'effect', 'human', 
                      'being', 'human being', 'humans','child', 'age', 'young', 'old', 'male', 'female', 'males', 
                      'females', 'child', 'children', 'condition', 'clinic', 'result', 'results', 'clinicaltrials.gov',
                      'gov', 'maintain', 'maintained','treatment', 'treatments', 'care', 'cares', 'self','people', 
                      'patient', 'patients', 'man', 'woman', 'men', 'women', 'use', 'used', 'individual', 'individuals', 
                      'person', 'persons', 'author','authors', 'outpatient', 'outpatients', 'conclusion', 'conclusions',
                      'introduction', 'method', 'methods', 'background', 'result', 'results', 'outcome', 'outcomes',
                      'clinical', 'control', 'controlled', 'conditions', 'study', 'compare', 'compared', 'comparison',
                      'objective', 'objectives', 'aim', 'aims', 'enrol', 'enrolled', 'different', 'difference', 'need',
                      'needed', 'with', 'without', 'randomisation', 'randomization', 'compare', 'compared' 'approach', 
                      'approached', 'strategy', 'strategies']
        stop_words = set(stopwords.words('english') + common_terms + list(string.punctuation))
        for token in word_tokenize(doc):
            if token.lower().isdigit() or token.lower() in stop_words or token.isnumeric() or len(token) <= 3:
                continue
            yield token.lower()
            
    def __vectorize(self, corpus):
        corpus = [list(self.__tokenize(doc)) for doc in corpus]
        
        texts  = TextCollection(corpus) 
    
        for doc in corpus:
            yield{
                    term: texts.tf_idf(term, doc)
                    for term in doc
                    }
        
    
    def tfidf_terms_rank(self, publications):
        temp_list = []
        for item in self.__vectorize(publications):
            temp_list.append(item)
        corpus_df = pd.DataFrame(temp_list)
        return corpus_df
    
if __name__=='__main__':
    pub = pd.read_csv(Path('.','tests', 'publication.csv'))
    
    si = searchIntervention()
    c = si.tfidf_terms_rank(pub['Abstract'])
    #iterate over each row to sort and selct top per document e.g 
    #h = c.iloc[5]
    #h = h.sort_values(ascending=False).reset_index(name='scores').rename(columns={'index':'drugs'})
    # interest  = h.iloc[...]



    
    