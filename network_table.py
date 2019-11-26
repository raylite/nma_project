#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 17:07:44 2019

@author: ja18581
"""
import pandas as pd
from pathlib import Path

from intervention_search import searchIntervention 
##create a dataframe of drugs used in a study against each other
#create on drug basis
#create of study basis


class buildTable():
    siObject = searchIntervention()
    def __init__(self, data_dir, filename):
        self.__rxfile_name = filename #name of file containing drug names
        self.__data_dir = data_dir
        self.__rx_based_drugs = self.siObject.load_drugs_list('data', 'rxnav_terms.csv')
        self.__ct_based_drugs = self.siObject.load_drugs_list('data', 'ctgov_resolved.csv')
        self.__outpath = Path('.', self.__data_dir, 'drug_article_compare.csv')
        self.__drugs_list = self.siObject.load_drugs_list(self.__data_dir, self.__rxfile_name) #a dataframe of drug names
    
    def term_driven_search(self):
        self.__drugs_list.fillna('', inplace=True)
        
        self.__drugs_list = self.__drugs_list.groupby('RXCUI', as_index=False).agg({'Terms':lambda x: ', '.join(x.unique()),
                                  'SCDC': lambda x: ', '.join(x.unique()),
                                  'SCD': lambda x: ', '.join(x.unique()),
                                  'Resolved': lambda x: ', '.join(x.unique())})
    
        for index, row in self.__drugs_list.iterrows():
            adrug = row['Terms']#name to search
            main_term = row['Resolved']#name for filtering. There is need to harmonize both names
                  
            related_publications = self.siObject.search_pubmed(adrug)
            
            compared_drugs = self.siObject.search_codrug(main_term, self.__drugs_list, related_publications)
            
            compared_drugs.to_csv(self.__outpath, index = False, mode = 'a')
        return True
            
    
    def article_driven_search(self):
        '''
        Take a dataframe of articles, iterate over one by one.
        for each one, check which of the drug terms appear therein
        '''
        #iterate publications dir and process
        self.__drugs_list = self.__drugs_list.groupby('RXCUI', as_index=False).agg({'Terms':lambda x: ', '.join(x.unique()),
                                  'SCDC': lambda x: ', '.join(x.unique()),
                                  'SCD': lambda x: ', '.join(x.unique()),
                                  'Resolved': lambda x: ', '.join(x.unique())})
    
        data_dir = Path('.', 'data', 'pubmed').glob('**/*.csv')
        for afile in data_dir:
            publications = pd.read_csv(afile)
            publications_w_drugs = self.siObject.search_codrug(None, self.__drugs_list, publications)
            publications_w_drugs.to_csv(afile, index=False, mode='o')
        return True


if __name__=='__main__':
    data_dir = 'data' 
    filename = 'drug.csv'#to be changed as appropriate
    interv_search = buildTable(data_dir, filename)
    
      
    #publications = interv_search.term_driven_search()
    
    publications = interv_search.article_driven_search()#for the direct pubmed all articles
    
    ##proceed with term  secondary terms search across the abstracts to know which other drugs is used alongside
    #the main term...also get version of the main term used.
    