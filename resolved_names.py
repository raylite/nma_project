#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:30:10 2019

@author: ja18581
"""

import requests
from requests.exceptions import HTTPError
from lxml import etree as ET
from bs4 import BeautifulSoup as bs
import pandas as pd
import numpy as np
import time

from drugs_processor import trial_reg_intervention_processor
from user_agent import userAgent


_BASE_URL = 'https://rxnav.nlm.nih.gov/REST/'

class rxnorm(userAgent):
    def __init__(self, drug_term = None):
        self.drug_name = None
        self.__get_related_info(self.drug_name)
        self.__user_agent = userAgent()
    
    def _get_approximate_terms(self):#uses the drud name
        '''
        takes a drug name as input
        
        returns a list of rxcuis of related terms
        '''
        url_params = {'term':self.drug_name, 'maxEntries':4}
        
        retry = 1
        success = False
        while not success and retry < 5:
            user_agent = self.__user_agent.get_random_ua()
            headers = {'user-agent': user_agent.strip(),
                   'referer': 'https://google.co.uk'}
            try:
                response = requests.get(_BASE_URL+'approximateTerm', 
                                        headers=headers, 
                                        params = url_params, timeout=15)
                success = True
                response.raise_for_status()
            except HTTPError as http_err:
                print(f"HTTP error occured: {http_err}")
                time.sleep(retry * 30)
                retry += 1
            except Exception as err:
                print (f'Uknown request error occured: {err}')
                time.sleep(retry *30)
                retry += 1
            else:
                return response
            
    def get_related_rxcuis(self, drug_term):
        self.drug_name = drug_term
        
        print(f"Processing for drug: {self.drug_name}")
        
        xml = self._get_approximate_terms()##get approximate terms for self.drug_name
        
        if xml.content:
            rxcui_list = self._etxract_unique_rxcuis(xml)
        else:
            time.sleep(120)
            xml = self._get_approximate_terms()
            rxcui_list = self._etxract_unique_rxcuis(xml)
        return rxcui_list

    def get_related_info(self, rxcui):
        '''
        takes a list of rxcuis
        returns list of associated drug names and dosages
        '''
        related_names_list = []
        pack_names = []
        main_contents = []            
                    
        IN = 'IN+SCDC+SCD'
        MIN = 'MIN+SCDC+SCD' #'multi-ingredient'
    
        if self.drug_name is not None:
            if '/' in self.drug_name or 'and' in self.drug_name: #its multi-ingredient
                url = _BASE_URL+'rxcui/{}/related?tty={}'.format(rxcui, MIN)
            else:
                url = _BASE_URL+'rxcui/{}/related?tty={}'.format(rxcui, IN)
            retry = 1
            success = False
            while not success and retry < 5:    
                user_agent = self.__user_agent.get_random_ua()
                headers = {'user-agent': user_agent.strip(),
                   'referer': 'https://google.co.uk'}
                try:
                    response = requests.get(url, headers=headers, timeout=15)
                    success = True
                    response.raise_for_status()
                except HTTPError as http_err:
                    print(f"HTTP error occured: {http_err}")
                    time.sleep(retry*30)
                    retry += 1
                except Exception as err:
                    print (f'Uknown request error occured {err}')
                    time.sleep(retry*30)
                    retry += 1
                else:
                    main_content, all_names, pack_doses = self._extract_related_drugs(response)
                    #print(f'Main content from processing {rxcui} are {main_content}')
                    related_names_list = related_names_list + all_names
                    pack_names = pack_names + pack_doses
                    main_contents = main_contents + main_content
                
            #print (f'All Names {related_names_list}, Pack Names: {pack_names}, main content: {main_contents}')
        
            return (list(np.unique(np.array(main_contents))), 
                    list(np.unique(np.array(related_names_list))), 
                    list(np.unique(np.array(pack_names)))
                    )
            
    __get_related_info = get_related_info   
    
    def all_display_names(self):
        url = _BASE_URL + 'displaynames'
        
        retry = 1
        success = False
        
        while not success and retry < 5:
            user_agent = self.__user_agent.get_random_ua()
            headers = {'user-agent': user_agent.strip(),
                   'referer': 'https://google.co.uk'}
            try:
                response = requests.get(url, headers=headers, timeout=15)
                success = True
                response.raise_for_status()
            except HTTPError as http_err:
                print(f"HTTP error occured: {http_err}")
                time.sleep(retry*30)
                retry += 1
            except Exception as err:
                print (f'Uknown request error occured {err}')
                time.sleep(retry*30)
                retry += 1
            else:
                all_display_names = self._extract_dispaly_names(response)
        
        return all_display_names
            
    
    def _extract_dispaly_names(self, xml):
        soup = bs(xml.content, 'xml')
        names = soup.findAll('term')
        names_list = [name.getText() for name in names]
        
        return names_list
        
    
    def _extract_related_drugs(self, xml):
        soup = bs(xml.content, 'xml')
        concepts = soup.findAll('conceptGroup')
        
        components = [name.getText() for concept in concepts if concept.find('tty').getText()=='SCDC' for name in concept.findAll('name')]
        c_synonyms = [synonym.getText() for concept in concepts if concept.find('tty').getText()=='SCDC' for synonym in concept.findAll('synonym')]
        drugs = [synonym.getText() for concept in concepts if concept.find('tty').getText()=='SCD' for synonym in concept.findAll('name')]
        d_synonyms = [synonym.getText() for concept in concepts if concept.find('tty').getText()=='SCD' for synonym in concept.findAll('synonym')]
        drug_name = [name.getText() for concept in concepts if concept.find('tty').getText()=='IN' or concept.find('tty').getText()=='MIN' for name in concept.findAll('name')]
        #dn_synonyms = [name.getText() for concept in concepts if concept.find('tty').getText()=='IN' or concept.find('tty').getText()=='MIN' for name in concept.findAll('synonym')]
        
        SCDC = components + c_synonyms #standard clinical drug components
        SCD = drugs + d_synonyms #clinical drugs/Pack
        return drug_name, SCDC, SCD
        
    def _etxract_unique_rxcuis(self, xml):
        tree = ET.XML(xml.content)
        rxcuis = None
        for tag in tree:
            rxcuis = tag.xpath(".//rxcui/text()")#make this list unique
        
        return list(np.unique(np.array(rxcuis)))
    
    def clean_terms(self, list_of_drug_terms, source = None):
        all_drug_related_terms = []
        
        for concept in list_of_drug_terms:
            if source == 'ct':
                concept = concept.split(',')
            elif source == 'rxnav':
                concept = concept.split('_')
            for term in concept:
                cui_list = self.get_related_rxcuis(term)
                print(f'CUI LIST: {cui_list}')
                for cui in cui_list:
                    if cui is None:
                        continue
                    print (f'Processing {cui} from the list {cui_list}')
                    main_content, components, drugs = self.get_related_info(cui)
                    if components or drugs:
                         all_drug_related_terms.append({'RXCUI': cui, 
                                                        'Terms': term,
                                                        'Resolved': ', '.join(main_content).strip(', '),
                                                        'SCDC': ', '.join(components).strip(', '),
                                                        'SCD': ', '.join(drugs).strip(', '),
                                                        'Source_name': concept
                                                    })
        all_drug_related_terms = pd.DataFrame(all_drug_related_terms)
    
        return all_drug_related_terms
    
    
class trial_register_handler(rxnorm):
    def __init__(self, terms_list):
        rxnorm.__init__(self)
        self.list_of_drug_terms = terms_list
        self.source = 'ct'
    
    def resolve_terms(self):
        terms_df = self.clean_terms(self.list_of_drug_terms, self.source)
        return terms_df
    
    
if __name__=='__main__':  
    dir_path ='ctgov_sep_interventions'
    intervention_type = 'drugs'
    
    #try rxnav display names as against trial reg.
# =============================================================================
#     print ("Starting on the RxNorm list")
#     rx = rxnorm()
#     
#     print ("Getting all dispay names")
#     all_names = rx.all_display_names()
#     
#     print ("Commencing the name resolution process on the RxNorm list")
#     all_terms = rx.clean_terms(all_names, 'rxnav')
#     all_terms.to_csv('/home/ja18581/nma_project/rxnav_terms.csv', index = False)
#     all_terms.to_excel('/home/ja18581/nma_project/rxnav_terms.xlsx', index = False)
#     
#     print("RxNav Completed....")
# =============================================================================
    
    trip_drugs = trial_reg_intervention_processor()
    
    drug_names_df = trip_drugs.load_drugs_list(dir_path, intervention_type)
    drug_names = trip_drugs.clean_drug_details(drug_names_df)
    drug_names.to_csv('/home/ja18581/nma_project/cleaned_ct_drug_names.csv', index = False)
    print("File saved successfully")
     
    terms = drug_names['cleaned_name']
    print ("Getting the rxcuis for approximate rxnorm concepts for the trial register phrases")
    
    drug_names = trial_register_handler(terms)
     
    print ("Commencing the name resolution process on the trial register list")
    resolved_names = drug_names.resolve_terms()
    
    resolved_names.to_csv('/home/ja18581/nma_project/ctgov_resolved.csv', index = False)
    resolved_names.to_excel('/home/ja18581/nma_project/ctgov_resolved.xlsx', index = False)
    
    print ("Completed the name resolution process on the trial register list")

    print("All process completed")
    
   