#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 13:41:15 2019

@author: ja18581
"""

##search pubmed for articles containing certain intervention of interest
from Bio import Entrez
import pandas as pd
from urllib.error import HTTPError
import time
from lxml import etree as ET

API_KEY = 'a2e2a3e33502aa03aa40735fb24c80de9008'
tool = "CiREx"

class pubmed_term_search(object):
    def __init__(self):
        Entrez.email = "b.k.olorisade@bristol.ac.uk"
        Entrez.api_key = API_KEY
        Entrez.tool = tool
        self.__list_of_ids = None
        self.__webenv = None
        self.__query_key = None
        self.__result_count = None
        
    def terms_search(self, terms):
        handle = Entrez.esearch(db="pubmed", term=terms, retmax = '1000', usehistory='y')
        record = Entrez.read(handle)
        self.__list_of_ids = record['IdList']
        
        self.__webenv = record["WebEnv"]
        self.__query_key = record["QueryKey"]
        self.__result_count = int(record['Count'])#len(self.__list_of_ids)
        
    def fetch_details(self):
        batch_size = 10
        abstract = None
        all_docs = []
           
        
        for start in range(0, self.__result_count, batch_size):
            end = min(self.__result_count, start+batch_size)
            if start%100 == 0 and start != 0:
                print(f"Going to download record {start+1} to {min(self.__result_count, end+90)}")
            elif start == 0: 
                print(f"Going to download record {start+1} to {min(self.__result_count, end+90)}")
                
            attempt = 0
            success = False
            while not success and attempt < 3:
                try:
                    success = True
                    fetch_handle = Entrez.efetch(db="pubmed",
                                                 rettype="abstract", retmode="xml",
                                                 retstart=start, retmax=batch_size,
                                                 webenv=self.__webenv, query_key=self.__query_key)
                    #try rettype=abstract, retmode = text
                except HTTPError as err:
                    attempt += 1
                    if 500 <= err.code <= 599:
                        print("Received error from server %s" % err)
                        print("Attempt %i of 3" % attempt)
                        time.sleep(15)
                    else:
                        raise
                
            data = fetch_handle.read()
            fetch_handle.close()
            tree = ET.XML(data)
            
            for doc in tree:
                
                PMID = doc.xpath(".//MedlineCitation/PMID[@Version='1']/text()") 
                title = doc.xpath(".//ArticleTitle/text()")                            
                authorslastnames = doc.xpath(".//Author/LastName/text()")
                authorsfirstnames = doc.xpath(".//Author/ForeName/text()")#zip both, extract and append
                nctid = doc.xpath(".//AccessionNumberList/AccessionNumber/text()")
                abstract = doc.xpath(".//AbstractText/text()")
                
                background = []
                method = []
                objective = []
                
                try:
                    if doc.xpath(".//Abstract/AbstractText[@NlmCategory='BACKGROUND']"):
                        background = doc.xpath(".//Abstract/AbstractText[@NlmCategory='BACKGROUND']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='BACKGROUND']"):
                        background = doc.xpath(".//Abstract/AbstractText[@Label='BACKGROUND']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='BACKGROUND:']"):
                        background = doc.xpath(".//Abstract/AbstractText[@Label='BACKGROUND:']/text()")
                    
                except:
                    background = ''
                    
                try:
                    if doc.xpath(".//Abstract/AbstractText[@NlmCategory='OBJECTIVE']"):
                        objective = doc.xpath(".//Abstract/AbstractText[@NlmCategory='OBJECTIVE']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='OBJECTIVES']"):
                        objective = doc.xpath(".//Abstract/AbstractText[@Label='OBJECTIVES']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='AIMS']"):
                        objective = doc.xpath(".//Abstract/AbstractText[@Label='AIMS']/text()")
                except:
                    objective = ''
                
                try:
                    if doc.xpath(".//Abstract/AbstractText[@NlmCategory='METHODS']"):
                        method = doc.xpath(".//Abstract/AbstractText[@NlmCategory='METHODS']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='METHODS']"):
                        method = doc.xpath(".//Abstract/AbstractText[@Label='METHODS']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='METHODS/DESIGN']"):
                        method = doc.xpath(".//Abstract/AbstractText[@Label='METHODS/DESIGN']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='METHODS AND RESULTS']"):
                        method = doc.xpath(".//Abstract/AbstractText[@Label='METHODS AND RESULTS']/text()")
                    elif doc.xpath(".//Abstract/AbstractText[@Label='MATERIALS AND METHODS:']"):
                        method = doc.xpath(".//Abstract/AbstractText[@Label='MATERIALS AND METHODS:']/text()")
                     
                except:
                    method = ''
                    
                docs = {
                        "PMID": ''.join(PMID) or None,
                        "Title": ' '.join(title) or None,
                        "Abstract": ' '.join(abstract) or None,
                        "Authors": ', '.join([" ".join([last, first]) for last, first in zip(authorslastnames, authorsfirstnames)]) or None,
                        "NCT_Number": ''.join(nctid) or None,
                        "Objective": ' '.join(objective) or None,
                        "Method": ' '.join(method) or None,
                        "Background": ' '.join(background) or None
                        }
                all_docs.append(docs) 
            
        return pd.DataFrame(all_docs)
    

    

