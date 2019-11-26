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
from sqlalchemy import create_engine
import logging
from logging.config import fileConfig
from pathlib import Path
from tqdm import tqdm
import csv
from multiprocessing.pool import ThreadPool
from multiprocessing import Value


API_KEY = 'a2e2a3e33502aa03aa40735fb24c80de9008'
tool = "CiREx"

all_rcts = '''(randomised[All Fields] OR ("random allocation"[MeSH Terms] OR ("random"[All Fields] AND "allocation"
[All Fields]) OR "random allocation"[All Fields] OR "randomized"[All Fields])) AND ("clinical trials as topic"
[MeSH Terms] OR ("clinical"[All Fields] AND "trials"[All Fields] AND "topic"[All Fields]) OR "clinical trials as 
topic"[All Fields] OR "trial"[All Fields] OR "placebo"[All Fields])'''

ti_rcts = '''(randomised[Title] OR randomized[Title]) AND trial[Title]'''

tiabs_rcts = '''(("randomised"[Title/Abstract] OR "randomized"[Title/Abstract]) AND"trial"[Title/Abstract]) '''



class pubmed_term_search(object):    
    def __init__(self, partial_list = None):
        Entrez.email = "b.k.olorisade@bristol.ac.uk"
        Entrez.api_key = API_KEY
        Entrez.tool = tool
        self.__webenv = None
        self.__query_key = None
        self.__result_count = None
        self.__db_engine = create_engine('mysql://root:p@55w0rd@127.0.0.1/pubmed_db?charset=utf8mb4&binary_prefix=true')
        self.__rounds = Value('i', 0)
        self.__batch_size = 10000
        self.__process_list = []
        self.__out_file = Path('.', 'logs', 'processed.txt')
        self.__partial_list = partial_list
        fileConfig('logging_config.ini')
        self.logger = logging.getLogger()
        self.dir = Path('.', 'data', 'pubmed')
        if not self.dir.exists():
            self.dir.mkdir(parents=True)
        self.string = ti_rcts
        
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        conn = self.db_engine.connect()
        conn.close()
        self.logger.debug('DB closure ensured')
        
    def search_all(self):
        search_string = self.string #possibility to download almost all records. 28million+ as at 25/7/19
        handle = Entrez.esearch(db="pubmed", term=search_string, mindate='1700', maxdate='2019',
                                retmax = self.__batch_size, usehistory='y')
        record = Entrez.read(handle)
        
        self.__webenv = record["WebEnv"]
        self.__query_key = record["QueryKey"]
        self.__result_count = int(record['Count'])
        self.logger.info(f'\nCOUNT: {self.__result_count}')
        self.__start = 0
        ##call over and over adding 1 to retmax to set restart value
        if self.__partial_list:
            print (f'''Working in partial mode with {len(self.__partial_list)} indices already processed 
                    and {(self.__result_count//self.__batch_size) - len(self.__partial_list)} to be processed''')
            self.__rounds.value = len(self.__partial_list)
        else:
            print (f'Working in full mode with {self.__result_count//self.__batch_size} indices to process')
    
    def retrieve_batch(self, start_index):#candidate for multithreading
        end = min(self.__result_count, start_index+self.__batch_size)
        if start_index%self.__batch_size == 0 and start_index != 0:
            self.logger.info(f"Downloading records {start_index+1} to {min(self.__result_count, end)}")
        elif start_index == 0: 
            self.logger.info(f"Downloading records {start_index+1} to {min(self.__result_count, end)}")
            
        attempt = 1
        success = False
        while not success and attempt < 5:
            try:
                fetch_handle = Entrez.efetch(db="pubmed",
                                             rettype="abstract", retmode="xml",
                                             retstart=start_index, retmax=self.__batch_size,
                                             webenv=self.__webenv, query_key=self.__query_key)
                success = True
            except HTTPError as err:
                self.logger.debug(f"Received error from server {err}")
                self.logger.info(f"Attempt {attempt} of 5 for reading records with starting index: {start_index}")
                attempt += 1
                success = False
                time.sleep(attempt * 15)                
            else:
                data = fetch_handle.read()
                fetch_handle.close()
                tree = ET.XML(data)
                self.extract_batch_data(tree, start_index)
        #return tree
        
    def extract_batch_data(self, tree, index):#candidate for multiprocessing
        abstract = None
        all_docs = []
        for doc in tree:
            PMID = doc.xpath(".//MedlineCitation/PMID[@Version='1']/text()") 
            title = doc.xpath(".//ArticleTitle/text()")                            
            authorslastnames = doc.xpath(".//Author/LastName/text()")
            authorsfirstnames = doc.xpath(".//Author/ForeName/text()")#zip both, extract and append
            nctid = doc.xpath(".//AccessionNumberList/AccessionNumber/text()")
            abstract = doc.xpath(".//AbstractText/text()") #see if oyu can exclude <AbstractText Label="Trial Registration" NlmCategory="UNASSIGNED">ClinicalTrials.gov identifier: NCT00331773.</AbstractText>
            ptype = doc.xpath(".//PublicationType/text()")
            
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
                    "Background": ' '.join(background) or None,
                    "Ptype": ptype or None
                    }
            
            all_docs.append(docs) 
        pub_set = pd.DataFrame(all_docs)
        self.__save_details(pub_set, index)
        #return pub_set
        
        
    def fetch_details(self):
        cycle = range(self.__start, self.__result_count, self.__batch_size)
        if self.__partial_list is not None:
            cycle = list(set(cycle) - set(self.__partial_list))
        pool = ThreadPool(processes=5)
        pool.map(self.retrieve_batch, cycle)
        pool.close()
# =============================================================================
#         for start in tqdm(cycle):
#             #retrieve batch block
#             print (f'START: {start}')
#             tree = self.retrieve_batch(start)
#             
#             pub_set = self.extract_batch_data(tree)
#                 #extraxt batch data block
#             
#             self.__save_details(pub_set)
#         
#         return True
# =============================================================================
    
    def __save_details(self, publications, index):
        filename_sufix = index//self.__batch_size
        path = Path(self.dir, f'rct_{filename_sufix}_{self.__rounds.value}.csv')
        with self.__rounds.get_lock():
            try:
                #publications.to_sql('articles', self.__db_engine, if_exists='append',  index = False)
                publications.to_csv(path, index = False)
            except:
                print(f"Data save attempt for index {index} failed")
            else:
                self.__rounds.value += 1
                with open(self.__out_file, 'a+') as f:
                    if self.__rounds.value <= 1:
                        f.write(f'{index}')
                    else:
                        f.write(f', {index}')
            self.logger.info(f'Details for start index {index} processed')

    
if __name__=='__main__':
    print("Starting...")
    fileConfig('logging_config.ini')
    logger = logging.getLogger()
    try:
        file_content = csv.reader(open(Path('.', 'logs', 'processed.txt'), newline=''), delimiter=',')
    except:
        processed = None
    else:
        processed = [int(item.strip()) for row in file_content for item in row]
        
    pubmed = pubmed_term_search(processed)
    
    pubmed.search_all()
    pubmed.fetch_details()
    print('Completed')
    
    del(pubmed)
