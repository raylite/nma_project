# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 11:29:10 2018

@author: raylite
"""

from Bio import Entrez
import pandas as pd
from urllib.error import HTTPError
import time
#import xml.etree.ElementTree as ET
from lxml import etree as ET

API_KEY = 'a2e2a3e33502aa03aa40735fb24c80de9008'
tool = "CiREx"


def PMID_Search(id_list):
    Entrez.email = "b.k.olorisade@bristol.ac.uk"
    Entrez.api_key = API_KEY
    Entrez.tool = tool
    
    search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    count = len(id_list)
    batch_size = 10
    abstract = None
    all_docs = []
       
    
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        if start%100 == 0 and start != 0:
            print("Going to download record {} to {}".format(start+1, end+90))
        elif start == 0: 
            print("Going to download record {} to {}".format(start+1, end+90))
            
        attempt = 0
        success = False
        while not success and attempt < 3:
            try:
                success = True
                fetch_handle = Entrez.efetch(db="pubmed",
                                             rettype="abstract", retmode="xml",
                                             retstart=start, retmax=batch_size,
                                             webenv=webenv, query_key=query_key)
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
            abstract = doc.xpath(".//AbstractText/text()")
            authorslastnames = doc.xpath(".//Author/LastName/text()")
            authorsfirstnames = doc.xpath(".//Author/ForeName/text()")#zip both, extract and append
            nctid = doc.xpath(".//AccessionNumberList/AccessionNumber/text()")
                    
            docs = {
                    "PMID": ''.join(PMID) or None,
                    "Title": ' '.join(title) or None,
                    "Abstract": ' '.join(abstract) or None,
                    "Authors": ', '.join([" ".join([last, first]) for last, first in zip(authorslastnames, authorsfirstnames)]) or None,
                    "NCT Number": ''.join(nctid) or None
                    }
            all_docs.append(docs) 
        
    return pd.DataFrame(all_docs)
