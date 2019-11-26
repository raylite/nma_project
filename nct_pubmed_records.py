#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:23:48 2019

@author: ja18581
"""
#import pandas as pd
from pubmed_downloader import PMID_Search
import MySQLdb
from sqlalchemy import create_engine
import pandas as pd

class nct_pubmed_records():
    '''
    Use a list of PMIDS to retrieve records of studies with NCT IDs from Pubmed
    
    input: pubmed_result.txt
    output: nct:pubmed_records.csv
    Fields: Title, NCT ID, PMID, Authors, Abstract
    '''
    def __init__(self):
        self.records = None
        
    def download_records(self, infile):
        with open(self.infile, 'r') as file_handle:
            self.records = PMID_Search(file_handle.readlines())
            
    def save_records(self, outfile):
        self.records.to_csv(outfile, encoding='utf8')
        
    def _create_db_connection(self):
        host = '127.0.0.1'
        user = 'root'
        password = 'p@55w0rd'
        db = 'ctgov_db'
        self._db_connection = MySQLdb.connect(host = host,
                               user = user,
                               password = password,
                               db = db,
                               use_unicode=True,
                               charset = 'utf8mb4')
        self._alchemy_engine = create_engine("mysql://{}:{}@{}/{}?charset=utf8mb4".format(user,password,host,db))
        return self._db_connection
        
    def save_records_to_db(self, tablename):
        connection = None
        try:
            connection = self._create_db_connection()
            self.records.to_sql(name = tablename, con = connection, if_exists='replace', index=False)
        except:##pandas only supports sqlalchemy
            connection = self._alchemy_engine.connect()
            self.records.to_sql(name = tablename, con = connection, if_exists='replace', index=False)
        connection.close()
        
    def load_records_from_db(self, tablename):
        query = "SELECT 'NCT Number', PMID, Title, Abstract, Authors FROM {table}".format(table=tablename)
        connection = None
        try:
            connection = self._create_db_connection()
            self.records = pd.read_sql(sql=query, con=connection, columns = ['NCT Number', 'PMID', 'Title', 'Abstract', 'Authors'])
        except:
            connection = self._alchemy_engine.connect()
            self.records = pd.read_sql(sql=query, con=connection, columns = ['NCT Number', 'PMID', 'Title', 'Abstract', 'Authors'])
            
        connection.close()


if __name__=='__main__':
    tablename = 'nct_pubmed_recordset'
    pmid_record_file_path = 'pubmed_result.txt'##list of pubmed Ids of NCT studies
    outfile = 'nct_pubmed_records.csv'
    records = nct_pubmed_records()
    
    action = input("Please choose an option. (E)xtract from Pubmed and save or (L)oad from DB: ")
    if action.lower() == 'e':
        records.download_records(pmid_record_file_path)
        records.save_records(outfile)
        records.save_records_to_db(tablename)
    elif action.lower() == 'l':
        records.load_records_from_db(tablename)
        records_df = records.records
    print("Process completed")
    