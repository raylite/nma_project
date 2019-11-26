#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:26:58 2019

@author: ja18581
"""

'''
Module to load records save in the ctgov_db as a dataframe
input: dbname, tablename
output: records dataframe
'''
import MySQLdb
import pandas as pd

def connect_db(dbname):
    try:
        connection = MySQLdb.connect(host = '127.0.0.1', 
                           user='root', 
                           password='p@55w0rd', 
                           db = dbname,
                           charset = 'utf8mb4',
                           use_unicode = True)
    except MySQLdb.Error as err:
        print("Database {} does not exists. Error {} raised".format('ctgov_db', err))
        
    return connection


def load_records(connection, tablename):
    cursor = connection.cursor()
    query = "SELECT nct_id, study_title, description, study_type, study_design, conditions, \
    interventions, eligibility, age, gender, arms FROM {} LIMIT 20".format(tablename)

    cursor.execute(query)
    records = cursor.fetchall()
    return records



if __name__=='__main__':
    dbname = 'ctgov_db'
    table = 'trial_register'
    
    db_conn = connect_db(dbname)
    records = pd.DataFrame(list(load_records(db_conn, table)), 
                           columns=["ID", "Title", "Description", "Type","Design", "Conditions", "Interventions",
                                    'Eligibility Criteria', 'Age', 'Gender', 'Arms'])
    
    db_conn.close()