#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 12:29:49 2019

@author: ja18581
"""

'''
Extraxt the locally stored XML records for all studies on the clinicaltrials.gov and save them in a mysql 
database.
dbname: ctgov_db, user: root, password:---
table name: trial_register
columns: nct_id, study_title, description, study_type, study_design, conditions, interventions,
          eligibility, age, gender, arms
'''
    

import os
import xml.etree.ElementTree as ET
import MySQLdb
import json

root = '/home/ja18581/ctgov_records'
count = 0

def iterfolder(dirname):
    for _, _, files in os.walk(dirname):
        for afile in files:
            yield afile

def readxml(rootdir):
    for afolder in os.scandir(rootdir):#
        list_of_records = []
        if afolder.is_dir():
            #fileiterator = iterfolder(afolder.path)
            
            for xmlfile in iterfolder(afolder.path):
                print ("Processing: ", os.path.join(afolder.path, xmlfile))
                tree = ET.parse(os.path.join(afolder.path, xmlfile))
                root = tree.getroot()
                list_of_records.append(parsexml(root))
            
            write_to_db(list_of_records)
    return "Success"
                

def parsexml(tree_root):
    print("Parsing")
    if tree_root.findall('.//nct_id'):
        nctid = tree_root.find('.//nct_id').text
    else:
        nctid = None
    
    if tree_root.findall('.//brief_title'):
        title = tree_root.find('.//brief_title').text 
    else:
        title = None
    
    if tree_root.findall('./detailed_description/textblock'):
        descr = tree_root.find('./detailed_description/textblock').text.strip()
    else:
        descr = None
    
    if tree_root.findall('.//study_type'):
        stype = tree_root.find('.//study_type').text 
    else:
        stype = None
    
    if tree_root.findall('./study_design_info'):
        design = json.dumps([{sd.tag: sd.text for sd in design} for design in tree_root.findall('./study_design_info')])
    else:
        design = None
    
    if tree_root.findall('.//condition'):
        condition = tree_root.find('.//condition').text 
    else:
        condition = None
    
    if tree_root.findall('./intervention'):
        interventions = json.dumps([{inter.tag: inter.text for inter in inters} for inters in tree_root.findall('./intervention')]) 
    else:
        interventions = None
    
    if tree_root.findall('./eligibility/criteria/textblock'):
        eligibility = tree_root.find('./eligibility/criteria/textblock').text
    else:
        eligibility = None
    
    if tree_root.findall('.//minimum_age'):
        if tree_root.findall('.//maximum_age'):
            age = json.dumps([{tree_root.find('.//minimum_age').tag: tree_root.find('.//minimum_age').text,
                   tree_root.find('.//maximum_age').tag: tree_root.find('.//maximum_age').text
                  }])
        else:
            age = {tree_root.find('.//minimum_age').tag: tree_root.find('.//minimum_age').text}
    else:
        age = None
    
    if tree_root.findall('.//gender'):
        gender = tree_root.find('.//gender').text
    else:
        gender = None
    if tree_root.findall('./arm_group'):
        arms = json.dumps([{arm.tag:arm.text for arm in arms} for arms in tree_root.findall('./arm_group')]) 
    else:
        arms = None
    
    return tuple([nctid, title, descr, stype, design, condition, interventions, eligibility, 
                       age, gender, arms]) 
def create_database(dbname):
    print("Creating DB for the first time")
    try:
        connection = MySQLdb.connect(host = '127.0.0.1', 
                           user='root', 
                           password='p@55w0rd')
        cursor = connection.cursor()
        cursor.execute("CREATE DATABASE {} CHARSET utf8mb4 COLLATE utf8mb4_unicode_ci".format(dbname))
        return connection
    except MySQLdb.Error as err:
        print (f"Failed to create database: {err}")
        exit(1)
            
def write_to_db(parsedxmllist):
    connection = None
    try:
        connection = MySQLdb.connect(host = '127.0.0.1', 
                           user='root', 
                           password='p@55w0rd', 
                           db = 'ctgov_db',
                           use_unicode=True,
                           charset = 'utf8mb4')
    except MySQLdb.Error as err:
        print("Database {} does not exists. Error {} raised".format('ctgov_db', err))
        #if err.errno == errorcode.ER_BAD_DB_ERROR:
        connection = create_database('ctgov_db')
        cursor = connection.cursor()
        cursor.execute("USE ctgov_db")
    
    
    try:
        table_query = 'CREATE TABLE IF NOT EXISTS trial_register (id INT AUTO_INCREMENT PRIMARY KEY, \
        nct_id TEXT, study_title  MEDIUMTEXT, description LONGTEXT, study_type TEXT, study_design JSON,\
        conditions TEXT, interventions JSON, eligibility LONGTEXT, age JSON, gender TEXT, \
        arms JSON)'
        cursor = connection.cursor()
        cursor.execute(table_query)
    except:
        print ("Table already exist.")
    
    insert_query = 'INSERT INTO trial_register(nct_id, study_title, description, study_type,\
                                                 study_design, conditions, interventions,\
                                                 eligibility, age, gender, arms) \
                      VALUES (%s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s)'
    global count
    
    try:
        count += 1
        with connection.cursor() as cursor:
            for record in parsedxmllist:
                cursor.execute(insert_query, record)
            connection.commit()
        print(f"Insert {count} Successful")
    finally:
        connection.close()

if __name__=='__main__':
    print ("Begin Processing...")
    operation = readxml(root)
    print (f"Opreation concluded with {operation}")