#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 14:03:11 2019

@author: ja18581
"""


'''
Module to extract only the interventions from the ctgov_db
input: dbname, tablename
output: interventions.csv
'''
import pandas as pd
import MySQLdb
import json

class retrv_interventions():
    '''
    Module to load and reformat only the intervetions from the ctgov_db
    '''
    def __init__(self, dbname=None, mtable=None, itable=None):
        '''
        :param dbname: path to the database
        :type dbname: str
        :param table: name of the database table
        :type table: str
        :returns: instance of the class
        '''
        self.db = dbname
        self.main_table = mtable
        self.interventions_table = itable
        self.user = 'root'
        self.password = 'p@55w0rd'
        self._db_connection = self.__connect_db()
        self._db_cursor = self._db_connection.cursor()
        
    #def __del__(self):
        #self._db_connection.close()
    
    def __connect_db(self):
        '''
        Utility to connect to the database
        :raises :class: MySQLdb.Error: Database does not exist
        '''
        try:
            return MySQLdb.connect(host = '127.0.0.1', 
                           user=self.user, 
                           password=self.password, 
                           db = self.db,
                           charset = 'utf8mb4',
                           use_unicode = True)
        except MySQLdb.Error as err:
            print(f"Database {self.db} does not exists. Error {err} raised")
        

    def load_interventions(self):
        '''
        :returns: a dataframe of interventions. Column: [NCT ID, type, name, description]
        '''
        query = "SELECT nct_id, interventions FROM {}".format(self.main_table)  
        with self._db_cursor as cursor:
            cursor.execute(query)
            intervention_records = cursor.fetchall()
                
            interventions_list = self.__format_interventions(intervention_records)
            self.save_interventions(interventions_list)
        return ('Success')
        
        
    def __format_interventions(self, intervention_records):
        intervention_list = []
        for study_record in intervention_records:
            if study_record[1] is not None: #check intervention is not None
               study_intervention  = json.loads(study_record[1])
               for intervention_record in study_intervention:
                   interventiion_dict = {'description':intervention_record.get('description'),#using .get returns None if key unavailable instead raising keyerror
                                         'type':intervention_record.get("intervention_type"), 
                                         'name':intervention_record.get("intervention_name"),
                                         'nct_id':study_record[0]
                                         }
                   intervention_list.append(interventiion_dict)
        return intervention_list
    
    def save_interventions(self, records_list, table=None):#
        #set vale for tablename
        if not self.interventions_table and table is not None:
            self.interventions_table = table
        elif not self.interventions_table and table is None:
            self.interventions_table = input("Please provide a name for the interventions table to create: ")
            
        #if first time, create table
        query = 'SHOW TABLES LIKE "{}"'.format(self.interventions_table)
        print(query)
        self._db_cursor.execute(query)
        result = self._db_cursor.fetchone()
        if result:
            print(f"Table {self.interventions_table} already exists. Bypassing creation..")
            
        else:
            query = "CREATE TABLE IF NOT EXISTS {} (id INT AUTO_INCREMENT PRIMARY KEY, nct_id TEXT,\
            name TEXT, type TEXT, description LONGTEXT)".format(self.interventions_table)
            
            
            self._db_cursor.execute(query)
            print(f"Table {self.interventions_table} successfully created")
        
        #save records in db
        with self._db_cursor as cursor:
            columns_names = ['nct_id', 'name', 'type', 'description']
            columns_names_str = ', '.join(columns_names)
            values_str = ', '.join('%s' for _ in range(len(columns_names)))
            
            for record_dict in records_list:
                query = 'INSERT INTO {table_name} ({columns}) \
                      VALUES ({values_str})'.format(table_name = self.interventions_table,
                      columns = columns_names_str, values_str = values_str) 
                values = [record_dict[column_name] for column_name in columns_names]
                print(query, values)
                cursor.execute(query, values)
            self._db_connection.commit()
            print(f'Interventions successfuly inserted')
            
            
    def load_formatted_interventions(self, tablename):
        column_names = ['nct_id', 'name', 'type', 'description']
        column_names_str = ', '.join(column_names)
        self.interventions_table = tablename
        
        query = "SELECT {columns} FROM {table}".format(table=self.interventions_table, 
                        columns=column_names_str) 
        try:
            self._db_cursor.execute(query)
            intervention_records = self._db_cursor.fetchall()
        except:
            print(f'Make sure that table {self.interventions_table} exists')
        
        
        interventions_df = pd.DataFrame(list(intervention_records), columns=column_names)
        return interventions_df
                


if __name__=='__main__':
    db_name = 'ctgov_db'
    full_table = 'trial_register'
    interventions_table = 'interventions'
    
    interventions = retrv_interventions(db_name, full_table, interventions_table)
    
    option = input("Specify the first character to select an option of which action you want to take. (E)xtraxt from main record/(L)oad formatted records: ")
        
    if option.lower() == 'e':
        intervention_records = interventions.load_interventions()
        print(f'Data extracted and saved to db. Execution exited with code {intervention_records}..')
    elif option.lower() == 'l':
        nct_interventions = interventions.load_formatted_interventions(interventions_table)
        print(f'Records loaded successfully')
        
    
    interventions._db_connection.close()#close database connection

    
    
