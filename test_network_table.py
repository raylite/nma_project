#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 13:06:44 2019

@author: ja18581
"""

import unittest
import HtmlTestRunner
import pandas as pd
from pathlib import Path
import os


from network_table import buildTable

#test data file

def initialize_test_data():
    global test_data_dir, test_drug_file, test_article_file, drugs_dir, all_drugs
    
    test_data_dir = 'tests'
    test_drug_file = 'drug.csv'
    test_article_file = 'publication.csv'
    
    print("initialization of test data for file handling complete.")
    
def clear_test_data():
    global test_data_dir, test_drug_file, test_article_file
    if test_data_dir is not None:
        test_data_dir = None
        print("Test data dir set to None")
        
    if test_drug_file is not None:
        test_drug_file = None
        print("Drug file name set to None.")
        
    if test_article_file is not None:
        test_article_file = None
        print("Article file name set to None.")
        

class TestNetworkTable(unittest.TestCase):
    build_table =None
    
    @classmethod
    def setUpClass(cls):
        initialize_test_data()
        print (' ')
        print(f'TestNetworkTable setUpClass')
        
    @classmethod
    def tearDownClass(cls):
        clear_test_data()
        print (' ')
        print(f'TestNetworkTable tearDownClass')
        
    def setUp(self):
        self.build_table = buildTable(test_data_dir, test_drug_file)
        print (' ')
        print(f'TestNetworkTable setUp')
        
    def tearDown(self):
        if self.build_table is not None:
            self.build_table = None
        print(f'TestNetworkTable tearDown')
        
        
    def test_term_driven_search(self):
        print('')
        print("****test_term_driven_search******")
        
        result = self.build_table.term_driven_search()
        print(f'loading drug list from dir. Result after operation = {result} expected = {True}')
        self.assertEqual(result, True)
        print(f'ensuring output file exists = {os.path.exists(Path(".", test_data_dir, "drug_article_compare.csv"))}, expected = True')
        self.assertEqual(os.path.exists(Path('.', test_data_dir, 'drug_article_compare.csv')), True)
        
    def test_article_driven_search(self):
        pass
        
        

def build_test_suite():
    # create unittest.TestSuite object.
    test_suite = unittest.TestSuite()
    # add each test function to the test suite object.
    test_suite.addTest(TestNetworkTable('test_term_driven_search'))
    test_suite.addTest(TestNetworkTable('test_article_driven_search'))
    
    return test_suite

def build_text_report():
    test_suite = build_test_suite()
    # create unittest.TextTestRunner() object.
    test_runner = unittest.TextTestRunner()
    # run the test suite.
    test_runner.run(test_suite)
    # run below code to run all test function
    # unittest.main()

        
    # generate html report.
def build_html_report():
    test_suite = build_test_suite()
    test_runner = HtmlTestRunner.HTMLTestRunner(output=Path('.', 'tests', 'html_report'))
    test_runner.run(test_suite)

if __name__ == "__main__":
    build_text_report()
    build_html_report()
    
    