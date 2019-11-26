#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 09:58:59 2019

@author: ja18581
"""

import unittest
import HtmlTestRunner
import pandas as pd
from pathlib import Path


from intervention_search import searchIntervention

#test data file
test_data_path = Path('.', 'tests', 'publication.csv')

def load_test_data():
    global test_data_file_object, test_data
    
    test_data_file_object = open(test_data_path, 'r')
    test_data = pd.read_csv(test_data_file_object)
    print("open and load test data as dataframe complete.")
    
def close_test_data_file():
    global test_data_file_object
    if test_data_file_object is not None:
        test_data_file_object.close()
        test_data_file_object = None
        print("Close file publication.csv complete.")

class TestInterventionSearch(unittest.TestCase):
    search_intervention =None
    
    @classmethod
    def setUpClass(cls):
        load_test_data()
        print ('')
        print(f'TestInterventionSearch setUpClass')
        
    @classmethod
    def tearDownClass(cls):
        close_test_data_file()
        print(f'TestInterventionSearch tearDownClass')
        
    def setUp(self):
        self.search_intervention = searchIntervention()
        print(f'TestInterventionSearch setUp')
        
    def tearDown(self):
        if self.search_intervention is not None:
            self.search_intervention = None
        print(f'TestInterventionSearch tearDown')
        
    
    def test_load_drugs_list(self):
        print("****test_load_drugs_list******")
        result = self.search_intervention.load_drugs_list('data', 'rxnav_terms.csv')
        print(f'load existing drug terms from rxnav of size = {result.shape} expect {(25482, 4)}')
        self.assertEqual(result.shape, (25482,4))
        
    def test_search_pubmed(self):
        print("****test_search_pubmed******")
        test_term = 'warfarin'
        result = self.search_intervention.search_pubmed(test_term)
        print(f'warfarin as term searched on pubmed = {isinstance(result, pd.DataFrame)} expect {True}')
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIsNotNone(result)
        
    def test_check_variant(self):
        print("****test_check_variant******")
        result = self.search_intervention.check_variant('warfarin', pd.read_csv(Path('.','data','rxnav_terms.csv')), test_data)
        print(f'warfarin variant from rxnav file = {result["main"].item()} expect warfarin' )
        self.assertEqual(result['main'].item(), 'warfarin')
    
    def test_search_codrug(self):
        print("****test_search_codrug******")
        result = self.search_intervention.search_codrug('warfarin', pd.read_csv(Path('.','data','rxnav_terms.csv')), test_data)
        print(f'warfarin codrug search from rxnav file = {result["co_drug"].item()} expect  Aspirin')
        self.assertIn('Aspirin', result['co_drug'].item())
        
        
def build_test_suite():
    # create unittest.TestSuite object.
    test_suite = unittest.TestSuite()
    # add each test function to the test suite object.
    test_suite.addTest(TestInterventionSearch('test_load_drugs_list'))
    test_suite.addTest(TestInterventionSearch('test_search_pubmed'))
    test_suite.addTest(TestInterventionSearch('test_check_variant'))
    test_suite.addTest(TestInterventionSearch('test_search_codrug'))
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