#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:18:34 2019

@author: ja18581
"""

import numpy as np
from pathlib import Path

class userAgent:
    def __init__(self):
        self.__random_agent = None
        self.__ua_file = 'ua_file.txt'
        
    def get_random_ua(self):
        try:
            with open(Path('.', self.__ua_file), 'rb') as f:
                lines = f.readlines()
            if len(lines) > 0:
                prng = np.random.RandomState()
                index = prng.permutation(len(lines) - 1)
                idx = np.asarray(index, dtype=np.integer)[0]
                self.__random_agent = lines[int(idx)]
        except Exception as ex:
            print('Exception in random_ua')
            print(str(ex))
        finally:
            return self.__random_agent
        
# =============================================================================
# if __name__=='__main__':
#     ua = userAgent()
#     user_agent = ua.get_random_ua()
#     print(user_agent.strip())
# =============================================================================
