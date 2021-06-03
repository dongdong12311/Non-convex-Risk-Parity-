#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 23:07:09 2021

@author: dongdong
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 15:45:24 2020

@author: dongdong
"""

import struct
import numpy as np
import pandas as pd
import os

for s in range(1,6):
    
    basedir = '../data/logportfoliodata/sourcedata/'+str(s)+'/'
    ns = [10, 50, 100, 150, 200 ,250 ,300, 350, 400 ,450 ,
          500 ,600, 700, 800, 900 ,1000,2000,3000,4000,5000 ]
    for n in ns:
        
     
        filename = basedir + str(n) + '.mat'
        target_dir =  '../data/logportfoliodata/'+str(s)+'/'
        try:
            os.mkdir(target_dir)
        except:
            pass
        target_dir = os.path.join(target_dir,str(n))
        try:
            os.mkdir(target_dir)
        except:
            pass
         
        all_item = ['Aeq','Sigma','beq','lb']
        types = [np.float64,np.float64,np.float64,np.float64,np.float64,
                 np.int,np.int,np.int,np.int,
                 np.int,np.int,np.float64,np.float64]
        pack_map = {np.float64:'d',np.int:'i'}
        import scipy.io
        data = scipy.io.loadmat(filename)
        for ind in range(len(all_item)):
            item = all_item[ind]
            temp = data[item]
         
            target_file = os.path.join(target_dir,all_item[ind]+'.mat') 
            pack_type = pack_map[types[ind]]
            
            with open(target_file,'wb+') as fd:
                byte_len = 8
                for i in temp:
                    for j in i:    
                        bd =  struct.pack(pack_type,j)
        
                        fd.write(bd)
