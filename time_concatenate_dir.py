#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os

lps = ['/scratch-shared/janssens/bomex-james/bomex100_e12/ppagg',
       '/scratch-shared/janssens/bomex-james/bomex100_ldelta_from100_12hr/ppagg']
ranges = [[0,48], [0,48]]
savePath = '/scratch-shared/janssens/bomex-james/bomex100_ldelta_from100_12hr/ppagg_merged'

def find_files(path):
    '''
    Find files that contain data fields.

    '''
    # Find all files that contain scenes
    files  = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if not isinstance(file, str):
                file = file.decode('utf-8')
            if '.npy' in file:
                fname = os.path.join(r, file)
                files.append(fname)
    
    files = np.sort(files)
    return files

def process(lps,savePath,ranges):
    for i in range(len(lps)):
        files = find_files(lps[i])
        for j in range(len(files)):
            arr = np.load(files[j])
            fname = files[j].split('/')[-1]
            if i == 0:
                if len(arr.shape) < 2:
                    if 'time' in fname:
                        np.save(savePath+'/'+fname, arr[ranges[i][0]:ranges[i][1]])
                    else:
                        np.save(savePath+'/'+fname, arr)
                else:
                    np.save(savePath+'/'+fname, arr[ranges[i][0]:ranges[i][1],:])
            else:
                try:
                    arr_stem = np.load(savePath+'/'+fname)
                except:
                    print('Could not load', savePath+'/'+fname, 'continuing...')
                    continue
                if len(arr.shape) < 2 and 'time' in fname:
                    arr_concat = np.concatenate([arr_stem, arr[ranges[i][0]:ranges[i][1]]])
                elif len(arr.shape) >= 2:
                    arr_concat = np.concatenate([arr_stem, arr[ranges[i][0]:ranges[i][1],:]], axis=0)
                else:
                    continue
                np.save(savePath+'/'+fname, arr_concat)
    return

process(lps,savePath,ranges)
