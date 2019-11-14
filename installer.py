#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:34:57 2018

@author: hannes
"""

import sys
import getopt
import os

#instructions
form='\npython3 ./installer.py\n\n\t\
                '


def main(argv):
    try:
       opts, args = getopt.getopt(argv,"h",[])
    except getopt.GetoptError:
       print ('{}'.format(form))
       sys.exit()
    for opt, arg in opts:
       if opt == '-h' or opt == '-help' or opt == '--help':
          print ('{}'.format(form))
          sys.exit()

    os.system('chmod a+x *.py && echo \"export PATH=\$PATH:{}\" >> ~/.bashrc'.format(os.getcwd()))


    sys.exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
