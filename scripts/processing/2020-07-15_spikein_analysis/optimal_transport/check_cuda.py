#!/usr/bin/env python
'''
DESCRIPTION

    <+DESCRIPT+>

FOR HELP

    python check_cuda.py --help

AUTHOR:      Jake Yeung (j.yeung@hubrecht.eu)
LAB:         Quantitative Biology Lab (https://www.hubrecht.eu/research-groups/van-oudenaarden-group/)
CREATED ON:  2020-11-23
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime
import torch

def main():
    print("Checking if available")
    print(torch.cuda.is_available())
    dev = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    print("Dev:")
    print(dev)

if __name__ == '__main__':
    main()
