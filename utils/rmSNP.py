#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re



def main():
    argvs = sys.argv
    for line in sys.stdin:
        line = line.replace('\n','')
        line = line.replace('\r','')
        lineCols = re.split('\t',line)
        if lineCols[5] == '-1' and lineCols[6] == '-1':
            continue
        else :
            print(line)

if __name__ == '__main__':
    main()
