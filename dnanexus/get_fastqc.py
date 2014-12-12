#!/usr/bin/env python
import argparse
import os
import sys
import subprocess

import dxpy
import requests
from dxencode import dxencode as dxencode

SERVER = 'https://www.encodeproject.org'
ASSAY_TYPE = 'whole genome bisulfite sequencing'
ASSAY_TERM_ID = 'OBI:0001863'
HEADERS = {'content-type': 'application/json'}
PROJECT_NAME = 'dna-me-pipeline'
APPLETS = {}


def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Get and parse fastqc results')


    ap.add_argument('-u', '--authid',
                    help='encodeD User id key (AUTHID)',
                    required=False)

    ap.add_argument('-p', '--authpw',
                    help='encodeD User id key (AUTHPW)',
                    required=False)


    return ap.parse_args()


def parse_summary(filename, summaries):

    for line in open(filename).readlines():
        (valid, test, fn) = line.strip('\n').split('\t')
        if valid not in ['PASS', 'FAIL', 'WARN']:
            print "Unknown status: %s" % line
        if not summaries.get(test):
            summaries[test] = { valid: 1 }
        elif not summaries[test].get(valid):
            summaries[test].update({valid: 1})
        else:
            summaries[test][valid] +=1

    return summaries

def main():
    cmnd = get_args()

    ## resolve projects
    project = dxencode.resolve_project(PROJECT_NAME)
    print 'Project: ' + project.describe()['name']
    pid =  project.get_id()

    counts = {}
    n = 0
    summaries = dxpy.find_data_objects(classname='file', folder='/runs', name='*_summary.txt', recurse=True, name_mode='glob', project=pid, return_handler=False)
    while summaries:
        try:
            flink = dxpy.dxlink(summaries.next())
            n = n+1
        except StopIteration:
            break
        fd = dxpy.describe(flink)
        fn = "fastqc/%s" % fd['name']
        if not os.path.isfile(fn):
            print 'Downloading: %s from %s' % (fn, fd['folder'])
            try:
                dxpy.download_dxfile(flink, fn)
            except Exception, e:
                print "Error %s" % e

        parse_summary(fn, counts)

    print "%s total summary files" % n
    for test in counts.keys():
        print "%s\t %s passed\t%s failed\t%s warnings" % (test, counts[test].get('PASS', 0), counts[test].get('FAIL', 0), counts[test].get('WARN', 0))
if __name__ == '__main__':
    main()
