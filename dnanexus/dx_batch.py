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


def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Set up DNA Methylation runs on DNA Nexus')

    ap.add_argument('-t', '--test',
                    help='Use test input folder',
                    action='store_true',
                    required=False)

    ap.add_argument('-l', '--maplambda',
                    help='Run against lambda genome (QC)',
                    action='store_true',
                    default=False,
                    required=False)

    ap.add_argument('-n', '--numberjobs',
                    help='Maximum Number of jobs to run',
                    type=int,
                    default=9999,
                    required=False)


    return ap.parse_args()


def main():
    cmnd = get_args()

    ## resolve projects
    (AUTHID, AUTHPW, SERVER) = dxencode.processkey('www')

    query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&files.file_format=fastq&frame=embedded&replicates.library.biosample.donor.organism.name=mouse' % ASSAY_TERM_ID
    res = requests.get(SERVER+query, headers=HEADERS, auth=(AUTHID, AUTHPW),allow_redirects=True, stream=True)

    exps = res.json()['@graph']

    n=0
    pid = os.getpid()
    if cmnd.maplambda:
        lambdaqc = '--maplambda'
    else:
        lambdaqc = ''

    for exp in exps:
        acc = exp['accession']
        if n >= cmnd.numberjobs:
            print "Stopping at %s replicates" % n
            break
        for rep in exp.get('replicates', []):
            try:
                runcmd = "./launchDnaMe.py %s --gzip -e %s --br %s --tr %s > runs/launch%s-%s-%s.%s%s.out" % (lambdaqc, acc, rep['biological_replicate_number'], rep['technical_replicate_number'],acc, rep['biological_replicate_number'], rep['technical_replicate_number'],pid,lambdaqc)
                print runcmd
                if not cmnd.test:
                    os.system(runcmd)
                n+=1
            except KeyError, e:
                print "%s failed: %s" % (acc, e)

if __name__ == '__main__':
    main()
