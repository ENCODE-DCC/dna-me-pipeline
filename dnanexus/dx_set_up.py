#!/usr/bin/env python
import argparse
import os
import sys
import subprocess

import dxpy
import requests

SERVER = 'https://www.encodeproject.org/'
ASSAY_TYPE = 'whole genome bisulfite sequencing'
HEADERS = {'content-type': 'application/json'}
PROJECT_NAME = 'dna-me-pipeline'
APPLETS = {}


def find_applet_by_name(applet_name, applets_project_id):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
    cached = '*'
    if (applet_name, applets_project_id) not in APPLETS:
        found = dxpy.find_one_data_object(classname="applet", name=applet_name,
                                          project=applets_project_id,
                                          zero_ok=False, more_ok=False, return_handler=True)
        APPLETS[(applet_name, applets_project_id)] = found
        cached = ''

    print cached + "Resolved %s to %s" % (applet_name, APPLETS[(applet_name, applets_project_id)].get_id())
    return APPLETS[(applet_name, applets_project_id)]

def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Set up DNA Methylation runs on DNA Nexus')


    ap.add_argument('-u', '--authid',
                    help='encodeD User id key (AUTHID)',
                    required=True)

    ap.add_argument('-p', '--authpw',
                    help='encodeD User id key (AUTHPW)',
                    required=True)

    return ap.parse_args()

def resolve_project(project_name, level=None):
    try:
        project = dxpy.find_one_project(name=project_name, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print 'Could not find 1 and only 1 project named %s; ' % format(project_name)
        sys.exit(1)

    return dxpy.DXProject(project['id'])

def main():
    args = get_args()

    ## resolve projects
    project = resolve_project(PROJECT_NAME)
    print 'Project: ' + project.describe()['name']
    pid =  project.get_id()

    applet = find_applet_by_name('fastqc-exp', pid )

    query = '/search/?type=experiment&assay_term_name=%s&award.rfa=ENCODE3&limit=all&frame=embedded' % ASSAY_TYPE
    res = requests.get(SERVER+query, headers=HEADERS, auth=(args.authid, args.authpw),allow_redirects=True, stream=True)

    exps = res.json()['@graph']

    #for exp in exps:
    acc = exps[0]['accession']
    applet.run({ "accession": acc}, project=pid)

if __name__ == '__main__':
    main()
