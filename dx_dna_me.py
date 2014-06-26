import argparse
import os
import sys
import subprocess

import dxpy


ENCODE_DNA_ME_PROJECT_NAME = 'Methylation Prototype'
''' This DNA Nexus project holds all the created applets and folders'''

REPLICATES_FOLDER = '/replicates'

GENOME_REFERENCES = {
# Note this should be referred to by: biosample.donor.organism.name for any dataset
    'mouse':  'mm10.fa.gz',
    'human':  'hg19.fa.gz'
}

REFERENCE_FILES = {}
APPLETS = {}

# TODO - load from pipeline object or .json text mockups
ANALYSIS_STEPS = [
    'index',
    'trim',
    'map',
    'extract'
]

def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Generate DNAnexus workflow for the ENCODE methylation pipeline.')

    ap.add_argument('-r', '--replicates',
                    help='Replicate fastq files.',
                    nargs='+',
                    required=False)

    ap.add_argument('-e', '--experiment',
                    help='ENCODED experiment accession',
                    nargs='+',
                    required=True)

    ap.add_argument('-g', '--gender',
                    help='Gender of sample',
                    choices=['m', 'f'],
                    default='m',
                    required=False)

    ap.add_argument('-p', '--paired',
                    help='Force use of paired-end pipeline',
                    action='store_true',
                    required=False)

    return ap.parse_args()

def find_reference_file_by_name(reference_name, applets_project_id):
    '''Looks up a reference file by name in the project that holds common tools. From Joe Dale's code.'''
    cached = '*'
    if (reference_name, applets_project_id) not in REFERENCE_FILES:
        found = dxpy.find_one_data_object(classname="file", name=reference_name,
                                          project=applets_project_id,
                                          folder='/Reference Data',
                                          recurse=True,
                                          zero_ok=False, more_ok=False, return_handler=True)
        REFERENCE_FILES[(reference_name, applets_project_id)] = found
        cached = ''

    print cached + "Resolved %s to %s" % (reference_name, REFERENCE_FILES[(reference_name, applets_project_id)].get_id())
    return dxpy.dxlink(REFERENCE_FILES[(reference_name, applets_project_id)])

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

def get_project(project_name):
    project = dxpy.find_projects(name=project_name, name_mode='glob', return_handler=True)

    project = [p for p in project]
    if len(project) < 1:
        project = dxpy.DXProject(dxpy.api.project_new({'name': project_name, 'summary': 'ChIP-Seq Pipeline'})['id'])
    elif len(project) > 1:
        print 'Found more than 1 project matching ' + project_name + '.'
        print 'Please provide a unique project!'
        sys.exit(1)
    else:
        project = project[0]

    return project

def populate_workflow(wf, replicates, experiment, paired, gender, organism, applets_project_id):
    '''This function will populate the workflow for the methyl-seq Pipeline.'''

    genome = find_reference_file_by_name(GENOME_REFERENCES[organism])
    # TODO somethink like loop over analysis_steps in pipeline objects
    ### INDEX
    index_input = {
        'genome': genome
    }
    stage_id = wf.add_stage(find_applet_by_name('index'), applets_project_id, stage_input=index_input, folder=experiment)
    index_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'meIndex'
    })
    ### TRIM
    trim_input = {
        'reads': replicates
    }
    stage_id = wf.add_stage(find_applet_by_name('trim'), applets_project_id, stage_input=trim_input, folder=experiment)
    trim_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'trimmed_reads'
    })

    ### MAP
    map_input = {
        'genome': genome,
        'trimmed_reads': trim_output,
        'meIndex': index_output
    }
    stage_id = wf.add_stage(find_applet_by_name('map'), applets_project_id, stage_input=map_input, folder=experiment)
    map_output = dxpy.dxlink({
        'stage': stage_id,
        'outputField': 'mapped_files'
    })

    ### EXTRACT
    extract_input = {
        'genome': genome,
        'mapped_files': map_output
    }
    stage_id = wf.add_stage(find_applet_by_name('extract'), applets_project_id, stage_input=extract_input, folder=experiment)



def copy_files(fids, project, folder):
    new_fids = []
    for fid in fids:
        (pid, fid) = fid.split(':')
        f = dxpy.DXFile(dxid=fid,project=pid)
        fn = f.describe()['name']
        found_file = dxpy.find_one_data_object(classname='file', project=project.get_id(), folder=folder, zero_ok=True, name=fn)

        if found_file is None:
            new_fids += [dxpy.dxlink(f.clone(project.get_id(), folder))]
        else:
            new_fids += [dxpy.dxlink(found_file)]

    return new_fids

def project_has_folder(project, folder):
    folders = project.list_folder()['folders']

    return folder in folders

def resolve_applets_project():
    try:
        project = dxpy.find_one_project(name=ENCODE_DNA_ME_PROJECT_NAME, name_mode='exact', return_handler=False)
    except:
        print 'Could not find 1 and only 1 project named {0}.'.format(ENCODE_DNA_ME_PROJECT_NAME)
        exit(0)

    return dxpy.DXProject(project['id]'])

def main():
    args = get_args()

    project = resolve_applets_project()
    #project = get_project(args.project_name)
    #print project.keys()
    print 'Project: ' + project.describe()['name']
    #print project.keys()
    project_folder = project_has_folder(project, args.experiment)
    if not project_folder:
        project_folder = project.new_folder(args.experiment)

    #TODO get all replicate ids from encoded DB from ENCSR (args.experiment)
    #TODO error out if ENCSR not found, status not complete etc.

    replicates_folder = project_folder+REPLICATES_FOLDER
    if not project_has_folder(project, replicates_folder):
        replicates = dxpy.find_data_objects(classname='file', name='*.fq.gz',
                                            name_mode='glob', project=project.get_id(),
                                            folder=replicates_folder, return_handler=False)
        replicates = [dxpy.dxlink(r) for r in replicates]
    else:
        if (len(args.replicates) < 1):
            sys.exit('Need to have at least 1 replicate file.')
        project.new_folder(replicates_folder, True)
        replicates = copy_files(args.replicates, project, replicates_folder)

    paired = args.paired
    gender = args.gender
    organism = 'human'
    #TODO determine paired or gender from ENCSR metadata
    # Now create a new workflow ()
    wf = dxpy.new_dxworkflow(title='dx_dna_me_wgsbs_'+args.experiment,
                             name='ENCODE  Whole-Genome Shotgun Bisulfite Analysis Pipeline',
                             description='The ENCODE Bismark pipeline for WGBS shotgun methylation analysis',
                             project=project.get_id())

    populate_workflow(wf, replicates, args.experiment, paired, gender, organism, project['id'])
    #TODO - run the workflow automatically
    #TODO - export template workflows

if __name__ == '__main__':
    main()

