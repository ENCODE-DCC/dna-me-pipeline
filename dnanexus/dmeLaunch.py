#!/usr/bin/env python
# dmeLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
#from dxencode import dxencode as dxencode
import dxencode as dxencode
from launch import Launch

class DmeLaunch(Launch):
    '''Descendent from Launch class with 'rampage' methods'''

    PIPELINE_NAME = "dna-me"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''
    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline " + \
                    "analysis for one replicate or combined replicates. "
    ''' This help title should name pipline and whether combined replicates are supported.'''
                    
    RESULT_FOLDER_DEFAULT = '/dname/'
    ''' This the default location to place results folders for each experiment.'''
    
    PIPELINE_BRANCH_ORDER = [ "TECH_REP", "BIO_REP" ] #, "COMBINED_REPS" ]
    '''A pipeline is frequently made of branches that flow into each other, such as replicate level to combined replicate.'''
    
    PIPELINE_BRANCHES = {
    #'''Each branch must define the 'steps' and their (artificially) linear order.'''
        "TECH_REP": {
                "ORDER": { "se": [ "dme-align-se" ],
                           "pe": [ "dme-align-pe" ] },
                "STEPS": {
                            "dme-align-se": {
                                "inputs": { "reads": "reads", "dme_ix": "dme_ix" },
                                "app": "dme-align-se", 
                                "params": { }, 
                                "results": {
                                    "bam_techrep":      "bam_techrep", 
                                    #"bam_techrep_qc":   "bam_techrep_qc", # Don't include this: old alignments didn't make it
                                    "map_techrep":      "map_techrep",
                                },
                            },
                            "dme-align-pe": {
                                "inputs": { "reads1": "reads1", "reads2": "reads2", "dme_ix": "dme_ix" }, 
                                "app": "dme-align-pe", 
                                "params": { "min_insert": "min_insert", "min_insert": "min_insert" }, 
                                "results": {
                                    "bam_techrep":      "bam_techrep", 
                                    "bam_techrep_qc":   "bam_techrep_qc",
                                    "map_techrep":      "map_techrep",
                                },
                            }, 
                }
        },
        "BIO_REP":  {
                "ORDER": { "se": [  "dme-extract-se" ],
                           "pe": [  "dme-extract-pe" ] },
                "STEPS": {
                            "dme-extract-pe": {
                                "inputs": { "bam_ABC": "bam_set", "map_report_ABC": "map_report_set", "dme_ix": "dme_ix" }, 
                                "app": "dme-extract-pe", 
                                "params": { "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_biorep":   "bam_biorep", 
                                    "bam_biorep_qc":"bam_biorep_qc",
                                    "map_biorep":   "map_biorep",
                                    "CpG_bed":      "CpG_bed", 
                                    "CHG_bed":      "CHG_bed", 
                                    "CHH_bed":      "CHH_bed", 
                                    "CpG_bb":       "CpG_bb", 
                                    "CHG_bb":       "CHG_bb", 
                                    "CHH_bb":       "CHH_bb",
                                    "mbias_report": "mbias_report", 
                                },
                            },
                            "dme-extract-se": {
                                "inputs": { "bam_ABC": "bam_set", "map_report_ABC": "map_report_set", "dme_ix": "dme_ix" }, 
                                "app": "dme-extract-se", 
                                "params": { "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_biorep":   "bam_biorep", 
                                    "bam_biorep_qc":"bam_biorep_qc",
                                    "map_biorep":   "map_biorep",
                                    "CpG_bed":      "CpG_bed", 
                                    "CHG_bed":      "CHG_bed", 
                                    "CHH_bed":      "CHH_bed", 
                                    "CpG_bb":       "CpG_bb", 
                                    "CHG_bb":       "CHG_bb", 
                                    "CHH_bb":       "CHH_bb",
                                    "mbias_report": "mbias_report", 
                                },
                            }, 
                }
        },
    }

    # NOTE: dme-align produces *_techrep_bismark.bam and dme-extract merges 1+ techrep bams into a *_bismark_biorep.bam.
    #       The reason for the name order is so thal older *_bismark.bam alignments are recognizable as techrep bams
    FILE_GLOBS = {
        #"reads":                    "/*.fq.gz",
        #"reads1":                   "/*.fq.gz",
        #"reads2":                   "/*.fq.gz",
        # dme-align-pe/se results:
        "bam_techrep":              "/*_bismark.bam", 
        "map_techrep":              "/*_techrep_bismark_map_report.txt",
        "bam_techrep_qc":           "/*_bismark_techrep_qc.txt",
        # dme-extract-pe/se inp/results:
        "bam_ABC":                  "/*_bismark.bam", 
        "map_report_ABC":           "/*_bismark_map_report.txt",
        "bam_biorep":               "/*_bismark_biorep.bam", 
        "bam_biorep_qc":            "/*_bismark_biorep_qc.txt", 
        "map_biorep":               "/*_bismark_techrep_map_report.txt",
        "CpG_bed":                  "/*_bismark_CpG.bed.gz", 
        "CHG_bed":                  "/*_bismark_CHG.bed.gz", 
        "CHH_bed":                  "/*_bismark_CHH.bed.gz", 
        "CpG_bb":                   "/*_bismark_CpG.bb", 
        "CHG_bb":                   "/*_bismark_CHG.bb", 
        "CHH_bb":                   "/*_bismark_CHH.beb", 
        "mbias_report":             "/*_bismark_mbias_report.txt",
    }

    REF_PROJECT_DEFAULT = "dna-me-pipeline"  # TODO: move all ref files to ref project!
    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should use ACCESSION based fileNames
        "dme_ix":   {
                        "hg19": {
                                "female":   "hg19_female_bismark_bowtie1_index.tgz",
                                "male":     "hg19_male_bismark_bowtie1_index.tgz"
                                },
                        "mm10": {
                                "female":   "mm10_male_bismark_bowtie1_index.tgz",
                                "male":     "mm10_male_bismark_bowtie1_index.tgz"
                                }
                        },
        }


    def __init__(self):
        Launch.__init__(self)
        
    def get_args(self):
        '''Parse the input arguments.'''
        ap = Launch.get_args(self,parse=False)
        
        # NOTE: Could override get_args() to have this non-generic control message
        #ap.add_argument('-c', '--control',
        #                help='The control bam for peak calling.',
        #                required=False)

        return ap.parse_args()

    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        psv = Launch.pipeline_specific_vars(self,args)
        
        # Some specific settings
        psv['nthreads']    = 8
        psv['min_insert']  = 0
        psv['max_insert']  = 500
        
        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        # TODO:  move all ref files to ref project and replace "/ref/" and self.REF_PROJECT_DEFAULT
        #bwaIx = self.psv['refLoc']+self.REFERENCE_FILES['bwa_index'][self.psv['genome']][self.psv['gender']]
        base_dir = "/ref/" + self.psv['genome'] + "/dna-me/"
        dmeIx = base_dir+self.REFERENCE_FILES["dme_ix"][self.psv['genome']][self.psv['gender']]
        #dmeIxFid = dxencode.find_file(dmeIx,dxencode.REF_PROJECT_DEFAULT)
        dmeIxFid = dxencode.find_file(dmeIx,self.REF_PROJECT_DEFAULT)
        if dmeIxFid == None:
            sys.exit("ERROR: Unable to locate Bismark index file '" + dmeIx + "'")
        else:
            priors['dme_ix'] = dmeIxFid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
    

    def add_combining_reps(self, psv):
        '''Defines how replicated are combined.'''
        # OVERRIDING parent because DNase-seq pipeline doesn't follow the standard replicate combination model
        debug=False
        
        reps = psv['reps']
        # In the 'standard combining model' PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
        # and all replicates are in psv['reps'] keyed as 'a','b',... and having rep['rep_tech'] = 'rep1_1'
        # All these simple reps will have rep['branch_id'] = "REP"
        
        # First, each tech_rep is processed individually
        bio_reps = []
        for rep_id in sorted( reps.keys() ):
            if len(rep_id) == 1: # single letter: simple replicate
                rep = reps[rep_id]
                rep['branch_id'] = "TECH_REP"
                if rep['br'] not in bio_reps:
                    bio_reps.append(rep['br'])
                else:
                    self.combined_reps = True  # More than one tech_rep per bio_rep so combining will be done!
                if debug:
                    print "DEBUG: rep: " + rep_id
                    
        # Next bio_reps have their technical replicates merged and processing continues
        for bio_rep in bio_reps:
            river = {}
            river['branch_id'] = "BIO_REP"
            river['tributaries'] = []
            river['rep_tech'] = 'reps' + str(bio_rep) + '_'  # reps1_1.2.3 is rep1_1 + rep1_2 + rep1_3
            river['br'] = bio_rep
            for tributary_id in sorted( reps.keys() ): 
                if len(tributary_id) == 1:
                    tributary = reps[tributary_id]
                    if tributary['br'] == bio_rep:
                        if len(river['tributaries']) > 0:
                            river['rep_tech'] += '.'
                        river['rep_tech'] += tributary['rep_tech'][5:]
                        river['tributaries'].append(tributary_id)
            assert len(river['tributaries']) >= 1  # It could be the case that there is one tech_rep for a bio_rep!
            # river_id for ['a','b'] = 'b-bio_rep1'
            river_id = river['tributaries'][-1] + '-bio_rep' + str(bio_rep)
            reps[river_id] = river
            # Special case of 2 allows for designating sisters
            if len(river['tributaries']) == 2:
                reps[river['tributaries'][0]]['sister'] = river['tributaries'][1]
                reps[river['tributaries'][1]]['sister'] = river['tributaries'][0]
            if debug:
                print "DEBUG: biorep: " + river_id + " tributaries: " + str(len(river['tributaries']))

        # Finally a pair of bio_reps are merged and processing finishes up
        if 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER:
            if len(bio_reps) == 2:
                self.combined_reps = True  # More than one bio_rep so combining will be done!
                sea = {} # SEA is the final branch into which all tributaries flow
                sea['branch_id'] = 'COMBINED_REPS'
                sea['tributaries'] = []
                sea['rep_tech'] = 'reps'
                for tributary_id in sorted( reps.keys() ):
                    if len(tributary_id) == 1:  # ignore the simple reps
                        continue 
                    tributary = reps[tributary_id]
                    if len(sea['tributaries']) > 0:
                        sea['rep_tech'] += '-'
                    sea['rep_tech'] += tributary['rep_tech'][4:]
                    sea['tributaries'].append(tributary_id)
            
                psv['rep_tech'] = sea['rep_tech']
                reps[self.SEA_ID] = sea
                # Special case of 2 allows for designating sisters
                reps[sea['tributaries'][0]]['sister'] = sea['tributaries'][1]
                reps[sea['tributaries'][1]]['sister'] = sea['tributaries'][0]
            #else:
            #    print "Found " + str(len(bio_reps)) + " bio_reps.  If exactly two, they would be combined."
            #print json.dumps(reps,indent=4,sort_keys=True)

    #######################


if __name__ == '__main__':
    '''Run from the command line.'''
    dmeLaunch = DmeLaunch()
    dmeLaunch.run()

