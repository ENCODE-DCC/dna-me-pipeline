#!/usr/bin/env python
# dmeLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
from launch import Launch
#from template import Launch # (does not use dxencode at all)

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
                                    "bam_techrep_pe":    "bam_techrep", 
                                    "bam_techrep_pe_qc": "bam_techrep_qc",
                                    "map_techrep_pe":    "map_techrep",
                                },
                            }, 
                }
        },
        "BIO_REP":  {
                "ORDER": { "se": [  "dme-extract-se", "dme-cx-to-bed-alt", "dme-bg-to-signal-alt" ],
                           "pe": [  "dme-extract-pe", "dme-cx-to-bed",     "dme-bg-to-signal" ] },
                "STEPS": {
                            "dme-extract-pe": {
                                "inputs": { "bam_pe_ABC": "bam_set", "map_report_pe_ABC": "map_report_set", "dme_ix": "dme_ix" }, 
                                "app": "dme-extract-pe", 
                                "params": { }, 
                                "results": {
                                    ###"bam_biorep":   "bam_biorep", ### Not holding on to biorep bam 
                                    "bam_biorep_qc":"bam_biorep_qc",
                                    "map_biorep":   "map_biorep",
                                    "mbias_report": "mbias_report", 
                                    "cx_report":    "cx_report", 
                                    "bg_gz":        "bg_gz", 
                                    "CpG_bed":      "CpG_bed",     "CpG_bb":       "CpG_bb", 
                                    "CHG_bed":      "CHG_bed",     "CHG_bb":       "CHG_bb", 
                                    "CHH_bed":      "CHH_bed",     "CHH_bb":       "CHH_bb",
                                },
                            },
                            "dme-extract-se": {
                                "inputs": { "bam_ABC": "bam_set", "map_report_ABC": "map_report_set", "dme_ix": "dme_ix" }, 
                                "app": "dme-extract-se", 
                                "results": {
                                    ###"bam_biorep":   "bam_biorep", ### Not holding on to biorep bam 
                                    "bam_biorep_qc":"bam_biorep_qc",
                                    "map_biorep":   "map_biorep",
                                    "mbias_report": "mbias_report", 
                                    "cx_report":    "cx_report", 
                                    "bg_gz":        "bg_gz", 
                                    "CpG_bed":      "CpG_bed",     "CpG_bb":       "CpG_bb", 
                                    "CHG_bed":      "CHG_bed",     "CHG_bb":       "CHG_bb", 
                                    "CHH_bed":      "CHH_bed",     "CHH_bb":       "CHH_bb",
                                },
                            }, 
                            "dme-cx-to-bed": {
                                "inputs": { "cx_report": "cx_report", "chrom_sizes": "chrom_sizes" }, 
                                "app": "dme-cx-to-bed", 
                                "results": {
                                    "CpG_bed":      "CpG_bed",     "CpG_bb":       "CpG_bb", 
                                    "CHG_bed":      "CHG_bed",     "CHG_bb":       "CHG_bb", 
                                    "CHH_bed":      "CHH_bed",     "CHH_bb":       "CHH_bb",
                                },
                            }, 
                            "dme-bg-to-signal": {
                                "inputs": { "bg_gz": "bg_gz", "chrom_sizes": "chrom_sizes" }, 
                                "app": "dme-bg-to-signal", 
                                "results": {
                                    "signal":       "signal",
                                },
                            }, 
                            "dme-cx-to-bed-alt": {
                                "inputs": { "cx_report": "cx_report", "chrom_sizes": "chrom_sizes" }, 
                                "app": "dme-cx-to-bed-alt", 
                                "results": {
                                    "CpG_bed":      "CpG_bed",     "CpG_bb":       "CpG_bb", 
                                    "CHG_bed":      "CHG_bed",     "CHG_bb":       "CHG_bb", 
                                    "CHH_bed":      "CHH_bed",     "CHH_bb":       "CHH_bb",
                                },
                            }, 
                            "dme-bg-to-signal-alt": {
                                "inputs": { "bg_gz": "bg_gz", "chrom_sizes": "chrom_sizes" }, 
                                "app": "dme-bg-to-signal-alt", 
                                "results": {
                                    "signal":       "signal",
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
        "bam_techrep_qc":           "/*_techrep_bismark_qc.txt",
        "bam_techrep_pe":           "/*_techrep_bismark_pe.bam", 
        "map_techrep_pe":           "/*_techrep_bismark_pe_map_report.txt",
        "bam_techrep_pe_qc":        "/*_techrep_bismark_pe_qc.txt",
        # dme-extract-pe/se inp/results:
        "bam_ABC":                  "/*_bismark.bam", 
        "map_report_ABC":           "/*_techrep_bismark_map_report.txt",
        "bam_pe_ABC":               "/*_techrep_bismark_pe.bam", 
        "map_report_pe_ABC":        "/*_techrep_bismark_pe_map_report.txt",
        ###"bam_biorep":               "/*_bismark_biorep.bam", ### Not holding on to biorep bam
        "bam_biorep_qc":            "/*_bismark_biorep_qc.txt", 
        "map_biorep":               "/*_bismark_biorep_map_report.txt",
        "signal":                   "/*_bismark_biorep.bw", 
        "CpG_bed":                  "/*_bismark_biorep_CpG.bed.gz", 
        "CpG_bb":                   "/*_bismark_biorep_CpG.bb", 
        "CHG_bed":                  "/*_bismark_biorep_CHG.bed.gz", 
        "CHG_bb":                   "/*_bismark_biorep_CHG.bb", 
        "CHH_bed":                  "/*_bismark_biorep_CHH.bed.gz", 
        "CHH_bb":                   "/*_bismark_biorep_CHH.bb", 
        "mbias_report":             "/*_bismark_biorep_mbias_report.txt",
        "cx_report":                "/*_bismark_biorep.CX_report.txt.gz",
        "bg_gz":                    "/*_bismark_biorep.bedGraph.gz",
    }

    GENOMES_SUPPORTED = ['GRCh38', 'hg19', 'mm10']
    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should use ACCESSION based fileNames
        "dme_ix":   {
                        "GRCh38": "GRCh38_XY_bismark_bowtie1_index.tgz",
                        "hg19":   "hg19_male_bismark_bowtie1_index.tgz",
                        "mm10":   "mm10_male_bismark_bowtie1_index.tgz"
                        },
        "chrom_sizes":  {
                        "GRCh38": "GRCh38_full_male.chrom.sizes",
                        "hg19":   "male.hg19.chrom.sizes",
                        "mm10":   "male.mm10.chrom.sizes"
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
        #bwaIx = self.psv['refLoc']+self.REFERENCE_FILES['bwa_index'][self.psv['genome']][self.psv['gender']]
        base_dir = '/' + self.psv['genome'] + "/dna-me/"
        dme_path = base_dir+self.REFERENCE_FILES["dme_ix"][self.psv['genome']]
        dme_fid = self.find_file(dme_path,self.REF_PROJECT_DEFAULT)
        if dme_fid == None:
            sys.exit("ERROR: Unable to locate Bismark index file '" + dme_path + "'")
        else:
            priors['dme_ix'] = dme_fid

        chrom_sizes = '/' + self.psv['genome'] + "/"+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']]
        chrom_sizes_fid = self.find_file(chrom_sizes,self.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        return priors
    

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
            if len(river['tributaries']) > 1:
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

