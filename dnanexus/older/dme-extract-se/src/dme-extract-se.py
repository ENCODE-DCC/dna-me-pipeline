#!/usr/bin/env python
# dme-extract-se.py  extract methylation using bismark
#
# Methyl extraction using bismark is slow and needs a great deal of storage.
# A dx execution instance mem1_hdd2_x32 has 32 cpus and 3360GB of storage.
# Only part of the bismark extraction can capitalize on multiple cpus, but
# unfortunately, the single threaded bismark2bedGraph and coverage2cytosine 
# have been running out of storage when run in an instance with "only"
# 1680GB.  So this will run the full bismark_methylation_extractor (with
# coverage reporting), on the large, expensive mem1_hdd2_x32, then continue
# with signal and bed/bb creation on a less expensive instance.

# The encodeD files graph relating steps and dependencies suggests that a 
# single step going from tech_rep bams all the way to bed/bb results is most
# appropriate.  Therefore:
# 1) This applet is run on a cheaper cpu with enough storage, but runs a sub-job
#    for techrep bam merge/sort and full bismark_methylation_extraction 
#    (--cytosine_report --CX_contexton) on a mem1_hdd2_x32.  Only the samtools sort,
#    extraction and pigz compression can benefit from 32 cpus.  However, the large
#    storage of this instance is proving necessary for the single threaded cytosine
#    report gneration.
# 2) Context reports to bed/bb and bedGraph to bw are run on the main instance with
#    hopefuuly smaller storage footprint.
# 3) Results of this applet are 3 bed and 3 bb files, plus qc, mbias report and a bw.
# NOTE: at this time, this version only differs from the PE version by a single 
#       parameter in the bismark_methylation_extractor comand. Nevertheless,
#       encodeD requires distinct pe and se analysis_steps so these are separate applets.

import sys, os, subprocess, json, requests, shlex, urlparse
import os.path
#from datetime import datetime
import dxpy

def run_cmd(cmd, out="stdout", append=False, silent=False):
    fh = None
    mode = 'w'
    if out != "stdout":
        if append:
            mode = 'a'
        fh = open(out, mode)
    if not silent:
        line = '> ' + cmd
        if fh != None:
            if mode == 'a':
                line += ' >> ' + out
            else:
                line += ' > ' + out
        print line
    if fh != None:
        subprocess.check_call(shlex.split(cmd),stdout=fh)
        fh.close()
    else:
        subprocess.check_call(shlex.split(cmd))
    
def append_line(line,filename):
    #mode = 'a'
    #if not os.path.isfile(filename):
    #    mode = 'w' 
    try:
        fh = open(filename, 'a')
        fh.write(line + '\n')
        fh.close()
    except:
        print >> sys.stderr, "* FAILED to append to " + filename
        
def bam_or_sam(biorep_bam, uncompress_bam, target_root, bam_size=0):
    '''Decides whether to extract from bam or uncompressed sam and returns the file name.'''
    # NOTE: Better to use sam and let extractor use more threads, but this takes up precious storage
    alignments=biorep_bam
    ncores=32
    # TODO: What bam_size constitutes too large for sam?  93,953,130,496 is fine!
    if bam_size == 0:
        bam_size = os.path.getsize(biorep_bam)
    if uncompress_bam and bam_size < 400000000000:
        alignments = target_root + ".sam"
        print "* extract(): Decompressing biorep bam (size: %d)..." % bam_size
        run_cmd('samtools view ' + biorep_bam, out=alignments)
        run_cmd('rm -f ' + biorep_bam) # STORAGE IS LIMITED
    else:
        ncores = ncores / 2
        if uncompress_bam:
            print "* Using compressed bam and %d cores because bam_size: %d exceeds limit." % (ncores,bam_size)
        else:
            print "* Using compressed biorep bam (size: %d) and %d cores..." % (bam_size, ncores)
    return (alignments, ncores)

###### merging bams
def merge_bams(bam_set, ncores):
    '''Merges techrep bams into biorep bam.'''
    # NOTE: dme-align produces *_techrep_bismark.bam and dme-extract merges 1+ techrep bams into a *_bismark_biorep.bam.
    #       The reason for the name 'word' order is so thal older *_bismark.bam alignments are recognizable as techrep bams

    target_root = ""
    merged = ""
    tech_reps = ""
    exp_id = ""
    rep_tech = ""
    
    for techrep_bam_dlink in reversed(bam_set):
        file_desc = dxpy.describe(techrep_bam_dlink)
        file_root = file_desc['name']
        print "* Working on '" + str(techrep_bam_dlink) + "' " + file_root
        file_root = file_root.replace('_techrep_bismark.bam','') 
        file_root = file_root.replace('_bismark.bam','') 
        if len(target_root) == 0:
            target_root=file_root
        else:
            target_root = file_root + '_' + target_root
            if len(merged) == 0:
                target_root += '_bismark_biorep' 
                merged = 's merged as'

        # Try to simplify the names
        if os.path.isfile('/usr/bin/parse_property.py'): 
            if len(exp_id) == 0:
                file_path = file_desc['folder'] + '/' + file_desc['name']
                exp_id = subprocess.check_output(shlex.split('parse_property.py -f %s --project %s --exp_id -q' \
                                                                        % (file_desc['id'], file_desc['project']) ))
                exp_id = ''.join(exp_id.split()) # Remove \n, etc.
                if len(exp_id) > 0:
                    print "* Discovered exp_id: '%s'" % exp_id
            if len(exp_id) > 0:
                rep_tech = subprocess.check_output(shlex.split('parse_property.py -f %s --project %s --rep_tech -q' \
                                                                        % (file_desc['id'], file_desc['project']) ))
                rep_tech = ''.join(rep_tech.split()) # Remove \n, etc.
        if len(rep_tech) > 0:
            print "* Discovered rep_tech: '%s'" % rep_tech
            if len(tech_reps) > 0:
                tech_reps = tech_reps + '_' + rep_tech
            else:
                tech_reps = rep_tech
                
        print "* Downloading %s_techrep_bismark.bam file..." % file_root
        dxpy.download_dxfile(techrep_bam_dlink, file_root + '_techrep_bismark.bam')
    
        if not os.path.isfile("sofar.bam"): 
            run_cmd('mv %s_techrep_bismark.bam sofar.bam' % file_root)
        else:
            print "* Merging in %s_techrep_bismark.bam..." % file_root
            # NOTE: keeps the first header
            run_cmd('samtools cat sofar.bam %s_techrep_bismark.bam' % file_root, out='merging.bam')
            run_cmd('mv merging.bam sofar.bam')
            run_cmd('rm %s_techrep_bismark.bam' % file_root) # STORAGE IS LIMITED
                

    if len(exp_id) > 0 and len(tech_reps) > 0:
        target_root = '%s_%s_bismark_biorep' % (exp_id, tech_reps)
        print "* Name biorep bam as: %s.bam" % target_root
    else:
        print "* Long biorep bam to be named: %s.bam" % target_root

    # At this point there is a 'sofar.bam' with one or more input bams
    if len(merged) == 0:
        target_root = file_root + "_bismark_biorep"
        run_cmd('mv sofar.bam %s.bam' % target_root)
        print "* Only one input file '%s.bam', no merging required." % target_root
    else:
        # sorting needed due to samtools cat
        print "* Sorting merged bam..."
        run_cmd('samtools sort -@ %d -m 1600M -f sofar.bam %s.bam' % (ncores,target_root) )
        run_cmd('rm sofar.bam') # STORAGE IS LIMITED
        print "* Files merged into '%s.bam'" % target_root
    
    return (target_root,target_root+'.bam')

def merge_map_reports(map_report_set, target_root):
    '''Merges techrep map_reports.'''

    # Working on map_reports now
    all_reports=""
    biorep_map_report = target_root + '_map_report.txt'
    append_line("### Combined Bismark map report for several technical replicates ###\n",biorep_map_report)
    for techrep_map_report_dlink in map_report_set:
        file_desc = dxpy.describe(techrep_map_report_dlink)
        file_root = file_desc['name']
        file_root = file_root.replace('_techrep_bismark_map_report.txt','') 
        file_root = file_root.replace('_bismark_map_report.txt','') 
        file_root = file_root.replace('_map_report.txt','')
        techrep_map_report = file_root + '_techrep_map_report.txt' 
        append_line("###################################",biorep_map_report)
        append_line("### Map report for ${file_root} ###",biorep_map_report)
        print "* Downloading %s_techrep_bismark_map_report.txt file..." % file_root
        dxpy.download_dxfile(techrep_map_report_dlink, techrep_map_report)
        run_cmd('cat ' + techrep_map_report, out=biorep_map_report,append=True)
        if len(all_reports) == 0:
            all_reports = techrep_map_report
        else:
            all_reports += ',' + techrep_map_report
        
    if all_reports == techrep_map_report: # only one
        run_cmd('mv %s %s' % (techrep_map_report,biorep_map_report) )
        all_reports = biorep_map_report
        
    return (biorep_map_report,all_reports)
    
def biorep_bam_qc_metrics(target_root, all_reports):
    '''Initial QC stats based upon map_report and flagstats.'''

    print "* Collect bam stats..."
    run_cmd('samtools flagstat %s.bam' % target_root, out='%s_flagstat.txt' % target_root )
    # NOTE: samtools stats may take longer than it is worth 
    #run_cmd('samtools stats %s.bam' % target_root, out='%s_samstats.txt' % target_root )
    #run_cmd('head -3 %s_samstats.txt' % target_root)
    #run_cmd('grep ^SN %s_samstats.txt | cut -f 2-' % target_root,out='%s_samstats_summary.txt' % target_root )

    print "* Prepare metadata..."
    qc_stats=''
    reads=0
    #read_len=0
    if os.path.isfile('/usr/bin/qc_metrics.py'): 
        qc_stats = subprocess.check_output(shlex.split('qc_metrics.py -n bismark_map -f ' + all_reports))
        meta = subprocess.check_output(shlex.split('qc_metrics.py -n samtools_flagstats -f %s_flagstat.txt' % target_root))
        qc_stats += ', ' + meta

        reads = int(subprocess.check_output(shlex.split('qc_metrics.py -n samtools_flagstats -f %s_flagstat.txt -k total' % target_root)))
        #meta= = subprocess.check_output(shlex.split('qc_metrics.py -n samtools_stats -d ':' -f %s_samstats_summary.txt' % target_root))
        #read_len=int(subprocess.check_output(shlex.split('qc_metrics.py -n samtools_stats -d ':' -f %s_samstats_summary.txt -k "average length"' % target_root)))
        #qc_stats += ', ' + meta
    qc_metrics = json.loads('{'+qc_stats+'}')

    # All qc to one file per target file:
    append_line("===== samtools flagstat =====", target_root+'_qc.txt')
    run_cmd('cat %s_flagstat.txt' % target_root, out=target_root+'_qc.txt',append=True,silent=True)
    #append_line("\n===== samtools stats =====", target_root+'_qc.txt')
    #run_cmd('cat %s_samstats.txt' % target_root,out=target_root+'_qc.txt',append=True,silent=True)
    
    return (qc_metrics, reads, target_root+'_qc.txt')

###### Bismark extraction
def bismark_full_extract(target_root, alignments, ncores):
    '''bismark_methylation_extractor'''
     
    print "* bismark_full_extract(): Analyse methylation in %s and using %d threads..." % (alignments, ncores)
    # TODO: missing: --no_overlap/--include_overlap --ignore_XXX
    # NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
    cmd = 'bismark_methylation_extractor --multicore %d --single-end -s --comprehensive --cytosine_report ' % (ncores)
    cmd += '--CX_context --ample_mem --output output/ --zero_based --genome_folder input ' + alignments
    run_cmd('mkdir -p output/')
    run_cmd(cmd)
    run_cmd('ls -l output/')
    if os.path.isfile('output/%s_splitting_report.txt' % target_root): 
        run_cmd('mv output/%s_splitting_report.txt %s_splitting_report.txt' % (target_root,target_root)) 
    run_cmd('mv output/%s.M-bias.txt %s_mbias_report.txt' % (target_root, target_root))

def bismark_qc_metrics(target_root, qc_metrics):
    # bismark_methylation_extractor
     
    if os.path.isfile('/usr/bin/qc_metrics.py') and os.path.isfile('%s_splitting_report.txt' % target_root): 
        print "* bismark_full_extract(): Extract metadata..."
        meta = subprocess.check_output(shlex.split('qc_metrics.py -n bismark_extract -f %s_splitting_report.txt' % target_root))
        qc_metrics.update(json.loads('{'+meta+'}'))
    return qc_metrics

@dxpy.entry_point("merge_extract_full")
def merge_extract_full(bam_set, map_report_set, dme_ix_dxlink, uncompress_bam, props):
    '''subjob runs bismark_methylation_extractor on mem1_hdd2_x32'''

    (target_root,biorep_bam) = merge_bams(bam_set, 32)
    (biorep_map,all_reports) = merge_map_reports(map_report_set, target_root)
    (qc_metrics, reads, biorep_bam_qc) = biorep_bam_qc_metrics(target_root, all_reports)
    
    print "* merge_extract_full(): Retrieve and uncompress index..."
    dme_ix = "dme_index.tar.gz"
    dxpy.download_dxfile(dme_ix_dxlink, dme_ix)
    run_cmd('tar -zxf ' + dme_ix)

    # NOTE: Better to use sam and let extractor use more threads, but this takes up precious storage
    (alignments, ncores) = bam_or_sam(biorep_bam, uncompress_bam, target_root)

    bismark_full_extract(target_root, alignments, ncores)
    qc_metrics = bismark_qc_metrics(target_root, qc_metrics)

    print "* merge_extract_full(): Retrieve split report..."
    append_line("\n===== bismark_methylation_extractor: splitting_report =====",biorep_bam_qc)
    run_cmd('cat %s_splitting_report.txt' % target_root,out=biorep_bam_qc,append=True,silent=True)

    # TODO: Is this even needed?  Currently we do to get the size!
    #if len(bam_set) > 1:  # Wouldn't need to do this unless there is a merge
    #    print "* merge_extract(): Storing biorep bam..."
    #    props_ex = props.copy()
    #    props_ex.update({ 'reads': str(reads) })
    #    biorep_bam_dxlink = dxpy.dxlink(dxpy.upload_local_file(biorep_bam,properties=props_ex,details=qc_metrics,wait_on_close=True))
    #else:
    #    biorep_bam_dxlink = bam_set[0]

    print "* merge_extract_full(): Storing extraction results..."
    biorep_bam_qc_dxfile = dxpy.upload_local_file(biorep_bam_qc,properties=props,details=qc_metrics)
    biorep_map_dxfile    = dxpy.upload_local_file(biorep_map,   properties=props,details=qc_metrics)
    run_cmd('pigz output/%s.CX_report.txt' % target_root)
    cx_report_dxfile     = dxpy.upload_local_file('output/%s.CX_report.txt.gz' % target_root)
    bedgraph_gz_dxfile   = dxpy.upload_local_file('output/%s.bedGraph.gz' % target_root)
    chrom_sizes_dxfile   = dxpy.upload_local_file('input/chrom.sizes')
    split_report_dxfile  = dxpy.upload_local_file(target_root+'_splitting_report.txt')
    mbias_report_dxfile  = dxpy.upload_local_file(target_root+'_mbias_report.txt',properties=props,details=qc_metrics)

    print "* merge_extract_full(): Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    return {
        #"biorep_bam_dxlink":    biorep_bam_dxfile,
        "biorep_bam_qc_dxlink": dxpy.dxlink(biorep_bam_qc_dxfile),
        "biorep_map_dxlink":    dxpy.dxlink(biorep_map_dxfile),
        "split_report_dxlink":  dxpy.dxlink(split_report_dxfile),
        "cx_report_dxlink":     dxpy.dxlink(cx_report_dxfile),
        "bedgraph_gz_dxlink":   dxpy.dxlink(bedgraph_gz_dxfile),
        "chrom_sizes_dxlink":   dxpy.dxlink(chrom_sizes_dxfile),
        "mbias_report_dxlink":  dxpy.dxlink(mbias_report_dxfile),
        "target_root":          target_root,
        "qc_metrics":           qc_metrics
    }

###### Post extraction routines
def bedmethyl(target_root, cx_report, chrom_sizes, cleanup=False):
    '''runs cxrepo-bed.py and bedToBigBed.'''
     
    print "* bedmethyl(): Run cxrepo-bed.py..."
    run_cmd('cxrepo-bed.py -o output/ output/' + cx_report)
    CpG_root = '%s_CpG' % target_root
    CHG_root = '%s_CHG' % target_root
    CHH_root = '%s_CHH' % target_root
    run_cmd('mv output/CG_%s  %s.bed' % (cx_report,CpG_root))
    run_cmd('mv output/CHG_%s %s.bed' % (cx_report,CHG_root))
    run_cmd('mv output/CHH_%s %s.bed' % (cx_report,CHH_root))
    if cleanup:
        run_cmd("rm -rf output/")

    print "* bedmethyl(): Convert to BigBed..."
    run_cmd('bedToBigBed %s.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 %s %s.bb' % (CpG_root,chrom_sizes,CpG_root))
    run_cmd('bedToBigBed %s.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 %s %s.bb' % (CHG_root,chrom_sizes,CHG_root))
    run_cmd('bedToBigBed %s.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 %s %s.bb' % (CHH_root,chrom_sizes,CHH_root))
    run_cmd('pigz %s.bed' % CpG_root)
    run_cmd('pigz %s.bed' % CHG_root)
    run_cmd('pigz %s.bed' % CHH_root)
    return (CpG_root + '.bed.gz',CHG_root + '.bed.gz',CHH_root + '.bed.gz',CpG_root + '.bb',CHG_root + '.bb',CHH_root + '.bb')

#@dxpy.entry_point("bedmethyl_io")
def bedmethyl_io(cx_report_dxlink, chrom_sizes_dxlink, target_root, qc_metrics, props):
    '''subjob runs cxrepo-bed.py, bedToBigBed on mem3_hdd2_x8'''

    print "* bedmethyl_io(): Retrieve CX report and chrom.sizes..."
    run_cmd('mkdir -p output/')
    cx_report = target_root + ".CX_report.txt"
    chrom_sizes = "chrom.sizes"
    dxpy.download_dxfile(chrom_sizes_dxlink, chrom_sizes)
    dxpy.download_dxfile(cx_report_dxlink, 'output/%s.gz' % cx_report)
    run_cmd('gunzip output/%s.gz' % cx_report)

    (CpG_bed,CHG_bed,CHH_bed,CpG_bb,CHG_bb,CHH_bb) = bedmethyl(target_root, cx_report, chrom_sizes) #,cleanup=True) # TODO
    
    print "* bedmethyl_io(): Storing bedmethyl results..."
    CpG_bed_dxfile = dxpy.upload_local_file(CpG_bed,properties=props,details=qc_metrics)
    CHG_bed_dxfile = dxpy.upload_local_file(CHG_bed,properties=props,details=qc_metrics)
    CHH_bed_dxfile = dxpy.upload_local_file(CHH_bed,properties=props,details=qc_metrics)

    CpG_bb_dxfile = dxpy.upload_local_file(CpG_bb,properties=props,details=qc_metrics)
    CHG_bb_dxfile = dxpy.upload_local_file(CHG_bb,properties=props,details=qc_metrics)
    CHH_bb_dxfile = dxpy.upload_local_file(CHH_bb,properties=props,details=qc_metrics)

    print "* bedmethyl_io(): Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    return {
        "CpG_bed_dxlink": dxpy.dxlink(CpG_bed_dxfile),
        "CHG_bed_dxlink": dxpy.dxlink(CHG_bed_dxfile),
        "CHH_bed_dxlink": dxpy.dxlink(CHH_bed_dxfile),
        "CpG_bb_dxlink":  dxpy.dxlink(CpG_bb_dxfile),
        "CHG_bb_dxlink":  dxpy.dxlink(CHG_bb_dxfile),
        "CHH_bb_dxlink":  dxpy.dxlink(CHH_bb_dxfile)
    }


def signal(target_root, bedGraph, chrom_sizes, cleanup=False):
    '''subjob runs bedGraphToBigWig on mem3_hdd2_x8'''

    print "* signal(): Convert to signal bedGraph to bigWig..."
    if bedGraph.endswith('.gz'):
        run_cmd('gunzip ' + bedGraph)
        bedGraph = bedGraph[0:-3]
    bigWig = target_root + ".bw"
    run_cmd('bedGraphToBigWig %s %s %s' % (bedGraph, chrom_sizes, bigWig))
    if cleanup:
        run_cmd('rm ' + bedGraph)
    
    return bigWig

def signal_io(bedgraph_gz_dxlink, chrom_sizes_dxlink, target_root, qc_metrics, props):
    '''subjob runs bedGraphToBigWig on mem3_hdd2_x8'''

    print "* signal_io(): Retrieve bedgraph and chrom.sizes..."
    bedGraph = target_root + ".bedGraph"
    bedGraph_gz = bedGraph + ".gz"
    chrom_sizes = "chrom.sizes"
    dxpy.download_dxfile(bedgraph_gz_dxlink, bedGraph_gz)
    dxpy.download_dxfile(chrom_sizes_dxlink, chrom_sizes)

    bigWig = signal(target_root, bedGraph_gz, chrom_sizes)

    print "* signal_io(): Storing signal results..."
    bigWig_dxfile = dxpy.upload_local_file(bigWig,properties=props,details=qc_metrics) #,cleanup=True) # TODO

    print "* signal_io(): Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    return {
        "bigWig_dxlink": dxpy.dxlink(bigWig_dxfile)
    }

@dxpy.entry_point("main")
def main(bam_set, map_report_set, dme_ix, uncompress_bam=True):

    # tool_versions.py --applet $script_name --appver $script_ver
    props = {}
    if os.path.isfile('/usr/bin/tool_versions.py'): 
        sw_versions = subprocess.check_output(['tool_versions.py', '--dxjson', 'dnanexus-executable.json'])
        props["SW"] = sw_versions
    
    print "* Value of bam_set:        '" + str(bam_set) + "'"
    print "* Value of map_report_set: '" + str(map_report_set) + "'"
    print "* Value of dme_ix:         '" + str(dme_ix) + "'"
    print "* Value of uncompress_bam: '" + str(uncompress_bam) + "'"

    print "* Calling merge_extract_full()..."
    inp = {
        'bam_set':        bam_set,
        'map_report_set': map_report_set, 
        'dme_ix_dxlink':  dme_ix,
        'uncompress_bam': uncompress_bam,
        'props':          props
    }
    extract_job = dxpy.new_dxjob(inp, "merge_extract_full")
    print "* Kicked off extract() and waiting..."
    extract_job.wait_on_done() # Wait because we want the qc_metrics to pass to other jobs.
    extract_out = extract_job.describe()['output']
    target_root = extract_out['target_root']
    qc_metrics = extract_out['qc_metrics']


    print "* Calling bedmethyl()..."
    # What is cheaper?  bedmethyl and signal in main or farm one out to a separate process?
    bedmethyl_out = bedmethyl_io(extract_out["cx_report_dxlink"], extract_out["chrom_sizes_dxlink"], target_root, qc_metrics, props)
    #inp = {
    #    'cx_report_dxlink':   extract_out["cx_report_dxlink"],
    #    'chrom_sizes_dxlink': extract_out["chrom_sizes_dxlink"],
    #    'target_root':        target_root,
    #    'qc_metrics':         qc_metrics,
    #    'props':              props
    #}
    #bedmethyl_job = dxpy.new_dxjob(inp, "bedmethyl_io")
    #print "* Kicked off bedmethyl() but not waiting waiting..."

    print "* Calling signal()..."
    signal_out = signal_io(extract_out["bedgraph_gz_dxlink"],extract_out["chrom_sizes_dxlink"],target_root,qc_metrics,props)

    print "* Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    #bedmethyl_job.wait_on_done() # Wait because we want the qc_metrics to pass to other jobs.
    #bedmethyl_out = bedmethyl_job.describe()['output']
    print "* Finished."

    return {
        # from extract() 
        #"bam_biorep":    extract_out['biorep_bam_dxlink'], 
        "bam_biorep_qc": extract_out['biorep_bam_qc_dxlink'], 
        "map_biorep":    extract_out['biorep_map_dxlink'],
        "mbias_report":  extract_out["mbias_report_dxlink"],
                
        # from signal() 
        "signal": signal_out["bigWig_dxlink"],
        
        # from bedmethyl() 
        "CpG_bed": bedmethyl_out["CpG_bed_dxlink"],
        "CHG_bed": bedmethyl_out["CHG_bed_dxlink"],
        "CHH_bed": bedmethyl_out["CHH_bed_dxlink"],
        "CpG_bb":  bedmethyl_out["CpG_bb_dxlink"],
        "CHG_bb":  bedmethyl_out["CHG_bb_dxlink"],
        "CHH_bb":  bedmethyl_out["CHH_bb_dxlink"],

        "metadata": json.dumps(qc_metrics) 
        }

dxpy.run()
