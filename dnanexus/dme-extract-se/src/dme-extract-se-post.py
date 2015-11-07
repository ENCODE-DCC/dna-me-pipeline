#!/usr/bin/env python
# dme-extract-se.py  extract methylation using bismark
#
# 1) Merges techrep bams and map_reports and generated qc metrics in main()
# 2) Extracts methylation... either extract() or:  
#    a) bismark_methylation_extractor in extractor()
#    b) bismark2bedGraph in biz_to_bg()
#    c) coverage2cytosine in coverage()
# 3) generates bigWig in signal()
# 4) generates bedmethyl and bigbeds in bedmethyl()
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

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
        
def bam_or_sam(biorep_bam, biorep_bam_dxlink, uncompress_bam, target_root, caller="extract"):
    '''Decides whether to extract from bam or uncompressed sam and returns the file name.'''
    # NOTE: Better to use sam and let extractor use more threads, but this takes up precious storage
    alignments=biorep_bam
    ncores=32
    # TODO: What bam_size constitutes too large for sam?  93,953,130,496 is fine!
    bam_size = dxpy.describe(biorep_bam_dxlink)['size']
    if uncompress_bam and bam_size < 400000000000:
        alignments = target_root + ".sam"
        print "* extract(): Decompressing biorep bam (size: %d)..." % bam_size
        run_cmd('samtools view ' + biorep_bam, out=alignments)
        run_cmd('rm -f ' + biorep_bam) # STORAGE IS LIMITED
    else:
        ncores = ncores / 2
        if uncompress_bam:
            print "* %s(): Using compressed bam and %d cores because bam_size: %d exceeds limit." % (caller,ncores,bam_size)
        else:
            print "* %s(): Using compressed biorep bam (size: %d) and %d cores..." % (caller,bam_size, ncores)
    return (alignments, ncores)

###### Bismark extraction
def bismark_simple_extract(target_root, alignments, ncores):
    '''bismark_methylation_extractor without coverage.'''
     
    print "* extractor(): Analyse methylation in %s and using %d threads..." % (alignments, ncores)
    # TODO: missing: --no_overlap/--include_overlap --ignore_XXX
    # NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
    cmd = 'bismark_methylation_extractor --multicore %d --single-end -s --comprehensive -output output/ %s' % (ncores, alignments)
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

@dxpy.entry_point("extractor")
def extractor(biorep_bam_dxlink, dme_ix_dxlink, uncompress_bam, target_root, qc_metrics, props):
    '''subjob runs bismark_methylation_extractor on mem1_hdd2_x32'''

    print "* extractor(): Retrieve merged bam and index..."
    biorep_bam = target_root + ".bam"
    dme_ix = "dme_index.tar.gz"
    dxpy.download_dxfile(biorep_bam_dxlink, biorep_bam)
    dxpy.download_dxfile(dme_ix_dxlink, dme_ix)

    print "* extractor(): Uncompress index..."
    run_cmd('tar -zxf ' + dme_ix)

    # NOTE: Better to use sam and let extractor use more threads, but this takes up precious storage
    (alignments, ncores) = bam_or_sam(biorep_bam, biorep_bam_dxlink, uncompress_bam, target_root, "extractor")

    bismark_simple_extract(target_root, alignments, ncores)
    qc_metrics = bismark_qc_metrics(target_root, qc_metrics)

    print "* extractor(): Storing extraction results..."
    split_report_dxfile = dxpy.upload_local_file(target_root+'_splitting_report.txt')
    chrom_sizes_dxfile = dxpy.upload_local_file('input/chrom.sizes')
    mbias_report_dxfile = dxpy.upload_local_file(target_root+'_mbias_report.txt',properties=props,details=qc_metrics)
    CpG_context_dxfile = dxpy.upload_local_file('output/CpG_context_%s.txt' % (target_root))
    CHG_context_dxfile = dxpy.upload_local_file('output/CHG_context_%s.txt' % (target_root))
    CHH_context_dxfile = dxpy.upload_local_file('output/CHH_context_%s.txt' % (target_root))

    print "* extractor(): Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    return {
        "CpG_context_dxlink":   dxpy.dxlink(CpG_context_dxfile),
        "CHG_context_dxlink":   dxpy.dxlink(CHG_context_dxfile),
        "CHH_context_dxlink":   dxpy.dxlink(CHH_context_dxfile),
        "split_report_dxlink":  dxpy.dxlink(split_report_dxfile),
        "chrom_sizes_dxlink":   dxpy.dxlink(chrom_sizes_dxfile),
        "mbias_report_dxlink":  dxpy.dxlink(mbias_report_dxfile),
        "qc_metrics":           qc_metrics
    }

###### all but extract
def bismark_coverage(target_root, CpG_context, CHG_context, CHH_context, gzip=True, cleanup=False):
    '''bismark2bedGraph.'''
     
    print "* bismark_coverage(): Generating bedGraph from context files..."
    bedGraph = '%s.bedGraph' % target_root
    cmd = 'bismark2bedGraph --CX_context --ample_mem  --zero --dir output/ -output '
    cmd += '%s %s %s %s' % (bedGraph, CpG_context, CHG_context, CHH_context)
    run_cmd(cmd)
    run_cmd('ls -l output/')
    if os.path.isfile('output/'+bedGraph+'.gz'):
        run_cmd("mv output/%s.gz %s.gz" % (bedGraph,bedGraph))
        if gzip:
            bedGraph = bedGraph + '.gz'
        else:
            run_cmd('gunzip %s.gz' % bedGraph)
    elif os.path.isfile('output/'+bedGraph):
        run_cmd("mv output/%s %s" % (bedGraph,bedGraph))
        if gzip:
            run_cmd('gzip ' + bedGraph)
            bedGraph = bedGraph + '.gz'

    print "* bismark_coverage(): Generating CX_context file..."
    cov_file = '%s.bismark.cov.gz' % (target_root)
    cx_report = '%s.CX_report.txt' % (target_root)
    cmd = 'coverage2cytosine --genome_folder input --CX_context --zero_based --dir output/ --output %s %s' \
                                                                                                % (cx_report, cov_file)
    run_cmd(cmd)
    run_cmd("mv output/%s %s" % (cx_report,cx_report))
    run_cmd('ls -l output/')
    if cleanup:
        run_cmd("rm -rf input/")
        run_cmd("rm -rf output/")

    return (bedGraph, cx_report)

def bedmethyl(target_root, cx_report, chrom_sizes, cleanup=False):
    '''runs cxrepo-bed.py and bedToBigBed.'''
     
    print "* bedmethyl(): Run cxrepo-bed.py..."
    run_cmd('cxrepo-bed.py -o output/ output/' + cx_report)
    run_cmd('mv output/CG_%s  %s_CpG.bed' % (cx_report,target_root))
    run_cmd('mv output/CHG_%s %s_CHG.bed' % (cx_report,target_root))
    run_cmd('mv output/CHH_%s %s_CHH.bed' % (cx_report,target_root))
    if cleanup:
        run_cmd("rm -rf output/")

    print "* bedmethyl(): Convert to BigBed..."
    run_cmd('bedToBigBed %s_CpG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 %s %s_CpG.bb' % (target_root,chrom_sizes,target_root))
    run_cmd('bedToBigBed %s_CHG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 %s %s_CHG.bb' % (target_root,chrom_sizes,target_root))
    run_cmd('bedToBigBed %s_CHH.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 %s %s_CHH.bb' % (target_root,chrom_sizes,target_root))
    CpG_root = '%s_CpG' % target_root
    CHG_root = '%s_CHG' % target_root
    CHH_root = '%s_CHH' % target_root
    run_cmd('gzip %s.bed' % CpG_root)
    run_cmd('gzip %s.bed' % CHG_root)
    run_cmd('gzip %s.bed' % CHH_root)
    return (CpG_root + '.bed.gz',CHG_root + '.bed.gz',CHH_root + '.bed.gz',CpG_root + '.bb',CHG_root + '.bb',CHH_root + '.bb')

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

def post_extraction(CpG_context_dxlink, CHG_context_dxlink, CHH_context_dxlink, dme_ix_dxlink, target_root, qc_metrics, props):
    '''runs everything after bismark simple extraction in the main instance'''

    print "* post_extraction(): Retrieve context files and index..."
    CpG_context = 'CpG_context_%s.txt' % target_root
    CHG_context = 'CHG_context_%s.txt' % target_root
    CHH_context = 'CHH_context_%s.txt' % target_root
    run_cmd('mkdir -p output/')
    dxpy.download_dxfile(CpG_context_dxlink, 'output/' + CpG_context)
    dxpy.download_dxfile(CHG_context_dxlink, 'output/' + CHG_context)
    dxpy.download_dxfile(CHH_context_dxlink, 'output/' + CHH_context)
    dme_ix = "dme_index.tar.gz"
    dxpy.download_dxfile(dme_ix_dxlink, dme_ix)

    print "* post_extraction(): Uncompress index..."
    run_cmd('tar -zxf ' + dme_ix)

    # First coverage:
    (bedGraph, cx_report) = bismark_coverage(target_root, CpG_context, CHG_context, CHH_context, gzip=False, cleanup=True)
    
    # Next beds
    chrom_sizes = "input/chrom.sizes"    
    (CpG_bed,CHG_bed,CHH_bed,CpG_bb,CHG_bb,CHH_bb) = bedmethyl(target_root, cx_report, chrom_sizes, cleanup=True)
    
    # Finally signal
    bigWig = signal(target_root, bedGraph, chrom_sizes, cleanup=True)

    print "* post_extraction(): Storing results..."
    CpG_bed_dxfile = dxpy.upload_local_file(CpG_bed,properties=props,details=qc_metrics)
    CHG_bed_dxfile = dxpy.upload_local_file(CHG_bed,properties=props,details=qc_metrics)
    CHH_bed_dxfile = dxpy.upload_local_file(CHH_bed,properties=props,details=qc_metrics)

    CpG_bb_dxfile = dxpy.upload_local_file(CpG_bb,properties=props,details=qc_metrics)
    CHG_bb_dxfile = dxpy.upload_local_file(CHG_bb,properties=props,details=qc_metrics)
    CHH_bb_dxfile = dxpy.upload_local_file(CHH_bb,properties=props,details=qc_metrics)

    bigWig_dxfile = dxpy.upload_local_file(bigWig,properties=props,details=qc_metrics)

    print "* post_extraction(): Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    return {
        "CpG_bed_dxlink": dxpy.dxlink(CpG_bed_dxfile),
        "CHG_bed_dxlink": dxpy.dxlink(CHG_bed_dxfile),
        "CHH_bed_dxlink": dxpy.dxlink(CHH_bed_dxfile),
        "CpG_bb_dxlink":  dxpy.dxlink(CpG_bb_dxfile),
        "CHG_bb_dxlink":  dxpy.dxlink(CHG_bb_dxfile),
        "CHH_bb_dxlink":  dxpy.dxlink(CHH_bb_dxfile),
        "bigWig_dxlink":  dxpy.dxlink(bigWig_dxfile)
    }


@dxpy.entry_point("main")
def main(bam_set, map_report_set, dme_ix, uncompress_bam=True, nthreads=8):

    # tool_versions.py --applet $script_name --appver $script_ver
    props = {}
    if os.path.isfile('/usr/bin/tool_versions.py'): 
        sw_versions = subprocess.check_output(['tool_versions.py', '--dxjson', 'dnanexus-executable.json'])
        props["SW"] = sw_versions
    
    print "* Value of bam_set:        '" + str(bam_set) + "'"
    print "* Value of map_report_set: '" + str(map_report_set) + "'"
    print "* Value of dme_ix:         '" + str(dme_ix) + "'"
    print "* Value of uncompress_bam: '" + str(uncompress_bam) + "'"
    print "* Value of nthreads:       '%d" % nthreads

    # NOTE: dme-align produces *_techrep_bismark.bam and dme-extract merges 1+ techrep bams into a *_bismark_biorep.bam.
    #       The reason for the name 'word' order is so thal older *_bismark.bam alignments are recognizable as techrep bams

    target_root = ""
    merged = ""
    tech_reps = ""
    exp_id = ""
    rep_tech = ""
    
    for techrep_bam_dlink in bam_set:
        file_desc = dxpy.describe(techrep_bam_dlink)
        file_root = file_desc['name']
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
            if len(exp_id) > 0:
                rep_tech = subprocess.check_output(shlex.split('parse_property.py -f %s --project %s --rep_tech -q' \
                                                                        % (file_desc['id'], file_desc['project']) ))
        if len(rep_tech) > 0:
            if len(tech_reps) > 0:
                tech_reps += '_' + rep_tech
            else:
                tech_reps = rep_tech
                
        print "* Downloading %s_techrep_bismark.bam file..." % file_root
        dxpy.download_dxfile(techrep_bam_dlink, file_root + '_techrep_bismark.bam')
    
        if not os.path.isfile("sofar.bam"): 
            run_cmd('mv %s_techrep_bismark.bam sofar.bam' % file_root)
        else:
            print "* Merging in %s_techrep_bismark.bam..." % file_root
            # NOTE: keeps the first header
            run_cmd('samtools cat sofar.bam %s_techrep_bismark.bam > merging.bam' % file_root)
            run_cmd('mv merging.bam sofar.bam')
            run_cmd('rm %s_techrep_bismark.bam' % file_root) # STORAGE IS LIMITED
                

    if len(exp_id) > 0 and len(tech_reps) > 0:
        target_root = '%s_%s_bismark_biorep' % (exp_id, tech_reps)

    # At this point there is a 'sofar.bam' with one or more input bams
    if len(merged) == 0:
        target_root = file_root + "_bismark_biorep"
        run_cmd('mv sofar.bam %s.bam' % target_root)
        print "* Only one input file '%s.bam', no merging required." % target_root
    else:
        # sorting needed due to samtools cat
        print "* Sorting merged bam..."
        run_cmd('samtools sort -@ %d 6G -f sofar.bam %.bam' % (nthreads,target_root) )
        run_cmd('rm sofar.bam') # STORAGE IS LIMITED
        print "* Files merged into '%s.bam'" % target_root

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
        os.system('cat %s >> %s' % (techrep_map_report,biorep_map_report))
        if len(all_reports) == 0:
            all_reports = techrep_map_report
        else:
            all_reports += ',' + techrep_map_report
        
    if all_reports == techrep_map_report: # only one
        run_cmd('cp %s %s' % (techrep_map_report,biorep_map_report) )
    
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

    print "* Storing biorep bam..."
    props_ex = props.copy()
    props_ex.update({ 'reads': str(reads) })
    biorep_bam_dxfile = dxpy.upload_local_file(target_root+'.bam',properties=props_ex,details=qc_metrics,wait_on_close=True)

    #print "* Calling extract()..."
    print "* Calling extractor()..."
    inp = {
        'biorep_bam_dxlink': dxpy.dxlink(biorep_bam_dxfile),
        'dme_ix_dxlink':     dme_ix,
        'uncompress_bam':    uncompress_bam,
        'target_root':       target_root,
        'qc_metrics':        qc_metrics,
        'props':             props
    }
    extract_job = dxpy.new_dxjob(inp, "extractor")
    print "* Kicked off extract() and waiting..."
    extract_job.wait_on_done() # Wait because we want the qc_metrics to pass to other jobs.

    print "* Additional metadata after extract()..."
    extract_out = extract_job.describe()['output']
    qc_metrics = extract_out['qc_metrics']
    
    print "* Retrieve split report..."
    dxpy.download_dxfile(extract_out["split_report_dxlink"], target_root + '_splitting_report.txt')
    append_line("\n===== bismark_methylation_extractor: splitting_report =====",target_root+'_qc.txt')
    run_cmd('cat %s_splitting_report.txt' % target_root,out=target_root+'_qc.txt',append=True,silent=True)

    print "* Calling post_extraction()..."
    post_extraction_out = post_extraction(extract_out["CpG_context_dxlink"], \
                                          extract_out["CHG_context_dxlink"], \
                                          extract_out["CHH_context_dxlink"], \
                                          dme_ix, target_root, qc_metrics, props)

    print "* Check storage..."
    run_cmd('ls -l')
    run_cmd('df -k .')

    print "* Uploading files..."
    biorep_bam_qc_dxlink     = dxpy.upload_local_file(target_root + '_qc.txt',        properties=props,details=qc_metrics)
    biorep_map_report_dxlink = dxpy.upload_local_file(target_root + '_map_report.txt',properties=props,details=qc_metrics)

    print "* Finished."

    return {
        # from main() 
        #"bam_biorep":    biorep_bam_dxlink, 
        "bam_biorep_qc": biorep_bam_qc_dxlink, 
        "map_biorep":    biorep_map_report_dxlink,
        
        # from extract() 
        "mbias_report": extract_job.get_output_ref("mbias_report_dxlink"),
        
        # from post_extraction() 
        "signal": post_extraction_out["bigWig_dxlink"],
        
        "CpG_bed": post_extraction_out["CpG_bed_dxlink"],
        "CHG_bed": post_extraction_out["CHG_bed_dxlink"],
        "CHH_bed": post_extraction_out["CHH_bed_dxlink"],
        
        "CpG_bb": post_extraction_out["CpG_bb_dxlink"],
        "CHG_bb": post_extraction_out["CHG_bb_dxlink"],
        "CHH_bb": post_extraction_out["CHH_bb_dxlink"],

        "metadata": json.dumps(qc_metrics) 
        }

dxpy.run()
