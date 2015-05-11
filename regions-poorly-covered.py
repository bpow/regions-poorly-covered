#!/usr/bin/env python2.7

import argparse
import sys
import gzip
import os
import subprocess

def run_depth_of_coverage(args):
    gatk_coverage_output = args.output_basename + '.coverage'
    fifo = None
    if args.fifo:
        os.mkfifo(gatk_coverage_output)
        fifo = subprocess.Popen('gzip < "%s" > "%s.gz"'%(gatk_coverage_output, gatk_coverage_output), shell=True)

    depth_of_coverage_command = [args.java, '-Xmx%s'%args.memory, '-jar', args.gatkjar,
       '-T', 'DepthOfCoverage', '-R', args.reference, '-I', args.bam_files,
       '-L', args.intervals, '-ct', str(args.coverage), '-omitIntervals',
       '--minBaseQuality', str(args.base_quality), '--minMappingQuality', str(args.mapping_quality),
       '-nt', str(args.threads), '-o', gatk_coverage_output]

    logfile = open(args.output_basename + '.coverage.log', 'w', buffering=1)
    doc = subprocess.Popen(depth_of_coverage_command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, bufsize=1)
    for line in iter(doc.stdout.readline, ''):
        if fifo is not None:
            if fifo.poll() is not None:
                raise RuntimeException('gzip process exited before DepthOfCoverage with returncode: %d'%fifo.returncode)
        sys.stdout.write(line)
        logfile.write(line)
    doc.wait()
    logfile.close()
    if args.fifo:
        fifo.wait()
        os.unlink(gatk_coverage_output)
        return gatk_coverage_output + '.gz'
    else:
        return gatk_coverage_output


def percent_coverage_below_threshold(coverage_file, coverage_threshold = 20, percent_threshold = 100):
    if coverage_file.endswith('.gz'):
        inf = gzip.open(coverage_file, 'r')
        coverage_file = coverage_file[:-3] # so we do not keep the '.gz' in the next filename
    else:
        inf = open(coverage_file, 'r')
    output_file = coverage_file + '.DoC.c%d.p%d.txt'%(coverage_threshold, percent_threshold)
    ouf = open(output_file, 'w')
    header = None
    for line in inf:
        line = line.rstrip()
        if header is None:
            header = line.split('\t')
            if header[0:3] != ['Locus', 'Total_Depth', 'Average_Depth_sample']:
                raise RuntimeError("Header is not what I expected:\n\t%s"%'\t'.join(header))
            number_of_samples = len(header) - 3
            header = ['Chrom', 'Position'] + header[1:] + ['N_below_threshold', 'Pct_below_threshold']
            ouf.write('\t'.join(header))
            ouf.write('\n')
        else:
            row = line.split('\t')
            row[0] = '\t'.join(row[0].split(':', 1)) # that is a weird way to do this...
            count_below = len([x for x in row[3:] if int(x) < coverage_threshold])
            pct_below = 100.0*count_below/number_of_samples
            if pct_below > percent_threshold:
                row += [str(count_below), "%.2f"%pct_below]
                ouf.write('\t'.join(row))
                ouf.write('\n')
    inf.close()
    ouf.close()
    return output_file

class BedLine:
    def __init__(self, chrom = None, start = -1, end = -1, total_avg_depth = 0, total_pct = 0):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.total_avg_depth = total_avg_depth
        self.total_pct = total_pct
    def __str__(self):
        return '\t'.join([str(x) for x in (self.chrom, self.start, self.end,
            self.total_avg_depth/(self.end-self.start), self.total_pct/(self.end-self.start))])

def rpc_to_bed(thresholded_file):
    inf = open(thresholded_file, 'r')
    output = thresholded_file + '.bed'
    ouf = open(output, 'w')

    working_output = BedLine()

    ouf.write('#CHROM\tSTART\tEND\tAVG_COV\tAVG_%_BELOW_THRESH\n')

    for line in inf:
        if line.startswith('Chrom\tPosition'):
            continue
        row = line.rstrip().split('\t')
        (chrom, pos) = row[0:2]
        pos = int(pos)
        site_avg = float(row[3])
        site_pct = float(row[-1])
        if chrom == working_output.chrom:
            if pos == working_output.end + 1:
                working_output.end = pos
                working_output.total_avg_depth += site_avg
                working_output.total_pct += site_pct
            else:
                ouf.write(str(working_output))
                ouf.write('\n')
                working_output = BedLine(chrom, pos-1, pos, site_avg, site_pct)
        else:
            if working_output.chrom is not None:
                ouf.write(str(working_output))
                ouf.write('\n')
            working_output = BedLine(chrom, pos-1, pos, site_avg, site_pct)

    ouf.write(str(working_output))
    ouf.write('\n')

    ouf.close()
    inf.close()
    return output

if '__main__' == __name__:
    parser = argparse.ArgumentParser(
        description='''Calculate regions of poor coverage''')
    parser.add_argument('-b', '--bam_files', help='file containing the list of bam files to process (one per line)', required=True)
    parser.add_argument('-i', '--intervals', help='intervals over which to calculate coverage', required=True)
    parser.add_argument('-o', '--output_basename', help='base name for output files', required=True)
    parser.add_argument('-r', '--reference', help='reference fasta file (must have associated ".dict")', required=True)
    parser.add_argument('-c', '--coverage', type=int, default=20,
        help='Cutoff for individual coverage')
    parser.add_argument('-p', '--percent', type=int, default=90,
        help='Only output lines where percentage of samples with coverage < "-c" is above this value')
    parser.add_argument('-q', '--base_quality', type=int, default=20,
        help='Minimum base quality to be considered for coverage calculations (--minBaseQuality in DoC)')
    parser.add_argument('-Q', '--mapping_quality', type=int, default=20,
        help='Minimum mapping quality to be considered for coverage calculations (--minMappingQuality in DoC)')
    parser.add_argument('-m', '--memory', default='18g',
        help='memory to allocate to java (number with suffix, like the java "-Xmx" option')
    parser.add_argument('-t', '--threads', default=1)
    parser.add_argument('-g', '--gatkjar', help='.jar file for the GenomeAnalysisToolkit', default='GenomeAnalysisTK.jar')
    parser.add_argument('-j', '--java', default='java', help='java executable')
    parser.add_argument('-s', '--skipDOC', help='Skip running DepthOfCoverage, and use this file as input for next steps')
    parser.add_argument('-f', '--fifo', action='store_true',
        help='use a FIFO to compress the DoC output as it is produced (this is probably a good idea in most cases)')

    args = parser.parse_args()

    if args.skipDOC is None:
        coverage_file = run_depth_of_coverage(args)
    else:
        coverage_file = args.skipDOC
    rpc_txt = percent_coverage_below_threshold(coverage_file, int(args.coverage), float(args.percent))

    rpc_to_bed(rpc_txt)
