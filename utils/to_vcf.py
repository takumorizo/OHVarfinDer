#!/usr/bin/env python
import sys
import re
import pysam

def prepare_vcf_header(ref_path, sample, score):
    fasta_obj =  pysam.FastaFile(ref_path)
    header = ''
    header += '##fileformat=VCFv4.2\n'
    header += '##FILTER=<ID=PASS,Description=\"All filters passed, used BF >= ' + str(score) +  ' as a filter.\">\n'
    header += '##reference=file://' + ref_path + '\n'
    for i in range(len(fasta_obj.lengths)):
        header += '##<ID=' + str(fasta_obj.references[i]) + ',length=' + str(fasta_obj.lengths[i]) + '>\n'
    header += '##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n'
    header += '##FORMAT=<ID=RCT,Number=1,Type=Integer,Description=\"Ref Count Tumor\">\n'
    header += '##FORMAT=<ID=ACT,Number=1,Type=Integer,Description=\"Alt Count Tumor\">\n'
    header += '##FORMAT=<ID=RCN,Number=1,Type=Integer,Description=\"Ref Count Normal\">\n'
    header += '##FORMAT=<ID=ACN,Number=1,Type=Integer,Description=\"Alt Count Normal\">\n'
    header += '##FORMAT=<ID=BF,Number=1,Type=Float,Description=\"Bayes factor\">\n'
    header += '##FORMAT=<ID=BF,Number=1,Type=Float,Description=\"Bayes factor\">\n'
    header += '##FORMAT=<ID=HAP,Number=1,Type=Integer,Description=\"Considered haplotype, 0(not used) or 1(used)\">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + str(sample) +  '\n'
    fasta_obj.close()
    return header

def convert_to_vcf_variant(Chr, start, end, ref, obs, fasta_obj):
    if ref != '-' and obs != '-': # SNV
        pos    = end
        vcf_ref = ref
        vcf_obs = obs
    elif ref == '-': # INS
        pos    = end
        vcf_ref = fasta_obj.fetch(Chr, start-1, start)
        vcf_obs = vcf_ref + obs
    elif obs == '-': # DEL
        pos    = start
        vcf_ref = fasta_obj.fetch(Chr, start-1, start) + ref
        vcf_obs = fasta_obj.fetch(Chr, start-1, start)
    else:
        raise Exception("unexpected mutation")
    return [Chr, pos, '.', vcf_ref, vcf_obs]

def main():
    argvs           = sys.argv
    ref_path        = argvs[1]
    ohvar_path      = argvs[2]
    output_vcf_path = argvs[3]
    sample_name     = argvs[4]
    min_score       = float(argvs[5]) if len(argvs) > 5 else 0.5

    is_first  = True
    fasta_obj =  pysam.FastaFile(ref_path)
    with open(ohvar_path, 'r') as ohvar, open(output_vcf_path, 'w') as output_vcf:
        header = prepare_vcf_header(ref_path, sample_name, min_score)
        output_vcf.writelines(header)
        for line in ohvar:
            output_list = []
            if is_first:
                is_first = False
                continue
            line = line.replace('\n','')
            line = line.replace('\r','')
            cols = re.split('\t',line)

            if cols[22] == '-':
                continue

            Chr, start, end, ref, obs = cols[0], int(cols[1]), int(cols[2]), cols[3], cols[4]
            ref_T, obs_T, ref_N, obs_N, BF = int(cols[5]), int(cols[6]), int(cols[7]), int(cols[8]), float(cols[22])
            HAP = 1 if cols[23] == 'considered_close_hetero_germline_variants' else 0

            QUAL = '.'
            FILTER  = 'PASS' if BF >= min_score else '.'
            INFO = '.'
            FORMAT_TAG = 'RCT:ACT:RCN:ACN:BF:HAP'
            FORMAT = ':'.join(map(str, [ref_T,obs_T,ref_N,obs_N,BF,HAP] ))

            output_list.extend(convert_to_vcf_variant(Chr, start, end, ref, obs, fasta_obj))
            output_list.extend([QUAL, \
                                FILTER, \
                                INFO, \
                                FORMAT_TAG, \
                                FORMAT ])
            outputString = '\t'.join(map(str, output_list))
            outputString = outputString.replace('\n', '')
            outputString = outputString + '\n'
            output_vcf.writelines(outputString)
    fasta_obj.close()

if __name__ == '__main__':
    main()
