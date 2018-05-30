#!/usr/bin/env python
import sys
import re
import pysam

def prepareBasicVCFHeader(refPath):
    fastaObj     =  pysam.FastaFile(refPath)
    outputHeader = ""
    outputHeader += "##fileformat=VCFv4.2\n"
    outputHeader += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    outputHeader += '##reference=file://' + refPath + '\n'
    for i in range(len(fastaObj.lengths)):
        outputHeader += '##<ID=' + str(fastaObj.references[i]) + ',length=' + str(fastaObj.lengths[i]) + '>\n'
    outputHeader += '##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n'
    outputHeader += '##INFO=<ID=RCT,Number=1,Type=Integer,Description=\"Ref Count Tumor\">\n'
    outputHeader += '##INFO=<ID=ACT,Number=1,Type=Integer,Description=\"Alt Count Tumor\">\n'
    outputHeader += '##INFO=<ID=RCN,Number=1,Type=Integer,Description=\"Ref Count Normal\">\n'
    outputHeader += '##INFO=<ID=ACN,Number=1,Type=Integer,Description=\"Alt Count Normal\">\n'
    outputHeader += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    fastaObj.close()
    return outputHeader

def convertToVCFVariant(Chr, start, end, ref, obs, fastaObj):
    if ref != '-' and obs != '-': # SNV
        pos    = end
        vcfRef = ref
        vcfObs = obs
    elif ref == '-': # INS
        pos    = end
        vcfRef = fastaObj.fetch(Chr, start-1, start)
        vcfObs = vcfRef + obs
    elif obs == '-': # DEL
        pos    = start
        vcfRef = fastaObj.fetch(Chr, start-1, start) + ref
        vcfObs = fastaObj.fetch(Chr, start-1, start)
    else:
        raise Exception("unexpected mutation")
    return [Chr, pos, '.', vcfRef, vcfObs]

def main():
    argvs         = sys.argv
    refPath       = argvs[1]
    ohvarPath     = argvs[2]
    outputVCFPath = argvs[3]

    isFirst = True
    fastaObj     =  pysam.FastaFile(refPath)
    with open(ohvarPath, 'r') as ohvar, open(outputVCFPath, 'w') as outputVCF:
        headerStr = prepareBasicVCFHeader(refPath)
        outputVCF.writelines(headerStr)

        for line in ohvar:
            outputList = []
            if isFirst:
                isFirst = False
                continue
            line = line.replace('\n','')
            line = line.replace('\r','')
            lineCols = re.split('\t',line)

            Chr, start, end, ref, obs = lineCols[0], int(lineCols[1]), int(lineCols[2]), lineCols[3], lineCols[4]
            refT, obsT, refN, obsN, BFScore = int(lineCols[5]), int(lineCols[6]), int(lineCols[7]), int(lineCols[8]), float(lineCols[22])
            outputList.extend(convertToVCFVariant(Chr, start, end, ref, obs, fastaObj))
            outputList.extend([BFScore, '.', 'RCT='+str(refT)+';ACT='+str(obsT)+';RCN='+str(refN)+';ACN='+str(obsN)])
            outputString = '\t'.join(map(str, outputList))
            outputString = outputString.replace('\n', '')
            outputString = outputString + '\n'
            outputVCF.writelines(outputString)
    fastaObj.close()

if __name__ == '__main__':
    main()


