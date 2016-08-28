# coding=utf-8
'''
    将文件
'''
import optparse
import sys
import time
import GFF3_decoding as gtf

__author__ = 'Guoliang Lin'
Softwarename = ''
version = '1.0.1'
bugfixs = ''
__date__ = '2016-08-10'

listnone = [-1, -1, -1]
dictstatistics = {}


def statistis(numbers):
    if dictstatistics.has_key(numbers):
        dictstatistics[numbers] += 1
    else:
        dictstatistics[numbers] = 1


def findpatten(genename):
    for xlist in gtf.genomeDict.keys():
        for subgene in gtf.genomeDict[xlist]:
            assert isinstance(subgene, gtf.GeneSubunit)
            if subgene.genename == genename:
                return [subgene.IsoformNum, subgene.maxExon, subgene.minExon]
    return listnone


filename = 'Homo_sapiens.GRCh38.85.gtf'
gtf.decodegff(filename)
intshuzi = 67385
for m in range(0, 6):
    locations = 'SRR0' + str(intshuzi) + '.fastq.sam.bcf_loaction.txt'
    with open(locations) as inputfile:
        with open('SRR0' + str(intshuzi) + "fullen", 'w') as fulleni:
            with open('SRR0' + str(intshuzi) + "has_fullintron", 'w') as hasfullintron:
                with open('SRR0' + str(intshuzi) + "has_intron", 'w') as hasintron:
                    with open('SRR0' + str(intshuzi) + "no_intron", 'w') as nointron:
                        for location in inputfile:
                            location = location.strip()

                            locationsc = location.split('|')[0]
                            # locationreverse = location[::-1]
                            listgener, liststatr = location.split('|')[1].split('@', 1)
                            liststat = liststatr.split('~')
                            listgene = listgener.split('~')
                            for x in range(0, len(liststat)):
                                listpatten = findpatten(listgene[x])
                                statistis(listpatten[1])
                                if int(liststat[x]) == 3:
                                    #                                if listgene[x]!='':
                                    fulleni.write(
                                        listgene[x] + '\t' + '\t'.join(map(str, listpatten)) + '\t' + locationsc + '\n')
                                elif int(liststat[x]) == 2:
                                    if listgene[x] != '':
                                        hasfullintron.write(listgene[x] + '\t' + '\t'.join(
                                            map(str, listpatten)) + '\t' + locationsc + '\n')
                                elif int(liststat[x]) == 1:
                                    if listgene[x] != '':
                                        hasintron.write(listgene[x] + '\t' + '\t'.join(
                                            map(str, listpatten)) + '\t' + locationsc + '\n')
                                else:
                                    if listgene[x] != "intergenic" and listgene[x] != '':
                                        nointron.write(listgene[x] + '\t' + '\t'.join(
                                            map(str, listpatten)) + '\t' + locationsc + '\n')
    intshuzi = intshuzi + 1
