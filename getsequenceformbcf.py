# coding=utf-8
__author__ = 'Guoliang Lin'
Softwarename = 'get sequence from bcf with pipeline'
version = '1.0.1'
bugfixs = ''
__date__ = '2016-03-02'
import optparse
import sys
import time
import base_merge
import GFF3_decoding


def printinformations():
    print("%s software version is %s in %s" % (Softwarename, version, __date__))
    print(bugfixs)
    print('Starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))


def programends():
    print('Ends at :' + time.strftime('%Y-%m-%d %H:%M:%S'))


class Loaction:
    """store loaction for fpkm intems"""

    def __init__(self, string):
        """init values"""
        assert isinstance(string, str)
        self.schaffold = string.split(':')[0]
        self.start = int(string.split(':')[1].split('-')[0])
        self.end = int(string.split(':')[1].split('-')[1])


def _parse_args():
    """Parse the command line for options."""
    usage = 'usage: %prog -i FILE.bcf -g FILE.gff3 -o OUTPREFIX'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i',
                      '--input', dest='input', type='string',
                      help='input segment.out file ')
    #    parser.add_option('-f','--fpkm',dest='fpkm_file',type='string',help='input fpkm file')
    #    parser.add_option('-v','--variation', dest='variation', type='string', help='input variation information file')
    parser.add_option('-g', '--gff3', dest='gff', help='gff3 file')
    parser.add_option('-o', '--output', dest='output', type='string', help='input variation information file')
    options, args = parser.parse_args()
    # positional arguments are ignored
    return options


def TrimHead(bcffile):  # 去掉带#号的部分
    for line in bcffile:
        if (line.find('#') == 0):
            pass
        else:
            break
    return line, bcffile


def writetodisk(file, liststring, segment, chromsome, outlocal, fulllenout):  # 写入磁盘
    """

    :rtype: object
    """
    genelist, fulllen = GFF3_decoding.IsFullLength(chromsome, int(segment[0]), int(segment[1]))
    if genelist != None:
        file.write(">" + chromsome + ':' + segment[0] + '-' + segment[1] + '|' + '_'.join(genelist) + '-' + '_'.join(
                [str(vluz) for vluz in fulllen]) + "\n")
        outlocal.write(
            ">" + chromsome + ':' + segment[0] + '-' + segment[1] + '|' + '_'.join(genelist) + '-' + '_'.join(
                    [str(vluz) for vluz in fulllen]) + "\n")
        file.write(''.join(liststring) + '\n')
        tmp = []
        for m in range(len(genelist)):
            tmp.append([genelist[m], fulllen[m]])
        if sum(fulllen) > 0:
            fulllenout.write('_'.join([k[0] for k in filter(lambda x: x[1] > 0, tmp)]) + '\n')


def iscontinuity(twoNumber):  # 判断是否连续 连续则返回True,
    """

    :type twoNumber: list
    """
    if len(twoNumber) == 1 or (
            (int(twoNumber[-1][-1]) - int(twoNumber[-2][-1])) == 1 and twoNumber[-1][0] == twoNumber[-2][0]):
        return True, True
    elif (int(twoNumber[-1][-1]) - int(twoNumber[-2][-1])) == 0 and twoNumber[-1][0] == twoNumber[-2][0]:

        return True, False
    else:
        return False, False


def getsinglebase(refalt, twoNumber):
    """
    :type refalt: list
    """
    m = ''
    # n=True
    # p=True
    assert isinstance(refalt[1], str)
    refalt[1] = refalt[1].replace('X', refalt[0])
    n, p = iscontinuity(twoNumber)
    if refalt[1].find(',') != -1:
        listforset = refalt[1].split(',')
        for element in listforset:
            m = m + element[0]
        # n,p =iscontinuity(twoNumber)
        m = m.replace('N', '')
        if p:
            return base_merge.merge(frozenset(m)), n
        else:
            return "", n
    else:
        if p:
            return refalt[1][0], n
        else:
            return '', n

        # if refalt[1].find('X')==-1:
        #     return refalt[1][0]
        # else:
        #     return refalt[0][0]


def addtosequence(item, twoNumber):
    """
    :type twoNumber: list
    :type sequence: str
    :type item: str
    """
    tmp = item.split()
    while len(twoNumber) >= 2:
        twoNumber.pop(0)
    twoNumber.append(tmp[0:2])
    return getsinglebase(tmp[3:5], twoNumber)


def getseqence(bcffile, line, ouputname):  # 处理数据
    outfst = open(ouputname + '.fasta', 'w')
    outlocal = open(ouputname + '_loaction.txt', 'w')
    fullenout = open(ouputname + '_fulllen.txt', 'w')
    sequence = []
    segment = []
    chromsome = ''
    twoNumber = []
    initbase, flag = addtosequence(line, twoNumber)
    sequence.append(initbase)
    segment.append(twoNumber[-1][1])
    for item in bcffile:
        initbase, flag = addtosequence(item, twoNumber)
        if flag:
            sequence.append(initbase)
        else:
            segment.append(twoNumber[-2][1])
            writetodisk(outfst, sequence, segment, twoNumber[-2][0], outlocal, fullenout)
            sequence = [initbase]
            segment = [twoNumber[-1][1]]
    fullenout.close()
    outlocal.close()
    outfst.close()


def openbcf(bcffilename, outputname):
    if bcffilename == None:
        bcffile = sys.stdin
    else:
        bcffile = open(bcffilename)
    line, bcffile = TrimHead(bcffile)  # 去除带#行数
    getseqence(bcffile, line, outputname)
    if bcffilename != "":
        bcffile.close()


def main():
    printinformations()
    fpkmdict = {}
    options = _parse_args()
    GFF3_decoding.decodegff(options.gff)
    openbcf(options.input, options.output)
    programends()


if __name__ == '__main__':
    sys.exit(main())
