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


def printinformations():
    print("%s software version is %s in %s" % (Softwarename, version, __date__))
    print(bugfixs)
    print('Starts at :' + time.strftime('%Y-%m-%d %H:%M:%S'))


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


def TrimHead(bcffile): #去掉带#号的部分
    for line in bcffile:
        if (line.find('#') == 0):
            pass
        else:
            break
    return line, bcffile


def writetodisk(file,liststring,segment,chromsome,outlocal): #写入磁盘
    """

    :rtype: object
    """
    file.write(">"+chromsome+':'+segment[0]+'-'+segment[1]+'\n')
    outlocal.write(">"+chromsome+':'+segment[0]+'-'+segment[1]+'\n')
    file.write(''.join(liststring)+'\n')

def iscontinuity(twoNumber):#判断是否连续
    """

    :type twoNumber: list
    """
    if len(twoNumber)==1 or ((int(twoNumber[-1][-1])-int(twoNumber[-2][-1]))==1 and twoNumber[-1][0]==twoNumber[-2][0]):
        return True
    else:
        return False

def getsinglebase(refalt):
    """
    :type refalt: list
    """
    m=[]
    n=''
    assert isinstance(refalt[1], str)
    refalt[1]=refalt[1].replace('X',refalt[0])
    if refalt[1].find(',')!=-1:
        listforset=refalt[1].split(',')
        for element in listforset:
            m.append(element[0])
        return base_merge.merge(frozenset(m))
    else:
        return refalt[1][0]

    # if refalt[1].find('X')==-1:
    #     return refalt[1][0]
    # else:
    #     return refalt[0][0]

def addtosequence(item,twoNumber):
    """
    :type twoNumber: list
    :type sequence: str
    :type item: str
    """
    tmp=item.split()
    while len(twoNumber)>=2:
        twoNumber.pop(0)
    twoNumber.append(tmp[0:2])
    return getsinglebase(tmp[3:5]),iscontinuity(twoNumber)



def getseqence(bcffile, line, ouputname): #处理数据
    outfst=open(ouputname+'.fasta','w')
    outlocal=open(ouputname+'_loaction.txt','w')
    sequence = []
    segment=[]
    chromsome=''
    twoNumber=[]
    initbase,flag=  addtosequence(line, twoNumber)
    sequence.append(initbase)
    segment.append(twoNumber[-1][1])
    for item in bcffile:
        initbase,flag =addtosequence(item, twoNumber)
        if flag :
            sequence.append(initbase)
        else:
            segment.append(twoNumber[-2][1])
            writetodisk(outfst,sequence,segment,twoNumber[-2][0],outlocal)
            sequence = [initbase]
            segment=[twoNumber[-1][1]]




def openbcf(bcffilename, outputname):
    if bcffilename == "":
        bcffile = sys.stdin
    else:
        bcffile = open(bcffilename)
    line, bcffile = TrimHead(bcffile)  # 去除带#行数
    getseqence(bcffile, line, outputname)
   # bcffile.closed()


def main():
    fpkmdict = {}
    options = _parse_args()
    openbcf(options.input,options.output)


if __name__ == '__main__':
    sys.exit(main())
