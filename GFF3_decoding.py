# coding=utf-8
"""
This modlue is for decoding  Gff3 file or gtf file
"""
import sys

geneDict = {}
genomeDict = {}


def Isoverlab(listone, listtwo):
    if int(listone[-1]) < int(listtwo[0]):
        listone.extend(listtwo)
    elif int(listone[-1]) > int(listtwo[-1]):
        pass
    else:
        listone[-1] = listtwo[-1]

    return listone


class Isoform:
    def __init__(self, listitem):
        self.IsoformDict = {}
        self.IsoformDict[listitem[2]] = [listitem[3:5]]
        self.Introns = []
        self.exonNum = 1

    def addmore(self, listitem):
        if self.IsoformDict.has_key(listitem[2]):
            self.IsoformDict[listitem[2]].append(listitem[3:5])
        else:
            self.IsoformDict[listitem[2]] = [listitem[3:5]]

    def builtIntron(self):
        if self.IsoformDict.has_key("UTR"):
            self.Introns.extend(self.IsoformDict["UTR"])
        self.Introns.extend(self.IsoformDict["exon"])
        self.Introns.sort(key=lambda x: int(x[0]))
        self.Introns = reduce(Isoverlab, self.Introns)
        self.Introns = [int(x) for x in self.Introns if True]
        #        self.Introns.insert(0, int(self.IsoformDict["transcript"][0][0]))
        self.Introns.insert(0, int(
            min(reduce(lambda x, y: [1000, min(int(x[0]), int(y[0]))], self.IsoformDict["transcript"]))))
        #        self.Introns.append(int(self.IsoformDict["transcript"][0][1]))
        self.Introns.append(max(reduce(lambda x, y: [0, max(int(x[1]), int(y[1]))], self.IsoformDict["transcript"])))
        self.Introns = [[self.Introns[2 * i], self.Introns[2 * i + 1]] for i in range(len(self.Introns) / 2)]

    def getexonNumber(self):
        if self.IsoformDict.has_key('exon'):
            self.exonNum = len(self.IsoformDict['exon'])



class Gene:
    """
    用来存储基因的信息，包括基因名，scaffold,以及起始和终止位点
    """

    def __init__(self, listitems, genename):
        """
        init values
        """
        assert isinstance(listitems, list)
        assert isinstance(genename, str)
        self.scaffold = listitems[0]
        self.start = int(listitems[3])
        self.end = int(listitems[4])
        self.genename = genename


class GeneSubunit(Gene):
    """
    用来存储gene的亚结构，外显子，内含子，utr等结构
    """

    def __init__(self, listitems, genename):
        """
        init values
        Exon Intron
        """
        Gene.__init__(self, listitems, genename)
        self.Isoforms = []
        self.IsoformNum = 0
        self.superIsoform = []
        self.CommonIntrons = []
        self.maxExon = 1
        self.minExon = 10000

    def AddIsoform(self, listitem):
        self.IsoformNum = self.IsoformNum + 1
        self.Isoforms.append(Isoform(listitem))

    def Additems(self, listitem):
        self.Isoforms[-1].addmore(listitem)

    def builtsuperIsoform(self):
        """
        构建出CommonIntrons
        """
        for x in self.Isoforms:
            assert isinstance(x, Isoform)
            x.getexonNumber()
            self.maxExon = max(self.maxExon, x.exonNum)
            self.minExon = min(self.minExon, x.exonNum)
            if x.IsoformDict.has_key('UTR'):
                self.superIsoform.extend(x.IsoformDict["UTR"])
            self.superIsoform.extend(x.IsoformDict['exon'])
        self.superIsoform.sort(key=lambda x: int(x[0]))
        self.superIsoform = map(lambda x: x[0:2], self.superIsoform)
        self.superIsoform = reduce(Isoverlab, self.superIsoform)
        self.superIsoform = [int(x) for x in self.superIsoform if True]
        self.superIsoform.insert(0, self.start)
        self.superIsoform.append(self.end)
        self.CommonIntrons = [[self.superIsoform[2 * i], self.superIsoform[2 * i + 1]] for i in
                              range(len(self.superIsoform) / 2) if True]





        # def AddExons(self, listitem):
        #     self.Exons.append(listitem)
        #
        # def AddFutrs(self, listitem):
        #     self.Futrs.append(listitem)
        #
        # def AddTutrs(self, listitem):
        #     self.Tutrs.append(listitem)
        #
        # def GenerateIntron(self):
        #     self.Exons.sort(key=lambda x: x[0])
        #     tmp=[self.Exons[0]]
        #     for element in self.Exons:
        #
        # def

        # def classifyitems(listitems, Is):


def B_Search(length, pos, scaffold):  # 二分搜索法
    start, end = 0, length
    cursor = (start + end) // 2
    while start != cursor:
        if genomeDict[scaffold][cursor].start > pos:
            end = cursor
        else:
            start = cursor
        cursor = (start + end) // 2
    return cursor


def IsCoverIntron(star, end, gene):
    assert isinstance(gene, GeneSubunit)
    for isoform in gene.Isoforms:
        for Intron in isoform.Introns:
            if star <= Intron[0] <= Intron[1] <= end:
                return True
    return False


def IsoverIntron(start, end, gene):
    assert isinstance(gene, GeneSubunit)
    for Intron in gene.CommonIntrons:
        if start <= Intron[0] <= Intron[1] <= end:
            return True
    return False


def GffPatternDet(start, end, gene):
    """
    IsOver的值为0 表示没有覆盖整个Intron,1表示有Intron,2表示覆盖一个全长的Intron,3表示覆盖全场
    """
    assert isinstance(gene, GeneSubunit)
    IsContinue = False
    IsOver = 0
    Overlab = False
    if end >= gene.end:
        IsContinue = True
        if start >= gene.end:  # gs<=ge<=as<=ae
            pass
        elif gene.start <= start < gene.end:  # gs<=as<=ge<=ae
            Overlab = True
            if IsoverIntron(start, end, gene):
                if IsCoverIntron(start, end, gene):
                    IsOver = 2
                else:
                    IsOver = 1
        else:  # as<=gs<=ge<=ae
            Overlab = True
            IsOver = 3
    elif end <= gene.start:  # as<=ae<=gs<=ge
        pass
    elif start <= gene.start:  # as<=gs<=ae<=ge
        Overlab = True
        if IsoverIntron(start, end, gene):
            if IsCoverIntron(start, end, gene):
                IsOver = 2
            else:
                IsOver = 1
    else:  # gs<=as<=ae<=ge
        Overlab = True
        if IsoverIntron(start, end, gene):
            if IsCoverIntron(start, end, gene):
                IsOver = 2
            else:
                IsOver = 1
    return IsContinue, Overlab, IsOver


# 如果是返回
def IsFullLength(scaffold, start, end, part=False):
    genelist = []
    fulllen = []
    IsContinue = True
    if genomeDict.has_key(scaffold):
        lens = len(genomeDict[scaffold])
        idx = B_Search(lens, start, scaffold)
        while IsContinue and idx < lens:
            IsContinue, IsOverlap, IsOverAll = GffPatternDet(start, end, genomeDict[scaffold][idx])
            if IsOverlap:
                genelist.append(genomeDict[scaffold][idx].genename)
                fulllen.append(IsOverAll)
            idx += 1
        if len(genelist) == 0:
            genelist.append("intergenic")
            fulllen.append(0)
        return genelist, fulllen
    else:
        print("%s scaffold information not in Gff file!", scaffold)
        return None, None


def decodegff(gtffilename):
    """
    根据GFF文件创建gene Isoform
    :rtype: str
    """
    with open(gtffilename) as gtffile:
        tmpgene = None
        for item in gtffile:
            if item.find("#") != 0:
                listitems = item.split("\t")
                # classifyitems(listitems)
                if listitems[2] == 'gene':
                    tmp = listitems[8].split(';')
                    genename = tmp[2].split(' ')[2].replace('"', '')
                    tmpgene = GeneSubunit(listitems[0:5], genename)
                break

        for item in gtffile:
            listitems = item.split("\t")
            # classifyitems(listitems)
            if listitems[2] == 'gene':
                tmpgene.builtsuperIsoform()
                if tmpgene != None and genomeDict.has_key(tmpgene.scaffold):
                    genomeDict[tmpgene.scaffold].append(tmpgene)
                elif tmpgene != None:
                    genomeDict[tmpgene.scaffold] = [tmpgene]
                tmp = listitems[8].split(';')
                genename = tmp[2].split(' ')[2].replace('"', '')
                tmpgene = GeneSubunit(listitems[0:5], genename)
            elif listitems[2] == 'transcript':
                if tmpgene.IsoformNum > 0:
                    tmpgene.Isoforms[-1].builtIntron()
                tmpgene.AddIsoform(listitems)
            else:
                tmpgene.Additems(listitems)
    for key, val in genomeDict.items():
        assert isinstance(val, list)
        val.sort(key=lambda x: x.start)  # 按照起始位点排序


def main():
    print("This is a Test!")


if __name__ == '__main__':
    sys.exit(main())
