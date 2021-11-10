from GTFclass import GTFinfo
import os, sys
import argparse
pythonVersion = eval('.'.join(sys.version.split('.')[:2]))

def get_args():
    ArgsParser = argparse.ArgumentParser(
        prog='GTFParser',
        usage='parse GTF file to get BED format files of different features',
        description=
        '''
 _______ _______ ______ ______         _______  ______ ________ ______
|  _____|__   __|  ____|  ___ |  /\   |_ ___  ||  ____| _____|| _____|
| |  ___   | |  | |___ |  >_> | /||\   | >_> _|| |_____ |_____| |__ _|
| | |_  |  | |  |  ___||  ___/ / /\ \  |  _ /  |____  |  ____||  _ |
| |___| |  | |  | |    | |    / /__\ \_| | \ \  ____| | |_____| | \ \ 
|_______|  |_|  |_|    |_|   /_/-__-\_\\\_|  \_\|______|_______|_|  \_\ ''',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help=True,
        epilog='This is a basic draft to run with specific gtf format. Welcome any useful advice! Let\'s strengthen it together'
    )
    ArgsParser.add_argument('-g','--gtf',
        help='provide a gtf file', required=True
    )
    ArgsParser.add_argument('--demand-features', type=str, default= 'gene exon intron', help='select the features you want. Use space as delimiter.\n\
            gene, exon, intron in default')
    ArgsParser.add_argument('-tss','--promoter-range', type=str,default='-2000,0', help='setting promoter range manually. Default -2000,0 (5`->3` direction).')
    ArgsParser.add_argument('-tts','--TTS-range', type=str,default='0,2000', help='setting terminal range manually. Default 0,2000 (5`->3` direction).')
    ArgsParser.add_argument('-5','--5SS-range', type=str,default='-3,8', help='5\' splice site range. Default -3,8 relative to 5\' exon-intron boundry (5`->3` direction).')
    ArgsParser.add_argument('-3','--3SS-range', type=str,default='-3,3', help='3\' splice site range. Default -3,3 relative to 3\'(5`->3` direction).')
    ArgsParser.add_argument('-ppt','--PPT-range', type=str,default='-30,-3', help='Polypyrimidine tract range. Default -20,-3 (5`->3` direction).')

    ArgsParser.add_argument('--split-output', action='store_true', help='Individual file for each feature.\n\
        In default output all feature records to stdout together')
    ArgsParser.add_argument('--split-prefix', type=str, default='GTF', help='prefix of individual file. Only function when spliting output')
    ArgsParser.add_argument('-u', '--uniq', action='store_true', help='Unique all individual files to minimize file size.\n\
        Only function when spliting output. In default not unique.')
    ArgsParser.add_argument('--no-overlapping',action='store_true',help='Except feature \"gene\", no overlapping region exists between any two features ')
    ArgsParser.add_argument('-v','--verbose', action='store_true',help='print time information')
    # show version
    ArgsParser.add_argument('--version', action='version', version='%(prog)s 1.0')
    ArgsDict = vars(ArgsParser.parse_args())

    return ArgsDict

def init_list():
    pass

################################
# [es1, ee1, es2, ee2, es3, ee3] (1th e: exon; 2th s: start; 2th e: end)
# splice site range definition: 1. 5'ss: 3bp in upstream exon and 8bp in dostream intron;
#                               2. 3'ss: 3bp in downstream exon, and 3bp in upstream intron
#                               3. PPT: -20 ~ -3 bp upsteam of 3'ss
#
################################
def classify_features(ExonPosList,strand='+',PromoterRange=(-2000,0),TTSRange=(0,+2000),SS5range=(-3,8),SS3range=(-3,3),PPTrange=(-20,-3)):
    LenExonList = len(ExonPosList)

    popbottom = ExonPosList[0]
    poptop = ExonPosList[-1]

    SSList = ExonPosList[1:-1]
    LenSSList = len(SSList)

    if strand == '+':
        # take care of promoter or TTS region out of gene body
        exonStart = max(PromoterRange[1] + popbottom, popbottom)
        exonEnd = min(TTSRange[0] + poptop, poptop)

        # here wanna implement no overlap function, but the files of exon and intron only contain several hundred records (unresoved problems). 
        # And no need to achieve this goal.

        #exonList = [exonStart] + \
        #    [ SSList[no] + SS5range[0] if no%2==0 else SSList[no] + SS3range[1] for no in range(LenSSList) ] \
        #        + [exonEnd]

        #intronList = [ SSList[no] + SS5range[1] if no%2==0 else SSList[no] + PPTrange[0] for no in range(LenSSList) ]

        promoter_pos = [PromoterRange[0] + popbottom, PromoterRange[1] + popbottom]
        TTS_pos = [TTSRange[0] + poptop, TTSRange[1] + poptop]

        intronNo = [ no+1 for no in range(int((LenSSList)/2)) ]
        exonNo = [ no+1 for no in range(int(LenExonList/2)) ]

        # ExonPosList here equals IntronPosList
        # much easier to understand basing on a certain splice site
        ss5 = [ (SSList[no]+SS5range[0], SSList[no]+SS5range[1]) for no in range(LenSSList) if no%2==0 ]
        ss3 = [ (SSList[no]+SS3range[0], SSList[no]+SS3range[1]) for no in range(LenSSList) if no%2==1 ]
        ppt = [ (SSList[no]+PPTrange[0], SSList[no]+PPTrange[1]) for no in range(LenSSList) if no%2==1 ]

    else:
        exonStart = min(poptop - PromoterRange[1], popbottom)
        exonEnd = max(popbottom - TTSRange[0], poptop)

        #exonList = [exonStart] + \
        #   [ SSList[no] - SS5range[0] if no%2==1 else SSList[no] - SS3range[1] for no in range(LenSSList) ]\
        #         + [exonEnd]

        #intronList = [ SSList[no] - SS5range[1] if no%2==1 else SSList[no] - PPTrange[0] for no in range(LenSSList) ]

        promoter_pos = [poptop - PromoterRange[1] , poptop - PromoterRange[0]]
        TTS_pos = [popbottom - TTSRange[1], popbottom - TTSRange[0]]

        intronNo = sorted([ no+1 for no in range(int(LenSSList/2)) ], reverse=True)
        exonNo = sorted([ no+1 for no in range(int(LenExonList/2)) ], reverse=True)

        ss5 = [ (SSList[no]-SS5range[1], SSList[no]-SS5range[0]) for no in range(LenSSList) if no%2==1 ]
        ss3 = [ (SSList[no]-SS3range[1], SSList[no]-SS3range[0]) for no in range(LenSSList) if no%2==0 ]
        ppt = [ (SSList[no]-PPTrange[1], SSList[no]-PPTrange[0]) for no in range(LenSSList) if no%2==0 ]

    if LenExonList == 2:
        intronNo = False
        ss5 = False
        ss3 = False
        ppt = False
        intronList = False
    else:
        intronList = SSList

    return promoter_pos, TTS_pos, intronNo, intronList, exonNo, ss5, ss3, ppt

def output_stdout(ExonPosList,record_model,promoter_up,promoter_down,TTS_up,TTS_down,SS5range=None,SS3range=None,PPTrange=None,strand='+'):
    promoter_pos,TTS_pos,intronNo,intronList,exonNo,SS5List,SS3List, PPTList = \
            classify_features(ExonPosList, strand, PromoterRange=(promoter_up, promoter_down), TTSRange=(TTS_up, TTS_down), \
                SS5range=SS5range, SS3range=SS3range, PPTrange=PPTrange)

    
    if promoter_pos[1] == 0:
        pass
    else:
        # output promoter, column 'number' means nothing
        if promoter_pos[0] < 0:
            promoter_pos[0] = 0
        record_model[1] = str(promoter_pos[0])
        record_model[2] = str(promoter_pos[1])
        record_model[6] = 'promoter-TSS'
        sys.stdout.write('\t'.join(record_model))

    # ouput exons, column 'number' represents exon number
    exonList = ExonPosList
    record_model[6] = 'exon'
    order = 0
    lenExonList = len(exonList)
    while order < lenExonList:
        Estart = str(exonList[order])
        Eend = str(exonList[order+1])
        exonno = str(exonNo[int(order/2)])

        # if any operation on exon position (e.g. no overlapping), 
        # using code below to prevernt from the condition that end pos is smaller than start one
        #if Estart<Eend:
        #    order += 2
        #    continue

        record_model[1] = Estart
        record_model[2] = Eend
        record_model[3] = exonno
        sys.stdout.write('\t'.join(record_model))
        order += 2
    # output introns, column 'number' represents intron number
    if intronNo:
        order = 0
        lenIntronList = len(intronList)
        record_model[6] = 'intron'
        while order < lenIntronList:
            Istart = str(intronList[order])
            Iend = str(intronList[order+1])
            intronno = str(intronNo[int(order/2)])
            
            #
            #if Istart<Iend:
            #    order += 2
            #    continue

            record_model[1] = Istart
            record_model[2] = Iend
            record_model[3] = intronno
            sys.stdout.write('\t'.join(record_model))
            order += 2

        # output 5' splice site, column 'number' represents intron number
        order = 0
        lenSS5List = len(SS5List)
        record_model[6] = '5ss'
        while order < lenSS5List:
            pos = SS5List[order]
            number = intronNo[order]
            ss5start = pos[0]
            ss5end = pos[1]

            record_model[1] = str(ss5start)
            record_model[2] = str(ss5end)
            record_model[3] = str(number)
            sys.stdout.write('\t'.join(record_model))
            order += 1

        # output 3' splice site, column 'number' represents intron number
        order = 0
        lenSS3List = len(SS3List)
        record_model[6] = '3ss'
        while order < lenSS3List:
            pos = SS3List[order]
            number = intronNo[order]
            ss3start = pos[0]
            ss3end = pos[1]

            record_model[1] = str(ss3start)
            record_model[2] = str(ss3end)
            record_model[3] = str(number)
            sys.stdout.write('\t'.join(record_model))
            order += 1

        # output ppt, column 'number' represents intron number
        order = 0
        lenPPTList = len(PPTList)
        record_model[6] = 'PPT'
        while order < lenPPTList:
            pos = PPTList[order]
            number = intronNo[order]
            PPTstart = pos[0]
            PPTend = pos[1]

            record_model[1] = str(PPTstart)
            record_model[2] = str(PPTend)
            record_model[3] = str(number)
            sys.stdout.write('\t'.join(record_model))
            order += 1

    # output TTS, column 'number' means nothing
    # drop out 0,0 record, which will cause bedtools intersect error.
    if TTS_pos[1] == 0:
        pass
    else:
        if TTS_pos[0] < 0:
            TTS_pos[0] = 0
        record_model[1] = str(TTS_pos[0])
        record_model[2] = str(TTS_pos[1])
        record_model[3] = '0'
        record_model[6] = 'TTS'
        sys.stdout.write('\t'.join(record_model))

def output_splitfiles(ExonPosList,record_model,promoter_up,promoter_down,TTS_up,TTS_down,SS5range=None,SS3range=None,PPTrange=None,strand='+'):
    promoter_pos,TTS_pos,intronNo,intronPosList,exonNo,SS5List,SS3List, PPTList = \
            classify_features(ExonPosList, strand, PromoterRange=(promoter_up, promoter_down), TTSRange=(TTS_up, TTS_down), \
                SS5range=SS5range, SS3range=SS3range, PPTrange=PPTrange)

    # output promoter, column 'number' means nothing
    # drop out 0,0 record, which will cause bedtools intersect error.
    if promoter_pos[1] == 0:
        pass
    else:
        if promoter_pos[0] < 0:
            promoter_pos[0] = 0
        record_model[1] = str(promoter_pos[0])
        record_model[2] = str(promoter_pos[1])
        record_model[6] = 'promoter-TSS'
        promoter_fo.write('\t'.join(record_model))

    # ouput exons, column 'number' represents exon number
    exonList = ExonPosList
    record_model[6] = 'exon'
    order = 0
    lenExonList = len(exonList)
    while order < lenExonList:
        Estart = exonList[order]
        Eend = exonList[order+1]
        exonno = exonNo[int(order/2)]

        #if Estart<Eend:
        #    order += 2
        #    continue

        record_model[1] = str(Estart)
        record_model[2] = str(Eend)
        record_model[3] = str(exonno)
        order += 2

        exon_fo.write('\t'.join(record_model))

    # output introns, column 'number' represents intron number
    if intronNo:
        order = 0
        lenIntronList = len(intronPosList)
        record_model[6] = 'intron'
        while order < lenIntronList:
            Istart = intronPosList[order]
            Iend = intronPosList[order+1]

            #if Istart<Iend:
            #    order += 2
            #    continue

            intronno = intronNo[int(order/2)]
            record_model[1] = str(Istart)
            record_model[2] = str(Iend)
            record_model[3] = str(intronno)
            order += 2

            intron_fo.write('\t'.join(record_model))

        # output 5' splice site, column 'number' represents intron number
        order = 0
        lenSS5List = len(SS5List)
        record_model[6] = '5ss'
        while order < lenSS5List:
            pos = SS5List[order]
            number = intronNo[order]
            ss5start = pos[0]
            ss5end = pos[1]

            record_model[1] = str(ss5start)
            record_model[2] = str(ss5end)
            record_model[3] = str(number)
            order += 1

            ss5_fo.write('\t'.join(record_model))

        # output 3' splice site, column 'number' represents intron number
        order = 0
        lenSS3List = len(SS3List)
        record_model[6] = '3ss'
        while order < lenSS3List:
            pos = SS3List[order]
            number = intronNo[order]
            ss3start = pos[0]
            ss3end = pos[1]

            record_model[1] = str(ss3start)
            record_model[2] = str(ss3end)
            record_model[3] = str(number)
            order += 1

            ss3_fo.write('\t'.join(record_model))

        # output ppt, column 'number' represents intron number
        order = 0
        lenPPTList = len(PPTList)
        record_model[6] = 'PPT'
        while order < lenPPTList:
            pos = PPTList[order]
            number = intronNo[order]
            PPTstart = pos[0]
            PPTend = pos[1]

            record_model[1] = str(PPTstart)
            record_model[2] = str(PPTend)
            record_model[3] = str(number)
            order += 1

            ppt_fo.write('\t'.join(record_model))

    # drop out 0,0 record, which will cause bedtools intersect error.
    if TTS_pos[1] == 0:
        pass
    else:
        # output TTS, column 'number' means nothing
        if TTS_pos[0] < 0:
            TTS_pos[0] = 0
        record_model[1] = str(TTS_pos[0])
        record_model[2] = str(TTS_pos[1])
        record_model[3] = '0'
        record_model[6] = 'TTS'
        TTS_fo.write('\t'.join(record_model))


#### parse parameters ####
ArgsDict = get_args()
gtf = ArgsDict['gtf']
features = ArgsDict['demand_features'].split(' ')

PromoterRange = ArgsDict['promoter_range']
TTSRange = ArgsDict['TTS_range']
SS5Range = ArgsDict['5SS_range']
SS3Range = ArgsDict['3SS_range']
PPTRange = ArgsDict['PPT_range']

SplitOutput = ArgsDict['split_output']
prefix = ArgsDict['split_prefix']
unique = ArgsDict['uniq']

verbose = ArgsDict['verbose']

promoter_up = int(PromoterRange.split(',')[0].strip(' '))
promoter_down = int(PromoterRange.split(',')[1].strip(' '))
TTS_up = int(TTSRange.split(',')[0].strip(' '))
TTS_down = int(TTSRange.split(',')[1].strip(' '))
SS5up = int(SS5Range.split(',')[0].strip(' '))
SS5down = int(SS5Range.split(',')[1].strip(' '))
SS3up = int(SS3Range.split(',')[0].strip(' '))
SS3down = int(SS3Range.split(',')[1].strip(' '))
PPTup = int(PPTRange.split(',')[0].strip(' '))
PPTdown = int(PPTRange.split(',')[1].strip(' '))

# split files to save result
if SplitOutput:
    exon_fo = open(prefix+'.exon.bed','w')
    intron_fo = open(prefix+'.intron.bed','w')
    ss5_fo = open(prefix+'.5ss.bed','w')
    ss3_fo = open(prefix+'.3ss.bed','w')
    gene_fo = open(prefix+'.gene.bed','w')
    ppt_fo = open(prefix+'.PPT.bed','w')
    promoter_fo = open(prefix+'.promoter.bed','w')
    TTS_fo = open(prefix+'.TTS.bed','w')
    if verbose:
        sys.stderr.write('''Starting to write results in following files:\n
{1}.gene.bed\n
{1}.exon.bed\n
{1}.intron.bed\n
{1}.5ss.bed\n
{1}.3ss.bed\n
{1}.PPT.bed\n
{1}.promoter.bed\n
{1}.TTS.bed\n''')

    FoDict = {
        'exon': exon_fo,
        'intron': intron_fo,
        'ss5':ss5_fo,
        'ss3':ss3_fo,
        'gene':gene_fo,
        'ppt':ppt_fo,
        'promoter':promoter_fo,
        'TTS':TTS_fo
    }
    FoList = [prefix+suffix+'.bed' for suffix in ['.exon','.intron','.5ss','.3ss','.gene','.PPT','.promoter','.TTS']]
else:pass

# init variables
ExonPosList = []
strand = '+'
chrom = ''
Gid = ''
Gbiotype = ''
Tid = ''

####################
# main function body
####################
with open(gtf, 'r') as gtffo:
    line = gtffo.readline()
    feature = ''
    last_line = False
    # remove file header
    while line.startswith('#'):
        line = gtffo.readline()

    # deal with first gene records
    while True:
        parser = GTFinfo(line)
        feature = parser.return_feature()
        if feature=='exon':
            break

        if feature == 'gene':
            Gid = parser.return_info('gene_id')
            Gbiotype = parser.return_info('gene_biotype')
            # output gene record
            Gstart, Gend = parser.exon_info()
            chrom, strand = parser.basic_info()
            gene_record = [chrom, Gstart, Gend, '0', '.', strand, 'gene', Gid, Gbiotype, '.\n']
            if SplitOutput:
                gene_fo.write('\t'.join(gene_record))
            else:
                sys.stdout.write('\t'.join(gene_record))
        elif feature == 'transcript':
            Tid = parser.return_info('transcript_id')
            chrom,strand = parser.basic_info()
        else:pass
        line = gtffo.readline()
    first_gene = 1

    while True:
        if line and not line.startswith('#'):
            # not reach bottom
            parser = GTFinfo(line)
            feature = parser.return_feature()

            if feature == 'exon':
                start,end = parser.exon_info()
                ExonPosList.append(start)
                ExonPosList.append(end)
                line = gtffo.readline();continue

            if feature == 'gene':
                Gid_n = parser.return_info('gene_id')
                Gbiotype_n = parser.return_info('gene_biotype')
                Gstart_n, Gend_n = parser.exon_info()
                chrom_n, strand_n = parser.basic_info()
                gene_record = [chrom_n, Gstart_n, Gend_n, '0', '.', strand_n, 'gene', Gid_n, Gbiotype_n, '.\n']
                if SplitOutput:
                    gene_fo.write('\t'.join(gene_record))
                else:
                    sys.stdout.write('\t'.join(gene_record))
                line = gtffo.readline();continue
            elif feature == 'transcript':
                Tid_n = parser.return_info('transcript_id')
                chrom_n,strand_n = parser.basic_info()
            elif feature != 'exon': ## skip CDS, stop codon and other features.
                line = gtffo.readline();continue

            if ExonPosList == []: # no exon in this gene
                line = gtffo.readline();continue
        # reach bottom
        else:
            last_line = True
            if ExonPosList == []:break

        # make sure pos in ascending order and modify coordinate different between GTF and BED
        ExonPosList = sorted([ eval(pos) - 1  if no%2==0 else eval(pos) for no,pos in zip(range(len(ExonPosList)),ExonPosList) ])

        # pattern: chrom, start, end, number,., strand, feature type, gene id, gene biotype, transcript id
        # intronNo_for5ss, intronNo_for3ss, intronNo_forPPT = intronNo, intronNo, intronNo
        record_model = [chrom,0,0,'0','.',strand,feature,Gid,Gbiotype,Tid+'\n']

        if SplitOutput:
            output_splitfiles(ExonPosList, record_model, promoter_up, promoter_down, TTS_up, TTS_down, \
                SS5range=(SS5up, SS5down), SS3range=(SS3up, SS3down), PPTrange=(PPTup, PPTdown),strand=strand)
        else:
            output_stdout(ExonPosList, record_model, promoter_up, promoter_down, TTS_up, TTS_down, \
                SS5range=(SS5up, SS5down), SS3range=(SS3up, SS3down), PPTrange=(PPTup, PPTdown), strand=strand)
        # clear exon position list
        ExonPosList = []
        # keep basic information
        if first_gene == 1:
            first_gene = 0
            pass
        else:
            strand, chrom, Gid, Gbiotype, Tid = strand_n, chrom_n, Gid_n, Gbiotype_n, Tid_n
        # exit when reaching bottom
        if last_line:break
        line = gtffo.readline()


    if SplitOutput:
        for fo in FoDict.values():
            fo.close()
        if unique:
            if verbose:
                sys.stderr.write('''Start unique files...\n''')
            for file in FoList:
                os.system('sort -k1,1 -k2,2n -k3,3n -k6,6 -k8,8 -u {0} > {0}.tmp;mv {0}.tmp {0}'.format(file))
            if verbose:
                sys.stderr.write('''All done''')