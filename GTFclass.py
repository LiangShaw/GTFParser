import os,sys

class GTFinfo:
    def __init__(self,gtfline):
        self.gtf_list = gtfline.split('\t')

        self.feature = self.gtf_list[2]

        extra_info = self.gtf_list[8]

        info_ls = extra_info.strip().split('"; ')
        #sys.stderr.write(str(info_ls)+'\n')
        self.info_dict = { f.split('"')[0].strip(' '):f.split('"')[1] for f in info_ls if '"' in f }
        # gene_id, transcript_id, gene_biotype, exon_number

    def return_feature(self):
        return self.feature

    def return_info(self,InfoName):
        return self.info_dict[InfoName]

    def basic_info(self):
        chrom = self.gtf_list[0]
        strand = self.gtf_list[6]
        return chrom, strand

    def exon_info(self):
        start = self.gtf_list[3]
        end = self.gtf_list[4]
        return start, end

if __name__ == '__main__':
    for gtf in sys.stdin:
        if gtf.startswith('#'):
            continue
        try:
            parser = GTFParser(gtf)
            print(parser.return_feature('gene_id'))
        except:
            print(gtf)
            os._exit(0)
