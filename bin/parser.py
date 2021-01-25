
import os
import sys



class SamtoolsDepth:
    def __init__(self, args):
        self.infn = args.infn
        self.min_depth = args.depth
        self.outfn = args.outfn

    def to_bed(self):
        outfh = open(self.outfn, 'w')
        infh = open(self.infn)
        line = infh.readline()
        items = line.rstrip('\n').split('\t')
        pre_chrom = items[0]
        pre_pos = items[1]

        for line in infh:
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            pos = items[1]
            depth = items[2]

            if int(depth) < int(self.min_depth):
                continue

            if chrom in [pre_chrom] and int(pos)-int(pre_pos) in [1]:
                pre_chrom = chrom
                pre_pos = pos
            elif chrome in [pre_chrom] and int(pos)-int(pre_pos) not in [1]:
                end = pos
                out_items = [chrome, start, end]
        #

class ExtractCoIntersect:
    def __init__(self, args):
        self.infns = args.inputs
        self.outfn  = args.output
        self.tempDir = args.temp
        if not os.path.isdir(self.tempDir):
            os.mkdir(self.tempDir)
        #
        for idx in range(1, len(self.infns)):
            iter_id = 'iter{0}'.format(idx)
            if idx in [1]:
                temp_fn = os.path.join(self.tempDir, 'ExtractCoIntersect.{0}.bed'.format(iter_id))
                _cmd_s = ['bedtools']
                _cmd_s.append('intersect')
                _cmd_s.append('-a')
                _cmd_s.append(self.infns[idx-1])
                _cmd_s.append('-b')
                _cmd_s.append(self.infns[idx])
                _cmd_s.append('>')
                _cmd_s.append(temp_fn)
                print(_cmd_s)
                os.system(' '.join(_cmd_s))
            else:
                pre_id = 'iter{0}'.format(idx-1)
                pre_fn = os.path.join(self.tempDir, 'ExtractCoIntersect.{0}.bed'.format(pre_id))
                temp_fn = os.path.join(self.tempDir, 'ExtractCoIntersect.{0}.bed'.format(iter_id))
                _cmd_s = ['bedtools']
                _cmd_s.append('intersect')
                _cmd_s.append('-a')
                _cmd_s.append(pre_fn)
                _cmd_s.append('-b')
                _cmd_s.append(self.infns[idx])
                _cmd_s.append('>')
                _cmd_s.append(temp_fn)
                print(_cmd_s)
                os.system(' '.join(_cmd_s))
        #
        _cmd_s = ['cp']
        _cmd_s.append(temp_fn)
        _cmd_s.append(self.outfn)
        print(_cmd_s)
        os.system(' '.join(_cmd_s))


class ParseCNVkit:
    def __init__(self, args):
        self.samples = args.samples
        #
        self.fdic = self.path_finder()
        #
        self.ref_dic = self.parse_refbed(args.refbed)
        self.seg_call_dic = self.parse_callcns()
        self.sed_anno_dic = self.parse_annobed()
        #
        self.write_ref_bed(args.outdir)
        self.write_comp_bed(args.outdir)

    def write_comp_bed(self, outdir):
        outfn = os.path.join(outdir, 'comp.ipsc57LBvsRef.bed')
        outfh = open(outfn, 'w')
        outfh.write('track name=compCNV description="CNV region of other samples in ipsc57_LB_vs_Ref"\n')
        outfh.write('#chrom chromStart chromEnd name\n')
        for region, info_dic in sorted(self.sed_anno_dic.items()):
            for sample, call_region_s in sorted(info_dic.items()):
                for call_region in call_region_s:
                    chrom, start, end = self.region_unpack(call_region)
                    #value = self.seg_call_dic[sample][call_region]['cn']
                    value = self.seg_call_dic[sample][call_region]['log2']
                    if abs(float(value)) >= 2.0:
                        items = [chrom, start, end, f'{sample}_{value}']
                        outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()

    def write_ref_bed(self, outdir):
        outfn = os.path.join(outdir, 'ref.ipsc57LBvsRef.bed')
        outfh = open(outfn, 'w')
        outfh.write('track name=ipsc57LBvsRef description="CNV region of ipsc57_LB_vs_Ref"\n')
        outfh.write('#chrom chromStart chromEnd name\n')
        for region, state in sorted(self.ref_dic.items()):
            chrom, start, end = self.region_unpack(region)
            items = [chrom, start, end, state.replace(' ', '_')]
            outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()

    def region_unpack(self, region):
        items = [region.split(':')[0]]
        items.extend(region.split(':')[1].split('_'))
        return items




    def parse_annobed(self):
        _dic = dict()
        for region in self.ref_dic:
            _dic.setdefault(region, {})
        for sample, info_dic in self.fdic.items():
            annobed = info_dic['annobed']
            for line in open(annobed):
                items = line.rstrip('\n').split('\t')
                call_chrom = items[0]
                call_start = items[1]
                call_end = items[2]
                call_region = f'{call_chrom}:{call_start}_{call_end}'
                call_genes = items[3].split(',')
                call_cn = items[4]
                anno_chrom = items[5]
                anno_start = items[6]
                anno_end = items[7]
                anno_region = f'{anno_chrom}:{anno_start}_{anno_end}'
                anno_state = items[8]
                _dic[anno_region].setdefault(sample, [])
                _dic[anno_region][sample].append(call_region)
        return _dic

    def parse_callcns(self):
        _dic = dict()
        for sample, info_dic in self.fdic.items():
            callcns = info_dic['callcns']
            for line in open(callcns):
                items = line.rstrip('\n').split('\t')
                if items[0] in ['chromosome']:
                    idx_dic = dict()
                    for idx, items in enumerate(items):
                        idx_dic.setdefault(items, idx)
                    continue
                chrom = items[idx_dic['chromosome']]
                start = items[idx_dic['start']]
                end = items[idx_dic['end']]
                genes = items[idx_dic['gene']].split(',')
                log2 = items[idx_dic['log2']]
                cn = items[idx_dic['cn']]
                depth = items[idx_dic['depth']]
                probes = items[idx_dic['probes']]
                weight = items[idx_dic['weight']]
                #
                seg_region = f'{chrom}:{start}_{end}'
                _dic.setdefault(sample, {}).setdefault(seg_region, {})
                _dic[sample][seg_region].setdefault('genes', genes)
                _dic[sample][seg_region].setdefault('log2', log2)
                _dic[sample][seg_region].setdefault('cn', cn)
                _dic[sample][seg_region].setdefault('depth', depth)
                _dic[sample][seg_region].setdefault('probes', probes)
                _dic[sample][seg_region].setdefault('weight', weight)
        return _dic

    def parse_refbed(self, refbed):
        state_dic = {'State_2': 'Deleted part in LB compared to ref',
                     'State_4': 'Amplified part in LB compared to ref',
                     'State_5': 'Amplified part in LB compared to ref'}
        ref_dic = dict()
        for line in open(refbed):
            items = line.rstrip('\n').split('\t')
            chrom = items[0]
            start = items[1]
            end = items[2]
            region = f'{chrom}:{start}_{end}'
            state = items[3]
            state_desc = state_dic[state]
            #
            ref_dic.setdefault(region, f'{state} {state_desc}')
        return ref_dic

    def path_finder(self):
        fdic = dict()
        for sample in self.samples:
            callcns = f'analysis/cnvkit/{sample}/{sample}.call.cns'
            if os.path.isfile(callcns):
                fdic.setdefault(sample, {}).setdefault('callcns', callcns)
            else:
                print(f'Error:Could not found {callcns}')
                sys.exit()
            annobed = f'analysis/cnvkit/{sample}/{sample}.call.anno.bed'
            if os.path.isfile(annobed):
                fdic.setdefault(sample, {}).setdefault('annobed', annobed)
            else:
                print(f'Error:Could not found {callcns}')
                sys.exit()
        return fdic




def main(args):

    if args.subparser_name in ['SamtoolsDepth']:
        samtools_depth = SamtoolsDepth(args)
        samtools_depth.to_bed()

    if args.subparser_name in ['ExtractCoIntersect']:
        co_intersect = ExtractCoIntersect(args)

    if args.subparser_name in ['ParseCNVkit']:
        cnvkit_obj = ParseCNVkit(args)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name')

    subparser = subparsers.add_parser('SamtoolsDepth')
    subparser.add_argument('--infn', help='outfn of samtools depth')
    subparser.add_argument('--depth', help='min depth for filter out')
    subparser.add_argument('--outfn', help='outfn bed fmt')

    subparser = subparsers.add_parser('ExtractCoIntersect')
    subparser.add_argument('--inputs', nargs='+', help='input file names bed fmt')
    subparser.add_argument('--output', help='out file name')
    subparser.add_argument('--temp', help='temp dir')

    subparser = subparsers.add_parser('ParseCNVkit')
    subparser.add_argument('--samples', nargs='+', help='sample names')
    subparser.add_argument('--refbed', help='reference bed')
    subparser.add_argument('--outdir', help='path for output')

    args = parser.parse_args()
    main(args)

