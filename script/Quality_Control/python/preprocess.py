from preprocessing import extract_umi_barcode, cut_adapt, MatchConfig, BarcodeConfig
import argparse
import os

class PreConfig:
    def __init__(self, locus, core, base_quality, sample, compression_level, cut, batch_size, gzip_output):
        self.locus = locus
        self.core = core
        self.base_quality = base_quality
        self.sample = sample
        self.compression_level = compression_level
        self.cut = cut
        self.batch_size = batch_size
        self.gzip_output = gzip_output

def str_to_bool(value):
    """Convert str to bool"""
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f'Boolean value expected, got: {value}')

def preprocess(preconfig, match_config, barcode_config, reads1, reads2, output, correct_barcode):

    output_name = 'fastq'

    # cutadapt
    if preconfig.cut:
        output_name = output_name + '_cut'
        output_path = os.path.join(output, output_name)
        try:
            os.mkdir(output_path)
        except:
            pass
        reads1, reads2 = cut_adapt(preconfig.locus, preconfig.core, reads1, reads2, output_path, preconfig.base_quality, preconfig.sample)
        print("\n")
    
    # extract UMI and barcode
    output_name = output_name + '_umi_barcode'
    output_path = os.path.join(output, output_name)
    try:
        os.mkdir(output_path)
    except:
        pass
    extract_umi_barcode(match_config, barcode_config, reads1, reads2, output_path, preconfig.sample, preconfig.compression_level, correct_barcode, preconfig.core, preconfig.batch_size, preconfig.gzip_output)
    print(output_path)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Preprocess fastq files')

    parser.add_argument('-r1', '--reads1', type=str, help='path to reads1 file')
    parser.add_argument('-r2', '--reads2', type=str, help='path to reads2 file')
    parser.add_argument('-o', '--output', type=str, help='path to output directory')

    parser.add_argument('-l', '--locus', type=str, default='CA', help='locus name')
    parser.add_argument('-c', '--core', type=int, default=8, help='number of cores to use')
    parser.add_argument('-q', '--base_quality', type=int, default=10, help='base quality score')
    parser.add_argument('-s', '--sample', type=str, help='sample name')
    parser.add_argument('-cl', '--compression_level', type=int, default=6, help='compression level for gzip output')
    parser.add_argument('-cut', '--cut', type=str_to_bool, default=False, help='perform cutadapt')
    parser.add_argument('-bs', '--batch_size', type=int, default=50000, help='number of read pairs per barcode extraction worker batch')
    parser.add_argument('-go', '--gzip_output', type=str_to_bool, default=True, help='whether to gzip barcode FASTQ output; false uses more disk space but can be faster')

    parser.add_argument('-l1', '--linker1', type = str, default = "GTGGCCGATGTTTCGCATCGGCGTACGACT", help = 'Linker 1 sequence')
    parser.add_argument('-l2', '--linker2', type = str, default = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG", help = 'Linker 2 sequence')
    parser.add_argument('-m', '--mm_rate', type = float, default = 0.05, help = 'Mismatch rate for linker sequences')

    parser.add_argument('-b1', '--barcodeA_whitelist', type = str, help = 'Path to barcode A whitelist file')
    parser.add_argument('-b2', '--barcodeB_whitelist', type = str, help = 'Path to barcode B whitelist file')
    parser.add_argument('-bmd', '--bc_max_dist', type = int, default = 1, help = 'Maximum distance for barcode correction')
    parser.add_argument('-cb', '--correct_barcode', type=str_to_bool, default=False, help='Whether to perform barcode correction')

    args = parser.parse_args()

    reads1 = args.reads1
    reads2 = args.reads2
    output = args.output

    locus = args.locus
    core = args.core
    base_quality = args.base_quality
    sample = args.sample
    compression_level = args.compression_level
    cut = args.cut
    batch_size = args.batch_size
    gzip_output = args.gzip_output

    linker1 = args.linker1
    linker2 = args.linker2
    mm_rate = args.mm_rate

    barcodeA_whitelist = args.barcodeA_whitelist
    barcodeB_whitelist = args.barcodeB_whitelist
    bc_max_dist = args.bc_max_dist
    correct_barcode = args.correct_barcode
    preconfig = PreConfig(locus, core, base_quality, sample, compression_level, cut, batch_size, gzip_output)

    match_config = MatchConfig(linker1, linker2, mm_rate)
    barcode_config = BarcodeConfig(barcodeA_whitelist, barcodeB_whitelist, bc_max_dist)

    preprocess(preconfig, match_config, barcode_config, reads1, reads2, output, correct_barcode)
