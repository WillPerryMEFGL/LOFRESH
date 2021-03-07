#!/bin/perl


unless ($#ARGV == 0)

{

   print "Usage: 2_cutadapt_WP2_CO1.pl FastqList_WP2\n";

die;
}


open (INLIST, "<$ARGV[0]") || die;


$indir = "/nfshome/store02/users/b.bsp81d/Lofresh_raw_data/WP2/";
$outdir = "/nfshome/store02/users/b.bsp81d/Lofresh_raw_data/WP2_gene_demultiplex/WP2_CO1/";


while (<INLIST>) {
$lib = $_;
chomp($lib);

$read1 = $lib . "_1.fq.gz";
$read2 = $lib . "_2.fq.gz";

$cutout1 = $lib . "_cutoutCO1.1.fastq.gz";
$cutout2 = $lib . "_cutoutCO1.2.fastq.gz";


# use Trimmomatic to clip adaptor sequences & trim low-quality bases
system("cutadapt -g GGNACNGGNTGAACNGTNTANCCNCC -a TANACNTCNGGNTGNCCNAANAANCA --discard-untrimmed -o $outdir/$cutout1 -p $outdir/$cutout2 $indir/$read1 $indir/$read2");


}