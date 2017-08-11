use strict;
use FileHandle;
use Env;
use PerlIO::gzip;
use Getopt::Long;
my $inputfile="EMC_Cartagene.GLM.CCPHighPosvsNeg.withSmoking.chrAll.qvalue.BH.bed";
my $outputfile="EMC_Cartagene.GLM.CCPHighPosvsNeg.withSmoking.chrAll.q0.01.dist100.bed";
my $help="";
my $qvalue=0.01;
my $dist=200;
my $qvcolind=11; ## qvalue column index
my $betaind=4; ## beta effect size column index

GetOptions ('help' =>\$help,
        'i|inputfile=s' => \$inputfile,
        'o|outputfile=s'  => \$outputfile,
        'd|dist=i' => \$dist,
        'q|qvalue=f' => \$qvalue,
        );
if ($help) {
 print_helpfile();
 exit();
}       

open INFILE, "<$inputfile" || die "Cannot open $inputfile for reading.\n";
open OUTFILE, ">$outputfile" || die "Cannot open $outputfile for writing.\n";

my $temp="";
my @temparray=();
my $chr="chr1";
my $start="";
my $end="";
my $qv="";
my $DMRqvmin = 0;
my $DMRqvmean = 0;
my $beta="";
my $DMCnum=0;
print OUTFILE "#chr\tstart\tend\tDMCnum\tDMRqvmin\tDMRqvmean\tmeanbeta\n";
$temp=<INFILE>; # filter the first header row
while ($temp=<INFILE>) {
           chomp($temp);
           @temparray = split("\t",$temp);
           $chr=$temparray[0];
           $qv=$temparray[$qvcolind];
           if($qv < $qvalue) {
             if($DMCnum==0){
                $start=$temparray[1];
                $end = $temparray[2];
                $beta=$temparray[$betaind];
                $DMCnum=1;
                $DMRqvmin = $qv;
		$DMRqvmean = $qv;
             }elsif( (($temparray[1] - $end) <= $dist) && $beta * $temparray[$betaind] > 0 ){
                $end = $temparray[2];
                $beta += $temparray[$betaind];
                $DMRqvmean += $qv;
		$DMCnum++;
                if($qv < $DMRqvmin) {$DMRqvmin = $qv;} 
             }else{
                #close a DMR if either distance >= cutoff or different beta effect direction.
                #$beta =$beta/$DMCnum;
                $beta = sprintf "%4.3g", $beta/$DMCnum;
                $DMRqvmean = sprintf "%4.3e", $DMRqvmean / $DMCnum;
		print OUTFILE "$chr\t$start\t$end\t$DMCnum\t$DMRqvmin\t$DMRqvmean\t$beta\n";
                ## start a new DMR region
                $DMCnum=0;
             }
	  }elsif($DMCnum>0) {
   		$beta = sprintf "%4.3g", $beta/$DMCnum;
   		$DMRqvmean = sprintf "%4.3e", $DMRqvmean / $DMCnum;
		print OUTFILE "$chr\t$start\t$end\t$DMCnum\t$DMRqvmin\t$DMRqvmean\t$beta\n";
   		$DMCnum=0;	  
	   }
}           
if($DMCnum>0) {
   $beta = sprintf "%4.3g", $beta/$DMCnum;
   $DMRqvmean = sprintf "%4.3e", $DMRqvmean / $DMCnum;	
   print OUTFILE "$chr\t$start\t$end\t$DMCnum\t$DMRqvmin\t$DMRqvmean\t$beta\n";
   $DMCnum=0;
} 
close(OUTFILE);

