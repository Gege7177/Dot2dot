#include "utils.h"

void print_version (void) {
  printf ("Dot2dot version 1.0 patch 3 build 22-06-18\n");
}

void print_usage (void) {
  printf ("dot [options] in () the same option in the config file\n");
  printf("   -s, --sequence (Sequence): \n");
  printf("   sequence file name (accept fasta/multifasta/fastq).\n");
  printf("   The file format is detected inspecting the sequences’ header.\n");
  printf("   -c, --config:   [Command line only]\n");
  printf("   configuration file name (required)\n");
  printf("   -o, --output (Outfile):\n");
  printf("   output file name. (bed/dot)\n");
  printf("   -l, --minmotif (MinMotifLen):   [default=2]\n");
  printf("   minimum size in bp of the motif sequence \n");
  printf("   -L, --maxmotif (MaxMotifLen): [default=30]\n"); 
  printf("   maximum size in bp of the motif sequence\n");
  printf("   -m, --minmatch (MinMatch): [range (0,1) default=1]\n");
  printf("    Minimum overall matching score normalized in the range [0,1]\n");
  printf("   -G, --maxgaps (MaxGaps): [default=0]\n");
  printf("   maximum number of mismatches in a motif (expressed in bp).\n");
  printf("   -I, --maxinsert (MaxInsert): [default=0]\n");
  printf("   maximum insert size expressed in bp.\n");
  printf("   -t, --threads (Threads): [default=1]\n");
  printf("   number of threads. Should be lower or equal to the number of sequences.\n");
  printf("   -v, --verbose: \n");
  printf("   print on the standard error the id of a sequence once loaded\n");
  printf("   -V, --version: \n");
  printf("   print dot-to-dot version and exit\n");
  printf("   -h, --help: \n");
  printf("   print a help page and exit\n");
  printf("\nConfig file options\n");
  printf("   FilterType: [default=NONE]\n");
  printf("   sets the level of filtering\n");
  printf("   * NONE: disable filtering\n");
  printf("   * THRESHOLD: filters-out only: too small TRs\n");
  printf("   * LIGHT: apply light filtering (see manual for details\n");
  printf("   * FAIR: apply fair filtering (see manual for details)\n");
  printf("   * HEAVY: apply heavy filtering (see manual for details)\n");
  printf("   MinTRLen: [default=12]\n");
  printf("   minimum length of a tandem repeat to be included in the output.\n");
  printf("   MinPurity: [range (0,1) default=0.1]\n");
  printf("   minimum purity of a tandem repeat to be included in the output.\n");
  printf("   Tolerance: [range (0,1) default=0]\n");
  printf("   Set a degree of tolerance to be accepted using heavy filtering\n");
  printf("   AllowOverlap: [range (Y/N) default=Y]\n");
  printf("   Enable/Disable  overlapping results in the final output.\n\n");
}
