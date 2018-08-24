# Dot2dot: Accurate Whole-Genome Tandem Repeats Discovery

<b>Motivation:</b> Large-scale sequencing projects have confirmed the hypothesis that eukaryotic DNA is
rich in repetitions whose functional role needs to be elucidated. In particular, tandem repeats (TRs) (i.e.
short, almost identical sequences that lie adjacent to each other) have been associated to many cellular
processes and, indeed, are also involved in several genetic disorders. The need of comprehensive lists
of TRs for association studies and the absence of a computational model able to capture their variability,
have revived research on discovery algorithms.

<b>Results:</b> Building upon the idea that sequence similarities can be easily displayed using graphical
methods, we formalized the structure that TRs induce in dot plot matrices where a sequence is compared
with itself. Leveraging on the observation that a compact representation of these matrices can be built
and searched in linear time, we developed Dot2dot: an accurate algorithm fast enough to be suitable
for whole-genome discovery of tandem repeats. Experiments on five manually-curated collections of
TRs have shown that Dot2dot is more accurate than other state-of-the-art methods, and completes the
analysis of the biggest known reference genome in about one day on a standard workstation.
