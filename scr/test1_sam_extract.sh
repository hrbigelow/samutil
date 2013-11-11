#!/bin/bash

# script to test samutil extract

# Synopsis:
# We want to confirm a few things:

# Note that an 'orphan' read is defined as a record in a BAM file that
# is marked as being 'paired in sequencing' but for which no other
# record in the BAM file represents its mate pair.  Note that this
# does NOT mean a case when only one read in the pair is mapped.  In
# that case, the unmapped mate pair will indeed appear as a BAM
# record, but just be unmapped.  And, so those records will be
# extracted as r1 and r2, just the same as if they were mapped.

# First, confirm that we get the same set of paired reads after
# extraction as the set that were initially aligned. This is achieved
# in steps 1-6. For technical reasons, it is easier to compare two
# fastq.gz files that are sorted by read id. The easiest way to sort a
# fastq.gz file by read id is to first align it, then extract the
# reads. (samutil extract outputs the extracted fastq.gz reads in read
# id order. So, steps 1 and 2 first generate the read-id-sorted fastq
# files.  Step 3 confirms the absence of orphans in the generated bam
# file. This is a basic integrity check that the bam file doesn't
# contain any orphans. (It is not invalid BAM syntax to contain
# orphans, but it makes no sense to have such a situation).  Steps 4-6
# merely confirm that you get out exactly what you put in, when
# starting with a set of reads, aligning them, then re-extracting
# them.

# Second, in the edge case that a bam file does contain orphans, are
# they treated correctly by samutil extract? To test this, start with
# a non-orphan-containing BAM file, and 'deal' the records 50 at a
# time into two separate BAM files.  These two separate BAM files will
# now contain the same number of distinct read ids that are orphans,
# and because of this, we should expect that samutil extract will find
# the same read-id-set of orphan reads in each of the two files
# (though, whether each one is read1 or read2 is random). To confirm
# the read-id-set equality, take both extracted orphan files and align
# them as quasi-paired end. This will cause bwa mem to use all four
# orientations (FF, FR, RF, RR) in alignment, but it should do so.

FQ1=$1
FQ2=$2
PREFIX=$3

# REF=/apps/gau/genomes/m_musculus/Mouse.B38/Mouse.B38.fasta
REF=/gau-intermediate2/hbigelow/sugar/CEUref.hg19.fasta
THREADS=12

#cat<<EOF
echo 'Step 1: Generate initial alignment'
bwa mem -t $THREADS $REF $FQ1 $FQ2 > $PREFIX.1.sam 2> $PREFIX.1.log

echo 'Step 2: Extracting reads to get ID-sorted fastq'
samutil sort -t $THREADS FRAGMENT $PREFIX.1.sam /dev/stdout 2>/dev/null | \
  samutil extract -t $THREADS /dev/stdin $PREFIX.2.{r1,r2,orphans}.fq.gz 2>/dev/null

echo 'Step 3: Confirm absence of orphans.'
echo 'Size of orphans file: ' $(stat -c '%s' $PREFIX.2.orphans.fq.gz)

echo 'Step 4: Realign extracted reads to same genome'
bwa mem -t $THREADS $REF $PREFIX.2.{r1,r2}.fq.gz > $PREFIX.4.sam 2> $PREFIX.4.log

echo 'Step 5: Extract as-is'
samutil sort -t $THREADS FRAGMENT $PREFIX.4.sam /dev/stdout 2>/dev/null | \
  samutil extract -t $THREADS /dev/stdin $PREFIX.5.{r1,r2,orphans}.fq.gz 2>/dev/null

echo 'Step 6: Confirm the output is identical to the output of 2'
echo 'Running diff -q'
zdiff -q $PREFIX.{2,5}.r1.fq.gz
zdiff -q $PREFIX.{2,5}.r2.fq.gz
echo 'End of running diff -q'

echo 'Step 7: Split the SAM file produced in step 4, into two, in order to simulate orphan-containing BAMs'
(samtools view -SH $PREFIX.4.sam 2>/dev/null; \
  samtools view -S $PREFIX.4.sam 2>/dev/null | awk '{ if (NR % 100 > 50) { print $0 } }') \
  > $PREFIX.7.chunk1.sam

(samtools view -SH $PREFIX.4.sam 2>/dev/null; \
  samtools view -S $PREFIX.4.sam 2>/dev/null | awk '{ if (NR % 100 <= 50) { print $0 } }') \
  > $PREFIX.7.chunk2.sam

echo 'Step 8: Extract reads from both chunks from step 7'
samutil sort -t $THREADS FRAGMENT $PREFIX.7.chunk1.sam /dev/stdout 2>/dev/null | \
  samutil extract -t $THREADS /dev/stdin $PREFIX.8.chunk1.{r1,r2,orphans}.fq.gz 2>/dev/null

samutil sort -t $THREADS FRAGMENT $PREFIX.7.chunk2.sam /dev/stdout 2>/dev/null | \
  samutil extract -t $THREADS /dev/stdin $PREFIX.8.chunk2.{r1,r2,orphans}.fq.gz 2>/dev/null

echo 'Step 9: Get union of orphans'
~/cc/samutil/merge_orphans.plx \
<(zcat $PREFIX.8.chunk1.orphans.fq.gz) \
<(zcat $PREFIX.8.chunk2.orphans.fq.gz) \
>(gzip -c - > $PREFIX.9.chunk1.mo.fq.gz) \
>(gzip -c - > $PREFIX.9.chunk2.mo.fq.gz)


echo 'Step 10: Realign reads'
bwa mem -t $THREADS $REF $PREFIX.8.chunk1.{r1,r2}.fq.gz > $PREFIX.10.chunk1.sam 2> $PREFIX.10.chunk1.log
bwa mem -t $THREADS $REF $PREFIX.8.chunk2.{r1,r2}.fq.gz > $PREFIX.10.chunk2.sam 2> $PREFIX.10.chunk2.log
bwa mem -t $THREADS $REF $PREFIX.9.chunk{1,2}.mo.fq.gz > $PREFIX.10.orphans.sam 2> $PREFIX.10.orphans.log

echo 'Step 11: Extract aligned reads again'
(samtools view -Sh $PREFIX.10.chunk1.sam 2>/dev/null; \
  samtools view -S $PREFIX.10.chunk2.sam 2>/dev/null; \
  samtools view -S $PREFIX.10.orphans.sam 2>/dev/null) | \
samutil sort -t $THREADS FRAGMENT /dev/stdin /dev/stdout 2>/dev/null | \
  samutil extract -t $THREADS /dev/stdin $PREFIX.11.{r1,r2,orphans}.fq.gz 2>/dev/null

echo 'Step 12: Confirm at least that the read sets are identical with step 2'
echo 'Running diff -q'
diff -q \
  <(zcat $PREFIX.11.r1.fq.gz | awk '{ if (NR % 4 == 1) { print $0 } }') \
  <(zcat $PREFIX.2.r1.fq.gz | awk '{ if (NR % 4 == 1) { print $0 } }')

diff -q \
  <(zcat $PREFIX.11.r2.fq.gz | awk '{ if (NR % 4 == 1) { print $0 } }') \
  <(zcat $PREFIX.2.r2.fq.gz | awk '{ if (NR % 4 == 1) { print $0 } }')
echo 'End of diff -q'

echo 'Size of orphans files (should be zero) : ' $(stat -c '%s' $PREFIX.11.orphans.fq.gz)

echo 'Done.'

#EOF
