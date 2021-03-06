# NWPA
Implementation of Needleman-Wunsch pairwise alignment of DNA sequences in R language. It supports sequences with ambigious base counts and gapped sequences. Moreover, this library provides two function to read and write FASTA files.

Installation: 

<code>devtools::install_github("Takheer/NWPA")</code>

Usage of <code>align</code>: 

<code>align("example_1.fas", "alignedSequences.fas")</code>

Will align two first sequences from <code>example_1.fas</code> file and will write them aligned into <code>alignedSequences.fas</code>
Default settings are <code>1</code> for <code>match</code>, <code>0</code> for <code>mismatch</code> and <code>-1</code> for <code>gap</code> insertion/deletion

You can also set your own values:

<code>align("example_1.fas", "alignedSequences.fas", match=0, mismatch=-2, gap=-1)</code>
