"""EyeOnCNV - detection and review of small copy-number events.

Pipeline for targeted-panel sequencing:

0. ``scripts/bam_to_depth_jumps.sh``: BAM -> positions flanking abrupt depth changes;
1. :mod:`eyeoncnv.intervals`: positions -> events, recurrence across samples;
2. :mod:`eyeoncnv.genes`: gene/transcript annotation from the capture BED;
3. :mod:`eyeoncnv.hgvs_coords`: c. coordinates of the event boundaries;
4. :mod:`eyeoncnv.highlight`: flag the events worth reviewing.

"""

__version__ = "0.1.0"
