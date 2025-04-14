# Test data for fusion detection

The dataset consists of hg38 chromosome 20 along with gencode 22 annotation.
It was generated accroding to these instructions:
 https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-generate-the-reference-files-for-a-non-supported-genome

Chr20 was chosen as it's small and contains at least one known fusion in the MCF7 cell line where both partners orginate from the same chromosome. 

The *.bt2 bowtie files are not needed for the JAFFAL pipeline but are checked for their presence anyhow. So these are just empty files.
