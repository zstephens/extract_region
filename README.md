# extract_region
a quick script for extracting portions of reads that span specified reference coordinates

## example
here's an example usage for extracting repeat expansions in C9orf72:

`python3 extract_region.py \ `  
`    -i input.bam \ `  
`    -o output.fa \ `  
`    -c chr9:27573485-27573546`  
