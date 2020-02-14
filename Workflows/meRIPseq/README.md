To use https://github.com/igvteam/igv-jupyter

1. name bam files like this:
`<sample>.<treatment>.<input/IP>.bam`

2. Metadata: header-free matrix of two columns: 1) `<sample>` 2) `<treatment>`, respectively

3. List of comparisons between conditions: BvsA,CvsB... where A, B... are `<treatment>`s, respectively

4. Run jupyter on your server, tunnel that to your local machine and enjoy.
