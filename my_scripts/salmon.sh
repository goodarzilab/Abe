JOBDIR=$1
FASTQDIR=$2
INDEX=$3

cd $JOBDIR
mkdir -p ./quants/

pwd 

for f in $FASTQDIR/*fastq.gz; do 
	samp=`basename ${f}`; 
	samp=${samp/.fastq.gz/}; 
	echo "Processing sample ${samp}"; 
    echo "Write to " ./quants/$samp;
	salmon quant -i $INDEX\
	-l A -r $f -p 18 --validateMappings -o ./quants/$samp; 
done
