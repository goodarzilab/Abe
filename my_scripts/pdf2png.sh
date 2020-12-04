pdf=$1
png=${pdf/.pdf/.png}

echo $pdf '>' $png
magick -density 300 $pdf $png
echo 'done!'
