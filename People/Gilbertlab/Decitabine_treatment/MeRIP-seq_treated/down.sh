cd ~/People/Gilbertlab/Decitabine_treatment/MeRIP-seq_treated/radar
mkdir -p ipage
mkdir -p plots


for f in `ls -d /flash/bin/iPAGEv1.0/PAGE_DATA/ANNOTATIONS/human_*gs*`; do 
base=`basename "$f"`; 

perl $PAGEDIR/page.pl --expfile=d_mtyl_T_vs_U.c.txt \
--species=$base --exptype=continuous --ebins=11 --nodups=1; 
mv -v d_mtyl_T_vs_U.c.txt_PAGE/ ipage/d_mtyl_T_vs_U_${base}/
cp -v ipage/d_mtyl_T_vs_U_${base}/d_mtyl_T_vs_U.c.txt.summary.pdf plots/ipage.d_mtyl_T_vs_U_${base}.pdf
done

# magick -density 300 $i -flatten -quality 90 $o
