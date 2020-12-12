for f in `ls -d /flash/bin/iPAGEv1.0/PAGE_DATA/ANNOTATIONS/human_ensembl*`; do
    base=`basename "$f"`
    echo '<table>' >> README.md
    echo '  <tr>' >> README.md
    echo '  <h2>'$base'<h2>' >> README.md
    echo '  <img src=lnTE_T_vs_U/'${base}'.all.png style="width:600px">' >> README.md
    echo '  <tr>' >> README.md
    echo '<table>' >> README.md
done
