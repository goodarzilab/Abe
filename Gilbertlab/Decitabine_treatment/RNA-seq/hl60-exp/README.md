## Decitabine treatment RNA-Seq time-series experiments
In order to test for any differences over multiple time points, once can use a design including the time factor, and then test using the **likelihood ratio test (LRT)**. Here, as we have control (DMSO) and treatment (Decitabine) time series, design formula containing the condition factor, the time factor, and the interaction of the two. In this case, using the likelihood ratio test with a reduced model which does not contain the interaction terms will test whether the condition induces a change in gene expression at any time point after the reference level time point (time 0).
(see [DESeq2 Time-series-experiments](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments) for more details)
## PCA
<img src=plots/plots/PCA_filtered.png style="width:600px">
## Volcano plot
<img src=plots/plots/Volcano_plot.png style="width:600px">
## Heatmap clustering
<img src=plots/plots/Heatmap_clustering.pdf.png style="width:600px">
<h1>Enrichment analysis<h1>
<table>
  <tr>
  <h2>human_ensembl<h2>
    <td><img src=6h_delta_exp/human_ensembl.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_encode_tf<h2>
    <td><img src=6h_delta_exp/human_ensembl_encode_tf.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_encode_tf.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_encode_tf.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c1<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c1.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c1.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c1.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c2<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c2.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c2.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c2.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c3<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c3.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c3.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c3.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c4<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c4.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c4.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c4.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c5<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c5.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c5.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c5.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c6<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c6.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c6.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c6.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_c7<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_c7.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_c7.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_c7.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_full<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_full.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_full.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_full.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_msigdb_h<h2>
    <td><img src=6h_delta_exp/human_ensembl_msigdb_h.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_msigdb_h.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_msigdb_h.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_all_gene_ids<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_all_gene_ids.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_all_gene_ids.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_all_gene_ids.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_all_gene_names<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_all_gene_names.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_all_gene_names.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_all_gene_names.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_3UTR<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_3UTR.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_5UTR<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_5UTR.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_coding_exons<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_coding_exons.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_coding_gene_ids_by_introns<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_coding_gene_ids_by_introns.all.png style="width:600px">
  <tr>
<table>
<table>
  <tr>
  <h2>human_ensembl_RBPs_DeepBind<h2>
    <td><img src=6h_delta_exp/human_ensembl_RBPs_DeepBind.all.png style="width:600px">
    <td><img src=72h_delta_exp/human_ensembl_RBPs_DeepBind.all.png style="width:600px">
    <td><img src=120h_delta_exp/human_ensembl_RBPs_DeepBind.all.png style="width:600px">
  <tr>
<table>
