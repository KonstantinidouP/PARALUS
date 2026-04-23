#!/bin/bash
# masterchef makes a bed file to use for annotation with paralus

unset ALL_PROXY
prog=$0
scriptdir=$(dirname "${prog}")

if [ $# -eq 0 ]; then
    echo 'PARalus v1.4.3
Please give me something to do 
No input files provided
Please type "--help" for instructions'
    exit 0
fi

if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    cat << 'EOF'
#####################################################
Usage:    
  /path/to/masterchef.sh -g genes.gtf -t tRNA_DB.bed
  
  Required fields:
    
    -g		genecode.gtf file  
    -t		tRNA file from tRNA_DB (optional)
    
####################################################
EOF
    exit 0
fi

while getopts hg:t:L: OPT; do
    case "$OPT" in
        g)
        	gtf=$OPTARG
            ;;
        t)
		tRNA_DB=$OPTARG
            ;;

        h) usage
        exit 0;; 
    esac
    
done

gtfname=${gtf%.gtf}

file_dir=$(dirname "${gtf}")

if [ ! -f "${gtf}" ] ; then
echo "error:
please provide a genes.gtf file using the -g option"
 exit 0
else
echo "GTF file provided"
fi

echo "extracting all the features in a master bed file written in the same directory of the GTF file wile creat a Paralus.bed file to use with Paralus
" 
		cat $gtf | awk -F "\t" '$3 == "exon" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "CDS"' | grep -v "##" | sed 's/"//g' | sed 's/;//g' | sed 's/\s/\t/g' | sed 's/\t\t/\t/g' | sed 's/"="/\t/g' | awk '{
		OFS="\t"; for (I=1; I<NF;I++) { if ($I== "gene_name") n=$(I+1); if ($I == "gene_type") t=$(I+1); } print $1,$4-1,$5,n,$3,$7,t}' \
		|  sed 's/,/-/g' | bedtools sort | uniq | awk -F "\t" '$5 !~ "UTR"' | awk -F "\t" '$7 !~ "TEC"' | awk -F "\t" '$4 !~ "Y_RNA"' \
		| awk -F "\t" '$4 !~ "pRNA"' | grep -v transcript | grep -v "antisense\|Metazoa_SRP\|misc_RNA\|sense_intronic\|-AS1\|mascRNA-menRNA\|pseudogene" > "${gtfname}"_gtf.bed

		echo "extracting gene names..."
		awk -F "\t" '{print $4}' "${gtfname}"_gtf.bed | sort | uniq > "${gtfname}"_genes

		echo "grouping annotation to create full body bed"
		bedtools groupby -i "${gtfname}"_gtf.bed -g 4 -c 1,2,3,4,5,6,7 -o distinct,min,max,distinct,distinct,distinct,distinct | awk -F "\t" '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8}' > "${gtfname}"_genes_1.bed		
		echo "creating final full genbody bed"
		uniq "${gtfname}"_genes_1.bed |  awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"gene",$6,$7 }' | bedtools sort | uniq > "${gtfname}"_genes.bed
		echo "creating tss and tes"

		awk -F "\t" '{ if ( $6 == "+") {OFS="\t"; print $1,$2,$2+1,$4,"tss",$6,$7"\n"$1,$3-1,$3,$4,"tes",$6,$7} else if ( $6 == "-") {OFS="\t"; print $1,$2,$2+1,$4,"tes",$6,$7"\n"$1,$3-1,$3,$4,"tss",$6,$7} }' "${gtfname}"_genes.bed > "${gtfname}"_tss_tes.bed
		echo " mergin start and stop codons"
		grep "start_codon\|stop_codon" "${gtfname}"_gtf.bed | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${gtfname}"_start_stop.bed
		echo "creting UTRs"
		cat "${gtfname}"_start_stop.bed "${gtfname}"_tss_tes.bed | bedtools sort > "${gtfname}"_tss_start_stop_tes.bed
		echo "removinf small rnas from gtf.bed and genes.bed and creating introns"
		grep -v -w "miRNA\|snoRNA\|snRNA\|sRNA"  "${gtfname}"_gtf.bed > "${gtfname}"_gtf_nomirs.bed
		grep -v -w "miRNA\|snoRNA\|snRNA\|sRNA"  "${gtfname}"_genes.bed >  "${gtfname}"_genes_nomirs.bed
		bedtools subtract -s -a "${gtfname}"_genes_nomirs.bed -b "${gtfname}"_gtf_nomirs.bed |  awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"intron",$6,$7 }' | bedtools sort > "${gtfname}"_introns.bed
		echo " create final 5utr and 3utr coordinates"
                awk  -F "\t" '{$5 == "start_codon" || $5 tes}' "${gtfname}"_tss_start_stop_tes.bed
		awk -F "\t" '($5 == "stop_codon") || ($5 == "tes") ' "${gtfname}"_tss_start_stop_tes.bed | grep -w protein_coding | bedtools groupby -g 4 -c 1,2,3,4,5,6,7 -o distinct,min,max,distinct,distinct,distinct,distinct | awk -F "\t" '{
		OFS="\t"; print $2,$3,$4,$5,"3utr",$7,$8}' | awk -F "\t" '{ if ( $6 == "+") {
 		OFS="\t"; print $1,$2+2,$3,$4,$5,$6,$7} else if ( $6 == "-") {
 		OFS="\t"; print $1,$2,$3-2,$4,$5,$6,$7} }'> "${gtfname}"_3utr.bed
 		awk -F "\t" '($5 == "start_codon") || ($5 == "tss") ' "${gtfname}"_tss_start_stop_tes.bed | grep -w protein_coding | bedtools groupby -g 4 -c 1,2,3,4,5,6,7 -o distinct,min,max,distinct,distinct,distinct,distinct | awk -F "\t" '{
 		OFS="\t"; print $2,$3,$4,$5,"5utr",$7,$8}' | awk -F "\t" '{ if ( $6 == "+") {
 		OFS="\t"; print $1,$2,$3-2,$4,$5,$6,$7} else if ( $6 == "-") {
 		OFS="\t"; print $1,$2+2,$3,$4,$5,$6,$7} }'> "${gtfname}"_5utr.bed
		cat "${gtfname}"_5utr.bed "${gtfname}"_3utr.bed | awk -F "\t" '$3 - $2 > 1' | bedtools sort > "${gtfname}"_3utr_5utr_2.bed
		echo "creating final utr stuff"  
		echo "extracting small rna information"
		grep -w "miRNA\|snoRNA\|snRNA\|sRNA" "${gtfname}"_gtf.bed | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${gtfname}"_gtf_mirs.bed
		echo "removing introns from utrs"
		bedtools subtract -s -a "${gtfname}"_3utr_5utr_2.bed -b "${gtfname}"_introns.bed | bedtools sort > "${gtfname}"_3utr_5utr.bed
		echo "mergin other than protein coding"
		grep -v -w protein_coding "${gtfname}"_gtf.bed | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${gtfname}"_gtf_no_protein_coding.bed
		echo "extracting CDS"
		awk -F "\t" '$5 == "CDS"' "${gtfname}"_gtf.bed | bedtools sort | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct > "${gtfname}"_cds.bed
		echo "extracting splice sites"
		awk -F "\t" '{ if ( $6 == "+") {OFS="\t"; print $1,$2-1,$2+1,$4,"5ss",$6,$7"\n"$1,$3-1,$3+1,$4,"3ss",$6,$7} else if ( $6 == "-") {OFS="\t"; print $1,$2-1,$2+1,$4,"3ss",$6,$7"\n"$1,$3-1,$3+1,$4,"5ss",$6,$7} }' "${gtfname}"_introns.bed > "${gtfname}"_5ss_3ss.bed
		echo "pooling all annotations together"
 		echo "preparing tRNAs"
 		awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,"exon",$6,"tRNA"}' $tRNA_DB > tRNAs.bed
		cat tRNAs.bed  "${gtfname}"_gtf_mirs.bed > smallRNA.bed
		cat tRNAs.bed "${gtfname}"_gtf_no_protein_coding.bed  "${gtfname}"_introns.bed "${gtfname}"_3utr_5utr.bed "${gtfname}"_5ss_3ss.bed "${gtfname}"_cds.bed "${gtfname}"_start_stop.bed | bedtools sort | uniq > "${gtfname}"_combined_no_small_rRNAs.bed
		bedtools subtract -s -a "${gtfname}"_combined_no_small_rRNAs.bed  -b smallRNA.bed > "${gtfname}"_combined_small_rRNAs_subtracted.bed
                cat smallRNA.bed "${gtfname}"_combined_small_rRNAs_subtracted.bed | sed 's/,/|/g' | awk -F "\t" '$3 - $2 > 1' | bedtools sort | uniq > "${gtfname}"_combined_small_rRNAs_subtracted_with_smallRNA.bed
		
		cat $gtf | awk -F "\t" '$3 == "exon" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "CDS"' | grep -v "##" | sed 's/"//g' | sed 's/;//g' | sed 's/\s/\t/g' | sed 's/\t\t/\t/g' | sed 's/"="/\t/g' | awk '{
                OFS="\t"; for (I=1; I<NF;I++) { if ($I== "gene_name") n=$(I+1); if ($I == "gene_type") t=$(I+1); } print $1,$4-1,$5,n,$3,$7,t}' \
                |  sed 's/,/-/g' | bedtools sort | uniq | awk -F "\t" '$4 !~ "Y_RNA"' | grep -v transcript  > "${gtfname}"_gtf_all.bed
		bedtools intersect -wa -v -s -a "${gtfname}"_gtf_all.bed -b "${gtfname}"_combined_small_rRNAs_subtracted_with_smallRNA.bed > "${gtfname}"_gtf_all_no_combined_small_rRNAs_subtracted_with_smallRNA.bed
		                cat "${gtfname}"_gtf_all_no_combined_small_rRNAs_subtracted_with_smallRNA.bed "${gtfname}"_combined_small_rRNAs_subtracted_with_smallRNA.bed | bedtools sort | uniq > "${gtfname}"_Paralus.bed
 






		rm "${gtfname}"_gtf.bed
		rm "${gtfname}"_genes.bed
		rm "${gtfname}"_genes_1.bed
		rm "${gtfname}"_3utr_5utr.bed
		rm "${gtfname}"_gtf_no_protein_coding.bed
		rm "${gtfname}"_tss_tes.bed
		rm "${gtfname}"_start_stop.bed
		rm "${gtfname}"_tss_start_stop_tes.bed
		rm "${gtfname}"_5utr.bed 
		rm "${gtfname}"_3utr.bed 
		rm "${gtfname}"_gtf_nomirs.bed
		rm "${gtfname}"_introns.bed
		rm "${gtfname}"_gtf_mirs.bed
		rm "${gtfname}"_cds.bed
		rm "${gtfname}"_5ss_3ss.bed
		rm "${gtfname}"_genes
		rm "${gtfname}"_3utr_5utr_2.bed
		rm "${gtfname}"_genes_nomirs.bed
		rm tRNAs.bed
		rm smallRNA.bed
		rm "${gtfname}"_combined_no_small_rRNAs.bed
		rm "${gtfname}"_combined_small_rRNAs_subtracted.bed
		rm "${gtfname}"_combined_small_rRNAs_subtracted_with_smallRNA.bed
		rm "${gtfname}"_gtf_all.bed
		rm "${gtfname}"_gtf_all_no_combined_small_rRNAs_subtracted_with_smallRNA.bed
		rm "${gtfname}"_3utr.bed
		rm "${gtfname}"_3utr_5utr_2.bed
		rm "${gtfname}"_5utr.bed



echo 

#################################
#####         DONE        #######
#################################
