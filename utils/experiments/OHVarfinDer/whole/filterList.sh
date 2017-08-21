#! /bin/bash

declare -a filterArray=(\
"low_mapping_quality" \
"germline_variant_too_close" \
"germline_variants_overlapped" \
"germline_indel_too_close" \
"too_many_softclips_nearby" \
"too_many_indels_nearby" \
"too_many_reads_in_window" \
"too_few_reads" \
"something_wrong" \
"ignore_BF_Computation" \
)

filterExpression=""
for ((i = 0; i < ${#filterArray[@]}; i++)) {
	filterExpression=" | grep -v ${filterArray[i]}"
}