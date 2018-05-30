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
"ignore_BF_Computation"
)

filterExpression="\"${filterArray[0]}"
for ((i = 1; i < ${#filterArray[@]}; i++)) {
	filterExpression="${filterExpression}|${filterArray[i]}"
}
filterExpression="${filterExpression}\""