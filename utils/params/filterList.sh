#! /bin/bash

declare -a filterArray=(\
"low_mapping_quality" \
"too_many_softclips_nearby" \
"germline_indel_too_close" \
"ignore_BF_Computation" \
)

filterExpression="${filterArray[0]}"
for ((i = 1; i < ${#filterArray[@]}; i++)) {
	filterExpression="${filterExpression}|${filterArray[i]}"
}