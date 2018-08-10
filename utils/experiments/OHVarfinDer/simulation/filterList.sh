#! /bin/bash

declare -a filterArray=(\
"something_wrong" \
)

filterExpression="${filterArray[0]}"
for ((i = 1; i < ${#filterArray[@]}; i++)) {
	filterExpression="${filterExpression}|${filterArray[i]}"
}