#! /bin/bash

declare -a filterArray=()

filterExpression=""
for ((i = 0; i < ${#filterArray[@]}; i++)) {
	filterExpression=" | grep -v ${filterArray[i]}"
}