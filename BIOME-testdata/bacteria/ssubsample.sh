#!/bin/bash

for file in MI*
	do
	tail -100000 $file >subsample.$file
	done
