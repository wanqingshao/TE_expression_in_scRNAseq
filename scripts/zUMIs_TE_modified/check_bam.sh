#!/bin/bash

echo $1
samtools view $1 | head -n 1
