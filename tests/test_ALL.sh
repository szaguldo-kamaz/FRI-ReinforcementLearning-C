#!/bin/bash

for i in test_*.sh; do
    echo
    echo $i
    echo
    ./$i
done
