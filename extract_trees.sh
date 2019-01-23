#!/bin/bash

path="$1"

for f in $(ls $path); do
  grep "ENS" $path/$f > $path/$f.tree
done
