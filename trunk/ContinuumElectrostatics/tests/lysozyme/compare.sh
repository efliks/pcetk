paste pdynamo/gintr.dat qmpb/gintr.dat  | awk 'NR > 1 { print ($5 - $11); }'
