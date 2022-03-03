#!/bin/bash
awk  'NR>1{print $3}' $1 > cosa

while IFS="" read -r p || [ -n "$p" ]
do
    if [ -f "$p" ]; then
        continue
    else 
        echo "$FILE does not exist."
    fi
done < cosa

echo 'Check done!'

rm cosa
