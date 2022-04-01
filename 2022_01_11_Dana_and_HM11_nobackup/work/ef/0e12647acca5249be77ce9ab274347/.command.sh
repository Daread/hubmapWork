#!/bin/bash -ue
echo 'const log_data = {' > log_data.js
for file in HEK_read_metrics.log
do
    samp_name=$(basename $file | sed 's/_read_metrics.log//')
    echo "\"$samp_name\" :  \`" >> log_data.js
    cat $file >> log_data.js
    echo "\`," >> log_data.js
done
sed -i '$ s/,$//' log_data.js
echo '}' >> log_data.js

echo 'const full_log_data = {' >> log_data.js
for file in HEK_full.log 
do  
    samp_name=$(basename $file | sed 's/_full.log//')
    echo "\"$samp_name\" :  \`" >> log_data.js
    cat $file >> log_data.js
    echo "\`," >> log_data.js
done
sed -i '$ s/,$//' log_data.js
echo '}' >> log_data.js
