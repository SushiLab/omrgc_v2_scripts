# Zip all files

for f in ../data/processed/*.txt
do
 echo "Processing $f";
 gzip -f $f;
done
