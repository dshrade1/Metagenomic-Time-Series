# ./listReadFiles.sh

# this all works
DATA_DIRS=*/
for i in $DATA_DIRS; do
DATA_FILES="$i"*
for k in $DATA_FILES; do
echo "${k:0:11} ${k:12:23}" >> readFileList2.txt
done
done

# this works
#cd $1
#DATA_DIRS=*
#for i in $DATA_DIRS; do
#cd "$i"
#DATA_FILES="$i"/*
#for k in $DATA_FILES; do
#echo "$i $k" >> readFileList.txt
#done
#done


