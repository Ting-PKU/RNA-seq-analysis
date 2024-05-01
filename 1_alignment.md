## 去接头
```
cat config |while read line;do arr=($line); fq1=${arr[0]}; fq2=${arr[1]};trim_galore -q 25 --phred33  -e 0.1 --stringency 2 --paired $fq1 $fq2 ;done
```
