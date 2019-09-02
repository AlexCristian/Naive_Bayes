fq=$1

echo converting ${fq} to ${fq%.txt}.fna.
cat ${fq} | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > ${fq%.txt}.fna
