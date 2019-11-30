for id in $(parsec  histories show_dataset_collection b4b580c4c286c45f 60311dd1b2959818 | jq .elements[].object.id -r); do 
  parsec histories download_dataset b4b580c4c286c45f $id .
#   parsec  histories download_dataset b4b580c4c286c45f e629965666ce4ce5 .
  # process
  FASTQ_FILE="$(echo *fastq-dump* } | cut -d' ' -f1).fastq.gz"
  mv *fastq-dump* $(echo *fastq-dump* } | cut -d' ' -f1).fastq.gz
  scHicDemultiplex -f "$FASTQ_FILE" -s ../samples.txt -b ../GSE94489_README.txt --threads 20
  # upload whatever via ftp
  rm $FASTQ_FILE
  curl -T "{$(echo demultiplexed/* | tr ' ' ',')}" ftp://galaxy.uni-freiburg.de -u wolffj@informatik.uni-freiburg.de:password
  rm demultiplexed -rf
done
