for id in $(parsec  histories show_dataset_collection b4b580c4c286c45f 60311dd1b2959818 | jq .elements[].object.id -r); do 
  parsec histories download_dataset $id --file_path file.fastq --use_default_filename
#   parsec  histories download_dataset b4b580c4c286c45f e629965666ce4ce5 .
  # process
  scHicDemultiplex -f file.fastq -s ../samples.txt -b ../GSE94489_README.txt -t 20
  # upload whatever via ftp
  curl -T "{$(echo *.gz | tr ' ' ',')}" ftp://galaxy.uni-freiburg.de -u wolffj@informatik.uni-freiburg.de:Password
  rm file.fastq
  rm *.gz
done