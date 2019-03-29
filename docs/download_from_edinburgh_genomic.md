Program to download from Ed Genomics website is "Aspera Connect".

export ASPERA_SCP_PASS=<password>

nohup ascp -P 33001 -O 33001 <user>@transfer.genomics.ed.ac.uk:11465_Kind_Peter/raw_data/20190204 . &