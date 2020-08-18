#!/usr/bin/bash

wget http://vald.astro.uu.se/~vald/FTP/MingjieJian.$1.gz
gunzip MingjieJian.$1.gz

# python3 line_data.py MingjieJian.$1 files/linelist/vald_$2

# rm MingjieJian.$1