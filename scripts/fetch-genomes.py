#!/usr/bin/python

# import urllib
import os
import sys
import time
import urllib.request

if len(sys.argv) != 3:
    print ("USAGE: fetch_genome.py <genome_id_list> <out_dir>")
    sys.exit(1)

url_template = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text"

if not os.path.exists(sys.argv[2]):
    os.mkdir(sys.argv[2])

for id in open(sys.argv[1]):
    id = id.strip()
    if id == "":
        continue

    sys.stdout.write("Fetching %s..." % id)
    sys.stdout.flush()
    out_file = os.path.join(sys.argv[2], id + ".fasta")
    if os.path.exists(out_file):
        print("already fetched")

    url = url_template
    try:
        with urllib.request.urlopen(url_template % id) as response:
            content = response.read()
            with open(out_file, "wb") as file:
                file.write(content)
        print(f"Content from {url_template % id} written to {out_file}")
    except urllib.error.URLError as e:
        print(f"Error opening URL: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

    # open(out_file, "w").write(urlopen(url_template % id).read())
    # print("Done")
    time.sleep(1.0/3)

