import xml.etree.cElementTree as et
import gzip
import time

colnames = ["uid","common_taxa","num_members","members"]
print_interval = 1000000

uniref = {}
members = []
num_unirefs = 0

outfile = open("./uniref100_member_taxa_tbl.csv", "w")
outfile.write(",".join(colnames) + "\n")
logfile = open("./uniref100_member_taxa_tbl.log", "w")

xml_file = gzip.open('./uniref100.xml.gz', 'r')
# xml_file = open('./uniref100.mini.xml', 'r')
stime = time.time()
for event, elem in et.iterparse(xml_file, events=("start", "end")):
    if elem.tag == "{http://uniprot.org/uniref}entry":
        if event == "start":
            uniref = {c:"" for c in colnames}
            members = []
            uniref['uid'] = elem.attrib['id']
        elif event == "end":
            uniref['members'] = ":".join(members)
            num_unirefs += 1
            tmp = ",".join([uniref[c] for c in colnames])
            outfile.write(tmp + "\n")
            if num_unirefs%print_interval == 0:
                logfile.write("{},{},{:.2f}min\n".format(num_unirefs, tmp, (time.time()-stime)/60))
                print("{},{},{:.2f}min".format(num_unirefs, uniref, (time.time()-stime)/60))
    else:
        if event == "end":
            continue
    if "type" in elem.attrib:
        if (elem.tag == "{http://uniprot.org/uniref}property") & (elem.attrib['type']=="common taxon ID"):
            uniref['common_taxa'] = elem.attrib['value']
        if (elem.tag == "{http://uniprot.org/uniref}property") & (elem.attrib['type']=="member count"):
            uniref['num_members'] = elem.attrib['value']
        if (elem.tag == "{http://uniprot.org/uniref}property") & (elem.attrib['type']=="NCBI taxonomy"):
            members.append(elem.attrib['value'])
    elem.clear()
    
logfile.write("{},{},{:.2f}min\n".format(num_unirefs, tmp, (time.time()-stime)/60))    
print("{},{},{:.2f}min".format(num_unirefs, uniref, (time.time()-stime)/60))
