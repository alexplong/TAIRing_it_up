import sys

# Convert GFF3 file into GTF formatted file

inFile = open(sys.argv[1],'r')

for line in inFile:
  if line[0] != '#':
    # Split line into columns by tab
    data = line.strip().split('\t')

    ID = ''

    # If the feature is a gene
    if data[2] == "gene":
      #get the id
      ID = data[-1].split('ID=')[-1].split(';')[0]

    # If the feature is anything else
    else:
      # Get the parent as the ID
      ID = data[-1].split('Parent=')[-1].split(';')[0]

    # Modify the last column
    data[-1] = 'gene_id "' + ID + '"; transcript_id "' + ID

    # Print out this new GTF line
    print '\t'.join(data)
