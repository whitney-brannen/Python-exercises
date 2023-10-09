# Whitney Brannen
# wbranne1@uncc.edu

"""

Explanation:

To complete this assignment, I first created a fasta parser generater to sort through each header in the fasta file and
yeild its corresponding sequence together.  I then sent this into a generator that would yield the amino acid counts from
a user selected reading frame.  I chose to have the user input their gene of choice and display the counts for only that
gene becasue of the large size of the file.  Next, I used many functions to get to the codon usage bias table.  First, I 
created a function that would create a dictionary for the sum number of counts for each codon (64 total).  Then I calculated the
total counts for each amino acid and create a corresponding list for the total counts of codons (the codons are in a list
in the order they appear in the aa_dict, the first codon = Met, the second and third= Phe and so on. I duplicated the total
number for the amount of codons the amino acid has, so the lists are the same length and each value at each index can be divided
by its cooresponing value in the other list).  After this, I divided the two lists to get the codon frequencies per amino
acid.  Next I created a dictionary for the single letter code per amino acid to include in the table to make it more informative.
Lastly, I formatted the table to be an 8 x 8 table seperated by tabs showing the codon, the amino acid it codes for, and
its frequency.  ***The table may take a few moments to load as it is parsing through the file, but will end up showing results.

"""

aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'], 
'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 
'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'], 
'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

single_letter_dict = {'M':['ATG'], 'F':['TTT', 'TTC'], 'L':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C':['TGT', 'TGC'], 'Y':['TAC', 'TAT'], 'W':['TGG'], 'P':['CCT', 'CCC', 'CCA', 'CCG'], 'H':['CAT', 'CAC'], 
	'Q':['CAA', 'CAG'], 'R':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'I':['ATT', 'ATC', 'ATA'], 'T':['ACT', 'ACC', 'ACA', 'ACG'], 
	'N':['AAT', 'AAC'], 'K':['AAA', 'AAG'], 'S':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'V':['GTT', 'GTC', 'GTA', 'GTG'], 
	'A':['GCT', 'GCC', 'GCA', 'GCG'], 'D':['GAT', 'GAC'], 'E':['GAA', 'GAG'], 'G':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}

# create dictionary (d) of just codons and their amino acid from the aa_dict
d= {}
for key, value in aa_dict.items(): 
	for values in value:
		d[values] = key


def fasta_parser(fh):
	# Source: https://stackoverflow.com/questions/29805642/learning-to-parse-a-fasta-file-with-python
	# Was having issues getting all lines of the sequence to join under their one header
	# and this was the only solution I could find that could also handle a large dataset.
    header, seq = None, [] # establishg tuple to put header and cooresponding seq into
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if header: yield (header, ''.join(seq))
            header, seq = line, [] # line as header and empty list for seq
        else:
            seq.append(line) # all other lines are appended to the list but seperated by line
    if header: yield (header, ''.join(seq)) # .join() joins the seq list together to be a complete string of amino acids


def count_codons(seq, frame):
	counts = {}
	for key in aa_dict:
		counts[key] = 0
	for i in range(0+frame, len(seq), 3):
		codon = seq[i: i+3]
		for aa, codons in aa_dict.items():
			if codon in codons:
				counts[aa] += 1
	yield aa, counts


def count_per_seq(seq):
	counts = {}
	for key in d:
		counts[key] = 0
	for i in range(0, len(seq), 3):
		codon = seq[i: i+3]
		for codons, aa in d.items():
			if codon in codons:
				counts[codons] += 1
	yield counts


def sum_per_codon(counts, lst):
	count = {}
	for index, value in enumerate(lst):
		if index in counts.keys():
			count[index] += value
		else:
			count[index] = value
	return count


def count_totals(sum_counts):
	total = {}
	number = list(d.keys())
	for index, codon in enumerate(number):
		for key, values in aa_dict.items():
			if codon in values:
				if key in total.keys():
					total[key] += sum_counts[index]
				else:
					total[key] = sum_counts[index]
	return total


def extend_totals(total):
	division = list("")
	for aa, total in total.items():
		if aa == "Met":
			division = [total]
		if aa == "Phe":
			division += [total] * 2
		if aa == "Leu":
			division += [total] * 6
		if aa == "Cys":
			division += [total] * 2
		if aa == "Tyr":
			division += [total] * 2
		if aa == "Trp":
			division += [total] 
		if aa == "Pro":
			division += [total] * 4
		if aa == "His":
			division += [total] * 2
		if aa == "Gln":
			division += [total] * 2
		if aa == "Arg":
			division += [total] * 6
		if aa == "Ile":
			division += [total] * 3
		if aa == "Thr":
			division += [total] * 4
		if aa == "Asn":
			division += [total] * 2
		if aa == "Lys":
			division += [total] * 2
		if aa == "Ser":
			division += [total] * 6
		if aa == "Val":
			division += [total] * 4
		if aa == "Ala":
			division += [total] * 4
		if aa == "Asp":
			division += [total] * 2
		if aa == "Glu":
			division += [total] * 2
		if aa == "Gly":
			division += [total] * 4
		if aa == "*":
			division += [total] * 3
	return division
	

def frequencies(sum_counts, total_vector):
	frequency = list()
	for i in range(len(sum_counts)):
		if total_vector[i] == 0:
			frequency += [0.00]
		else:
			frequency += [round((sum_counts[i]/total_vector[i]),3)]
	return frequency


def list_of_codons(): 
	codon = list("")
	for key, value in aa_dict.items():
		codon += value
	return codon


def codon_freq_tuples(codon, freqs):
	final_list = []
	for i in range(len(codon)):
		final_list += [(codon[i], freqs[i])]
	return final_list


def single_letter_aa(final_list):
	single_letter = []
	dictionary = []
	for key, values in single_letter_dict.items():
		for value in values:
			dictionary += (key, value)
	for codon, freq in final_list:
		if codon in dictionary:
			i = dictionary.index(codon)
			pos = dictionary[i-1]
			single_letter += [(codon, pos, freq)]
	return single_letter


def final_format(single_letter_format):
	final_format = "Codon Usage Bias Table for M. domestica: (Codon, Amino Acid, Frequency)\n\n"
	for index, (codon, letter, freq) in enumerate(single_letter_format):
		if index == 7 or index == 15 or index == 23 or index == 31 or index == 39 or index == 47 or index == 55:
			final_format += f"{codon} {letter} {freq} \n\n"
		else:
			final_format += f"{codon} {letter} {freq} \t"
	return final_format


def main():
	# Part 1
	gene = input("Enter the gene record you are interested in for amino acid count: \n")
	with open('practice.fa') as fh:
		for header, seq in fasta_parser(fh):
			if gene in header:
				print(f"Found Gene:\n{header}")
				frame = int(input("Please select a reading frame (0, +1, +2): "))
				for aa, counts in count_codons(seq, frame): # only counting codons for chosen gene because of large file size
					final = f"\nAA counts at reading frame {frame}:\n"
					for key, value in counts.items():
						final += f"{key}: {value}\n" # formatting for readable list
					print(final)

	# Part 2
		print("\nCalculating Codon Usage Bias Table...\n\n")

		for counts in count_per_seq(seq): # count each codon in all seqs
			lst = list(counts.values())
			sum_counts = list(sum_per_codon(counts, lst).values()) # sum all counts for total for document
	
		# total counts per amino acid
		totals = count_totals(sum_counts)

		# replicate totals per amino acid to be consistent with # of codons
		total_vector = extend_totals(totals)

		# divide total per codon by total per aa for frequencies
		freqs = list(frequencies(sum_counts, total_vector))

		# create list of codons
		codon = list_of_codons()
	
		# match codon to frequency in the form of tuples
		final_list = codon_freq_tuples(codon, freqs)

		# add in the single aa letter code for visualization purposes:
		# inspiration from images at: http://www.geneinfinity.org/sp/sp_codonusage.html 
		single_letter_format = single_letter_aa(final_list)

		print(final_format(single_letter_format) + "\n")
					

if __name__ == '__main__':
	main()


