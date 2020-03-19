import random

# Create a list of the 24 sample names (1-24)
eid_list = list(range(1, 25))
# Shuffle the list twice
random.shuffle(eid_list)
random.shuffle(eid_list)

# Print the formatted output names
for eid in eid_list:
	print('2020-GCRSR-{:04d}'.format(eid))

print()

# Same for the over-clustered samples (names will be 25-48)
ocid_list = list(range(25, 49))
random.shuffle(ocid_list)
random.shuffle(ocid_list)

for ocid in ocid_list:
	print('2020-GCRSR-{:04d}'.format(ocid))
