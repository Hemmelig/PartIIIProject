import pickle
import sys
file_name = sys.argv[1]

pkl_file = open(file_name, 'rb')
dict = pickle.load(pkl_file)
pkl_file.close()

print(dict)
print(dict[0].stats)
