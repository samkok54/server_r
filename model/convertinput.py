import sys
input_filename = sys.argv[1] #M1SN70SD777AL06.1.out
output_filename = sys.argv[2] #M1SN70SD777AL06.txt
SNP_number = int(sys.argv[3])

with open(input_filename,'r') as f:
    content = f.readlines()
f.close()
content = [x.split() for x in content] 

with open(output_filename,'w') as f2:
    for i in range(SNP_number):
        f2.write('SNP'+str(i+1)+'\t')
    f2.write('CLASS\r')
    for i in range(len(content)):
        for j in range(1,10+1):
            f2.write(content[i][j]+'\t')
        for j in range(21,(SNP_number-10)+21):
            f2.write(content[i][j]+'\t')
        f2.write(content[i][0]+'\r')
f2.close()
