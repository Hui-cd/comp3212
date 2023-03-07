import numpy as np
blosum50 =np.loadtxt('blosum50.txt', dtype=int)
blosum_index_map = dict(zip(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],[i for i in range(len(blosum50))]))

def get_blosum_score(a,b):
    
    return blosum50[blosum_index_map[a],blosum_index_map[b]]


def smith_waterman(protein1,protein2):
    d = -8
    height = len(protein2)+1
    width = len(protein1)+1
    
    matrix = np.zeros((height,width))
    # smith waterman algorithm achieve
    for i in range(1,height):
        for j in range(1,width):
            matrix[i,j] = max(0,matrix[i-1,j]+d,matrix[i,j-1]+d,matrix[i-1,j-1]+get_blosum_score(protein1[j-1],protein2[i-1]))
    # get max value and max value coordinate
    max_value = 0
    max_value_coordinate = (0,0)
    for i in range(height):
        for j in range(width):
            if matrix[i,j]>max_value:
                max_value = matrix[i,j]
                max_value_coordinate = (i,j)
    # get path
    path = []
    i = max_value_coordinate[0]
    j = max_value_coordinate[1]
    while i>0 or j>0:
        if i>0 and j>0 and matrix[i,j] == matrix[i-1,j-1]+get_blosum_score(protein1[j-1],protein2[i-1]):
            path.append((i-1,j-1,matrix[i,j],'diagonal'))
            i-=1
            j-=1
        elif i>0 and matrix[i,j] == matrix[i-1,j]+d:
            path.append((i-1,j,matrix[i,j],'up'))
            i-=1
        elif j>0 and matrix[i,j] == matrix[i,j-1]+d:
            path.append((i,j-1,matrix[i,j],'left'))
            j-=1
    # get protein1 and protein2 alignment
    protein1_alignment = ''
    protein2_alignment = ''
    for i in range(len(path)-1,-1,-1):
        if path[i][3] == 'diagonal':
            protein1_alignment += protein1[path[i][1]]
            protein2_alignment += protein2[path[i][0]]
        elif path[i][3] == 'up':
            protein1_alignment += '-'
            protein2_alignment += protein2[path[i][0]]
        elif path[i][3] == 'left':
            protein1_alignment += protein1[path[i][1]]
            protein2_alignment += '-'
    return path,matrix, protein1_alignment, protein2_alignment

print(smith_waterman('SALPQPTTPVSSFTSGSMLGRTDTALTNTYSAL','PSPTMEAVTSVEASTASHPHSTSSYFATTYYHLY'))


print(smith_waterman('HEAGAWGHEE','PAWHEAE'))