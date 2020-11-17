import random

MAX_VOLUMN = 10
MAX_WEIGHT = 10
MAX_COORD = 10

def generate_test(N, M, file_name):
    with open(file_name, 'w') as f:
        warehouse = random.choices(range(MAX_COORD), k=2)
        f.writelines(' '.join([str(warehouse[0]), str(warehouse[1])]) + '\n')
        f.writelines(' '.join([str(N), str(M)]) + '\n')
        for i in range(N):
            xcoord, ycoord = random.choices(range(MAX_COORD), k=2)
            volumn = random.randint(1, MAX_VOLUMN)
            weight = random.randint(1, MAX_WEIGHT)
            f.writelines(' '.join([str(xcoord), str(ycoord), str(volumn), str(weight)]) + '\n')