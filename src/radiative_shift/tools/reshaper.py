def reshape_to_blocks(matrix, block_m, block_n):
    B = matrix.reshape((-1, block_m, matrix.shape[1] // block_n, block_n))
    return B.transpose((0, 2, 1, 3))

# Quick test for function

# block_m = 3  # rows in block
# block_n = 3  # cols in block
# A = np.arange(6*6).reshape(6,6)  # sample data
# B = A.reshape((A.shape[0] // block_m, block_m, A.shape[1]//block_n, block_n))
# C = B.transpose((0,2,1,3))
#
# print(C)
