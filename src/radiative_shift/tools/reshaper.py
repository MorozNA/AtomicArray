def reshape_to_blocks(matrix, block_m, block_n):
    B = matrix.reshape((-1, block_m, matrix.shape[1] // block_n, block_n))
    return B.transpose((0, 2, 1, 3))


def reshape_to_matrix(matrix, t=False):
    if t:
        B = matrix.transpose((1, 2, 0, 3))
    else:
        B = matrix.transpose((0, 2, 1, 3))
    return B.reshape((B.shape[0] * B.shape[1], B.shape[2] * B.shape[3]))

