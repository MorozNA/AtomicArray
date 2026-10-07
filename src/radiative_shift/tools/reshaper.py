def matrix_to_blocks(matrix, block_m, block_n):
    """Return axes (block row, block column, row in block, column in block)

    block_m and block_n are the row and column sizes of each block
    """
    B = matrix.reshape((-1, block_m, matrix.shape[1] // block_n, block_n))
    return B.transpose((0, 2, 1, 3))


def blocks_to_matrix(matrix, t=False):
    """Flatten a block tensor with axes (block row, block column, row, column)

    t=True exchanges the block-grid axes and preserves the axes within blocks
    """
    if t:
        B = matrix.transpose((1, 2, 0, 3))
    else:
        B = matrix.transpose((0, 2, 1, 3))
    return B.reshape((
        B.shape[0] * B.shape[1],
        B.shape[2] * B.shape[3],
    ))
