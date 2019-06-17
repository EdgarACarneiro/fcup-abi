class NumMatrix:

    mat: [[float]]

    def __init__(self, rows, cols, val=0.0):
        """create a matrix with the given dimensions: 
        number of rows and columns filled with the given value (by omission zero) """
        self.mat = [[val] * cols for _ in range(0, rows)]

    def __getitem__(self, n):
        """Get the row at the given n index"""
        return self.mat[n]

    def num_rows(self):
        """Get the number of rows"""
        return len(self.mat)

    def num_cols(self):
        """Get the number of columns"""
        return len(self.mat[0])

    def get_value(self, i, j):
        """Get value at index i, j"""
        if i > j:
            return self.mat[i][j]
        else:
            return self.mat[j][i]

    def set_value(self, i, j, value):
        """Set matrix value at cell i, j"""
        if i > j:
            self.mat[i][j] = value
        else:
            self.mat[j][i] = value

    def print_mat(self):
        """Pretty print the matrix"""
        # Stringfying matrix
        s = [[str(e) for e in row] for row in self.mat]

        # Length of each matrix column
        len_s = [max(map(len, col)) for col in zip(*s)]

        # Cell formatation
        formatation = '\t'.join('{{:{}}}'.format(x) for x in len_s)

        # Apply cell formation to each matrix element
        pretty_mat = [formatation.format(*row) for row in s]

        # Print pretty matrix
        print('\n'.join(pretty_mat))

    def min_dist_indexes(self):
        """Returns the row and column of the minimum matrix value"""
        m = self.mat[1][0]
        res = (1, 0)
        for i in range(1, self.num_rows()):
            for j in range(i):
                if self.mat[i][j] < m:
                    m = self.mat[i][j]
                    res = (i, j)
        return res

    def add_row(self, newrow):
        """Add the given row to the matrix"""
        self.mat.append(newrow)

    def add_col(self, newcol):
        """Add the given col to the matrix"""
        for r in range(self.num_rows()):
            self.mat[r].append(newcol[r])

    def remove_row(self, ind):
        """Remove the given row"""
        del self.mat[ind]

    def remove_col(self, ind):
        """Remove the given column"""
        for r in range(self.num_rows()):
            del self.mat[r][ind]

    def copy(self):
        """Create a new NumMatrix instance that is
        a copy of the current instance"""
        new_matrix = NumMatrix(self.num_rows(), self.num_cols())
        for i in range(self.num_rows()):
            for j in range(self.num_cols()):
                new_matrix.mat[i][j] = self.mat[i][j]
        return new_matrix
