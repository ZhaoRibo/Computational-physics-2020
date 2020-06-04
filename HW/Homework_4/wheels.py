import copy

"""***************************
    * 矩阵类定义于矩阵运算 *
***************************"""


class Matrix(object):
    """
    简单定义矩阵类
    """
    def __init__(self, origin_list):
        assert isinstance(origin_list, list), "Input is not a list."

        if isinstance(origin_list[0], list):
            right_list = origin_list
        else:
            right_list = []
            right_list.append(origin_list)

        self.matrix = right_list
        self.shape = [len(self.matrix), len(self.matrix[0])]
        self.row = len(self.matrix)
        self.col = len(self.matrix[0])

    def __getitem__(self, index):
        """索引读取重载"""
        return self.matrix[index]

    def __setitem__(self, key, value):
        """索引赋值重载"""
        self.matrix[key] = value

    def get_row(self, n):
        """获得第n行"""
        return self.matrix[n]

    def get_col(self, n):
        """获得第n列"""
        n_col = [[self.matrix[i][n]] for i in range(self.row)]
        return n_col

    def pop_row(self, n):
        """删除第n行"""
        self.matrix.pop(n)
        self.shape[0] -= 1
        self.row -= 1

    def pop_col(self, n):
        """删除第n列"""
        for i in range(self.row):
            self.matrix[i].pop(n)
        self.shape[1] -= 1
        self.col -= 1

    def __add__(self, other):
        """矩阵加法"""
        assert isinstance(other, Matrix), "第二个实例不是Matrix"
        assert other.shape == self.shape, "两矩阵维度不匹配，不能相加"

        add_mat = Matrix(zero_mat(self.row, self.col))
        for i in range(self.row):
            for j in range(self.col):
                add_mat.matrix[i][j] = self.matrix[i][j] + other.matrix[i][j]
        return add_mat

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """矩阵减法"""
        assert isinstance(other, Matrix), "第二个实例不是Matrix"
        assert other.shape == self.shape, "两矩阵维度不匹配，不能相减"

        sub_mat = Matrix(zero_mat(self.row, self.col))
        for i in range(self.row):
            for j in range(self.col):
                sub_mat.matrix[i][j] = self.matrix[i][j] - other.matrix[i][j]
        return sub_mat

    def __rsub__(self, other):
        assert isinstance(other, Matrix), "第二个实例不是Matrix"
        assert other.shape == self.shape, "两矩阵维度不匹配，不能相减"

        sub_mat = Matrix(zero_mat(self.row, self.col))
        for i in range(self.row):
            for j in range(self.col):
                sub_mat.matrix[i][j] = other.matrix[i][j] - self.matrix[i][j]
        return sub_mat

    def __pos__(self):
        """正号"""
        return Matrix(self)

    def __neg__(self):
        """负号"""
        mat = Matrix(zero_mat(self.row, self.col))
        for i in range(self.row):
            for j in range(self.col):
                mat.matrix[i][j] = -self.matrix[i][j]
        return mat

    def __mul__(self, other):
        """矩阵乘法；右侧数乘"""
        if isinstance(other, Matrix):
            assert self.col == other.row, "两矩阵行列不匹配，不能相乘"

            if (self.row == 1) & (other.col == 1):
                mul_n = 0
                for i in range(self.col):
                    mul_n += self[0][i] * other[i][0]
                return mul_n
            else:
                mul_mat = Matrix(zero_mat(self.row, other.col))
                for i in range(self.row):
                    for j in range(other.col):
                        for k in range(self.col):
                            mul_mat.matrix[i][
                                j] += self.matrix[i][k] * other.matrix[k][j]
                return mul_mat
        elif isinstance(other, int) | isinstance(other, float):
            mul_mat = Matrix(zero_mat(self.row, self.col))
            for i in range(self.row):
                for j in range(self.col):
                    mul_mat.matrix[i][j] = self.matrix[i][j] * other
            return mul_mat
        else:
            assert False, "第二个实例不是Matrix或数"

    def __rmul__(self, other):
        return self * other

    def __matmul__(self, other):
        """矩阵内积"""
        assert isinstance(other, Matrix), "输入参数不是矩阵，不能点乘"
        assert other.shape == self.shape, "矩阵维数不匹配，不能点乘"

        mul_mat = Matrix(zero_mat(self.row, self.col))
        for i in range(self.row):
            for j in range(self.col):
                mul_mat[i][j] = self.matrix[i][j] * other.matrix[i][j]
        return mul_mat

    def __rmatmul__(self, other):
        return self @ other

    def __pow__(self, other):
        """矩阵求幂"""
        assert self.row == self.col, "矩阵不是方阵，不能求幂"
        assert isinstance(other, int), "指数不是整数"
        assert other > 0, "指数不为正"

        pow_mat = Matrix(identity_mat(self.row))
        for i in range(other):
            pow_mat = pow_mat * self
        return pow_mat

    def __eq__(self, other):
        """=="""
        if isinstance(other, Matrix):
            return (self.shape == other.shape
                    and all(a == b for a, b in zip(self, other)))
        else:
            return NotImplemented

    def T(self):
        """矩阵转置"""
        t_mat = Matrix(zero_mat(self.col, self.row))
        for i in range(self.col):
            for j in range(self.row):
                t_mat.matrix[i][j] = self.matrix[j][i]
        return t_mat

    def trace(self):
        """矩阵的迹"""
        assert self.row == self.col, "矩阵不是方阵，无法求迹"
        t = 0
        for i in range(self.row):
            t += self.matrix[i][i]
        return t

    def det(self):
        """矩阵行列式"""
        assert self.row == self.col, "矩阵不是方阵，无法计算行列式"

        d = 1
        mat = lu_decomposition(self)
        for i in range(self.row):
            d = d * mat.matrix[i][i]
        return d

    def conj(self):
        """共轭矩阵"""
        m = copy.deepcopy(self)
        for i in range(self.row):
            for j in range(self.col):
                if isinstance(self.matrix[i][j], complex):
                    m.matrix[i][j] = self.matrix[i][j].conjugate()
        return m


def zero_mat(row, col):
    """
    生成零矩阵，Matrix(zero_mat(int row,int col))
    """
    assert isinstance(row, int) & isinstance(col,
                                             int), "Type of input is not int."

    mat = [[0 for i in range(col)] for i in range(row)]
    return mat


def identity_mat(n):
    """
    生成单位矩阵，Matrix(identity_mat(int n))
    """
    assert isinstance(n, int), "Type of input is not int."

    mat = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        mat[i][i] = 1
    return mat


def diagonal_mat(dia_list):
    """
    生成对角矩阵
    """
    assert isinstance(dia_list, list), "Input must be a list"

    n = len(dia_list)
    mat = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        mat[i][i] = dia_list[i]
    return mat


"""
线性方程组求解
直接求解
"""


def gauss_method(M, b):
    '''
    高斯消元法->上三角矩阵
    '''
    n = M.row
    for k in range(n - 1):
        assert M[k][k] != 0, "系数矩阵主元存在‘0’项"
        for i in range(k + 1, n):
            c = -M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] = M[i][j] + c * M[k][j]
            b[i][0] = b[i][0] + c * b[k][0]


def col_gauss(M, b):
    '''
    列主元消元法->上三角矩阵
    '''
    n = M.row
    for k in range(n - 1):
        a = k
        for l in range(k, n):
            if abs(M[l][k]) > abs(M[a][k]):
                a = l
        if a != k:
            d = copy.deepcopy(M[k])
            M[k] = copy.deepcopy(M[a])
            M[a] = copy.deepcopy(d)
            d = b[k]
            b[k] = b[a]
            b[a] = d
        assert M[k][k] != 0, "系数矩阵主元存在‘0’项"
        for i in range(k + 1, n):
            c = -M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] = M[i][j] + c * M[k][j]
            b[i][0] = b[i][0] + c * b[k][0]


def lu_decomposition(M):
    '''基于高斯消元法的直接LU分解算法'''
    n = M.row
    lu = copy.deepcopy(M)
    for k in range(n - 1):
        assert lu[k][k] != 0, "系数矩阵主元 存在‘0’项"
        for i in range(k + 1, n):
            lu[i][k] = lu[i][k] / lu[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                lu[i][j] = lu[i][j] - lu[i][k] * lu[k][j]
    return lu


def col_lu_decomposition(M):
    '''列主元消元法的直接LU分解算法'''
    n = M.row
    p = [i for i in range(n)]
    lu = copy.deepcopy(M)
    for k in range(n - 1):
        a = k
        for l in range(k, n):
            if abs(M[l][k]) > abs(M[a][k]):
                a = l
        if a != k:
            d = copy.deepcopy(M[k])
            M[k] = copy.deepcopy(M[a])
            M[a] = copy.deepcopy(d)
            d = p[k]
            p[k] = p[a]
            p[a] = d
        assert lu[k][k] != 0, "系数矩阵主元 存在‘0’项"
        for i in range(k + 1, n):
            lu[i][k] = lu[i][k] / lu[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                lu[i][j] = lu[i][j] - lu[i][k] * lu[k][j]
    return lu, p


def diagonal_msolve(M, b):
    '''对角矩阵的解法'''
    n = M.row
    x = copy.deepcopy(b)
    for i in range(n):
        x[i][0] = b[i][0] / M[i][i]
    return x


def up_tri_msolve(M, b):
    '''上三角形矩阵的回代算法'''
    n = M.row
    x = copy.deepcopy(b)
    for i in reversed(range(n)):
        assert M[i][i] != 0, "系数矩阵主元存在‘0’项"
        for j in range(n - 1, i, -1):
            x[i][0] = x[i][0] - M[i][j] * x[j][0]
        x[i][0] = x[i][0] / M[i][i]
    return x


def low_tri_msolve(M, b):
    '''下三角形矩阵的前代算法'''
    n = M.row
    x = copy.deepcopy(b)
    for i in range(n):
        assert M[i][i] != 0, "系数矩阵主元存在‘0’项"
        for j in range(i):
            x[i][0] = x[i][0] - M[i][j] * x[j][0]
        x[i][0] = x[i][0] / M[i][i]
    return x


def sym_tri_msolve(M, b):
    '''对称矩阵的解法'''
    # n = M.row
    pass


def Cholesky(M):
    '''
    对称正定矩阵的Cholesky分解（平方根法）
    A = LL^T
    '''
    n = M.row
    for j in range(n):
        for k in range(j - 1):
            M[j][j] = M[j][j] - M[j][k]**2
        M[j][j] = M[j][j]**0.5
        for i in range(j + 1, n):
            for k in range(j - 1):
                M[i][j] = M[i][j] - M[i][k] * M[j][k]
            M[i][j] = M[i][j] / M[j][j]

