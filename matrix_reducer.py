# Only one public function, call using find_dependencies(matrix)
from itertools import chain
import numpy as np


def mod_2_representation(mat):
    """
    returns mod 2 representation of the matrix mat
    """
    for row in mat:
        for i in range(len(row)):
            row[i] = 0 if row[i] % 2 == 0 else 1
    return mat


def add_mod_2(r1, r2):
    """
    adds two rows of the matrix mod 2
    """
    # assert (len(r1) == len(r2), "lengths of two rows are not equal!")
    for i in range(len(r1)):
        r1[i] += r2[i]
        r1[i] %= 2
    return r1


def eliminate(mat):
    """
    returns RREF of mat after gaussian elim, given mat is a mod 2 matrix. Also returns boolean list of
    independent/dependent rows.
    """
    mat = mat.transpose()  # transposing so we can work on rows and then transpose at the end
    marks = [False] * mat.shape[1]
    for index, row in enumerate(mat):
        for i in range(len(row)):
            if row[i] == 1:  # found 1 in this row, pivot
                marks[i] = True  # mark column as having a 1
                for j in chain(range(index), range(index+1, mat.shape[0])):
                    if mat[j][i] == 1:
                        mat[j] = add_mod_2(mat[j], row)  # add the two rows together
                break
    mat = mat.transpose()
    return mat, marks


def get_solution_rows(mat, marks):
    """

    :param mat:
    :param marks:
    :return: left nullspace of mat mod 2
    """
    solutions = []
    for index, mark in enumerate(marks):
        if not mark:  # found a row s.t. mark is False, i.e. a dependent row.
            soln = []
            for i, num in enumerate(mat[index]):
                if num == 1:
                    for idx, row in enumerate(mat):
                        if row[i] == 1 and marks[index]:
                            soln.append(idx)
                    # soln.append(i)
            soln.append(index)
            solutions.append(soln)
    ret = []
    for soln in solutions:
        row = np.zeros(mat.shape[1])
        for i in soln:
            row[i] = 1
        ret.append(row)
    return np.array(ret)


def find_dependencies(mat):
    """
    :param mat: 2d np array where each row is the exponent vector of the prime factorization.
    :return: ret: 2d np array where each row is a vector in the mod 2 left nullspace of mat.
    """
    ret = mod_2_representation(mat).transpose()
    ret, marks = eliminate(ret)
    return get_solution_rows(ret, marks)
