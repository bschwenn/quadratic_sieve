# Only one public function, call using find_dependencies(matrix)
from itertools import chain
import numpy as np


def mod_2_representation(mat):
    """
    returns mod 2 representation of the matrix mat
    """
    return mat % 2

def add_mod_2(r1, r2):
    """
    adds two rows of the matrix mod 2
    """
    return r1 ^ r2


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
    free_row_indices = [idx for idx,truth in enumerate(marks) if not truth]
    solutions = []
    for free_row_index in free_row_indices:
        solution = []
        one_columns = [index for index, value in enumerate(mat[free_row_index]) if mat[free_row_index][index] == 1]
        for row_index, row in enumerate(mat):
            for column_index in one_columns:
                if row[column_index] == 1 and marks[row_index]:
                    solution.append(row_index)
                    break
        if solution:
            solution.append(free_row_index)
            solutions.append([1 if i in solution else 0 for i in range(mat.shape[0])])
    return solutions

def find_dependencies(mat):
    """
    :param mat: 2d np array where each row is the exponent vector of the prime factorization.
    :return: ret: 2d np array where each row is a vector in the mod 2 left nullspace of mat.
    """
    ret = mod_2_representation(mat)
    ret, marks = eliminate(ret)
    return get_solution_rows(ret, marks)
