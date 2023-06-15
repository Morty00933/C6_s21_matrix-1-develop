
#include "./s21_matrix.h"

#define M A->matrix
#define L A.matrix
#define BAD_MATRIX 1

int s21_check_matrix(matrix_t *A) {
  int error = 0;
  if (A->matrix != NULL && A->columns >= 1 && A->rows >= 1 && A != NULL)
    error = 1;
  return error;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (rows <= 0 || columns <= 0) {
    return BAD_MATRIX;
  }

  result->rows = rows;
  result->columns = columns;
  result->matrix = (double **)malloc(rows * sizeof(double *));

  if (result->matrix == NULL) {
    return ERROR_CAL;
  }

  for (int i = 0; i < rows; ++i) {
    result->matrix[i] = (double *)malloc(columns * sizeof(double));
    if (result->matrix[i] == NULL) {
      for (int j = 0; j < i; ++j) {
        free(result->matrix[j]);
      }
      free(result->matrix);
      return ERROR_CAL;
    }
  }

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      result->matrix[i][j] = 0;
    }
  }

  return OK;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) free(A->matrix[i]);
    free(A->matrix);
    A->columns = 0;
    A->rows = 0;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error = 1;

  // если обе матрицы имеют размер 1х1, проверим их значение без циклов
  if (A->rows == 1 && B->rows == 1 && A->columns == 1 && B->columns == 1) {
    if (fabs(A->matrix[0][0] - B->matrix[0][0]) > EPS) error = 0;
  } else if ((A->columns == B->columns && A->rows == B->rows) &&
             (A->columns + B->columns > 2 && A->rows + B->rows > 2)) {
    for (int i = 0; i < A->rows; i++)
      for (int a = 0; a < A->columns; a++)
        if (fabs(A->matrix[i][a] - B->matrix[i][a]) > EPS) error = 0;
  } else {
    error = 0;
  }
  return error;
}
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) {
    return ERROR;
  }
  if (A->rows != B->rows || A->columns != B->columns) {
    return ERROR_CAL;
  }
  int error = s21_create_matrix(A->rows, A->columns, result);
  if (error != OK) {
    return error;
  }
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
  return OK;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) && s21_check_matrix(B)) {
    if (A->rows == B->rows && A->columns == B->columns) {
      error = s21_create_matrix(A->rows, A->columns, result);
      if (!error) {
        for (int i = 0; i < A->rows; i++)
          for (int a = 0; a < A->columns; a++)
            result->matrix[i][a] = A->matrix[i][a] - B->matrix[i][a];
      }
    } else {
      error = ERROR_CAL;
    }
  } else {
    error = ERROR;
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = OK;
  if (!s21_check_matrix(A)) {
    error = ERROR;
  } else {
    s21_create_matrix(A->rows, A->columns, result);
    int rows = A->rows, columns = A->columns;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < columns; j++) {
        result->matrix[i][j] = number * A->matrix[i][j];
      }
    }
  }
  return error;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = (s21_check_matrix(A) && s21_check_matrix(B))
                ? (B->columns == A->rows || B->rows == A->columns) ? 0 : 2
                : 1;
  if (err == 0) {
    s21_create_matrix(A->rows, B->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += (A->matrix[i][k] * B->matrix[k][j]);
        }
      }
    }
  }
  return err;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = OK;

  if (s21_check_matrix(A)) {
    if (A->rows == 1 ||
        A->columns == 1) {  // проверяем, если одна из размерностей равна 1
      error = s21_create_matrix(A->columns, A->rows, result);
      if (!error) {
        for (int i = 0; i < A->rows; i++)
          for (int a = 0; a < A->columns; a++)
            result->matrix[a][i] = A->matrix[i][a];  // меняем индексы
      }
    } else {
      error = s21_create_matrix(A->columns, A->rows, result);
      if (!error) {
        for (int i = 0; i < A->rows; i++)
          for (int a = 0; a < A->columns; a++)
            result->matrix[a][i] = A->matrix[i][a];  // меняем индексы
      }
    }
  } else {
    error = ERROR;
  }
  return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A)) {
    if (A->columns == A->rows) {
      error = s21_create_matrix(A->rows, A->rows, result);
      if (A->columns == 1) {
        result->matrix[0][0] = A->matrix[0][0];  // Заменяем M на A
      } else {
        if (error == OK) {
          for (int m = 0; m < result->rows; m++) {
            for (int n = 0; n < result->columns; n++) {
              matrix_t tmp;
              s21_minor_creat(
                  m, n, A,
                  &tmp);  // Исправляем опечатку в функции s21_minor_create
              result->matrix[m][n] = s21_minor_calc(&tmp) * pow(-1, m + n);
              s21_remove_matrix(&tmp);
            }
          }
        }
      }
    } else {
      error = ERROR_CAL;
    }
  } else {
    error = ERROR;
  }
  return error;
}

void s21_minor_creat(int rows, int columns, matrix_t *A, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A)) {
    error = s21_create_matrix(A->rows - 1, A->columns - 1, result);
    int a = 0, b = 0;
    if (error == OK) {
      for (int m = 0; m < A->rows; m++)
        for (int n = 0; n < A->columns; n++)
          if (m != rows && n != columns) {
            result->matrix[a][b++] = A->matrix[m][n];  // заменяем M на A
            if (b == result->columns) a++, b = 0;
          }
    }
  }
}

double s21_minor_calc(matrix_t *A) {
  double result = 0;
  if (A->columns == A->rows) {
    if (A->columns == 1) {
      result = M[0][0];
    } else if (A->columns == 2) {
      result = M[0][0] * M[1][1] - M[0][1] * M[1][0];
    } else {
      for (int j = 0; j < A->columns; ++j) {
        matrix_t matrix_tmp = {NULL, 0, 0};
        s21_minor_creat(0, j, A, &matrix_tmp);
        result += M[0][j] * pow(-1, j) * s21_minor_calc(&matrix_tmp);
        s21_remove_matrix(&matrix_tmp);
      }
    }
  }
  return result;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = OK;
  if (s21_check_matrix(A)) {
    if (A->columns == A->rows) {
      *result = s21_minor_calc(A);
    } else {
      error = ERROR_CAL;
    }
  } else {
    error = ERROR;
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = 0;
  if (!s21_check_matrix(A)) {
    error = ERROR;
  } else {
    if (A->rows != A->columns) {
      error = ERROR_CAL;
    } else {
      double determinate = 0;
      error = s21_determinant(A, &determinate);
      if (determinate != 0) {
        if (!error) {
          matrix_t tmp_complements;
          matrix_t tmp_complements_trans;
          error = s21_calc_complements(A, &tmp_complements);
          error = s21_transpose(&tmp_complements, &tmp_complements_trans);
          error =
              s21_mult_number(&tmp_complements_trans, 1 / determinate, result);
          s21_remove_matrix(&tmp_complements);
          s21_remove_matrix(&tmp_complements_trans);
        }
      } else {
        error = ERROR_CAL;
      }
    }
  }
  return error;
}
