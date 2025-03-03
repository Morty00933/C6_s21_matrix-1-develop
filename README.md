# s21_matrix Library

Этот репозиторий содержит реализацию библиотеки для работы с матрицами на языке C. Библиотека включает в себя функции для создания, удаления, сравнения, арифметических операций, транспонирования, вычисления определителя, матрицы дополнений и нахождения обратной матрицы.

---

## Функциональные возможности

### Основные операции с матрицами
- **s21_create_matrix** – создание матрицы заданного размера и инициализация всех элементов нулями.
- **s21_remove_matrix** – освобождение памяти, занятой матрицей, и обнуление её размеров.
- **s21_check_matrix** – проверка корректности матрицы (наличие ненулевого указателя, положительных размеров).

### Арифметические операции
- **s21_eq_matrix** – проверка равенства двух матриц (с учетом заданной точности EPS).
- **s21_sum_matrix** – сложение двух матриц одинаковых размеров.
- **s21_sub_matrix** – вычитание одной матрицы из другой.

### Операции умножения
- **s21_mult_number** – умножение матрицы на скалярное число.
- **s21_mult_matrix** – умножение двух матриц (при корректном соответствии размерностей).

### Трансформация матрицы
- **s21_transpose** – транспонирование матрицы.

### Вычисление определителя и обратной матрицы
- **s21_determinant** – вычисление определителя матрицы.
- **s21_calc_complements** – вычисление матрицы дополнений (матрица миноров с учетом знаковых коэффициентов).
- **s21_minor_creat** и **s21_minor_calc** – создание подматрицы (минора) и рекурсивное вычисление её определителя.
- **s21_inverse_matrix** – нахождение обратной матрицы при ненулевом определителе.

---

## Структура проекта

- **s21_matrix.h**  
  Заголовочный файл, содержащий определения констант, типов (структура `matrix_t`) и прототипы всех функций работы с матрицами.

- **s21_matrix.c**  
  Исходный файл с реализациями функций для создания, обработки и арифметических операций над матрицами.

- **s21_matrixtest.c**  
  Набор модульных тестов, написанных с использованием фреймворка [Check](https://libcheck.github.io/check/), для проверки корректности работы реализованных функций.

- **Makefile**  
  Скрипт для сборки библиотеки и тестового модуля. (Файл Makefile должен быть настроен для компиляции проекта и запуска тестов.)

---

## Сборка и запуск тестов

### Требования
- Компилятор C (например, gcc)
- Библиотека [Check](https://libcheck.github.io/check/) для запуска модульных тестов

### Инструкции по сборке
1. Откройте терминал в корневой директории проекта.
2. Для сборки проекта выполните команду:
   ```bash
   make
   ```
   Это скомпилирует библиотеку и тестовый модуль.

### Запуск тестов
Запустить тесты можно с помощью:
```bash
make test
```
или запустив скомпилированный тестовый бинарный файл напрямую:
```bash
./s21_matrixtest
```

---

## Использование библиотеки

Чтобы использовать функции библиотеки в своём проекте:
1. Подключите заголовочный файл:
   ```c
   #include "s21_matrix.h"
   ```
2. При компиляции подключите исходный файл `s21_matrix.c` или соберите статическую/динамическую библиотеку.

Пример использования:
```c
#include "s21_matrix.h"
#include <stdio.h>

int main(void) {
    matrix_t A, B, R;
    // Создаем две матрицы 3x3
    if (s21_create_matrix(3, 3, &A) == OK && s21_create_matrix(3, 3, &B) == OK) {
        // Инициализируем матрицы
        A.matrix[0][0] = 1; A.matrix[0][1] = 2; A.matrix[0][2] = 3;
        A.matrix[1][0] = 4; A.matrix[1][1] = 5; A.matrix[1][2] = 6;
        A.matrix[2][0] = 7; A.matrix[2][1] = 8; A.matrix[2][2] = 9;
        B.matrix[0][0] = 9; B.matrix[0][1] = 8; B.matrix[0][2] = 7;
        B.matrix[1][0] = 6; B.matrix[1][1] = 5; B.matrix[1][2] = 4;
        B.matrix[2][0] = 3; B.matrix[2][1] = 2; B.matrix[2][2] = 1;
        
        // Складываем матрицы
        if (s21_sum_matrix(&A, &B, &R) == OK) {
            // Выводим результат (например, R.matrix[0][0])
            printf("R[0][0] = %lf\n", R.matrix[0][0]);
            s21_remove_matrix(&R);
        }
        s21_remove_matrix(&A);
        s21_remove_matrix(&B);
    }
    return 0;
}
```

---

## Авторские заметки

- Библиотека предназначена для образовательных целей и демонстрирует собственную реализацию базовых алгоритмов работы с матрицами.
- Тестовый модуль покрывает широкий спектр функционала: от проверки корректности операций сложения и вычитания до вычисления определителя и обратной матрицы.
- При использовании функций арифметических операций необходимо учитывать, что размеры матриц должны соответствовать требованиям (например, для умножения число столбцов первой матрицы должно равняться числу строк второй).

---

## Лицензия

Проект распространяется под лицензией MIT.

