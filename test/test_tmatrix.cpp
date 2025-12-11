#include "tmatrix.h"
#include <gtest.h>

TEST(TDynamicMatrix, can_create_matrix_with_positive_length)
{
    ASSERT_NO_THROW(TDynamicMatrix<int> m(5));
}

TEST(TDynamicMatrix, cant_create_too_large_matrix)
{
    ASSERT_ANY_THROW(TDynamicMatrix<int> m(MAX_MATRIX_SIZE + 1));
}

TEST(TDynamicMatrix, throws_when_create_matrix_with_negative_length)
{
    ASSERT_ANY_THROW(TDynamicMatrix<int> m(-5));
}

TEST(TDynamicMatrix, can_create_copied_matrix)
{
    TDynamicMatrix<int> m(5);
    ASSERT_NO_THROW(TDynamicMatrix<int> m1(m));
}

TEST(TDynamicMatrix, copied_matrix_is_equal_to_source_one)
{
    TDynamicMatrix<int> m(3);
    m.at(1).at(1) = 5; // безопасно
    TDynamicMatrix<int> m1(m);
    EXPECT_EQ(m, m1);
}

TEST(TDynamicMatrix, copied_matrix_has_its_own_memory)
{
    TDynamicMatrix<int> m(3);
    m.at(1).at(1) = 5;
    TDynamicMatrix<int> m1(m);
    m1.at(1).at(1) = 7;
    EXPECT_NE(m.at(1).at(1), m1.at(1).at(1));
}

TEST(TDynamicMatrix, can_get_size)
{
    TDynamicMatrix<int> m(4);
    EXPECT_EQ(4, m.size());
}

TEST(TDynamicMatrix, can_set_and_get_element)
{
    TDynamicMatrix<int> m(3);
    m.at(1).at(2) = 7; // безопасно
    EXPECT_EQ(m.at(1).at(2), 7);
}

TEST(TDynamicMatrix, throws_when_set_element_with_negative_index)
{
    TDynamicMatrix<int> m(3);
    EXPECT_ANY_THROW(m.at(-1).at(0));
}

TEST(TDynamicMatrix, throws_when_set_element_with_too_large_index)
{
    TDynamicMatrix<int> m(3);
    EXPECT_ANY_THROW(m.at(3).at(0));
}

TEST(TDynamicMatrix, can_assign_matrix_to_itself)
{
    TDynamicMatrix<int> m(3);
    m.at(1).at(1) = 9;
    m = m;
    EXPECT_EQ(m.at(1).at(1), 9);
}

TEST(TDynamicMatrix, can_assign_matrices_of_equal_size)
{
    TDynamicMatrix<int> m1(3);
    m1.at(1).at(1) = 4;
    TDynamicMatrix<int> m2(3);
    m2 = m1;
    EXPECT_EQ(m2.at(1).at(1), 4);
}

TEST(TDynamicMatrix, assign_operator_change_matrix_size)
{
    TDynamicMatrix<int> m1(3), m2(5);
    m2 = m1;
    EXPECT_EQ(m2.size(), 3);
}

TEST(TDynamicMatrix, can_assign_matrices_of_different_size)
{
    TDynamicMatrix<int> m1(3);
    m1.at(1).at(1) = 4;
    TDynamicMatrix<int> m2(5);
    m2 = m1;
    EXPECT_EQ(m2.size(), 3);
}
TEST(TDynamicMatrix, can_assign_matrices_of_different_size_content)
{
    TDynamicMatrix<int> m1(3);
    m1.at(1).at(1) = 4;
    TDynamicMatrix<int> m2(5);
    m2 = m1;
    EXPECT_EQ(m2.at(1).at(1), 4);
}

TEST(TDynamicMatrix, compare_equal_matrices_return_true)
{
    TDynamicMatrix<int> m1(3), m2(3);
    m1.at(0).at(0) = m2.at(0).at(0) = 5;
    EXPECT_TRUE(m1 == m2);
}

TEST(TDynamicMatrix, compare_matrix_with_itself_return_true)
{
    TDynamicMatrix<int> m(3);
    EXPECT_TRUE(m == m);
}

TEST(TDynamicMatrix, matrices_with_different_size_are_not_equal)
{
    TDynamicMatrix<int> m1(3), m2(4);
    EXPECT_NE(m1, m2);
}

TEST(TDynamicMatrix, can_add_matrices_with_equal_size)
{
    TDynamicMatrix<int> m1(2), m2(2);
    m1.at(0).at(0) = 1; m2.at(0).at(0) = 2;
    TDynamicMatrix<int> r = m1 + m2;
    EXPECT_EQ(r.at(0).at(0), 3);
}

TEST(TDynamicMatrix, cant_add_matrices_with_not_equal_size)
{
    TDynamicMatrix<int> m1(2), m2(3);
    EXPECT_ANY_THROW(m1 + m2);
}

TEST(TDynamicMatrix, can_subtract_matrices_with_equal_size)
{
    TDynamicMatrix<int> m1(2), m2(2);
    m1.at(0).at(0) = 3; m2.at(0).at(0) = 1;
    TDynamicMatrix<int> r = m1 - m2;
    EXPECT_EQ(r.at(0).at(0), 2);
}

TEST(TDynamicMatrix, cant_subtract_matrixes_with_not_equal_size)
{
    TDynamicMatrix<int> m1(2), m2(3);
    EXPECT_ANY_THROW(m1 - m2);
}