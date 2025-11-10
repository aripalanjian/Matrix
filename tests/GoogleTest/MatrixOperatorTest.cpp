//
// Created by Ari Palanjian on 11/5/25.
//

#include "../../matrix.hpp"
#include "gtest/gtest.h"

using Neo::Matrix;

class MatrixTest : public testing::Test {
    protected:
        MatrixTest() = default;
        ~MatrixTest() override = default;
        //https://stackoverflow.com/questions/4118025/how-can-i-make-a-constructor-which-lets-me-construct-with-a-braced-init-list
        Matrix<int> lhs{{{1, 2, 3, 4}, {4, 3, 2, 1}, {2, 2, 2, 2}}, 3, 4};
        Matrix<int> addSubRhs{{{1, 2, 3, 4}, {4, 3, 2, 1}, {2, 2, 2, 2}}, 3, 4};
        Matrix<int> mulDivRhs{{{1, 2, 3, 4}, {4, 3, 2, 1}, {2, 2, 2, 2},{1,1,1,1}}, 4, 4};

        Matrix<int> resAdd{{{2, 4, 6, 8}, {8, 6, 4, 2}, {4,4,4,4}}, 3, 4};
        Matrix<int> resSub{{{0,0,0,0}, {0,0,0,0}, {0,0,0,0}}, 3, 4};
        Matrix<int> resMul{{{19,18,17,16},{21,22,23,24},{16,16,16,16}}, 3, 4};
        // Matrix<int> resDiv{};
};

TEST_F(MatrixTest, Operators) {
    EXPECT_TRUE(*(lhs+addSubRhs) == resAdd);
    EXPECT_TRUE(*(lhs-addSubRhs) == resSub);
    EXPECT_TRUE(*(lhs*mulDivRhs) == resMul);
    // EXPECT_TRUE(*(lhs/mulDivRhs) == resDiv);
}


