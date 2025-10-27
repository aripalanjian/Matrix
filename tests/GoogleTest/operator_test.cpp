//
// Created by Ari Palanjian on 10/20/25.
//

#include <gtest/gtest.h>
#include "../../matrix.hpp"

TEST(Fundamental, IdentityMultiplication) {
    constexpr size_t m = 4;
    constexpr size_t n = 4;
    const int data[4][4] = {{1, 2, 3, 4}, {4, 3, 2, 1}, {2, 2, 2, 2}, {3, 3, 2, 1}};

    const auto data2 = new int *[m];

    for(size_t i = 0; i < m; i++) {
        data2[i] = new int[n];
        for(size_t j = 0; j < n; j++) {
            data2[i][j] = data[i][j];
        }
    }
    //Do Stuff
    Neo::Matrix<int> testMatrix{data2, m, n};
    const Neo::Matrix<int> IdentityMatrix{"I", n};
    EXPECT_TRUE(testMatrix.isSquare());
    EXPECT_TRUE(IdentityMatrix.isSquare());
    EXPECT_TRUE(testMatrix.getDims() == IdentityMatrix.getDims());
    auto tmp = testMatrix*IdentityMatrix;
    EXPECT_TRUE(tmp != std::nullopt);
    EXPECT_TRUE(*tmp == testMatrix);
    //Cleanup
    for (int i = 0; i < m; i++) {
        delete[] data2[i];
    }
    delete[] data2;
}

TEST(Fundamental, UnequalDims) {
    constexpr size_t m = 4;
    constexpr size_t n = 4;
    const int data[4][4] = {{1, 2, 3, 4}, {4, 3, 2, 1}, {2, 2, 2, 2}, {3, 3, 2, 1}};

    const auto data2 = new int *[m];

    for(size_t i = 0; i < m; i++) {
        data2[i] = new int[n];
        for(size_t j = 0; j < n; j++) {
            data2[i][j] = data[i][j];
        }
    }
    //Do Stuff
    Neo::Matrix<int> testMatrix{data2, m, n};
    Neo::Matrix<int> UnequalMatrix{3,n};
    EXPECT_FALSE(testMatrix.getDims() == UnequalMatrix.getDims());
    EXPECT_FALSE(testMatrix == UnequalMatrix);
    //Cleanup
    for (int i = 0; i < m; i++) {
        delete[] data2[i];
    }
    delete[] data2;
}
// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}


