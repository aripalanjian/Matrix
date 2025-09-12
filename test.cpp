#include "matrix.hpp"
#include <iostream>
#include <time.h>

void printDuration(timespec start, timespec finish){
    time_t computationTimeS;
    long computationTimeNS;

    if (start.tv_nsec > finish.tv_nsec){
        computationTimeS = finish.tv_sec - 1 - start.tv_sec;
        computationTimeNS = (long)1e9 + finish.tv_nsec - start.tv_nsec;
    } else {
        computationTimeS = finish.tv_sec - start.tv_sec;
        computationTimeNS = finish.tv_nsec - start.tv_nsec;
    }
    std::cout << "Computation time: " << computationTimeS << "." << computationTimeNS << "s\n";
}

void multiplicationTest(){
    std::cout << "Multiplication Test:\n";
    Matrix m1(3);
    m1.setDataRand();
    m1.print();

    Matrix m2(3,4);
    m2.setDataRand();
    m2.print();

    Matrix* res = m1 * m2;
    if(res){
        res->print();
        delete res;
    }

    Matrix m3(2,4);
    m3.setDataRand();
    m3.print();

    res = m1 * m3;
    if(res){
        res->print();
        delete res;
    }

}

void additionTest(){
    std::cout << "Addition Test:\n";
    Matrix m1(3);
    m1.setDataRand();
    m1.print();
    if(auto result = m1 + m1){
        result.value().print();
    }
}

void runTests(){
    additionTest();
    subtractionTest();
    multiplicationTest();
    divisionTest();
    determinantTest();
    transposeTest();

}

int main(){
    struct timespec start;
    struct timespec finish;
    
    clock_gettime(CLOCK_REALTIME, &start);
    multiplicationTest();
    additionTest();
    clock_gettime(CLOCK_REALTIME, &finish);

    printDuration(start,finish);

    return 0;
}

/*
TODO:
    incorporate setting default to rand or all 0s*/