#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <cstdio>
#include "DFN.hpp"
#include "Utils.hpp"

using namespace DFN_Library;

//Test con 2 fratture, 1 traccia passante

TEST(ImportFractures, ValidFile) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_1_file.txt";

    EXPECT_TRUE(ImportFractures(file_name, dfn, piano));
}

TEST(ImportFractures, NumberTraces) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_1_file.txt";
    ImportFractures(file_name, dfn, piano);

    EXPECT_EQ(dfn.NumberTraces, 1);
}

TEST(ImportFractures, NumberFractures) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_1_file.txt";
    ImportFractures(file_name, dfn, piano);

    EXPECT_EQ(dfn.NumberFractures, 2);
}

TEST(ImportFractures, Tips) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_1_file.txt";
    ImportFractures(file_name, dfn, piano);

    for (const auto& tip : dfn.Tips) {
        EXPECT_FALSE(tip.second);
    }

}

TEST(ImportFractures, Plane) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_1_file.txt";
    ImportFractures(file_name, dfn, piano);

    EXPECT_EQ(piano.Plane[0][0], 0);
    EXPECT_EQ(piano.Plane[0][1], 0);
    EXPECT_NE(piano.Plane[0][2], 0);
    EXPECT_EQ(piano.Plane[0][3], 0);
}

TEST(ImportFractures, Retta) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_1_file.txt";
    ImportFractures(file_name, dfn, piano);
    array<unsigned int, 2> coppia = {0,0};

    EXPECT_EQ(dfn.Retta[coppia][0][0], 0);
    EXPECT_NE(dfn.Retta[coppia][0][1], 0);
    EXPECT_EQ(dfn.Retta[coppia][0][2], 0);
    EXPECT_NE(dfn.Retta[coppia][1][0], 0);
    EXPECT_EQ(dfn.Retta[coppia][1][1], 0);
    EXPECT_EQ(dfn.Retta[coppia][1][2], 0);
}

//test 2 fratture, 1 traccia (passante per una sola frattura)

TEST(ImportFractures, Tips2) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_2_file.txt";
    ImportFractures(file_name, dfn, piano);
    array<unsigned int, 2> coppia1 = {0,0};
    array<unsigned int, 2> coppia2 = {0,1};

    EXPECT_TRUE(dfn.Tips[coppia1]);
    EXPECT_FALSE(dfn.Tips[coppia2]);
}

//test nessuna traccia

TEST(ImportFractures, ValidFile2) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_3_file.txt";

    EXPECT_TRUE(ImportFractures(file_name, dfn, piano));
}

TEST(ImportFractures, NumberTraces2) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_3_file.txt";
    ImportFractures(file_name, dfn, piano);

    EXPECT_EQ(dfn.NumberTraces, 0);
}

TEST(ImportFractures, NumberFractures2) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/test_3_file.txt";
    ImportFractures(file_name, dfn, piano);

    EXPECT_EQ(dfn.NumberFractures, 81);
}

//test file non presente

TEST(ImportFractures, FileDoesNotExist) {
    DFN dfn;
    Piano piano;

    std::string file_name = "Google_test/non_existing_file.txt";

    EXPECT_FALSE(ImportFractures(file_name, dfn, piano));
}




