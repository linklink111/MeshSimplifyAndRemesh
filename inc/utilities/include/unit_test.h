#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#define TEST_SUITE_BEGIN \
    int main( int argc, char **argv ){ \
        bool TEST_SUITE_OK = true; \
        int TEST_SUITE_PASSED=0, TEST_SUITE_FAILED=0; \
        std::vector< bool (*)(void) > TEST_CASES;

#define TEST_SUITE_END \
    std::cout << "Running test suite with " << TEST_CASES.size() << " cases..." << std::endl; \
    for( int i=0; i<TEST_CASES.size(); i++ ){ \
        if( TEST_CASES[i]() ){ \
            TEST_SUITE_PASSED++; \
        } else { \
            TEST_SUITE_FAILED++; \
            TEST_SUITE_OK=false; \
        } \
    } \
    if( TEST_SUITE_OK ){ \
        std::cout << "[" << TEST_SUITE_PASSED << "/" << TEST_CASES.size() << "] test cases passed!" << std::endl; \
        return 0; \
    } else { \
        std::cout << "[" << TEST_SUITE_FAILED << "/" << TEST_CASES.size() << "] test cases failed!" << std::endl; \
        return 1; \
    } \
}

#define TEST_CASE_BEGIN( name ) { \
    auto test_case = []( void ) -> bool { \
        const char *TEST_CASE_NAME = #name; \
        bool TEST_CASE_OK = true; \
        std::cout << "  Running test case '" << TEST_CASE_NAME << "'...";

#define TEST_CASE_END \
        if( TEST_CASE_OK ) std::cout << "passed!" << std::endl; \
        return TEST_CASE_OK; \
    }; \
    TEST_CASES.push_back( test_case ); \
}

#define CHECK( A ) { \
    bool CHECK_TMP = (A); \
    if( !CHECK_TMP ){ \
        std::cout << std::endl << "    check(" << #A << ") " << (CHECK_TMP ? "passed!" : "failed!") << std::endl; \
        std::cout << "      " << __FILE__ << ", line: " << __LINE__ << std::endl; \
        std::cout << "      " << "Argument " << #A << " evaluated to false!" << std::endl; \
    } \
    TEST_CASE_OK &= CHECK_TMP; \
}

#define CHECK_EQUAL( A, B ) { \
    bool CHECK_EQUAL_TMP = (A) == (B); \
    if( !CHECK_EQUAL_TMP ){ \
        std::cout << std::endl << "    check_equal(" << #A << ", " << #B << ") " << (CHECK_EQUAL_TMP ? "passed!" : "failed!") << std::endl; \
        std::cout << "      " << __FILE__ << ", line: " << __LINE__ << std::endl; \
        std::cout << "      " << "Called values:"; \
        std::cout << " " << #A << "=" << (A); \
        std::cout << " " << #B << "=" << (B) << std::endl; \
    } \
    TEST_CASE_OK &= CHECK_EQUAL_TMP; \
}

#define CHECK_CLOSE( A, B, TOL ) {\
    bool CHECK_CLOSE_TMP = (A)-(B) <= (TOL) && (B)-(A) <= (TOL); \
    if( !CHECK_CLOSE_TMP ){ \
        std::cout << std::endl << "    check_close(" << #A << ", " << #B << ") " << (CHECK_CLOSE_TMP ? "passed!" : "failed!") << std::endl; \
        std::cout << "      " << __FILE__ << ", line: " << __LINE__ << std::endl; \
        std::cout << "      " << "Called values:"; \
        std::cout << " " << #A << "=" << (A); \
        std::cout << " " << #B << "=" << (B); \
        std::cout << " " << "tol=" << TOL << std::endl; \
    } \
    TEST_CASE_OK &= CHECK_CLOSE_TMP; \
}

#endif
