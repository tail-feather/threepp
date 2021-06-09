#pragma once
#ifdef TEST_RUN_SINGLE_TEST
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif  // STELLA_TEST_RUN_SINGLE_TEST
