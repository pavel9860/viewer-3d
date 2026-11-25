#include "harness.hpp"

// cadcore_tests [NAME_PART]: runs the test cases whose name contains NAME_PART (all without it).
int main(int argc, char** argv) { return harness::run_all(argc > 1 ? argv[1] : nullptr); }
