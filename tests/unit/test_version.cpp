#include <catch2/catch_test_macros.hpp>

#include "core/Version.h"

TEST_CASE("version is reported", "[core]") {
    REQUIRE(shipsim::version() == "2.0.0");
}
